from random import sample

from pandas.core import base
from . import interstellar_setup
from . import settingImporter
import subprocess
import os
import re
import glob
import sys
import copy

class UnknownError(Exception):
    pass


def genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key):
    today_now=param_dict[os.path.basename(sampledir)]["today_now"]
    qoption=qcfg["QOPTION"]
    qoption=qoption.replace("<mem>",qcfg[mem_key])
    jid_now=cmd+today_now
    qcmd_base=["qsub","-e",sampledir+"/qlog","-o",sampledir+"/qlog","-cwd","-N",jid_now,qoption]
    return qcmd_base


def checkRequiredFile(key,flist):
    for f in flist:
        if re.search(key,f):
            return True
    return False


def run(sampledir_list,cfg_raw,qcfg,is_qsub,is_multisample,param_dict,proj_dir,cfgpath):
    cfg=settingImporter.config_extract_value_demulti(cfg_raw)
    cfg_ext=settingImporter.config_extract_value_ext(cfg_raw)
    
    for key in cfg:
        if "READ1_STRUCTURE" in key or "READ2_STRUCTURE" in key or "INDEX1_STRUCTURE" in key or "INDEX2_STRUCTURE" in key:
            dest_segments=cfg[key].split("+")
            if cfg[key]:
                for i in dest_segments:
                    if not i in cfg["available_seg"]:
                        raise UnknownError("The segment "+i+" is not availbale for the read structure configuration.")

    for segment in cfg["TARGET"].split(","):
        if not segment=="" and not segment in cfg_ext["value_segment"]:
            raise UnknownError("Segment "+segment+" is not defined in the config file.")
        
    cmd="demultiplex"
    njobdict=dict()
    for sampledir in sampledir_list:
        file_endfix="_correct_result.tsv.gz"
        file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/mk_sval/*") if re.search(file_endfix,i)]
        is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
        if is_qsub:
            mem_key="mem_"+cmd
            print("Running qsub jobs...: "+cmd,flush=True)
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key)
                
            for f in file_pool:
                if is_qc:
                    raw_qual=sampledir+"/value_extraction/_work/qc/"+re.sub(file_endfix,"_srcQual.QC.tsv.gz",os.path.basename(f))
                else:
                    raw_qual=sampledir+"/value_extraction/_work/import/"+re.sub(file_endfix,"_srcQual.tsv.gz",os.path.basename(f))
                correct_qual=os.path.dirname(f)+"/"+re.sub(file_endfix,"_correct_srcQual.tsv.gz",os.path.basename(f))
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                qcmd_now=qcmd_base+[sampledir+"/sh/"+cmd+".sh",outname_now,f,correct_qual,raw_qual]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Demultiplex', file=sys.stderr")
                    sys.exit(1)
            njobdict[sampledir]=len(file_pool)
            
        else:
            for f in file_pool:
                if is_qc:
                    raw_qual=sampledir+"/value_extraction/_work/qc/"+re.sub(file_endfix,"_srcQual.QC.tsv.gz",os.path.basename(f))
                else:
                    raw_qual=sampledir+"/value_extraction/_work/import/"+re.sub(file_endfix,"_srcQual.tsv.gz",os.path.basename(f))
                correct_qual=os.path.dirname(f)+"/"+re.sub(file_endfix,"_correct_srcQual.tsv.gz",os.path.basename(f))
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                cmd_now=[sampledir+"/sh/"+cmd+".sh",outname_now,f,correct_qual,raw_qual]
                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                if s.returncode != 0:
                    print("Job failed: Demultiplex', file=sys.stderr")
                    sys.exit(1)
   
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Demultiplex",jid_now,sampledir+"/qlog",njobdict[sampledir])
    
    #merging files
    for sampledir in sampledir_list:
        demulti_dirs=glob.glob(sampledir+"/demultiplex/_work/*")
        key_list=[]
        for demulti_chunk_dir in demulti_dirs:
            target_file=glob.glob(demulti_chunk_dir+"/*_demulti_key_list.txt")[0]
            with open(target_file, mode="r") as r:
                for line in r:
                    key_list.append(line.replace("\n",""))
        key_set=set(key_list)
        demulti_filename_list=[f for f in glob.glob(sampledir+"/demultiplex/_work/*/*")]
        for key in key_set:
            target_files=[t for t in demulti_filename_list if key in os.path.basename(t)]
            cmd=["cat"]+target_files+[">",sampledir+"/demultiplex/out/demultiplex"+key]
            cmd=" ".join(cmd)
            s=subprocess.run(cmd,shell=True)
            if s.returncode != 0:
                print("Job failed: Demultiplex file merge', file=sys.stderr")
                sys.exit(1)

