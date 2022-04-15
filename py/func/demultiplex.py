from random import sample

from pandas.core import base
from . import interstellar_setup
from . import settingImporter
import subprocess
import os
import re
import glob
import sys
import shutil

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
    shell_template=cfg_raw["general"]["SET_SHELL_ENV"]
    cfg=settingImporter.config_extract_value_demulti(cfg_raw)
    cfg_ext,dict_to_terminal=settingImporter.config_extract_value_ext(cfg_raw)
    
    for key in cfg:
        if "READ1_STRUCTURE" in key or "READ2_STRUCTURE" in key or "INDEX1_STRUCTURE" in key or "INDEX2_STRUCTURE" in key:
            dest_segments=cfg[key].split("+")
            if cfg[key]:
                for i in dest_segments:
                    if not i in cfg["available_seg"] and not re.search(r'^\"[^\"]+\"',i):
                        raise UnknownError("The segment "+i+" is not availbale for the read structure configuration.")
                    elif re.search(r'^\"[^\"]+\"',i):
                        seq=re.sub('\"',"",i)
                        for char in seq:
                            if char not in "ATGCN":
                                raise UnknownError("Only the nucleotides ATGCN can be used for constant sequences.")


    for segment in cfg["TARGET"].split(","):
        if not segment=="" and not segment in cfg_ext["value_segment"]:
            raise UnknownError("Segment "+segment+" is not defined in the config file.")
        
    cmd="demultiplex"
    njobdict=dict()
    for sampledir in sampledir_list:
        L_tmp=glob.glob(sampledir+"/demultiplex/_work/*")
        for l in L_tmp:
            shutil.rmtree(l) #clearance

        file_endfix="_correct_result.tsv.gz"
        file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/mk_sval/*") if re.search(file_endfix,i)]
        is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
        if is_qsub:
            mem_key="mem_"+cmd
            print("Running qsub jobs...: Demultiplex",flush=True)
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
    print("Merging files...",flush=True)
    njobdict=dict()
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

        njobdict[sampledir]=len(key_set)
        if is_qsub:
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,"demulti_filemerge",mem_key)
        for key in key_set:
            target_files=[t for t in demulti_filename_list if key in os.path.basename(t)]
            if cfg["FORMAT"]=="tsv":
                key_gunzip=key.replace(".gz","")
                cmd1=["gunzip -c"]+[target_files[0]]+["| head -n1 >",sampledir+"/demultiplex/out/demultiplex.header"+key_gunzip]
                cmd1=" ".join(cmd1)
                cmd2=["echo"]+target_files+["| xargs cat | zgrep -v Header >",sampledir+"/demultiplex/out/demultiplex.content"+key_gunzip]
                cmd2=" ".join(cmd2)
                cmd3=["cat",sampledir+"/demultiplex/out/demultiplex.header"+key_gunzip,sampledir+"/demultiplex/out/demultiplex.content"+key_gunzip+" | gzip -c > ",sampledir+"/demultiplex/out/demultiplex"+key]
                cmd3=" ".join(cmd3)
                cmd4=["rm",sampledir+"/demultiplex/out/demultiplex.header"+key_gunzip,sampledir+"/demultiplex/out/demultiplex.content"+key_gunzip]
                cmd4=" ".join(cmd4)
                if not is_qsub:
                    for cmd in [cmd1,cmd2,cmd3,cmd4]:
                        s=subprocess.run(cmd,shell=True)
                        if s.returncode != 0:
                            print("Job failed: Demultiuplexed file merge', file=sys.stderr")
                            sys.exit(1)
                else:
                    shell_lines=[]
                    with open(shell_template,mode="rt") as r:
                        for line in r:
                            line=line.replace("\n","")
                            shell_lines.append(line)
                    shell_lines+=[cmd1,cmd2,cmd3,cmd4]
                    shell_lines="\n".join(shell_lines)+"\n"
                    with open(sampledir+"/sh/demulti_filemerge.sh",mode="wt") as w:
                        w.write(shell_lines)

                    qcmd_now=qcmd_base+[sampledir+"/sh/demulti_filemerge.sh"]
                    qcmd_now=" ".join(qcmd_now)
                    s=subprocess.run(qcmd_now,shell=True)
                    if s.returncode != 0:
                        print("qsub failed: Demultiuplexed file merge', file=sys.stderr")
                        sys.exit(1)
                
            else:
                cmd=["echo"]+target_files+["| xargs cat >",sampledir+"/demultiplex/out/demultiplex"+key]
                cmd=" ".join(cmd)
                
                if not is_qsub:
                    s=subprocess.run(cmd,shell=True)
                    if s.returncode != 0:
                        print("Job failed: Demultiuplexed file merge', file=sys.stderr")
                        sys.exit(1)
                else:
                    shell_lines=[]
                    with open(shell_template,mode="rt") as r:
                        for line in r:
                            line=line.replace("\n","")
                            shell_lines.append(line)
                    shell_lines+=[cmd]
                    shell_lines="\n".join(shell_lines)+"\n"
                    with open(sampledir+"/sh/demulti_filemerge.sh",mode="wt") as w:
                        w.write(shell_lines)

                    qcmd_now=qcmd_base+[sampledir+"/sh/demulti_filemerge.sh"]
                    qcmd_now=" ".join(qcmd_now)
                    s=subprocess.run(qcmd_now,shell=True)
                    if s.returncode != 0:
                        print("qsub failed: Demultiuplexed file merge', file=sys.stderr")
                        sys.exit(1)
                
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now="demulti_filemerge"+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Demultiplex",jid_now,sampledir+"/qlog",njobdict[sampledir])