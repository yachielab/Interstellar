from posixpath import basename
from random import sample

from pandas.core import base
# from . import exe_import
# from . import exe_qc
# from . import exe_to_bt
# from . import exe_correct
# from . import exe_mk_sval
from . import interstellar_setup
from . import settingImporter
import subprocess
import os
import re
import glob
import sys
import time


def genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,n_core):
    today_now=param_dict[os.path.basename(sampledir)]["today_now"]
    qoption=qcfg["QOPTION"]
    qoption=qoption.replace("<mem>",qcfg[mem_key])
    qoption=qoption.replace("<num_cores>",str(n_core))
    jid_now=cmd+today_now
    qcmd_base=["qsub","-e",sampledir+"/qlog","-o",sampledir+"/qlog","-cwd","-N",jid_now,qoption]
    return qcmd_base


def run(sampledir_list,cfg_raw,qcfg,is_qsub,is_multisample,param_dict,proj_dir):
    used_commands=[]
    cfg,dict_to_terminal=settingImporter.config_extract_value_ext(cfg_raw)
    d=settingImporter.func_check(cfg)

    
    # import
    t = time.time()
    cmd="import"
    njobdict=dict()
    for sampledir in sampledir_list:
        endfix_input=param_dict[os.path.basename(sampledir)]["file_suffix"] 
        if is_qsub:
            mem_key="mem_import"
            # qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"])
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,1) # single core if qsub
            
            #collect input files
            file_pool=[]
            for r in ['read1','read2','index1','index2']:
                if r in param_dict[os.path.basename(sampledir)]["read_valid"] and not param_dict[os.path.basename(sampledir)]["read_valid"][r]=="":
                    file_pool.append([])

            all_files=glob.glob(sampledir+"/filesplit/*/*")
            for fileprefix in param_dict[os.path.basename(sampledir)]["target_prefix_list"]:
                idx=0
                for r in ['read1','read2','index1','index2']:
                    if r in param_dict[os.path.basename(sampledir)]["read_valid"] and not param_dict[os.path.basename(sampledir)]["read_valid"][r]=="":
                        target_files=[i for i in all_files if re.search("/"+r+"_"+fileprefix+r"_[0-9]+/.*"+endfix_input,i)]
                        file_pool[idx]+=sorted(target_files)
                        idx+=1
            # input_file_list=[self.input_read_files[i] for i in ['read1','read2','index1','index2'] if i in self.input_read_files]
            # basename_list=[os.path.basename(i) for i in input_file_list]
            
            # file_pool=[]
            # chunk_identifier_pool=os.path.basenaem(glob.glob(sampledir+"/filesplit/*"))
            # for f_name in basename_list:
            #     file_pool.append([sampledir+"/filesplit/"+f_name.replace(endfix_input,"")+i+"."+endfix_input for i in chunk_identifier_pool])
            # else:
            #     file_pool=[]
            #     if flg_file==1:
            #         for f_name in input_file_list:
            #             file_pool.append([f_name])
            #     else:
            #         prefix_pool=[re.sub(r"_[^_]+\."+endfix_input,"",os.path.basename(i)) for i in glob.glob(input_file_list[0]+"/*")]
            #         for f_name in input_file_list:
            #             file_pool.append([glob.glob(f_name+"/"+p+"*")[0] for p in prefix_pool])
            for infile in zip(*file_pool):
                # infile_now=[infile[i] for i in range(len(param_dict[os.path.basename(sampledir)]["read_valid"]))]
                outname_now=os.path.basename(infile[0].replace("."+endfix_input,""))
                qcmd_now=qcmd_base+[sampledir+"/sh/import.sh",outname_now]+list(infile)
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                used_commands.append(qcmd_now)
                if s.returncode != 0:
                    print("qsub failed: import", file=sys.stderr)
                    sys.exit(1)
            njobdict[sampledir]=len(file_pool[0])
        else:
            #collect input files
            file_pool=[]
            for r in ['read1','read2','index1','index2']:
                if r in param_dict[os.path.basename(sampledir)]["read_valid"] and not param_dict[os.path.basename(sampledir)]["read_valid"][r]=="":
                    file_pool.append([])

            for fileprefix in param_dict[os.path.basename(sampledir)]["target_prefix_list"]:
                idx=0
                for r in ['read1','read2','index1','index2']:
                    if r in param_dict[os.path.basename(sampledir)]["read_valid"] and not param_dict[os.path.basename(sampledir)]["read_valid"][r]=="":
                        tartget_file_dict=param_dict[os.path.basename(sampledir)]["tartget_file_dict"]
                        target_files=tartget_file_dict[r]
                        target_files=[i for i in target_files if re.search(fileprefix+r".*"+endfix_input,os.path.basename(i))]
                        file_pool[idx]+=sorted(target_files)
                        idx+=1

            for infile in zip(*file_pool):
                outname_now=os.path.basename(infile[0].replace("."+endfix_input,""))
                cmd_now=[sampledir+"/sh/import.sh",outname_now]+list(infile)
                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                used_commands.append(cmd_now)
                if s.returncode != 0:
                    print("script failed: import", file=sys.stderr)
                    sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Read segmentation",jid_now,sampledir+"/qlog",njobdict[sampledir])
    print("Elapsed time for Read segmentation",round(round(time.time() - t)/60,2),"minutes\n")


    #QC
    if "QUALITY_FILTER" in cfg["functions_used"]:
        t = time.time()
        cmd="qc"
        njobdict=dict()
        for sampledir in sampledir_list:
            if is_qsub:
                mem_key="mem_qc"
                # qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"])
                qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,1) # single core if qsub
                file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/import/*") if re.search("_srcSeq.tsv.gz",i)]
                for f in file_pool:
                    outname_now=re.sub(r"_srcSeq\.tsv\.gz$","",os.path.basename(f))
                    qcmd_now=qcmd_base+[sampledir+"/sh/qc.sh",outname_now,f,re.sub(r"_srcSeq\.","_srcQual.",f)]
                    qcmd_now=" ".join(qcmd_now)
                    s=subprocess.run(qcmd_now,shell=True)
                    used_commands.append(qcmd_now)
                    if s.returncode != 0:
                        print("qsub failed: Quality filtering", file=sys.stderr)
                        sys.exit(1)
                njobdict[sampledir]=len(file_pool)
            else:
                seq_file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/import/*") if re.search("_srcSeq.tsv.gz",i)]
                seq_file_pool_concat = ",".join(seq_file_pool)
                qual_file_pool_concat = re.sub(r"_srcSeq\.","_srcQual.",seq_file_pool_concat)
                outname_now=",".join([re.sub(r"_srcSeq\.tsv\.gz$","",os.path.basename(f)) for f in seq_file_pool])

                cmd_now=[sampledir+"/sh/qc.sh", outname_now, seq_file_pool_concat, qual_file_pool_concat]
                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                used_commands.append(cmd_now)
                if s.returncode != 0:
                    print("Job failed: Quality filtering", file=sys.stderr)
                    sys.exit(1)
                
        if is_qsub:
            for sampledir in sampledir_list:
                jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
                interstellar_setup.job_wait("Quality filtering",jid_now,sampledir+"/qlog",njobdict[sampledir])
        print("Elapsed time for quality filtering",round(round(time.time() - t)/60,2),"minutes\n")


    #to_bt
    if "BARTENDER_CORRECTION" in cfg["functions_used"]:
        t = time.time()
        cmd="to_bt"
        outname_now="to_bt"
        for sampledir in sampledir_list:
            if is_qsub:
                mem_key="mem_to_bt"
                # qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"])
                qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,1) # single core if qsub
                is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.pkl",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
                if is_qc:
                    qcmd_now=qcmd_base+[sampledir+"/sh/to_bt.sh",outname_now,'"'+sampledir+"/value_extraction/_work/qc/*_srcSeq.QC.pkl"+'"']                        
                else:
                    qcmd_now=qcmd_base+[sampledir+"/sh/to_bt.sh",outname_now,'"'+sampledir+"/value_extraction/_work/import/*_srcSeq.pkl"+'"']

                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                used_commands.append(qcmd_now)
                if s.returncode != 0:
                    print("qsub failed: Bartender", file=sys.stderr)
                    sys.exit(1)
            else:
                is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.pkl",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
                if is_qc:
                    cmd_now=[sampledir+"/sh/to_bt.sh",outname_now,'"'+sampledir+"/value_extraction/_work/qc/*_srcSeq.QC.pkl"+'"']                        
                else:
                    cmd_now=[sampledir+"/sh/to_bt.sh",outname_now,'"'+sampledir+"/value_extraction/_work/import/*_srcSeq.pkl"+'"']

                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                used_commands.append(cmd_now)
                if s.returncode != 0:
                    print("Job: Bartender", file=sys.stderr)
                    sys.exit(1)
        if is_qsub:    
            njobs=1
            for sampledir in sampledir_list:
                jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
                interstellar_setup.job_wait("Bartender",jid_now,sampledir+"/qlog",njobs)
        print("Elapsed time for Bartender",round(round(time.time() - t)/60,2),"minutes\n")


    #correct
    t = time.time()
    cmd="correct"
    outname_now="Interstellar"
    for sampledir in sampledir_list:
        if is_qsub:
            mem_key="mem_correct"
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"]) # This process uses maximum amount of CPUs even in qsub mode. Config is fresh since it's improted from the raw file.
            is_qc=interstellar_setup.checkRequiredFile("_srcCount.QC.pkl.gz",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
            if is_qc:
                qcmd_now=qcmd_base+[sampledir+"/sh/correct.sh",outname_now,'"'+sampledir+"/value_extraction/_work/qc/*_srcCount.QC.pkl.gz"+'"']
            else:
                qcmd_now=qcmd_base+[sampledir+"/sh/correct.sh",outname_now,'"'+sampledir+"/value_extraction/_work/import/*_srcCount.pkl.gz"+'"']

            qcmd_now=" ".join(qcmd_now)
            s=subprocess.run(qcmd_now,shell=True)
            used_commands.append(qcmd_now)
            if s.returncode != 0:
                print("qsub failed: Sequence error correction", file=sys.stderr)
                sys.exit(1)
        else:
            is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
            if is_qc:
                cmd_now=[sampledir+"/sh/correct.sh",outname_now,'"'+sampledir+"/value_extraction/_work/qc/*_srcCount.QC.pkl.gz"+'"']                        
            else:
                cmd_now=[sampledir+"/sh/correct.sh",outname_now,'"'+sampledir+"/value_extraction/_work/import/*_srcCount.pkl.gz"+'"']

            cmd_now=" ".join(cmd_now)
            s=subprocess.run(cmd_now,shell=True)
            used_commands.append(cmd_now)
            if s.returncode != 0:
                print("Job failed: Sequence error correction", file=sys.stderr)
                sys.exit(1)
    if is_qsub:    
        njobs=1
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Sequence error correction",jid_now,sampledir,njobs)
    print("Elapsed time for sequence correction",round(round(time.time() - t)/60,2),"minutes\n")


    
    #mk_sval
    t = time.time()
    cmd="mk_sval"
    njobdict=dict()
    for sampledir in sampledir_list:
        if is_qsub:
            mem_key="mem_mk_sval"
            # qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"])
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,1) # single core if qsub
            qcmd_base=qcmd_base+[sampledir+"/sh/mk_sval.sh"]

            is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.pkl",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
            if is_qc:
                file_endfix="_srcSeq.QC.pkl"
                file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/qc/*") if re.search("_srcSeq.QC.pkl",i)]              
            else:
                file_endfix="_srcSeq.pkl"
                file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/import/*") if re.search("_srcSeq.pkl",i)]
                
            for f in file_pool:
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                qual_file_path = sampledir+"/value_extraction/_work/import/"+outname_now+"_srcQual.pkl"
                qcmd_now=qcmd_base+[outname_now,f,qual_file_path,sampledir+"/value_extraction/_work/correct/Interstellar_srcCorrect.pkl.gz"]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                used_commands.append(qcmd_now)
                if s.returncode != 0:
                    print("qsub failed: Sequence to value conversion", file=sys.stderr)
                    sys.exit(1)

            njobdict[sampledir]=len(file_pool)
        else:
            is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.pkl",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
            if is_qc:
                file_endfix="_srcSeq.QC.pkl"
                file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/qc/*") if re.search("_srcSeq.QC.pkl",i)]              
            else:
                file_endfix="_srcSeq.pkl"
                file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/import/*") if re.search("_srcSeq.pkl",i)]

            seq_file_pool_concat = ",".join(file_pool)
            qual_file_pool_concat = ",".join([sampledir+"/value_extraction/_work/import/"+re.sub(file_endfix,"",os.path.basename(f))+"_srcQual.pkl" for f in file_pool])
            outname_now = ",".join([re.sub(file_endfix,"",os.path.basename(f)) for f in file_pool])
            cmd_now = [sampledir+"/sh/mk_sval.sh", outname_now, seq_file_pool_concat, qual_file_pool_concat, sampledir+"/value_extraction/_work/correct/Interstellar_srcCorrect.pkl.gz"]
            cmd_now=" ".join(cmd_now)
            s=subprocess.run(cmd_now,shell=True)
            used_commands.append(cmd_now)
            if s.returncode != 0:
                print("Job failed: Sequence to value conversion", file=sys.stderr)
                sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Sequence to value conversion",jid_now,sampledir+"/qlog",njobdict[sampledir])
    print("Elapsed time for value table generation",round(round(time.time() - t)/60,2),"minutes\n")


    #Generating a samplesheet
    if is_multisample:
        valuefile_list=[]
        for sampledir in sampledir_list:
            valuefile_list_tmp=[[i,os.path.basename(sampledir)] for i in glob.glob(sampledir+"/value_extraction/_work/mk_sval/*_correct_srcValue.pkl")]
            valuefile_list+=valuefile_list_tmp
        valuefile_list=["\t".join(i) for i in valuefile_list]
        samplesheet_table="\n".join(valuefile_list)+"\n"
        with open(proj_dir+"/_multisample/samplesheet/samplesheet.tsv",mode="wt") as w:
            w.write(samplesheet_table)
        
    
    for sampledir in sampledir_list:
        out_files=glob.glob(sampledir+"/value_extraction/_work/mk_sval/*")
        key_list=["_correct_result.tsv.gz"]
        
        for key in key_list:
            target_files=[t for t in out_files if re.search(key+r"$",os.path.basename(t))]
            if len(target_files)>0:
                cmd1=["cat"]+[target_files[0]]+[" | gunzip -c | head -n1 >",sampledir+"/value_extraction/out/corrected_table.header"]
                cmd1=" ".join(cmd1)
                cmd2=["echo"]+target_files+["| xargs cat | zgrep -v Header >",sampledir+"/value_extraction/out/corrected_table.content"]
                cmd2=" ".join(cmd2)
                cmd3=["cat",sampledir+"/value_extraction/out/corrected_table.header",sampledir+"/value_extraction/out/corrected_table.content | gzip -c > ",sampledir+"/value_extraction/out/corrected_table.tsv.gz"]
                cmd3=" ".join(cmd3)
                cmd4=["rm",sampledir+"/value_extraction/out/corrected_table.header",sampledir+"/value_extraction/out/corrected_table.content"]
                cmd4=" ".join(cmd4)
                for cmd in [cmd1,cmd2,cmd3,cmd4]:
                    s=subprocess.run(cmd,shell=True)
                    if s.returncode != 0:
                        print("Job failed: Tagged file merge", file=sys.stderr)
                        sys.exit(1)
