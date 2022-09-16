from random import sample

from pandas.core import base
# from . import exe_import
# from . import exe_qc
# from . import exe_to_bt
# from . import exe_correct
# from . import exe_mk_sval
from . import interstellar_setup
from . import settingImporter
from . import segmentImporter
import subprocess
import os
import re
import glob
import sys
import time
import pandas as pd

class UnknownError(Exception):
    pass


def genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,n_core):
    today_now=param_dict[os.path.basename(sampledir)]["today_now"]
    qoption=qcfg["QOPTION"]
    qoption=qoption.replace("<mem>",qcfg[mem_key])
    qoption=qoption.replace("<num_cores>",str(n_core))
    jid_now=cmd+today_now
    qcmd_base=["qsub","-e",sampledir+"/qlog","-o",sampledir+"/qlog","-cwd","-N",jid_now,qoption]
    return qcmd_base


def checkRequiredFile(key,flist):
    for f in flist:
        if re.search(key,f):
            return True
    return False


def writeStringInput(input_str,filepath):
    with open(filepath,mode="wt") as w:
        w.write(input_str)


def configRewrite(cfgpath,outdir,outnamedict):
    with open(cfgpath,mode="rt") as r:
        cfg_orig=[re.sub("\\t|\\n| ","",i) for i in r]
    cfg_new=[]
    
    for line in cfg_orig:
        line_split=line.split("=")
        outlist=[]
        if line_split[0] in outnamedict:
            if not re.search(r"^WHITELIST_ASSIGNMENT\(",line_split[1]):
                raise UnknownError("Barcode correspondence can only by applied to WHITELIST_ASSIGNMENT.")

            outnameprefix=outnamedict[line_split[0]] #-> bc_sort_1th
            destname=line_split[0].split(",")        #-> [dest1,dest2]

            for dest_each in destname:
                targetwhitelistpath=glob.glob(outdir+"/"+outnameprefix+"*"+dest_each+"*_sorted_whitelist.tsv")
                outlist+=targetwhitelistpath
                if len(targetwhitelistpath)>1:
                    raise UnknownError("More than 1 file were found in "+outdir)
            
            outlist=",".join(outlist)
            line=re.sub(r"(path:)[^,\)]+", "\\1"+outlist, line)
            line=re.sub(r"(correspondence_table:)[^,\)]+", "\\1"+"TRUE", line)
            if not "path:" in line:
                line=re.sub(r"\)$", ",path:"+outlist+")", line)
        cfg_new.append(line)

    with open(outdir+"/sorted.conf",mode="wt") as w:
        cfg_new="\n".join(cfg_new)+"\n"
        w.write(cfg_new)


def run(sampledir_list,cfg_raw,qcfg,is_qsub,is_multisample,param_dict,proj_dir,cfgpath):
    cfg=settingImporter.config_extract_value_trans(cfg_raw)
    cfg["bc_sort"]=False
    cfg_ext,dict_to_terminal=settingImporter.config_extract_value_ext(cfg_raw)
    func_dict=settingImporter.func_check_trans(cfg,dict_to_terminal)
    func_dict_ext=settingImporter.func_check(cfg_ext)

    d=settingImporter.getAllocation(func_dict_ext,cfg,cfg_ext)
    if is_multisample and cfg["bc_sort"]:
        raise UnknownError("You cannot execute multi-sample mode with barcode correspondence tables.")
    for key in cfg:
        if "READ1_STRUCTURE" in key or "READ2_STRUCTURE" in key or "INDEX1_STRUCTURE" in key or "INDEX2_STRUCTURE" in key:
            dest_segments=cfg[key].split("+")
            if cfg[key]:
                for i in dest_segments:
                    if re.search(r'^\"[^\"]+\"',i):
                        seq=re.sub('\"',"",i)
                        for nuc in seq:
                            if nuc not in ["A","T","G","C","N"]:
                                raise UnknownError("Only the nucleotides ATGCN can be used for constant sequences.")
                    else:
                        if not i in cfg["available_seg"]:
                            raise UnknownError("The segment "+i+" is not availbale for the read structure configuration.")

    cmd="buildTree"
    t = time.time()
    njobdict=dict()
    for sampledir in sampledir_list:
        file_endfix="_correct_srcValue.pkl"
        file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/mk_sval/*") if re.search(file_endfix,i)]
        if is_qsub:
            mem_key="mem_buildTree"
            # print("Running qsub jobs...: buildTree",flush=True)
            # qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"])
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,1) # single core if qsub
                
            for f in file_pool:
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                qcmd_now=qcmd_base+[sampledir+"/sh/buildTree.sh",outname_now,f]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Building tree", file=sys.stderr)
                    sys.exit(1)
            njobdict[sampledir]=len(file_pool)
        else:
            file_pool_concat = "\n".join(file_pool)
            outname_now="\n".join([re.sub(file_endfix,"",os.path.basename(f)) for f in file_pool])

            writeStringInput(file_pool_concat,sampledir+"/sh/srcValues_buildTree.txt")
            writeStringInput(outname_now,sampledir+"/sh/outnames_buildTree.txt")

            cmd_now=[sampledir+"/sh/buildTree.sh",sampledir+"/sh/outnames_buildTree.txt",sampledir+"/sh/srcValues_buildTree.txt"]
            cmd_now=" ".join(cmd_now)
            s=subprocess.run(cmd_now,shell=True)
            if s.returncode != 0:
                print("Job failed: Building tree", file=sys.stderr)
                sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Building tree",jid_now,sampledir+"/qlog",njobdict[sampledir])
    print("Elapsed time for building value trees",round(round(time.time() - t)/60,2),"minutes\n")


    cmd="mergeTree"
    t = time.time()
    file_pool = glob.glob(proj_dir+"/*/value_translation/_work/buildTree/*_Tree.pkl.gz")
    file_pool_concat = "\n".join(file_pool)
    writeStringInput(file_pool_concat, sampledir_list[0]+"/sh/valueTrees_mergeTree.txt")

    if is_qsub:
        mem_key="mem_mergeTree"
        # print("Running qsub jobs...: Merging trees",flush=True)
        qoption=qcfg["QOPTION"]
        qoption=qoption.replace("<mem>",qcfg[mem_key])
        qoption=qoption.replace("<num_cores>","1")
        # qoption=qoption.replace("<num_cores>",str(cfg_raw["general"]["NUM_CORES"])) # This process uses maximum amount of CPUs even in qsub mode. Config is fresh since it's improted from the raw file.
        qcmd_base=["qsub","-e",sampledir_list[0]+"/qlog","-o",sampledir_list[0]+"/qlog","-cwd","-N",cmd+param_dict[os.path.basename(sampledir_list[0])]["today_now"],qoption]
            
        outname_now="merge"
        qcmd_now=qcmd_base+[sampledir_list[0]+"/sh/mergeTree.sh",outname_now,sampledir_list[0]+"/sh/valueTrees_mergeTree.txt"]
        qcmd_now=" ".join(qcmd_now)
        s=subprocess.run(qcmd_now,shell=True)
        if s.returncode != 0:
            print("qsub failed: Merging tree", file=sys.stderr)
            sys.exit(1)
        njob=1
    else:
        outname_now="merge"
        cmd_now=[sampledir_list[0]+"/sh/mergeTree.sh",outname_now,sampledir_list[0]+"/sh/valueTrees_mergeTree.txt"]
        cmd_now=" ".join(cmd_now)
        s=subprocess.run(cmd_now,shell=True)
        if s.returncode != 0:
            print("Job failed: Merging tree", file=sys.stderr)
            sys.exit(1)
    if is_qsub:
        jid_now=cmd+param_dict[os.path.basename(sampledir_list[0])]["today_now"]
        interstellar_setup.job_wait("Merging tree",jid_now,sampledir_list[0]+"/qlog",njob)
    print("Elapsed time for merging value trees",round(round(time.time() - t)/60,2),"minutes\n")



    cmd="convert"
    t = time.time()
    njobdict=dict()
    if is_multisample:
        mergetree_dir=proj_dir+"/_multisample/mergeTree"

        # for i in sampledir_list:
        #     print(i,flush=True)
        #     print(mergetree_dir+"/merge_"+os.path.basename(i)+"_mergeTree.pkl.gz",flush=True)
        #     print(glob.glob(mergetree_dir+"/merge_"+os.path.basename(i)+"_mergeTree.pkl.gz"),flush=True)
        #     print("##",flush=True)


        sampledir_to_mergetree={k:glob.glob(mergetree_dir+"/merge_"+os.path.basename(k)+"_mergeTree.pkl.gz")[0] for k in sampledir_list}
    else:
        sampledir_to_mergetree={k:k+"/value_translation/_work/mergeTree/merge_mergeTree.pkl.gz" for k in sampledir_list}

    for sampledir in sampledir_list:
        file_endfix="_correct_srcValue.pkl"
        file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/mk_sval/*") if re.search(file_endfix,i)]
        mergetree=sampledir_to_mergetree[sampledir]
        
        if is_qsub:
            mem_key="mem_convert"
            # print("Running qsub jobs...: Value optimization",flush=True)
            # qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"])
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,1) # single core if qsub
                
            for f in file_pool:
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                qcmd_now=qcmd_base+[sampledir+"/sh/convert.sh",outname_now,mergetree,f,re.sub(file_endfix,"_correct_srcQual.pkl",f)]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Value optimization", file=sys.stderr)
                    sys.exit(1)
            njobdict[sampledir]=len(file_pool)
        else:
            file_pool_concat = "\n".join(file_pool)
            qual_file_pool_concat = "\n".join([sampledir+"/value_extraction/_work/mk_sval/"+re.sub(file_endfix,"",os.path.basename(f))+"_correct_srcQual.pkl" for f in file_pool])
            outname_now="\n".join([re.sub(file_endfix,"",os.path.basename(f)) for f in file_pool])

            writeStringInput(file_pool_concat, sampledir+"/sh/srcValue_convert.txt")
            writeStringInput(qual_file_pool_concat, sampledir+"/sh/srcQual_convert.txt")
            writeStringInput(outname_now, sampledir+"/sh/outname_convert.txt")

            cmd_now=[sampledir+"/sh/convert.sh",sampledir+"/sh/outname_convert.txt",mergetree,sampledir+"/sh/srcValue_convert.txt",sampledir+"/sh/srcQual_convert.txt"]
            cmd_now=" ".join(cmd_now)
            s=subprocess.run(cmd_now,shell=True)
            if s.returncode != 0:
                print("Job failed: Value optimization", file=sys.stderr)
                sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Value optimization",jid_now,sampledir+"/qlog",njobdict[sampledir])
    print("Elapsed time for generating optimized value tables",round(round(time.time() - t)/60,2),"minutes\n")


    # if cfg["bc_sort"]:
    #     cmd="bc_sort"
    #     # print("Runnning qsub jobs...: bc_sort",flush=True)
    #     mergetree=sampledir_list[0]+"/value_translation/_work/mergeTree/merge_mergeTree.pkl.gz"
    #     s2v=glob.glob(sampledir_list[0]+"/value_extraction/_work/mk_sval/*_sseq_to_svalue.pkl.gz")[0]
    #     njobs=0
    #     outnamedict=dict()
        
    #     for n,dest_seg in enumerate(cfg["val2table"]):
    #         tbl_now=cfg["val2table"][dest_seg]
    #         outname_now="bcsort_"+str(n)+"th"
    #         outnamedict[dest_seg]=outname_now
    #         if is_qsub:
    #             njobs+=1
    #             mem_key="mem_bc_sort"
    #             print("Running qsub jobs...: Barcode correspondence",flush=True)
    #             qoption=qcfg["QOPTION"]
    #             qoption=qoption.replace("<mem>",qcfg[mem_key])
    #             qcmd_base=["qsub","-e",sampledir_list[0]+"/qlog","-o",sampledir_list[0]+"/qlog","-cwd","-N",cmd+param_dict[os.path.basename(sampledir_list[0])]["today_now"],qoption]
                    
    #             qcmd_now=qcmd_base+[sampledir_list[0]+"/sh/bc_sort.sh",outname_now,mergetree,s2v,tbl_now]
    #             qcmd_now=" ".join(qcmd_now)
    #             s=subprocess.run(qcmd_now,shell=True)
    #             if s.returncode != 0:
    #                 print("qsub failed: Barcode correspondence", file=sys.stderr)
    #                 sys.exit(1)
    #         else:
    #             cmd_now=[sampledir_list[0]+"/sh/bc_sort.sh",outname_now,mergetree,s2v,tbl_now]
    #             cmd_now=" ".join(cmd_now)
    #             s=subprocess.run(cmd_now,shell=True)
    #             if s.returncode != 0:
    #                 print("Job failed: Barcode correspondence", file=sys.stderr)
    #                 sys.exit(1)
    #     if is_qsub:
    #         for sampledir in sampledir_list:
    #             jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
    #             interstellar_setup.job_wait("Barcode correspondence",jid_now,sampledir+"/qlog",njob)

    #     configRewrite(cfgpath,sampledir_list[0]+"/value_translation/_work/bc_sort",outnamedict)





    cmd="export"
    t = time.time()
    mem_key="mem_export"
    njobdict=dict()
    for sampledir in sampledir_list:
        file_endfix="_correct_srcValue.pkl"
        is_qc=checkRequiredFile("_srcSeq.QC.pkl",glob.glob(sampledir+"/value_extraction/_work/qc/*"))
        if is_qc:
            file_endfix="_srcSeq.QC.pkl"
            file_prefix=[os.path.basename(i).replace(file_endfix,"") for i in glob.glob(sampledir+"/value_extraction/_work/qc/*") if re.search(file_endfix,i)]
        else:
            file_endfix="_srcSeq.pkl"
            file_prefix=[os.path.basename(i.replace(file_endfix,"")) for i in glob.glob(sampledir+"/value_extraction/_work/import/*") if re.search(file_endfix,i)]
        njobdict[sampledir]=len(file_prefix)

        if is_qsub:
            for outname_now in file_prefix:
                dval =sampledir+"/value_translation/_work/convert/"+outname_now+"_converted_value.pkl"
                dqual=sampledir+"/value_translation/_work/convert/"+outname_now+"_converted_qual.pkl"
                sseq =sampledir+"/value_extraction/_work/mk_sval/"+outname_now+"_correct_result.pkl"
                squal=sampledir+"/value_extraction/_work/import/"+outname_now+"_srcQual.pkl"
                
                if is_multisample:
                    sizedict=proj_dir+"/_multisample/mergeTree/merge_"+os.path.basename(sampledir)+"_size_info.pkl.gz"
                else:
                    sizedict=sampledir+"/value_translation/_work/mergeTree/merge_size_info.pkl.gz"
                
                # qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,cfg_raw["general"]["NUM_CORES"])
                qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key,1) # single core if qsub
                if cfg["bc_sort"]:
                    qcmd_now=qcmd_base+[sampledir+"/sh/export_bc_sort.sh",outname_now,dval,dqual,sseq,squal,sizedict]
                else:
                    qcmd_now=qcmd_base+[sampledir+"/sh/export.sh",outname_now,dval,dqual,sseq,squal,sizedict]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Export", file=sys.stderr)
                    sys.exit(1)
        else:
            dval = "\n".join([sampledir+"/value_translation/_work/convert/"+p+"_converted_value.pkl" for p in file_prefix])
            dqual= "\n".join([sampledir+"/value_translation/_work/convert/"+p+"_converted_qual.pkl" for p in file_prefix])
            sseq = "\n".join([sampledir+"/value_extraction/_work/mk_sval/"+p+"_correct_result.pkl" for p in file_prefix])
            squal = "\n".join([sampledir+"/value_extraction/_work/import/"+p+"_srcQual.pkl" for p in file_prefix])
            outname_now = "\n".join(file_prefix)

            writeStringInput(dval, sampledir+"/sh/optVal_export.txt")
            writeStringInput(dqual, sampledir+"/sh/qualVal_export.txt")
            writeStringInput(sseq, sampledir+"/sh/srcSeq_export.txt")
            writeStringInput(squal, sampledir+"/sh/srcQual_export.txt")
            writeStringInput(outname_now, sampledir+"/sh/outnames_export.txt")
            
            if is_multisample:
                sizedict=proj_dir+"/_multisample/mergeTree/merge_"+os.path.basename(sampledir)+"_size_info.pkl.gz"
            else:
                sizedict=sampledir+"/value_translation/_work/mergeTree/merge_size_info.pkl.gz"

            if cfg["bc_sort"]:
                cmd_now=[sampledir+"/sh/export_bc_sort.sh",sampledir+"/sh/outnames_export.txt",sampledir+"/sh/optVal_export.txt",sampledir+"/sh/qualVal_export.txt",sampledir+"/sh/srcSeq_export.txt",sampledir+"/sh/srcQual_export.txt",sizedict]
            else:
                cmd_now=[sampledir+"/sh/export.sh",sampledir+"/sh/outnames_export.txt",sampledir+"/sh/optVal_export.txt",sampledir+"/sh/qualVal_export.txt",sampledir+"/sh/srcSeq_export.txt",sampledir+"/sh/srcQual_export.txt",sizedict]
            cmd_now=" ".join(cmd_now)
            s=subprocess.run(cmd_now,shell=True)
            if s.returncode != 0:
                print("Job failed: Export", file=sys.stderr)
                sys.exit(1)

    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Sequence export",jid_now,sampledir+"/qlog",njobdict[sampledir])

    for sampledir in sampledir_list:
        out_files=glob.glob(sampledir+"/value_translation/_work/export/*")
        key_list=["_R1.fastq.pkl","_R2.fastq.pkl","_I1.fastq.pkl","_I2.fastq.pkl"]
        
        for key in key_list:
            # target_files=[t for t in out_files if re.search(key+r"$",os.path.basename(t))]
            # target_files=sorted(target_files)
            # if len(target_files)>0:
            #     cmd=["echo"]+target_files+[" | xargs cat >",sampledir+"/value_translation/out/translated"+key]
            #     cmd=" ".join(cmd)
            #     s=subprocess.run(cmd,shell=True)
            #     if s.returncode != 0:
            #         print("Job failed: Sequence export", file=sys.stderr)
            #         sys.exit(1)
        
            target_files=[t for t in out_files if re.search(key+r"$",os.path.basename(t))]

            for idx,file_4set in enumerate(segmentImporter.split_yield_fastq_per_lines(target_files, n=4)):
                Output = pd.concat([pd.read_pickle(x) for x in file_4set])
                if idx == 0:
                    Output.to_csv(sampledir+"/value_translation/out/translated"+key.replace("pkl","gz"),mode="w",compression="gzip",sep="\t",index=False,header=False)
                else:
                    Output.to_csv(sampledir+"/value_translation/out/translated"+key.replace("pkl","gz"),mode="a",compression="gzip",sep="\t",index=False,header=False)
    print("Elapsed time for generating FASTQ files",round(round(time.time() - t)/60,2),"minutes\n")
