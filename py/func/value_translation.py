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


def configRewrite(cfgpath,outdir,outnamedict,func_dict):
    with open(cfgpath,mode="rt") as r:
        cfg_orig=[re.sub("\\t|\\n| ","",i) for i in r]
    cfg_new=[]
    
    for line in enumerate(cfg_orig):
        line_split=line.split("=")
        outlist=[]
        if line_split[0] in outnamedict:
            if not re.search("^WHITELIST_ASSIGNMENT(",line_split[1]):
                raise UnknownError("Barcode correspondence can only by applied to WHITELIST_ASSIGNMENT.")

            outnameprefix=outnamedict[line_split[0]] #-> bc_sort_1th
            destname=line_split[0].split(",")        #-> [dest1,dest2]
            for dest_each in enumerate(destname):
                targetwhitelistpath=glob.glob(outdir+"/"+outnameprefix+"*"+dest_each+"*_sorted_whitelist.tsv")
                outlist+=targetwhitelistpath
                if len(targetwhitelistpath)>1:
                    raise UnknownError("More than 1 file were found in "+outdir)
            
            outlist=",".join(outlist)
            line=re.sub("(path:)[^,\)]+", "\\1"+outlist, line)
            line=re.sub("(correspondence_table:)[^,\)]+", "\\1"+"", line)
        cfg_new.append(line)

    with open(outdir+"/sorted.conf",mode="wt") as w:
        cfg_new="\n".join(cfg_new)+"\n"
        w.write(cfg_new)


def run(sampledir_list,cfg_raw,qcfg,is_qsub,is_multisample,param_dict,proj_dir,cfgpath):
    cfg=settingImporter.config_extract_value_trans(cfg_raw)
    cfg_ext=settingImporter.config_extract_value_ext(cfg_raw)
    func_dict=settingImporter.func_check_trans(cfg)
    func_dict_ext=settingImporter.func_check(cfg_ext)

    d=settingImporter.getAllocation(func_dict_ext,cfg,cfg_ext)
    if is_multisample and cfg["bc_sort"]:
        raise UnknownError("You cannot execute multi-sample mode with barcode correspondence tables.")
    for key in cfg:
        if "READ1_STRUCTURE" in key or "READ2_STRUCTURE" in key or "INDEX1_STRUCTURE" in key or "INDEX2_STRUCTURE" in key:
            dest_segments=cfg[key].split("+")
            if cfg[key]:
                for i in dest_segments:
                    if not i in cfg["available_seg"]:
                        raise UnknownError("The segment "+i+" is not availbale for the read structure configuration.")

    

    cmd="buildTree"
    njobdict=dict()
    for sampledir in sampledir_list:
        file_endfix="_correct_srcValue.tsv.gz"
        file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/mk_sval/*") if re.search(file_endfix,i)]
        if is_qsub:
            mem_key="mem_buildTree"
            print("Running qsub jobs...: buildTree",flush=True)
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key)
                
            for f in file_pool:
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                qcmd_now=qcmd_base+[sampledir+"/sh/buildTree.sh",outname_now,f]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Building tree', file=sys.stderr")
                    sys.exit(1)
            njobdict[sampledir]=len(file_pool)
        else:
            for f in file_pool:
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                cmd_now=[sampledir+"/sh/buildTree.sh",outname_now,f]
                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                if s.returncode != 0:
                    print("Job failed: Building tree', file=sys.stderr")
                    sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Building tree",jid_now,sampledir+"/qlog",njobdict[sampledir])



    cmd="mergeTree"    
    if is_qsub:
        mem_key="mem_mergeTree"
        print("Running qsub jobs...: Merging trees",flush=True)
        qoption=qcfg["QOPTION"]
        qoption=qoption.replace("<mem>",qcfg[mem_key])
        qcmd_base=["qsub","-e",sampledir_list[0]+"/qlog","-o",sampledir_list[0]+"/qlog","-cwd","-N",cmd,qoption]
            
        outname_now="merge"
        qcmd_now=qcmd_base+[sampledir_list[0]+"/sh/mergeTree.sh",outname_now,'"'+proj_dir+"/*/value_translation/_work/buildTree/*_Tree.pkl.gz"+'"']
        qcmd_now=" ".join(qcmd_now)
        s=subprocess.run(qcmd_now,shell=True)
        if s.returncode != 0:
            print("qsub failed: Building tree', file=sys.stderr")
            sys.exit(1)
        njob=1
    else:
        outname_now="merge"
        cmd_now=[sampledir_list[0]+"/sh/mergeTree.sh",outname_now,'"'+proj_dir+"/*/value_translation/_work/buildTree/*_Tree.pkl.gz"+'"']
        cmd_now=" ".join(cmd_now)
        s=subprocess.run(cmd_now,shell=True)
        if s.returncode != 0:
            print("Job failed: Building tree', file=sys.stderr")
            sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Building tree",jid_now,sampledir+"/qlog",njob)



    cmd="convert"
    njobdict=dict()
    print("Runnning qsub jobs...: convert",flush=True)
    if is_multisample:
        mergetree_dir=proj_dir+"/_multisample/mergeTree"
        sampledir_to_mergetree={k:glob.glob(mergetree_dir+"/merge_"+os.path.basename(k)+"_mergeTree.pkl.gz")[0] for k in sampledir_list}
    else:
        sampledir_to_mergetree={k:k+"/value_translation/_work/mergeTree/merge_mergeTree.pkl.gz" for k in sampledir_list}

    for sampledir in sampledir_list:
        file_endfix="_correct_srcValue.tsv.gz"
        file_pool=[i for i in glob.glob(sampledir+"/value_extraction/_work/mk_sval/*") if re.search(file_endfix,i)]
        mergetree=sampledir_to_mergetree[sampledir]
        
        if is_qsub:
            mem_key="mem_convert"
            print("Running qsub jobs...: Value optimization",flush=True)
            qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key)
                
            for f in file_pool:
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                qcmd_now=qcmd_base+[sampledir+"/sh/convert.sh",outname_now,mergetree,f,re.sub(file_endfix,"_correct_srcQual.tsv.gz",f)]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Value optimization', file=sys.stderr")
                    sys.exit(1)
            njobdict[sampledir]=len(file_pool)
        else:
            for f in file_pool:
                outname_now=re.sub(file_endfix,"",os.path.basename(f))
                cmd_now=[sampledir+"/sh/convert.sh",outname_now,mergetree,f,re.sub(file_endfix,"_correct_srcQual.tsv.gz",f)]
                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                if s.returncode != 0:
                    print("Job failed: Value optimization', file=sys.stderr")
                    sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Value optimization",jid_now,sampledir+"/qlog",njobdict[sampledir])



    if cfg["bc_sort"]:
        cmd="bc_sort"
        print("Runnning qsub jobs...: bc_sort",flush=True)
        mergetree=sampledir_list[0]+"/value_translation/_work/mergeTree/merge_mergeTree.pkl.gz"
        s2v=glob.glob(sampledir_list[0]+"/value_translation/_work/mk_sval/*_sseq_to_svalue.pkl.gz")[0]
        njobs=0
        outnamedict=dict()
        
        for n,dest_seg in enumerate(cfg["val2table"]):
            tbl_now=cfg["val2table"][dest_seg]
            outname_now="bcsort_"+str(n)+"th"
            outnamedict[dest_seg]=outname_now
            if is_qsub:
                njobs+=1
                mem_key="mem_bc_sort"
                print("Running qsub jobs...: Barcode correspondence",flush=True)
                qoption=qcfg["QOPTION"]
                qoption=qoption.replace("<mem>",qcfg[mem_key])
                qcmd_base=["qsub","-e",sampledir_list[0]+"/qlog","-o",sampledir_list[0]+"/qlog","-cwd","-N",cmd,qoption]
                    
                qcmd_now=qcmd_base+[sampledir_list[0]+"/sh/bc_sort.sh",outname_now,mergetree,s2v,tbl_now]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Building tree', file=sys.stderr")
                    sys.exit(1)
            else:
                cmd_now=[sampledir_list[0]+"/sh/bc_sort.sh",outname_now,mergetree,s2v,tbl_now]
                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                if s.returncode != 0:
                    print("Job failed: Building tree', file=sys.stderr")
                    sys.exit(1)
        if is_qsub:
            for sampledir in sampledir_list:
                jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
                interstellar_setup.job_wait("Building tree",jid_now,sampledir+"/qlog",njob)

        configRewrite(cfgpath,sampledir_list[0]+"/value_translation/_work/bc_sort",outnamedict)





    cmd="export"
    mem_key="mem_export"
    njobdict=dict()
    for sampledir in sampledir_list:
        file_endfix="_correct_srcValue.tsv.gz"
        is_qc=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(sampledir+"/value_extraction/_work/*"))
        if is_qc:
            file_endfix="_srcSeq.QC.tsv.gz"
            file_prefix=[os.path.basename(i).replace(file_endfix,"") for i in glob.glob(sampledir+"/value_extraction/_work/qc/*") if re.search(file_endfix,i)]
        else:
            file_endfix="_srcSeq.tsv.gz"
            file_prefix=[os.path.basename(i.replace(file_endfix,"")) for i in glob.glob(sampledir+"/value_extraction/_work/import/*") if re.search(file_endfix,i)]
        njobdict[sampledir]=len(file_prefix)

        for outname_now in file_prefix:
            dval =sampledir+"/value_translation/_work/convert/"+outname_now+"_converted_value.tsv.gz"
            dqual=sampledir+"/value_translation/_work/convert/"+outname_now+"_converted_qual.tsv.gz"
            if is_qc:
                sseq =sampledir+"/value_extraction/_work/qc/"+outname_now+file_endfix
                squal=sampledir+"/value_extraction/_work/qc/"+outname_now+"_srcQual.QC.tsv.gz"
            else:
                sseq =sampledir+"/value_extraction/_work/import/"+outname_now+file_endfix
                squal=sampledir+"/value_extraction/_work/import/"+outname_now+"_srcQual.tsv.gz"
            
            if is_multisample:
                sizedict=proj_dir+"/_multisample/mergeTree/merge_"+os.path.basename(sampledir)+"_size_info.pkl.gz"
            else:
                sizedict=sampledir+"/value_translation/_work/mergeTree/merge_size_info.pkl.gz"
            
            if is_qsub:
                qcmd_base=genCmdBase(param_dict,sampledir,qcfg,cmd,mem_key)
                if cfg["bc_sort"]:
                    qcmd_now=qcmd_base+[sampledir+"/sh/export_bc_sort.sh",outname_now,dval,dqual,sseq,squal,sizedict]
                else:
                    qcmd_now=qcmd_base+[sampledir+"/sh/export.sh",outname_now,dval,dqual,sseq,squal,sizedict]
                qcmd_now=" ".join(qcmd_now)
                s=subprocess.run(qcmd_now,shell=True)
                if s.returncode != 0:
                    print("qsub failed: Export', file=sys.stderr")
                    sys.exit(1)
            else:
                if cfg["bc_sort"]:
                    cmd_now=[sampledir+"/sh/export_bc_sort.sh",outname_now,dval,dqual,sseq,squal,sizedict]
                else:
                    cmd_now=[sampledir+"/sh/export.sh",outname_now,dval,dqual,sseq,squal,sizedict]
                cmd_now=" ".join(cmd_now)
                s=subprocess.run(cmd_now,shell=True)
                if s.returncode != 0:
                    print("Job failed: Export', file=sys.stderr")
                    sys.exit(1)
    if is_qsub:
        for sampledir in sampledir_list:
            jid_now=cmd+param_dict[os.path.basename(sampledir)]["today_now"]
            interstellar_setup.job_wait("Value optimization",jid_now,sampledir+"/qlog",njobdict[sampledir])

    for sampledir in sampledir_list:
        out_files=glob.glob(sampledir+"/value_translation/_work/export/*")
        key_list=["_R1.fastq.gz","_R2.fastq.gz","_I1.fastq.gz","_I2.fastq.gz"]
        
        for key in key_list:
            target_files=[t for t in out_files if re.search(key+r"$",os.path.basename(t))]
            if len(target_files)>0:
                cmd=["cat"]+target_files+[">",sampledir+"/value_translation/out/translated"+key]
                cmd=" ".join(cmd)
                s=subprocess.run(cmd,shell=True)
                if s.returncode != 0:
                    print("Job failed: Tagged file merge', file=sys.stderr")
                    sys.exit(1)
