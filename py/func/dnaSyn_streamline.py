from . import settingImporter
from . import settingRequirementCheck
import subprocess
import glob
import os
import re
import sys
import datetime
import time
import pandas as pd

def setUpWorkDir(workdir,cmds):
    cmds=cmds.split(",")
    cmds=["qlog","filesplit","config","sh"]+cmds
    for i in cmds:
        shell_cmd="mkdir -p "+workdir+"/"+i
        subprocess.run(shell_cmd,shell=True)


def checkRequiredFile(key,flist):
    for f in flist:
        if re.search(key,f):
            return True
    return False
    # print(targetdir,"doesn't have required files for",cmd)
    # raise ArguementError()


def pipelineCheck(cmds,cfg,workdir):
    cmds=cmds.split(",")
    cmds_done=[]
    for cmd in cmds:
        if cmd in ["import","qc"]:
            if not "import" in cfg:
                raise KeyError("'import' section is required in the config file")
        elif cmd in ["correct","mk_sval"]:
            if not "correct" in cfg:
                raise KeyError("'correct' section is required in the config file")
        elif cmd in ["to_bt"]:
            if not "to_bt" in cfg:
                raise KeyError("'to_bt' section is required in the config file")
        elif cmd in ["buildTree","mergeTree","convert","bc_sort"]:
            if not "convert" in cfg:
                raise KeyError("'convert' section is required in the config file")
        elif cmd in ["export"]:
            if not "export" in cfg:
                raise KeyError("'export' section is required in the config file")
        elif cmd in ["tag"]:
            if not "tag" in cfg:
                raise KeyError("'tag' section is required in the config file")
        elif cmd in ["demultiplex"]:
            if not "demultiplex" in cfg:
                raise KeyError("'demultiplex' section is required in the config file")
        
        cmds_done.append(cmd)

        if cmd=="qc":
            if not "import" in cmds_done:
                flist=glob.glob(workdir+"/import/*")
                flg=checkRequiredFile("srcSeq.tsv.gz",flist)
                if not flg:
                    print(flist,workdir+"/import doesn't have required files for",cmd)
                    raise ArguementError()

        if cmd=="to_bt":
            if not "import" in cmds_done and not "qc" in cmds_done:
                flist=glob.glob(workdir+"/import/*")
                flg1=checkRequiredFile("srcSeq.tsv.gz",flist)
                flist=glob.glob(workdir+"/qc/*")
                flg2=checkRequiredFile("srcSeq.QC.tsv.gz",flist)
                if not flg1 and not flg2:
                    print(flist,workdir+"/import or qc doesn't have required files for",cmd)
                    raise ArguementError()

        if cmd=="correct":
            if not "import" in cmds_done and not "qc" in cmds_done:
                flist=glob.glob(workdir+"/import/*")
                flg1=checkRequiredFile("srcCount.pkl.gz",flist)
                flist=glob.glob(workdir+"/qc/*")
                flg2=checkRequiredFile("srcCount.QC.pkl.gz",flist)
                if not flg1 and not flg2:
                    print(flist,workdir+"/import or qc doesn't have required files for",cmd)
                    raise ArguementError()
        
        if cmd=="mk_sval":
            if not "import" in cmds_done and not "qc" in cmds_done:
                flist=glob.glob(workdir+"/import/*")
                flg1=checkRequiredFile("srcSeq.tsv.gz",flist)
                flist=glob.glob(workdir+"/qc/*")
                flg2=checkRequiredFile("srcSeq.QC.tsv.gz",flist)
                if not flg1 and not flg2:
                    print(flist,workdir+"/import or qc doesn't have required files for",cmd)
                    raise ArguementError()
            if not "correct" in cmds_done:
                flist=glob.glob(workdir+"/correct/*")
                flg=checkRequiredFile("srcCorrect.pkl.gz",flist)
                if not flg:
                    print(flist,workdir+"/correct doesn't have required files for",cmd)
                    raise ArguementError()

        if cmd=="buildTree":
            if not "mk_sval" in cmds_done:
                flist=glob.glob(workdir+"/mk_sval/*")
                flg=checkRequiredFile("srcValue.tsv.gz",flist)
                if not flg:
                    print(flist,workdir+"/mk_sval doesn't have required files for",cmd)
                    raise ArguementError()
        
        if cmd=="mergeTree":
            if not "buildTree" in cmds_done:
                flist=glob.glob(workdir+"/buildTree/*")
                flg=checkRequiredFile("Tree.pkl.gz",flist)
                if not flg:
                    print(flist,workdir+"/buildTree doesn't have required files for",cmd)
                    raise ArguementError()
        
        if cmd=="convert":
            if not "mk_sval" in cmds_done:
                flist=glob.glob(workdir+"/mk_sval/*")
                flg=checkRequiredFile("srcValue.tsv.gz",flist)
                if not flg:
                    print(flist,workdir+"/mk_sval doesn't have required files for",cmd)
                    raise ArguementError()
            if not "mergeTree" in cmds_done:
                flist=glob.glob(workdir+"/mergeTree/*")
                flg=checkRequiredFile("mergeTree.pkl.gz",flist)
                if not flg:
                    print(flist,workdir+"/mergeTree doesn't have required files for",cmd)
                    raise ArguementError()
        
        if cmd=="bc_sort":
            if not "mk_sval" in cmds_done:
                flist=glob.glob(workdir+"/mk_sval/*")
                flg=checkRequiredFile("sseq_to_svalue.pkl.gz",flist)
                if not flg:
                    print(flist,workdir+"/mk_sval doesn't have required files for",cmd)
                    raise ArguementError()
        
        if cmd=="export":
            if not "import" in cmds_done and not "qc" in cmds_done:
                flist=glob.glob(workdir+"/import/*")
                flg1=checkRequiredFile("srcSeq.tsv.gz",flist)
                flist=glob.glob(workdir+"/qc/*")
                flg2=checkRequiredFile("srcSeq.QC.tsv.gz",flist)
                if not flg1 and not flg2:
                    print(flist,workdir+"/import or qc doesn't have required files for",cmd)
                    raise ArguementError()
            if not "convert" in cmds_done:
                flist=glob.glob(workdir+"/convert/*")
                flg=checkRequiredFile("converted_value.tsv.gz",flist)
                if not flg:
                    print(flist,workdir+"/convert doesn't have required files for",cmd)
                    raise ArguementError()
            if not "mergeTree" in cmds_done:
                flist=glob.glob(workdir+"/mergeTree/*")
                flg=checkRequiredFile("size_info.pkl.gz",flist)
                if not flg:
                    print(flist,workdir+"/mergeTree doesn't have required files for",cmd)
                    raise ArguementError()

        if cmd=="tag" or cmd=="demultiplex":
            if not "import" in cmds_done and not "qc" in cmds_done:
                flist=glob.glob(workdir+"/import/*")
                flg1=checkRequiredFile("srcSeq.tsv.gz",flist)
                flist=glob.glob(workdir+"/qc/*")
                flg2=checkRequiredFile("srcSeq.QC.tsv.gz",flist)
                if not flg1 and not flg2:
                    print(flist,workdir+"/import or qc doesn't have required files for",cmd)
                    raise ArguementError()
            if not "mk_sval" in cmds_done:
                flist=glob.glob(workdir+"/mk_sval/*")
                flg=checkRequiredFile("correct_result.tsv.gz",flist)
                if not flg:
                    print(flist,workdir+"/mk_sval doesn't have required files for",cmd)
                    raise ArguementError()
    return cmds


def generateShellTemplate(template,cmdline,shelloutdir,shelloutname,):
    with open(template,mode="rt") as r, open(shelloutdir+"/"+shelloutname+".sh",mode="wt") as w:
        for line in r:
            w.write(line)
        w.write("\n")
        w.write(cmdline)
        w.write("\n")


def judgeEndFix(input_files):
    for n,f in enumerate(input_files):
        if n==0:
            if re.search(r"fastq$",f):
                flg_file=1
                endfix="fastq"
            elif re.search(r"fastq.gz$",f):
                flg_file=1
                endfix="fastq.gz"
            elif re.search(r"fq$",f):
                flg_file=1
                endfix="fq"
            elif re.search(r"fq.gz$",f):
                flg_file=1
                endfix="fq.gz"
            else:
                flg_file=0
                endfix="dir"
        else:
            if re.search(r"fastq$",f):
                flg_file_now=1
                endfix_now="fastq"
            elif re.search(r"fastq.gz$",f):
                flg_file=1
                endfix_now="fastq.gz"
            elif re.search(r"fq$",f):
                flg_file=1
                endfix_now="fq"
            elif re.search(r"fq.gz$",f):
                flg_file=1
                endfix_now="fq.gz"
            else:
                flg_file_now=0
            
            #All inputs shou;d have the same format
            if flg_file==flg_file_now:
                if flg_file==1 and not endfix_now==endfix:
                    raise ArguementError("Inconsistent input files.")
            else:
                raise ArguementError("Inconsistent input files.")
    return endfix,flg_file


def judgeNumFiles(pat,targetdir,expected_num):
    path_pat=targetdir+"/*"+pat+"*"
    if len(glob.glob(path_pat))==expected_num:
        return
    else:
        raise UnknownError("Previous job seems failed. qsub error?")


def extractOutnameFromDirectory(targetdir,endfix):
    L=glob.glob(targetdir+"/*")
    outnamelist=[re.sub(r"_[^_]+\."+endfix,"",os.path.basename(i)) for i in L]
    return outnamelist


def jobCheck(jid,outdir,n_jobs):
    s=subprocess.run("qacct -j "+jid+" | grep -E 'failed|exit_status' > "+outdir+"/qlog.tmp",shell=True)
    stat_table = pd.read_csv(outdir+"/qlog.tmp",delim_whitespace=True,header=None)
    if s.returncode != 0:
        print("streamline: job check failed.', file=sys.stderr")
        sys.exit(1)
    # os.remove(outdir+"/qlog.tmp")
    if sum(stat_table[1])>0:
        raise UnknownError("qsub seems failed.")
    elif stat_table.shape[0]==2*n_jobs:
        return True
    else:
        return False


def job_wait(cmd,jid_now,outdir,njobs):
    status_ok=False
    print("Waiting for the jobs...:",cmd,flush=True)
    while not status_ok:
        time.sleep(10)
        status_ok=jobCheck(jid_now,outdir+"/qlog",njobs)
    print("Jobs completed:",cmd,flush=True)
    return


def configRewrite(cfgpath,outdir):
    with open(cfgpath,mode="rt") as r:
        cfg_orig=[i.replace("\n","") for i in r]
    for wl in glob.glob(outdir+"/*_sorted_whitelist.tsv"):
        target_segment_name=os.path.basename(wl)
        target_segment_name=re.sub(r"_sorted_whitelist\.tsv$","",wl)
        target_segment_name=re.sub(r"^streamline_","",wl)
        flg=0
        flg_target=0
        for n,l in enumerate(cfg_orig):
            if re.search(r"\[export\]",l):
                flg=1
            elif re.search(r"^\[",l):
                flg=0

            if flg==1:
                if "d_val:"+target_segment_name+",":
                    flg_target=1
                if flg_target==1:
                    if "whitelist:" in l:
                        cfg_orig[n]=re.sub(r"whitelist:.+[^,]","whiteslist:"+wl,l)
                        flg_target=0
    with open(outdir+"/sorted.conf",mode="wt") as w:
        cfg_orig="\n".join(cfg_orig)+"\n"
        w.write(cfg_orig)


class UnknownError(Exception):
    pass


class ArguementError(Exception):
    pass


class STREAMLINE_SETTINGS(object):
    def __init__(self,opt):
        self.opt=opt
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_all=settingImporter.readconfig(cfgPath)
        cfg_stream=cfg_all["streamline"]
        cfg_stream=settingRequirementCheck.setDefaultValueInConfig("streamline",cfg_stream,self.opt.distribute)

        #setup workdir
        setUpWorkDir(self.opt.outdir,cfg_stream["pipeline"])
        
        #pipeline check
        self.pipeline=pipelineCheck(cfg_stream["pipeline"],cfg_all,self.opt.outdir)
        
        #streamline option
        self.orig_cfgpath=cfgPath
        self.config=cfg_stream
        self.distribute=self.opt.distribute
        self.split=self.opt.split
        self.testrun=self.opt.testrun
        self.outname=self.opt.outname
        self.outdir=self.opt.outdir

        #import option
        self.import_opt=dict()
        self.import_opt["read1"]=self.opt.read1
        self.import_opt["read2"]=self.opt.read2
        self.import_opt["index1"]=self.opt.index1
        self.import_opt["index2"]=self.opt.index2
        # self.import_opt["is_gz"]=self.opt.input_fastq_gzipped
        self.import_opt["is_simple"]=self.opt.simple
        self.import_opt["flash"]=self.opt.flash

        #correct option
        self.correct_opt=dict()
        self.correct_opt["no_show_summary"]=self.opt.no_show_summary
        self.correct_opt["yaxis_scale"]=self.opt.yaxis_scale

        #mk_sval option
        self.mk_sval_opt=dict()
        self.mk_sval_opt["resultonly"]=self.opt.resultonly

        #buildTree, mergeTree and convert option
        self.convert_opt=dict()
        self.convert_opt["samplemerge"]=self.opt.samplemerge
        self.convert_opt["samplesheet"]=self.opt.samplesheet

        #bc_sort option
        self.bc_sort_opt=dict()
        self.bc_sort_opt["table"]=self.opt.table

        #export option
        self.export_opt=dict()
        self.export_opt["export_bclist"]=self.opt.export_bclist

        #tag and demulti option
        self.demulti_opt=dict()
        self.demulti_opt["export_tsv"]=self.opt.export_tsv


class STREAMLINE_EXE(object):
    def __init__(self,settings):
        self.settings=settings

    def generateShellScripts(self):
        if self.settings.distribute and self.settings.split:
            #Generate shell scripts for file splitting by seqkit
            input_read_files=[i for i in [self.settings.import_opt['read1'],self.settings.import_opt['read2'],self.settings.import_opt['index1'],self.settings.import_opt['index2']] if not i==""]
            # if not input_read_files:
            #     raise ArguementError("Sequencing read files are not found")
            
            sh_cmd_list_template=["seqkit","split2","-s",str(self.settings.config["chunksize"]),"-O",self.settings.outdir+"/filesplit"]
            nfile_divmod=divmod(len(input_read_files),2)
            fileindex=0
            for n in range(nfile_divmod[0]):
                files_now=["-1",input_read_files[2*n],"-2",input_read_files[2*n+1]]
                sh_cmd_list=sh_cmd_list_template+files_now
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","seqkit_split"+str(fileindex))
                fileindex+=1
            if len(input_read_files)<4:
                for n in range(2*nfile_divmod[0],2*nfile_divmod[0]+nfile_divmod[1]):
                    sh_cmd_list=sh_cmd_list_template+[input_read_files[n]]
                    sh_cmd_line=" ".join(sh_cmd_list)
                    generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","seqkit_split"+str(fileindex))

        for cmd in self.settings.pipeline:
            #Generate shell scripts for specified commands
            if cmd=="import":
                sh_cmd_list=["dnaSynergizer","import","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/import","-o","$1"]
                # if self.settings.import_opt["is_gz"]:
                #     sh_cmd_list.append("-fq_gz")
                if self.settings.import_opt["is_simple"]:
                    sh_cmd_list.append("-simple")
                if not self.settings.import_opt["flash"]=="":
                    sh_cmd_list+=["-flash",self.settings.import_opt["flash"]]
                
                #input reads
                input_read_files={i:self.settings.import_opt[i] for i in ['read1','read2','index1','index2'] if not self.settings.import_opt[i]==""}
                self.input_read_files=input_read_files
                
                for n,i in enumerate(input_read_files):
                    if i=="read1":
                        sh_cmd_list+=["-R1","$"+str(n+2)]
                    elif i=="read2":
                        sh_cmd_list+=["-R2","$"+str(n+2)]
                    elif i=="index1":
                        sh_cmd_list+=["-I1","$"+str(n+2)]
                    elif i=="index2":
                        sh_cmd_list+=["-I2","$"+str(n+2)]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","import")
            
            elif cmd=="qc":
                sh_cmd_list=["dnaSynergizer","qc","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/qc","-o","$1","-rs","$2","-rq","$3"]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","qc")

            elif cmd=="to_bt":
                sh_cmd_list=["dnaSynergizer","to_bt","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/to_bt","-o","$1","-rs","$2"]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","to_bt")
            
            elif cmd=="correct":
                sh_cmd_list=["dnaSynergizer","correct","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/correct","-o","$1","-ip","$2","-yscale",self.settings.correct_opt["yaxis_scale"]]
                if self.settings.correct_opt["no_show_summary"]:
                    sh_cmd_list.append("-no_show_summary")
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","correct")
            
            elif cmd=="mk_sval":
                sh_cmd_list=["dnaSynergizer","mk_sval","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/mk_sval","-o","$1","-rs","$2","-rq","$3","-crp","$4"]
                if self.settings.mk_sval_opt["resultonly"]:
                    sh_cmd_list.append("-resultonly")
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","mk_sval")
            
            elif cmd=="buildTree":
                sh_cmd_list=["dnaSynergizer","buildTree","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/buildTree","-o","$1","-sv","$2"]
                if self.settings.convert_opt["samplemerge"]:
                    if self.settings.convert_opt["samplesheet"]=="":
                        raise KeyError("Samplesheet is required for samplemerge")
                    
                    sh_cmd_list.append("-samplemerge")
                    sh_cmd_list+=["-samplesheet",self.settings.convert_opt["samplesheet"]]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","buildTree")

            elif cmd=="mergeTree":
                sh_cmd_list=["dnaSynergizer","mergeTree","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/mergeTree","-o","$1","-lp","$2"]
                if self.settings.convert_opt["samplemerge"]:
                    if self.settings.convert_opt["samplesheet"]=="":
                        raise KeyError("Samplesheet is required for samplemerge")
                    
                    sh_cmd_list.append("-samplemerge")
                    sh_cmd_list+=["-samplesheet",self.settings.convert_opt["samplesheet"]]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","mergeTree")
            
            elif cmd=="convert":
                sh_cmd_list=["dnaSynergizer","convert","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/convert","-o","$1","-tree","$2","-sv","$3","-sq","$4"]
                if self.settings.convert_opt["samplemerge"]:
                    if self.settings.convert_opt["samplesheet"]=="":
                        raise KeyError("Samplesheet is required for samplemerge")
                    
                    sh_cmd_list.append("-samplemerge")
                    sh_cmd_list+=["-samplesheet",self.settings.convert_opt["samplesheet"]]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","convert")

            elif cmd=="bc_sort":
                sh_cmd_list=["dnaSynergizer","bc_sort","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/bc_sort","-o","$1","-tree","$2","-sseq_to_svalue","$3","-tbl",self.settings.bc_sort_opt["table"]]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","bc_sort")

            elif cmd=="export":
                if not "bc_sort" in self.settings.pipeline:
                    sh_cmd_list=["dnaSynergizer","export","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/export","-o","$1","-dv","$2","-dq","$3","-rs","$4","-rq","$5","-size","$6"]
                else:
                    sh_cmd_list=["dnaSynergizer","export","-conf",self.settings.outdir+"/bc_sort/sorted.conf","-d",self.settings.outdir+"/export","-o","$1","-dv","$2","-dq","$3","-rs","$4","-rq","$5","-size","$6"]
                if self.settings.export_opt["export_bclist"]:
                    sh_cmd_list.append("-export_bclist")
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","export")

            elif cmd=="demultiplex":
                sh_cmd_list=["dnaSynergizer","demultiplex","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/demultiplex","-o","$1","-cs","$2","-cq","$3","-rq","$4"]
                if self.settings.demultiplex_opt["export_tsv"]:
                    sh_cmd_list.append("-export_tsv")
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","demultiplex")

            elif cmd=="tag":
                sh_cmd_list=["dnaSynergizer","tag","-conf",self.settings.orig_cfgpath,"-d",self.settings.outdir+"/tag","-o","$1","-cs","$2","-cq","$3","-rq","$4"]
                sh_cmd_line=" ".join(sh_cmd_list)
                generateShellTemplate(self.settings.config["template_shellscript"],sh_cmd_line,self.settings.outdir+"/sh","tag")
        
        subprocess.run("chmod +x "+self.settings.outdir+"/sh/*",shell=True)


    def run(self):
        used_commands=[]
        input_file_list=[self.input_read_files[i] for i in ['read1','read2','index1','index2'] if i in self.input_read_files]
        basename_list=[os.path.basename(i) for i in input_file_list]

        #single CPU mode
        if not self.settings.distribute:
            for cmd in self.settings.pipeline:
                if cmd=="import":
                    cmdlist=[self.settings.outdir+"/sh/import.sh",self.settings.outname]
                    # input_files=[self.input_read_files[i] for i in ['read1','read2','index1','index2'] if i in self.input_read_files]
                    cmdlist+=input_file_list
                    
                elif cmd=="qc":
                    seqfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcSeq.tsv.gz"
                    qualfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcQual.tsv.gz"
                    cmdlist=[self.settings.outdir+"/sh/qc.sh",self.settings.outname,seqfile,qualfile]
                    
                elif cmd=="to_bt":
                    flg=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc"))
                    if flg:
                        seqfile=self.settings.outdir+"/qc/"+self.settings.outname+"_srcSeq.QC.tsv.gz"
                    else:
                        seqfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcSeq.tsv.gz"
                    cmdlist=[self.settings.outdir+"/sh/to_bt.sh",self.settings.outname,seqfile]
                    
                elif cmd=="correct":
                    flg=checkRequiredFile("_srcCount.QC.pkl.gz",glob.glob(self.settings.outdir+"/qc"))
                    if flg:
                        cntfiles='"'+self.settings.outdir+"/qc/"+self.settings.outname+"*_srcCount.QC.pkl.gz"+'"'
                    else:
                        cntfiles='"'+self.settings.outdir+"/import/"+self.settings.outname+"*_srcCount.pkl.gz"+'"'
                    cmdlist=[self.settings.outdir+"/sh/correct.sh",self.settings.outname,cntfiles]
                    
                elif cmd=="mk_sval":
                    flg=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc"))
                    if flg:
                        seqfile=self.settings.outdir+"/qc/"+self.settings.outname+"_srcSeq.QC.tsv.gz"
                        qualfile=self.settings.outdir+"/qc/"+self.settings.outname+"_srcQual.QC.tsv.gz"
                    else:
                        seqfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcSeq.tsv.gz"
                        qualfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcQual.tsv.gz"
                    correctionfile=self.settings.outdir+"/correct/"+self.settings.outname+"_srcCorrect.pkl.gz"
                    cmdlist=[self.settings.outdir+"/sh/mk_sval.sh",self.settings.outname,seqfile,qualfile,correctionfile]
                    
                elif cmd=="buildTree":
                    valuefile=self.settings.outdir+"/mk_sval/"+self.settings.outname+"_correct_srcValue.tsv.gz"
                    cmdlist=[self.settings.outdir+"/sh/buildTree.sh",self.settings.outname,valuefile]
                    
                elif cmd=="mergeTree":
                    trees='"'+self.settings.outdir+"/buildTree/"+self.settings.outname+"*_Tree.pkl.gz"+'"'
                    cmdlist=[self.settings.outdir+"/sh/mergeTree.sh",self.settings.outname,trees]
                    
                elif cmd=="convert":
                    tree=self.settings.outdir+"/mergeTree/"+self.settings.outname+"_mergeTree.pkl.gz"
                    valuefile=self.settings.outdir+"/mk_sval/"+self.settings.outname+"_correct_srcValue.tsv.gz"
                    qualfile =self.settings.outdir+"/mk_sval/"+self.settings.outname+"_correct_srcQual.tsv.gz"
                    cmdlist=[self.settings.outdir+"/sh/convert.sh",self.settings.outname,tree,valuefile,qualfile]

                elif cmd=="bc_sort":
                    tree=self.settings.outdir+"/mergeTree/"+self.settings.outname+"_mergeTree.pkl.gz"
                    s2v=self.settings.outdir+"/mk_sval/"+self.settings.outname+"_sseq_to_svalue.pkl.gz"
                    cmdlist=[self.settings.outdir+"/sh/bc_sort.sh",self.settings.outname,tree,s2v]

                elif cmd=="export":
                    flg=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc"))
                    if flg:
                        seqfile=self.settings.outdir+"/qc/"+self.settings.outname+"_srcSeq.QC.tsv.gz"
                        qualfile=self.settings.outdir+"/qc/"+self.settings.outname+"_srcQual.QC.tsv.gz"
                    else:
                        seqfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcSeq.tsv.gz"
                        qualfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcQual.tsv.gz"
                    dval =self.settings.outdir+"/convert/"+self.settings.outname+"_converted_value.tsv.gz"
                    dqual=self.settings.outdir+"/convert/"+self.settings.outname+"_converted_qual.tsv.gz"
                    sizefile=self.settings.outdir+"/mergeTree/"+self.settings.outname+"_size_info.pkl.gz"
                    cmdlist=[self.settings.outdir+"/sh/export.sh",self.settings.outname,dval,dqual,seqfile,qualfile,sizefile]

                elif cmd=="demultiplex":
                    flg=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc"))
                    if flg:
                        qualfile=self.settings.outdir+"/qc/"+self.settings.outname+"_srcQual.QC.tsv.gz"
                    else:
                        qualfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcQual.tsv.gz"
                    resfile=self.settings.outdir+"/mk_sval/"+self.settings.outname+"_correct_result.tsv.gz"
                    cqual=self.settings.outdir+"/mk_sval/"+self.settings.outname+"_correct_srcQual.tsv.gz"
                    cmdlist=[self.settings.outdir+"/sh/demultiplex.sh",self.settings.outname,resfile,cqual,qualfile]
                    
                elif cmd=="tag":
                    flg=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc"))
                    if flg:
                        qualfile=self.settings.outdir+"/qc/"+self.settings.outname+"_srcQual.QC.tsv.gz"
                    else:
                        qualfile=self.settings.outdir+"/import/"+self.settings.outname+"_srcQual.tsv.gz"
                    resfile=self.settings.outdir+"/mk_sval/"+self.settings.outname+"_correct_result.tsv.gz"
                    cqual=self.settings.outdir+"/mk_sval/"+self.settings.outname+"_correct_srcQual.tsv.gz"
                    cmdlist=[self.settings.outdir+"/sh/demultiplex.sh",self.settings.outname,resfile,cqual,qualfile]

                cmdline=" ".join(cmdlist)
                s=subprocess.run(cmdline,shell=True)
                used_commands.append(cmdline)
                if s.returncode != 0:
                    print("streamline:",cmd,'failed.', file=sys.stderr)
                    sys.exit(1)

        #qsub mode
        else:
            today_now=str(datetime.datetime.today()).replace(" ","")
            today_now=today_now.replace(":","")
            today_now=today_now.replace(".","")
            today_now=today_now.replace("-","")
            
            if input_file_list:
                #Judge input format: file or directory / endfix determination
                endfix,flg_file=judgeEndFix(input_file_list)

                #Prepare input files
                if flg_file==1: #split or not
                    if self.settings.split:
                        #fastq split
                        print("Runnning qsub jobs...: file split",flush=True)
                        qoption=self.settings.config["qoption"]
                        qoption=qoption.replace("<mem>",self.settings.config["mem_max"])
                        qcmd_base=["qsub",qoption,"-e",self.settings.outdir+"/qlog","-o",self.settings.outdir+"/qlog","-cwd"]
                        for i in glob.glob(self.settings.outdir+"/sh/seqkit*"):
                            qcmd_now=qcmd_base+["-N","seqkitsplit"+today_now,i]
                            qcmd_now=" ".join(qcmd_now)
                            print(qcmd_now)
                            s=subprocess.run(qcmd_now,shell=True)
                            used_commands.append(qcmd_now)
                            if s.returncode != 0:
                                print("streamline: qsub for seqkit split failed.', file=sys.stderr")
                                sys.exit(1)
                        time.sleep(30)
                        status_ok=False
                        print("Waiting for the jobs...: file split",flush=True)
                        while not status_ok:
                            time.sleep(10)
                            status_ok=jobCheck("seqkitsplit"+today_now,self.settings.outdir+"/qlog",len(glob.glob(self.settings.outdir+"/sh/seqkit*")))
                        print("Jobs completed: file split")
                        # filename_to_search=basename_list[0].replace("."+endfix,"")
                        # judgeNumFiles(filename_to_search,self.settings.outdir+"/filesplit",len(glob.glob(self.settings.outdir+"/filesplit"))/len(input_files))
                        # outname_list=[os.path.basename(i).replace("."+endfix) for i in glob.glob(self.settings.outdir+"/filesplit/"+filename_to_search)]
                    # else:
                        #no split
                        # outname_list=[basename_list[0]]
                
                else:
                    #Inputs are directories containing already-split files
                    for f in glob.glob(input_file_list[0]):
                        if re.search(r"fastq$",f):
                            endfix="fastq"
                        elif re.search(r"fastq.gz$",f):
                            endfix="fastq.gz"
                        elif re.search(r"fq$",f):
                            endfix="fq"
                        elif re.search(r"fq.gz$",f):
                            endfix="fq.gz"
                        else:
                            raise ArguementError("Input files should be fastq or fastq.gz")
                        break

                    # outname_list=extractOutnameFromDirectory(input_file_list[0],endfix)
            
            jid_prev=""
            for n_cmd,cmd in enumerate(self.settings.pipeline):
                # do_qhold=True if n_cmd>0 else False
                qoption=self.settings.config["qoption"]
                jid_now=cmd+today_now
                qcmd_base=["qsub","-e",self.settings.outdir+"/qlog","-o",self.settings.outdir+"/qlog","-cwd","-N",jid_now]

                if cmd=="import":
                    print("Runnning qsub jobs...: import",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_import"])
                    qcmd_now=qcmd_base+[qoption]
                    
                    if self.settings.split:
                        file_pool=[]
                        chunk_identifier_pool=list(set(re.sub(r"^.+\.","",i.replace("."+endfix,"")) for i in glob.glob(self.settings.outdir+"/filesplit/*")))
                        for f_name in basename_list:
                            file_pool.append([self.settings.outdir+"/filesplit/"+f_name.replace(endfix,"")+i+"."+endfix for i in chunk_identifier_pool])
                    else:
                        file_pool=[]
                        if flg_file==1:
                            for f_name in input_file_list:
                                file_pool.append([f_name])
                        else:
                            prefix_pool=[re.sub(r"_[^_]+\."+endfix,"",os.path.basename(i)) for i in glob.glob(input_file_list[0]+"/*")]
                            for f_name in input_file_list:
                                file_pool.append([glob.glob(f_name+"/"+p+"*")[0] for p in prefix_pool])

                    for infile in zip(*file_pool):
                        infile_now=[infile[i] for i in range(len(input_file_list))]
                        qcmd_now+=[self.settings.outdir+"/sh/import.sh",infile_now[0].replace("."+endfix,"")]+infile_now
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: import', file=sys.stderr")
                            sys.exit(1)

                    njobs=len(file_pool[0])
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now 
                    
                elif cmd=="qc":
                    print("Runnning qsub jobs...: qc",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_qc"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/qc.sh"]
                    file_pool=[i for i in glob.glob(self.settings.outdir+"/import/*") if re.search("_srcSeq.tsv.gz",i)]
                    for f in file_pool:
                        outname_now=re.sub("_srcSeq.tsv.gz","",os.path.basename(f))
                        qcmd_now=qcmd_base+[outname_now,f,re.sub(r"_srcSeq\.","_srcQual.",f)]
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: qc', file=sys.stderr")
                            sys.exit(1)
                    
                    njobs=len(file_pool)
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now

                elif cmd=="to_bt":
                    print("Runnning qsub jobs...: to_bt",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_to_bt"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/to_bt.sh"]
                    is_qc=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc/*"))
                    if is_qc:
                        outname_now="streamline"
                        qcmd_now=qcmd_base+[outname_now,'"'+self.settings.outdir+"/qc/*_srcSeq.QC.tsv.gz"+'"']                        
                    else:
                        outname_now="streamline"
                        qcmd_now=qcmd_base+[outname_now,'"'+self.settings.outdir+"/import/*_srcSeq.tsv.gz"+'"']

                    qcmd_now=" ".join(qcmd_now)
                    s=subprocess.run(qcmd_now,shell=True)
                    used_commands.append(qcmd_now)
                    if s.returncode != 0:
                        print("qsub failed: to_bt', file=sys.stderr")
                        sys.exit(1)

                    njobs=1
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now 

                elif cmd=="correct":
                    print("Runnning qsub jobs...: correct",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_correct"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/correct.sh"]
                    is_qc=checkRequiredFile("_srcCount.QC.pkl.gz",glob.glob(self.settings.outdir+"/qc/*"))
                    if is_qc:
                        outname_now="streamline"
                        qcmd_now=qcmd_base+[outname_now,'"'+self.settings.outdir+"/qc/*_srcCount.QC.pkl.gz"+'"']
                    else:
                        outname_now="streamline"
                        qcmd_now=qcmd_base+[outname_now,'"'+self.settings.outdir+"/import/*_srcCount.pkl.gz"+'"']

                    qcmd_now=" ".join(qcmd_now)
                    s=subprocess.run(qcmd_now,shell=True)
                    used_commands.append(qcmd_now)
                    if s.returncode != 0:
                        print("qsub failed: correct', file=sys.stderr")
                        sys.exit(1)

                    njobs=1
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now       
                         
                elif cmd=="mk_sval":
                    print("Runnning qsub jobs...: mk_sval",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_mk_sval"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/mk_sval.sh"]

                    is_qc=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc/*"))
                    if is_qc:
                        file_endfix="_srcSeq.QC.tsv.gz"
                        file_pool=[i for i in glob.glob(self.settings.outdir+"/qc/*") if re.search("_srcSeq.QC.tsv.gz",i)]              
                    else:
                        file_endfix="_srcSeq.tsv.gz"
                        file_pool=[i for i in glob.glob(self.settings.outdir+"/import/*") if re.search("_srcSeq.tsv.gz",i)]
                        
                    for f in file_pool:
                        outname_now=re.sub(file_endfix,"",os.path.basename(f))
                        qcmd_now=qcmd_base+[outname_now,f,re.sub(r"_srcSeq\.","_srcQual.",f),self.settings.outdir+"/correct/streamline_srcCorrect.pkl.gz"]
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: mk_sval', file=sys.stderr")
                            sys.exit(1)

                    njobs=len(file_pool)
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now   
                
                elif cmd=="buildTree":
                    print("Runnning qsub jobs...: buildTree",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_buildTree"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/buildTree.sh"]
                      
                    file_endfix="_srcValue.tsv.gz"
                    file_pool=[i for i in glob.glob(self.settings.outdir+"/mk_sval/*") if re.search(file_endfix,i)]    
                        
                    for f in file_pool:
                        outname_now=re.sub(file_endfix,"",os.path.basename(f))
                        qcmd_now=qcmd_base+[outname_now,f]
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: buildTree', file=sys.stderr")
                            sys.exit(1)

                    njobs=len(file_pool)
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now   

                elif cmd=="mergeTree":
                    print("Runnning qsub jobs...: mergeTree",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_mergeTree"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/mergeTree.sh"]
                    outname_now="streamline"
                    qcmd_now=qcmd_base+[outname_now,'"'+self.settings.outdir+"/buildTree/*_Tree.pkl.gz"+'"']

                    qcmd_now=" ".join(qcmd_now)
                    s=subprocess.run(qcmd_now,shell=True)
                    used_commands.append(qcmd_now)
                    if s.returncode != 0:
                        print("qsub failed: mergeTree', file=sys.stderr")
                        sys.exit(1)

                    njobs=1
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now     
           
                elif cmd=="convert":
                    print("Runnning qsub jobs...: convert",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_convert"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/convert.sh"]
                      
                    file_endfix="_srcValue.tsv.gz"
                    file_pool=[i for i in glob.glob(self.settings.outdir+"/mk_sval/*") if re.search(file_endfix,i)]
                    mergetree=glob.glob(self.settings.outdir+"/mergeTree/*_mergeTree.pkl.gz")[0]
                        
                    for f in file_pool:
                        outname_now=re.sub(file_endfix,"",os.path.basename(f))
                        qcmd_now=qcmd_base+[outname_now,mergetree,f,re.sub(file_endfix,"_srcQual.tsv.gz",f)]
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: convert', file=sys.stderr")
                            sys.exit(1)

                    njobs=len(file_pool)
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now   

                elif cmd=="bc_sort":
                    print("Runnning qsub jobs...: bc_sort",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_bc_sort"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/bc_sort.sh"]
                    outname_now="streamline"
                    mergetree=glob.glob(self.settings.outdir+"/mergeTree/*_mergeTree.pkl.gz")[0]
                    s2v=glob.glob(self.settings.outdir+"/mk_sval/*_sseq_to_svalue.pkl.gz")[0]
                    qcmd_now=qcmd_base+[outname_now,mergetree,s2v]

                    qcmd_now=" ".join(qcmd_now)
                    s=subprocess.run(qcmd_now,shell=True)
                    used_commands.append(qcmd_now)
                    if s.returncode != 0:
                        print("qsub failed: bc_sort', file=sys.stderr")
                        sys.exit(1)

                    njobs=1
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now
                    configRewrite(self.settings.orig_cfgpath,self.settings.outdir+"/bc_sort")

                elif cmd=="export":
                    print("Runnning qsub jobs...: export",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_export"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/export.sh"]

                    is_qc=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc/*"))
                    if is_qc:
                        file_endfix="_srcSeq.QC.tsv.gz"
                        file_prefix=[os.path.basename(i).replace(endfix,"") for i in glob.glob(self.settings.outdir+"/qc/*") if re.search(endfix,i)]
                    else:
                        file_endfix="_srcSeq.tsv.gz"
                        file_prefix=[os.path.basename(i.replace(endfix,"")) for i in glob.glob(self.settings.outdir+"/import/*") if re.search(endfix,i)]
                    
                    for outname_now in file_prefix:
                        dval=self.settings.outdir+"/convert/"+outname_now+"_converted_value.tsv.gz"
                        dqual=self.settings.outdir+"/convert/"+outname_now+"_converted_qual.tsv.gz"
                        if is_qc:
                            sseq =self.settings.outdir+"/qc/"+outname_now+file_endfix
                            squal=self.settings.outdir+"/qc/"+outname_now+"_srcQual.QC.tsv.gz"
                        else:
                            sseq =self.settings.outdir+"/import/"+outname_now+file_endfix
                            squal=self.settings.outdir+"/import/"+outname_now+"_srcQual.tsv.gz"
                        sizedict=self.settings.outdir+"/mergeTree/streamline_size_info.pkl.gz"
                        qcmd_now=qcmd_base+[outname_now,dval,dqual,sseq,squal,sizedict]
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: export', file=sys.stderr")
                            sys.exit(1)

                    njobs=len(file_pool)
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now   
                
                elif cmd=="demultiplex":
                    print("Runnning qsub jobs...: demultiplex",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_demultiplex"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/demultiplex.sh"]
                    is_qc=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc/*"))                    
                    file_prefix=[os.path.basename(i.replace("_correct_result.tsv.gz","")) for i in glob.glob(self.settings.outdir+"/mk_sval/*") if re.search("_correct_result.tsv.gz",i)]

                    for outname_now in file_prefix:
                        result=self.settings.outdir+"/mk_sval/"+outname_now+"_correct_result.tsv.gz"
                        cqual =self.settings.outdir+"/mk_sval/"+outname_now+"_correct_srcQual.tsv.gz"
                        if is_qc:
                            squal=self.settings.outdir+"/qc/"+outname_now+"_srcQual.QC.tsv.gz"
                        else:
                            squal=self.settings.outdir+"/import/"+outname_now+"_srcQual.tsv.gz"
                        qcmd_now=qcmd_base+[outname_now,result,cqual,squal]
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: demultiplex', file=sys.stderr")
                            sys.exit(1)

                    njobs=len(file_pool)
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now   

                elif cmd=="tag":
                    print("Runnning qsub jobs...: tag",flush=True)
                    qoption=qoption.replace("<mem>",self.settings.config["mem_tag"])
                    qcmd_base=qcmd_base+[qoption,"-hold_jid",jid_prev,self.settings.outdir+"/sh/tag.sh"]
                    is_qc=checkRequiredFile("_srcSeq.QC.tsv.gz",glob.glob(self.settings.outdir+"/qc/*"))                    
                    file_prefix=[os.path.basename(i.replace("_correct_result.tsv.gz","")) for i in glob.glob(self.settings.outdir+"/mk_sval/*") if re.search("_correct_result.tsv.gz",i)]

                    for outname_now in file_prefix:
                        result=self.settings.outdir+"/mk_sval/"+outname_now+"_correct_result.tsv.gz"
                        cqual =self.settings.outdir+"/mk_sval/"+outname_now+"_correct_srcQual.tsv.gz"
                        if is_qc:
                            squal=self.settings.outdir+"/qc/"+outname_now+"_srcQual.QC.tsv.gz"
                        else:
                            squal=self.settings.outdir+"/import/"+outname_now+"_srcQual.tsv.gz"
                        qcmd_now=qcmd_base+[outname_now,result,cqual,squal]
                        qcmd_now=" ".join(qcmd_now)
                        s=subprocess.run(qcmd_now,shell=True)
                        used_commands.append(qcmd_now)
                        if s.returncode != 0:
                            print("qsub failed: tag', file=sys.stderr")
                            sys.exit(1)

                    njobs=len(file_pool)
                    job_wait(cmd,jid_now,self.settings.outdir,njobs)
                    jid_prev=cmd+today_now   


        with open(self.settings.outdir+"/commandlog.txt",mode="wt") as w:
            for line in cmdlist:
                w.write(line+"\n")

            
            
