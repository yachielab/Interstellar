import subprocess
import re
import glob
import datetime
import sys
import time
import os
import pandas as pd

def checkRequiredFile(key,flist):
    for f in flist:
        if re.search(key,f):
            return True
    return False

    
def setUpSampleDir(sampledir,cmds_execute):
    cmds_execute=["qlog","filesplit","config","sh"]+cmds_execute
    for i in cmds_execute:
        shell_cmd="mkdir -p "+sampledir+"/"+i
        subprocess.run(shell_cmd,shell=True)

    for dir in ["value_extraction","value_translation","tag","demultiplex"]:
        if not dir in cmds_execute:
            continue

        shell_cmd="mkdir -p "+sampledir+"/"+dir+"/_work "+sampledir+"/"+dir+"/out"
        subprocess.run(shell_cmd,shell=True)
        if dir=="value_extraction":
            shell_cmd="mkdir -p "+sampledir+"/"+dir+"/_work/import"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/qc"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/mk_sval"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/correct"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/to_bt"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/bartender"
            subprocess.run(shell_cmd,shell=True)
        elif dir=="value_translation":
            shell_cmd="mkdir -p "+sampledir+"/"+dir+"/_work/buildTree"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/mergeTree"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/convert"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/bc_sort"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/export"
            subprocess.run(shell_cmd,shell=True)


def judgeEndFix(input_files):
    for n,f in enumerate(input_files):
        if n==0:
            if re.search(r"fastq$",f):
                endfix="fastq"
            elif re.search(r"fastq.gz$",f):
                endfix="fastq.gz"
            elif re.search(r"fq$",f):
                endfix="fq"
            elif re.search(r"fq.gz$",f):
                endfix="fq.gz"
            
        else:
            if re.search(r"fastq$",f):
                endfix_now="fastq"
            elif re.search(r"fastq.gz$",f):
                endfix_now="fastq.gz"
            elif re.search(r"fq$",f):
                endfix_now="fq"
            elif re.search(r"fq.gz$",f):
                endfix_now="fq.gz"
            
            #All inputs should have the same format
            if not endfix_now==endfix:
                raise ArguementError("File format should be same across all input files.")       
    return endfix


def fastqDirCheck(valid_dirs):
    if valid_dirs =="":
        raise EmptyError("At least 1 input fastq directory should be provided.")
    
    #collect all input files
    all_files=[]
    for d in valid_dirs:
        all_files+=glob.glob(d+"/*")
    
    #judge suffix
    file_suffix=judgeEndFix(all_files)
    return file_suffix


def generateShellTemplate(template,cmdline,shelloutdir,shelloutname):
    with open(os.path.expanduser(template),mode="rt") as r, open(os.path.expanduser(shelloutdir+"/"+shelloutname+".sh"),mode="wt") as w:
        for line in r:
            w.write(line)
        w.write("\n")
        w.write(cmdline)
        w.write("\n")


def jobCheck(jid,outdir,n_jobs):
    s=subprocess.run("qacct -j "+jid+" | grep -E 'failed|exit_status' > "+outdir+"/qlog.tmp",shell=True)
    if os.path.getsize(outdir+"/qlog.tmp")==0:
        return False

    stat_table = pd.read_csv(outdir+"/qlog.tmp",delim_whitespace=True,header=None)
    if s.returncode != 0:
        print("streamline: job check failed.', file=sys.stderr")
        sys.exit(1)
    # os.remove(outdir+"/qlog.tmp")
    if sum(stat_table[1])>0:
        raise UnknownError("qsub failed.")
    elif stat_table.shape[0]==2*n_jobs:
        return True
    else:
        return False


def job_wait(cmd,jid_now,outdir,njobs):
    status_ok=False
    print("Waiting for the jobs...:",cmd,flush=True)
    while not status_ok:
        time.sleep(10)
        status_ok=jobCheck(jid_now,outdir,njobs)
    print("Jobs completed:",cmd,flush=True)
    return


class EmptyError(Exception):
    pass


class UnknownError(Exception):
    pass


class ArguementError(Exception):
    pass


class SETUP_SETTINGS(object):
    def __init__(self,cfg,qcfg,cmds_execute,sampledir,samplename):
        self.cfg=cfg
        self.qcfg=qcfg
        self.cmds_execute=cmds_execute
        self.sampledir=sampledir
        self.samplename=samplename

    def settingGetter(self):
        #setup sampledir
        setUpSampleDir(self.sampledir,self.cmds_execute)

        if "value_extraction" in self.cmds_execute:

            #FASTQ file name parsing
            config_val_ext=self.cfg["value_extraction"]
            read1_dir=config_val_ext["READ1_DIR"]
            read2_dir=config_val_ext["READ2_DIR"]
            index1_dir=config_val_ext["INDEX1_DIR"]
            index2_dir=config_val_ext["INDEX2_DIR"]

            tmp_order=["read1","read2","index1","index2"]
            valid_dirs=[i for i in [read1_dir,read2_dir,index1_dir,index2_dir] if not i == ""]
            read_valid={k:"" for k in tmp_order}
            for n,i in enumerate([read1_dir,read2_dir,index1_dir,index2_dir]):
                if not i == "":
                    read_valid[tmp_order[n]]=i

            file_suffix=fastqDirCheck(valid_dirs)
            
            if not self.cfg["general"]["SAMPLESHEET"]=="":
                # multi-sample
                samplelist=[]
                with open(os.path.expanduser(self.cfg["general"]["SAMPLESHEET"]),mode="rt") as r:
                    for line in r:
                        line=line.replace("\n","")
                        samplelist.append(line.split("\t"))
                target_prefix_list=[l[0] for l in samplelist if l[1]==self.samplename]
                # if not self.qcfg == "":
                #     # qsub = True -> split files
                #     for t in target_prefix_list:
                #         shell_cmd="mkdir -p "+self.sampledir+"/filesplit/"+t
                #         subprocess.run(shell_cmd,shell=True)
            else:
                # single-sample
                all_inputs=[sorted(glob.glob(i+"/*")) for i in valid_dirs]
                target_prefix_list=[]
                for read_tup in zip(*all_inputs):
                    read_list=list(read_tup)
                    read_list=[os.path.basename(i) for i in read_list]
                    prefix = []
                    for x in zip(*read_list):
                        if len(set(x)) == 1:
                            prefix.append(x[0])
                        else:
                            break
                target_prefix_list.append("".join(prefix))
                # if not self.qcfg == "":
                #     # qsub = True -> split files
                #     shell_cmd="mkdir -p "+self.sampledir+"/filesplit/"+"".join(prefix)
                #     subprocess.run(shell_cmd,shell=True)
                        
            # Store values
            self.valid_dirs=valid_dirs
            self.read_valid=read_valid
            self.file_suffix=file_suffix
            self.target_prefix_list=target_prefix_list

class SETUP(object):
    def __init__(self,settings,is_qsub,is_multisample,cfgpath):
        self.settings=settings
        self.is_qsub=is_qsub
        self.is_multisample=is_multisample
        self.cfgpath=cfgpath

    def generateShellScripts(self):
        shelldir=self.settings.sampledir+"/sh"
        self.shelldir=shelldir
        read_iden_dict={}
        today_now=str(datetime.datetime.today()).replace(" ","")
        today_now=today_now.replace(":","")
        today_now=today_now.replace(".","")
        today_now=today_now.replace("-","")
        self.today_now=today_now

        if self.is_qsub:
            #Generate shell scripts for file splitting by seqkit
            for prefix in self.settings.target_prefix_list:
                input_read_files=[]

                ##
                print(self.settings.read_valid)
                print(self.settings.target_prefix_list)
                ##
                
                for r in ["read1","read2","index1","index2"]:
                    if r in self.settings.read_valid and not self.settings.read_valid[r]=="":
                        input_read_files.append(glob.glob(self.settings.read_valid[r]+"/"+prefix+"*")[0])
                
                sh_cmd_list_template=["seqkit","split2","-s",str(self.settings.cfg["general"]["CHUNKSIZE"])]
                nfile_divmod=divmod(len(input_read_files),2)
                fileindex=0
                for n in range(nfile_divmod[0]):
                    files_now=["-1",input_read_files[2*n],"-2",input_read_files[2*n+1]]
                    sh_cmd_list=sh_cmd_list_template+["-O",self.settings.sampledir+"/filesplit/"+prefix+"_"+str(fileindex)]+files_now
                    sh_cmd_line=" ".join(sh_cmd_list)
                    generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"seqkit_split_"+prefix+"_"+str(fileindex))
                    fileindex+=1
                if len(input_read_files)<4:
                    for n in range(2*nfile_divmod[0],2*nfile_divmod[0]+nfile_divmod[1]):
                        sh_cmd_list=sh_cmd_list_template+["-O",self.settings.sampledir+"/filesplit/"+prefix+"_"+str(fileindex)]+[input_read_files[n]]
                        sh_cmd_line=" ".join(sh_cmd_list)
                        generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"seqkit_split_"+prefix+"_"+str(fileindex))
                
                #get read identifier
                read_iden_dict[prefix]={}
                for n,r in enumerate(["read1","read2","index1","index2"]):
                    if r in self.settings.read_valid and not self.settings.read_valid[r]=="":
                        read_iden_dict[r]=input_read_files[n].replace(self.settings.file_suffix,"").replace(prefix,"")
        else:
            #Generate shell scripts for file splitting by seqkit
            for prefix in self.settings.target_prefix_list:
                input_read_files=[]
                for r in ["read1","read2","index1","index2"]:
                    if not self.settings.read_valid[r]=="":
                        input_read_files.append(glob.glob(self.settings.read_valid[r]+"/"+prefix+"*")[0])
                
                #get read identifier
                read_iden_dict[prefix]={}
                for n,r in enumerate(["read1","read2","index1","index2"]):
                    if not self.settings.read_valid[r]=="":
                        read_iden_dict[r]=input_read_files[n].replace(self.settings.file_suffix,"").replace(prefix,"")
        self.read_iden_dict=read_iden_dict

            
        
        #Generate shell scripts for specified commands
        if "value_extraction" in self.settings.cmds_execute:
            outdir=self.settings.sampledir+"/value_extraction/_work/"

            
            # import
            sh_cmd_list=["Interstellar-exec","import","-conf",self.cfgpath,"-d",outdir+"/import","-o","$1"]
            if not self.settings.cfg["value_extraction"]["FLASH"]=="":
                sh_cmd_list+=["-flash",self.settings.cfg["value_extraction"]["FLASH"]]
            for n,i in enumerate(self.settings.read_valid):
                if i=="read1" and not self.settings.read_valid[i]=="":
                    sh_cmd_list+=["-R1","$"+str(n+2)]
                elif i=="read2" and not self.settings.read_valid[i]=="":
                    sh_cmd_list+=["-R2","$"+str(n+2)]
                elif i=="index1" and not self.settings.read_valid[i]=="":
                    sh_cmd_list+=["-I1","$"+str(n+2)]
                elif i=="index2" and not self.settings.read_valid[i]=="":
                    sh_cmd_list+=["-I2","$"+str(n+2)]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"import")
            # qc
            sh_cmd_list=["Interstellar-exec","qc","-conf",self.cfgpath,"-d",outdir+"/qc","-o","$1","-rs","$2","-rq","$3"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"qc")
            # to_bt
            sh_cmd_list=["Interstellar-exec","to_bt","-conf",self.cfgpath,"-d",outdir+"/to_bt","-o","$1","-rs","$2"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"to_bt")
            # correct
            sh_cmd_list=["Interstellar-exec","correct","-conf",self.cfgpath,"-d",outdir+"/correct","-o","$1","-ip","$2"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"correct")
            # mk_sval
            sh_cmd_list=["Interstellar-exec","mk_sval","-conf",self.cfgpath,"-d",outdir+"/mk_sval","-o","$1","-rs","$2","-rq","$3","-crp","$4"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"mk_sval")
        

        if "value_translation" in self.settings.cmds_execute:
            outdir=self.settings.sampledir+"/value_translation/_work/"

            #buildTree
            sh_cmd_list=["Interstellar-exec","buildTree","-conf",self.cfgpath,"-d",outdir+"/buildTree","-o","$1","-sv","$2"]
            if not self.settings.cfg["general"]["SAMPLESHEET"]=="":
                sh_cmd_list.append("-samplemerge")
                sh_cmd_list+=["-samplesheet",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/samplesheet/samplesheet.tsv"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"buildTree")
            #mergeTree
            sh_cmd_list=["Interstellar-exec","mergeTree","-conf",self.cfgpath,"-o","$1","-lp","$2"]
            if not self.settings.cfg["general"]["SAMPLESHEET"]=="":
                sh_cmd_list.append("-samplemerge")
                sh_cmd_list+=["-samplesheet",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/samplesheet/samplesheet.tsv","-d",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/mergeTree"]
            else:
                sh_cmd_list+=["-d",outdir+"/mergeTree"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"mergeTree")
            #convert
            sh_cmd_list=["Interstellar-exec","convert","-conf",self.cfgpath,"-d",outdir+"/convert","-o","$1","-tree","$2","-sv","$3","-sq","$4"]
            if not self.settings.cfg["general"]["SAMPLESHEET"]=="":
                sh_cmd_list.append("-samplemerge")
                sh_cmd_list+=["-samplesheet",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/samplesheet/samplesheet.tsv"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"convert")
            # bc_sort (cannot executed in multi-sample mode.)
            sh_cmd_list=["Interstellar-exec","bc_sort","-conf",self.cfgpath,"-d",outdir+"/bc_sort","-o","$1","-tree","$2","-sseq_to_svalue","$3","-tbl","$4"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"bc_sort")
            # export: normal
            sh_cmd_list=["Interstellar-exec","export","-conf",self.cfgpath,"-d",outdir+"/export","-o","$1","-dv","$2","-dq","$3","-rs","$4","-rq","$5","-size","$6","-export_bclist"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"export")
            # export: bc_sort
            sh_cmd_list=["Interstellar-exec","export","-conf",outdir+"/bc_sort/sorted.conf","-d",outdir+"/export","-o","$1","-dv","$2","-dq","$3","-rs","$4","-rq","$5","-size","$6","-export_bclist"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"export_bc_sort")


        # demultiplex
        if "demultiplex" in self.settings.cmds_execute:
            outdir=self.settings.sampledir+"/demultiplex/_work"
            sh_cmd_list=["Interstellar-exec","demultiplex","-conf",self.cfgpath,"-d",outdir,"-o","$1","-cs","$2","-cq","$3","-rq","$4"]
            if self.settings.cfg["demultiplex"]["FORMAT"]=="tsv":
                sh_cmd_list.append("-export_tsv")
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"demultiplex")


        # tag
        if "tag" in self.settings.cmds_execute:
            outdir=self.settings.sampledir+"/tag/_work"
            sh_cmd_list=["Interstellar-exec","tag","-conf",self.cfgpath,"-d",outdir,"-o","$1","-cs","$2","-cq","$3","-rq","$4"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["TEMPLATE_SHELLSCRIPT"],sh_cmd_line,shelldir,"tag")
    
        subprocess.run("chmod u+x "+shelldir+"/*",shell=True)


    def fastq_split(self):
        used_commands=[]

        print("Runnning qsub jobs...: Split FASTQ files",flush=True)
        print(self.settings.sampledir)
        qoption=self.settings.qcfg["QOPTION"]
        qoption=qoption.replace("<mem>",self.settings.qcfg["MEM_MAX"])
        qcmd_base=["qsub",qoption,"-e",self.settings.sampledir+"/qlog","-o",self.settings.sampledir+"/qlog","-cwd"]
        for i in glob.glob(self.shelldir+"/seqkit*"):
            qcmd_now=qcmd_base+["-N","FASTQ_split"+self.today_now,i]
            qcmd_now=" ".join(qcmd_now)
            print(qcmd_now+"\n")
            s=subprocess.run(qcmd_now,shell=True)
            used_commands.append(qcmd_now)
            if s.returncode != 0:
                print("streamline: qsub for seqkit split failed.', file=sys.stderr")
                sys.exit(1)
        time.sleep(15)
        self.today_now=self.today_now
