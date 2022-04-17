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

    for dir in ["value_extraction","value_translation","annotate_header","demultiplex"]:
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
            # shell_cmd+=" "+sampledir+"/"+dir+"/_work/bc_sort"
            shell_cmd+=" "+sampledir+"/"+dir+"/_work/export"
            subprocess.run(shell_cmd,shell=True)


def judgeEndFix(input_files):
    for n,f in enumerate(input_files):
        if n==0:
            if re.search(r"\.fastq$",f):
                endfix=".fastq"
            elif re.search(r"\.fastq\.gz$",f):
                endfix=".fastq.gz"
            elif re.search(r"\.fq$",f):
                endfix=".fq"
            elif re.search(r"\.fq\.gz$",f):
                endfix=".fq.gz"
            
        else:
            if re.search(r"\.fastq$",f):
                endfix_now=".fastq"
            elif re.search(r"\.fastq\.gz$",f):
                endfix_now=".fastq.gz"
            elif re.search(r"\.fq$",f):
                endfix_now=".fq"
            elif re.search(r"\.fq.\gz$",f):
                endfix_now=".fq.gz"
            
            #All inputs should have the same format
            if not endfix_now==endfix:
                raise ArguementError("File format should be same across all input files.")       
    return endfix


def fastqDirCheck(valid_paths):
    if valid_paths =="":
        raise EmptyError("At least 1 input directory or FASTQ file paths should be provided.")
    
    #collect all input files if the path is directory
    all_files=[]
    for p in valid_paths:
        if os.path.isdir(p):
            all_files+=glob.glob(p+"/*")
        else:
            all_files+=glob.glob(p)
    
    #judge suffix
    file_suffix=judgeEndFix(all_files)
    return file_suffix


def getFiles(read_valid):
    file_dict={}
    for read_now in read_valid:
        path_now=read_valid[read_now]
        if os.path.isdir(path_now):
            file_dict[read_now]=glob.glob(path_now+"/*")
        else:
            file_dict[read_now]=glob.glob(path_now)
    
    #judge suffix
    return file_dict


def generateShellTemplate(template,cmdline,shelloutdir,shelloutname):
    with open(os.path.expanduser(template),mode="rt") as r, open(os.path.expanduser(shelloutdir+"/"+shelloutname+".sh"),mode="wt") as w:
        for line in r:
            w.write(line)
        w.write("\n")
        w.write(cmdline)
        w.write("\n")


def jobCheck(jid,outdir,n_jobs):
    print("qacct -d 7 -j "+jid+" 2> /dev/null | grep -E 'failed|exit_status' > "+outdir+"/qlog.tmp")

    s=subprocess.run("qacct -d 7 -j "+jid+" 2> /dev/null | grep -E 'failed|exit_status' > "+outdir+"/qlog.tmp",shell=True)
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
        status_ok=jobCheck(jid_now,outdir,njobs)
        if status_ok:
            break
        time.sleep(20)
        
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
            read1_path=config_val_ext["READ1_PATH"]
            read2_path=config_val_ext["READ2_PATH"]
            index1_path=config_val_ext["INDEX1_PATH"]
            index2_path=config_val_ext["INDEX2_PATH"]

            tmp_order=["read1","read2","index1","index2"]
            valid_paths=[i for i in [read1_path,read2_path,index1_path,index2_path] if not i == ""]
            read_valid={k:"" for k in tmp_order}
            for n,i in enumerate([read1_path,read2_path,index1_path,index2_path]):
                if not i == "":
                    read_valid[tmp_order[n]]=i

            file_suffix=fastqDirCheck(valid_paths)
            tartget_file_dict=getFiles(read_valid)
            
            if not self.cfg["general"]["SAMPLESHEET"]=="":
                # multi-sample
                samplelist=[]
                with open(self.cfg["general"]["SAMPLESHEET"],mode="rt") as r:
                    for line in r:
                        line=line.replace("\n","")
                        if line!="":
                            samplelist.append(line.split("\t"))
                target_prefix_list=[l[0] for l in samplelist if l[1]==self.samplename]

                # if not self.qcfg == "":
                #     # qsub = True -> split files
                #     for t in target_prefix_list:
                #         shell_cmd="mkdir -p "+self.sampledir+"/filesplit/"+t
                #         subprocess.run(shell_cmd,shell=True)
            # else:
                # single-sample
                # all_inputs=[sorted(glob.glob(i+"/*")) for i in valid_paths]
                # target_prefix_list=[]
                # for read_tup in zip(*all_inputs):
                #     read_list=list(read_tup)
                #     read_list=[os.path.basename(i) for i in read_list]
                #     prefix = []
                #     for x in zip(*read_list):
                #         if len(set(x)) == 1:
                #             prefix.append(x[0])
                #         else:
                #             break
                #     target_prefix_list.append("".join(prefix))
                # if not self.qcfg == "":
                #     # qsub = True -> split files
                #     shell_cmd="mkdir -p "+self.sampledir+"/filesplit/"+"".join(prefix)
                #     subprocess.run(shell_cmd,shell=True)
                        
            # Store values
            self.valid_paths=valid_paths
            self.read_valid=read_valid
            self.file_suffix=file_suffix
            self.target_prefix_list=target_prefix_list
            self.tartget_file_dict=tartget_file_dict

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
        tartget_file_dict=self.settings.tartget_file_dict

        if self.is_qsub:
            #Generate shell scripts for file splitting by seqkit
            for prefix in self.settings.target_prefix_list:
                input_read_files=[]

                reads_available=[]
                for r in ["read1","read2","index1","index2"]:
                    if r in self.settings.read_valid and not self.settings.read_valid[r]=="":
                        input_file_pool=tartget_file_dict[r]
                        target_files=[i for i in input_file_pool if re.search(prefix,i)]
                        input_read_files.append(sorted(target_files)) # like [[a_R1-001.fq,a_R1-002.fq,b_R1-001.fq], [a_R2-001.fq,a_R2-002.fq,b_R2-001.fq]]
                        reads_available.append(r)

                for idx,fileset_tup in enumerate(zip(*input_read_files)):
                    if idx==0:
                        nreads=len(fileset_tup)
                    else:
                        if nreads != len(fileset_tup):
                            raise UnknownError("Number of reads should be same across input files.")

                    sh_cmd_list_template=["seqkit","split2","-s",str(self.settings.qcfg["NUM_READS"])]
                    # nfile_divmod=divmod(nreads,2)
                    fileindex=0
                    for read_idx,f in enumerate(fileset_tup):
                        read_now=reads_available[read_idx]
                        sh_cmd_list=sh_cmd_list_template+["-O",self.settings.sampledir+"/filesplit/"+read_now+"_"+prefix+"_"+str(idx)]+["-1",f]
                        sh_cmd_line=" ".join(sh_cmd_list)
                        generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"seqkit_split_"+prefix+"_"+str(idx)+"_"+read_now)



                    # for n in range(nfile_divmod[0]):
                    #     files_now=["-1",fileset_tup[2*n],"-2",fileset_tup[2*n+1]]
                    #     sh_cmd_list=sh_cmd_list_template+["-O",self.settings.sampledir+"/filesplit/"+prefix+"_"+str(idx)+"_"+str(fileindex)]+files_now
                    #     sh_cmd_line=" ".join(sh_cmd_list)
                    #     generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"seqkit_split_"+prefix+"_"+str(idx)+"_"+str(fileindex))
                    #     fileindex+=1
                    # if len(fileset_tup)<4:
                    #     for n in range(2*nfile_divmod[0],2*nfile_divmod[0]+nfile_divmod[1]):
                    #         sh_cmd_list=sh_cmd_list_template+["-O",self.settings.sampledir+"/filesplit/"+prefix+"_"+str(idx)+"_"+str(fileindex)]+[fileset_tup[n]]
                    #         sh_cmd_line=" ".join(sh_cmd_list)
                    #         generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"seqkit_split_"+prefix+"_"+str(idx)+"_"+str(fileindex))
                
                #get read identifier
                # read_iden_dict[prefix]={}
                # for n,r in enumerate(["read1","read2","index1","index2"]):
                #     if r in self.settings.read_valid and not self.settings.read_valid[r]=="":
                #         read_iden_dict[r]=os.path.basename(input_read_files[n]).replace(self.settings.file_suffix,"").replace(prefix,"")
        # else:
        #     #Generate shell scripts for file splitting by seqkit
        #     for prefix in self.settings.target_prefix_list:
        #         input_read_files=[]
        #         for r in ["read1","read2","index1","index2"]:
        #             if not self.settings.read_valid[r]=="":
        #                 input_file_pool=glob.glob(r+"/*") if os.path.isdir(r) else glob.glob(r)
        #                 target_files=[i for i in input_file_pool if re.search(prefix,i)]
        #                 input_read_files.append(target_files[0]) #just grab a representative file name for the sample
                
                #get read identifier
        #         read_iden_dict[prefix]={}
        #         for n,r in enumerate(["read1","read2","index1","index2"]):
        #             if not self.settings.read_valid[r]=="":
        #                 read_iden_dict[r]=os.path.basename(input_read_files[n]).replace(self.settings.file_suffix,"").replace(prefix,"")
        # self.read_iden_dict=read_iden_dict
 
        
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
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"import")

            # qc
            sh_cmd_list=["Interstellar-exec","qc","-conf",self.cfgpath,"-d",outdir+"/qc","-o","$1","-rs","$2","-rq","$3"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"qc")
            
            # to_bt
            sh_cmd_list=["Interstellar-exec","to_bt","-conf",self.cfgpath,"-d",outdir+"/to_bt","-o","$1","-rs","$2"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"to_bt")
            
            # correct
            sh_cmd_list=["Interstellar-exec","correct","-conf",self.cfgpath,"-d",outdir+"/correct","-o","$1","-ip","$2"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"correct")
            
            # mk_sval
            sh_cmd_list=["Interstellar-exec","mk_sval","-conf",self.cfgpath,"-d",outdir+"/mk_sval","-o","$1","-rs","$2","-rq","$3","-crp","$4"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"mk_sval")
        

        if "value_translation" in self.settings.cmds_execute:
            outdir=self.settings.sampledir+"/value_translation/_work/"

            #buildTree
            sh_cmd_list=["Interstellar-exec","buildTree","-conf",self.cfgpath,"-d",outdir+"/buildTree","-o","$1","-sv","$2"]
            if not self.settings.cfg["general"]["SAMPLESHEET"]=="":
                sh_cmd_list.append("-samplemerge")
                sh_cmd_list+=["-samplesheet",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/samplesheet/samplesheet.tsv"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"buildTree")
            
            #mergeTree
            sh_cmd_list=["Interstellar-exec","mergeTree","-conf",self.cfgpath,"-o","$1","-lp","$2"]
            if not self.settings.cfg["general"]["SAMPLESHEET"]=="":
                sh_cmd_list.append("-samplemerge")
                sh_cmd_list+=["-samplesheet",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/samplesheet/samplesheet.tsv","-d",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/mergeTree"]
            else:
                sh_cmd_list+=["-d",outdir+"/mergeTree"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"mergeTree")
            
            #convert
            sh_cmd_list=["Interstellar-exec","convert","-conf",self.cfgpath,"-d",outdir+"/convert","-o","$1","-tree","$2","-sv","$3","-sq","$4"]
            if not self.settings.cfg["general"]["SAMPLESHEET"]=="":
                sh_cmd_list.append("-samplemerge")
                sh_cmd_list+=["-samplesheet",self.settings.cfg["general"]["PROJECT_DIR"]+"/_multisample/samplesheet/samplesheet.tsv"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"convert")
            # bc_sort (cannot executed in multi-sample mode.)
            # sh_cmd_list=["Interstellar-exec","bc_sort","-conf",self.cfgpath,"-d",outdir+"/bc_sort","-o","$1","-tree","$2","-sseq_to_svalue","$3","-tbl","$4"]
            # sh_cmd_line=" ".join(sh_cmd_list)
            # generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"bc_sort")
            
            # export: normal
            sh_cmd_list=["Interstellar-exec","export","-conf",self.cfgpath,"-d",outdir+"/export","-o","$1","-dv","$2","-dq","$3","-rs","$4","-rq","$5","-size","$6","-export_bclist"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"export")
            # export: bc_sort
            # sh_cmd_list=["Interstellar-exec","export","-conf",outdir+"/bc_sort/sorted.conf","-d",outdir+"/export","-o","$1","-dv","$2","-dq","$3","-rs","$4","-rq","$5","-size","$6","-export_bclist"]
            # sh_cmd_line=" ".join(sh_cmd_list)
            # generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"export_bc_sort")


        # demultiplex
        if "demultiplex" in self.settings.cmds_execute:
            outdir=self.settings.sampledir+"/demultiplex/_work"
            sh_cmd_list=["Interstellar-exec","demultiplex","-conf",self.cfgpath,"-d",outdir,"-o","$1","-cs","$2","-cq","$3","-rq","$4"]
            if self.settings.cfg["demultiplex"]["FORMAT"]=="tsv":
                sh_cmd_list.append("-export_tsv")
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"demultiplex")


        # annotate_header
        if "annotate_header" in self.settings.cmds_execute:
            outdir=self.settings.sampledir+"/annotate_header/_work"
            sh_cmd_list=["Interstellar-exec","annotate_header","-conf",self.cfgpath,"-d",outdir,"-o","$1","-cs","$2","-cq","$3","-rq","$4"]
            sh_cmd_line=" ".join(sh_cmd_list)
            generateShellTemplate(self.settings.cfg["general"]["SET_SHELL_ENV"],sh_cmd_line,shelldir,"annotate_header")
    
        subprocess.run("chmod u+x "+shelldir+"/*",shell=True)


    def fastq_split(self):
        used_commands=[]

        print("\nRunnning qsub jobs...: Split FASTQ files for sample "+self.settings.samplename,flush=True)
        qoption=self.settings.qcfg["QOPTION"]
        qoption=qoption.replace("<mem>",self.settings.qcfg["MEM_MIN"])
        qcmd_base=["qsub",qoption,"-e",self.settings.sampledir+"/qlog","-o",self.settings.sampledir+"/qlog","-cwd"]
        for i in glob.glob(self.shelldir+"/seqkit*"):
            qcmd_now=qcmd_base+["-N","FASTQ_split"+self.today_now,i]
            qcmd_now=" ".join(qcmd_now)
            s=subprocess.run(qcmd_now,shell=True)
            used_commands.append(qcmd_now)
            if s.returncode != 0:
                print("streamline: qsub for seqkit split failed.', file=sys.stderr")
                sys.exit(1)
        time.sleep(15)
        self.today_now=self.today_now
