from . import settingRequirementCheck
from . import settingImporter
from . import segmentImporter
from . import interstellar_setup
import regex
import datetime
import os
import gzip
import collections
import pickle
import subprocess
import pandas as pd
import glob
import csv
import shutil


class InputError(Exception):
    pass


class settings_import(object):
    def __init__(self,opt):
        self.opt=opt
    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext = settingImporter.config_extract_value_ext(cfg)
        
        input_files=[i for i in [self.opt.config,self.opt.read1,self.opt.read2,self.opt.index1,self.opt.index2,self.opt.outdir] if not i==""]
        settingRequirementCheck.pathExistCheck(input_files)

        flg_file=0
        for i in [self.opt.read1,self.opt.read2,self.opt.index1,self.opt.index2]:
            if not i=="":
                if flg_file==0:
                    flg_gz=1 if regex.search(r"gz$",i) else 0
                    flg_file=1
                else:
                    flg_gz_tmp=1 if regex.search(r"gz$",i) else 0
                    if not flg_gz==flg_gz_tmp:
                        raise InputError("All files should be commonly gzipped or ungzipped.")
                    
        self.input_fastq_gzipped=True if flg_gz==1 else False
        # self.simple=self.opt.simple
        self.components=cfg_value_ext["segments"]
        self.flash=cfg_value_ext["FLASH"]
        self.outFilePath_and_Prefix=self.opt.outdir+"/"+self.opt.outname
        today_now=str(datetime.datetime.today())
        today_now=regex.sub(r"\.|:| ","-",today_now)
        random_string=segmentImporter.randomname(15)
        tmpdir="_".join([self.outFilePath_and_Prefix,today_now,random_string])
        os.mkdir(tmpdir)
        self.tmpdir=tmpdir

        regexDict=dict(Read1_src=[],Read2_src=[],Index1_src=[],Index2_src=[],merge_src=[])
        for pat in cfg_value_ext:
            if regex.search("^READ1_STRUCTURE",pat):
                regexDict["Read1_src"].append(cfg_value_ext[pat])
            elif regex.search("^READ2_STRUCTURE",pat):
                regexDict["Read2_src"].append(cfg_value_ext[pat])
            elif regex.search("^INDEX1_STRUCTURE",pat):
                regexDict["Index1_src"].append(cfg_value_ext[pat])
            elif regex.search("^INDEX2_STRUCTURE",pat):
                regexDict["Index2_src"].append(cfg_value_ext[pat])
            elif regex.search("^READ_FLASH_STRUCTURE",pat):
                regexDict["merge_src"].append(cfg_value_ext[pat])
        regexDict={k:regexDict[k] for k in regexDict if regexDict[k]}
        self.regexDict=regexDict
        reads_default=["Read1_src","Read2_src","Index1_src","Index2_src"]
        srcReadKeys=[]
        for cnt,readPath in enumerate([self.opt.read1,self.opt.read2,self.opt.index1,self.opt.index2]):
            if readPath!="":
                srcReadKeys.append(reads_default[cnt])
        srcReadPaths=[i for i in [self.opt.read1,self.opt.read2,self.opt.index1,self.opt.index2] if not i==""]
        readPathDict={k:v for k,v in zip(srcReadKeys,srcReadPaths) if v!=""}
        self.flash_gzipped_reads=[]
        self.flash_gzipped_reads,self.src_readPathDict=segmentImporter.merge_reads_flash2(readPathDict,self.flash,self.input_fastq_gzipped,tmpdir,cfg_value_ext)
        
        # if self.simple:
        #     patternDict={}
        #     for readKey in regexDict:
        #         pat_now=regexDict[readKey][0]
        #         pat_now=regex.sub(r"><","_",pat_now)
        #         pat_now=regex.sub(r"\<|\>","",pat_now).split("_")
        #         patternDict[readKey]={k.split(":")[0]:int(k.split(":")[1]) for k in pat_now}
        #     self.patternDict=patternDict
        # else:
        regexDictCompiled={}
        for readKey in regexDict:
            regexDictCompiled[readKey]=[regex.compile(i) for i in regexDict[readKey]]
        self.regexDictCompiled=regexDictCompiled
        func_dict=settingImporter.func_check(cfg_value_ext)
        self.barcodes=func_dict["barcode_list"]
        
        
class BARISTA_IMPORT(object):
    parsedSeqDict=None
    parsedQualDict=None
    def __init__(self,settings):
        self.settings=settings

    def importFastq(self):
        fastqDict={}
        for readKey in self.settings.src_readPathDict:
            fq_path=self.settings.src_readPathDict[readKey]
            if readKey in self.settings.flash_gzipped_reads:
                fastqDict[readKey]=segmentImporter.sequenceGenerator(fq_path,self.settings,from_flash=True)
            else:
                fastqDict[readKey]=segmentImporter.sequenceGenerator(fq_path,self.settings)
        self.fastqDict=fastqDict


    # def simple_extract(self):
    #     fastq_tab_dict={}
    #     counterDict={}
    #     tmpdir=self.settings.tmpdir
    #     for nread,readKey in enumerate(self.settings.src_readPathDict):
    #         fqpath_now=self.settings.src_readPathDict[readKey]
    #         tmp_tsv_path="_".join([tmpdir+"/",readKey])+".tsv.gz"
    #         cmdlist=["seqkit fx2tab",fqpath_now,"-o",tmp_tsv_path]
    #         cmd=" ".join(cmdlist)
    #         subprocess.run(cmd,shell=True)

    #         fastq_now=pd.read_csv(tmp_tsv_path,sep="\t",header=None,chunksize=500000)
    #         fastq_tab_dict[readKey]=fastq_now
    #     fastq_key=list(fastq_tab_dict.keys())
    #     fastq_val=list(fastq_tab_dict.values())
    #     for chunkcnt,chunklist in enumerate(zip(*fastq_val)):
    #         print("chunk",chunkcnt,"started...")
    #         export_seq=pd.DataFrame()
    #         export_qual=pd.DataFrame()
    #         for readcnt,fq_chunk in enumerate(chunklist):
    #             fq_chunk=fq_chunk.replace("[ |\t].+$","",regex=True)
    #             if readcnt==0:
    #                 export_seq["Header"]="@"+fq_chunk[0]
    #                 export_qual["Header"]="@"+fq_chunk[0]
    #             readkey_now=fastq_key[readcnt]
    #             pos=0
    #             for component in self.settings.patternDict[readkey_now]:
    #                 if not self.settings.patternDict[readkey_now][component]==0:
    #                     export_seq[component]=fq_chunk[1].str[pos:pos+self.settings.patternDict[readkey_now][component]]
    #                     export_qual[component]=fq_chunk[2].str[pos:pos+self.settings.patternDict[readkey_now][component]]
    #                     pos+=self.settings.patternDict[readkey_now][component]
    #                 else:
    #                     export_seq[component]=fq_chunk[1].str[pos:]
    #                     export_qual[component]=fq_chunk[2].str[pos:]
    #                     break
    #         if chunkcnt==0:
    #             export_seq.to_csv(self.settings.outFilePath_and_Prefix+"_srcSeq.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
    #             export_qual.to_csv(self.settings.outFilePath_and_Prefix+"_srcQual.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
    #         else:    
    #             export_seq.to_csv(self.settings.outFilePath_and_Prefix+"_srcSeq.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
    #             export_qual.to_csv(self.settings.outFilePath_and_Prefix+"_srcQual.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
            
    #         for i in self.settings.barcodes:
    #             if chunkcnt==0:
    #                 counterDict[i]=collections.Counter(export_seq[i])
    #             else:
    #                 counter_tmp=collections.Counter(export_seq[i])
    #                 counterDict[i].update(counter_tmp)
    #     shutil.rmtree(tmpdir)
    #     with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCount.pkl.gz",mode="wb") as p:
    #         pickle.dump(counterDict,p)


    def extractComponents(self):
        parsedSeqDict=collections.defaultdict(list)
        parsedQualDict=collections.defaultdict(list)
        counterDict={}
        tmpdir=self.settings.tmpdir
        headerSplitRegex=regex.compile(r"[ |\t]+")

        for nread,readKey in enumerate(self.settings.src_readPathDict):
            print(readKey,flush=True)
            segment_parsed=segmentImporter.parseSegmentFromRegex(self.settings.regexDict[readKey])
            segment_parsed_set=set(segment_parsed)
            
            if self.settings.flash:
                merge_components=segmentImporter.parseSegmentFromRegex(self.settings.regexDict["merge_src"])
                merge_components_set=set(merge_components)
            for nrow,line in enumerate(self.fastqDict[readKey]):
                if (nrow+1)%4==1:
                    header=headerSplitRegex.split(line)[0]
                    parsedSeqDict["Header"].append(header)
                    parsedQualDict["Header"].append(header)
                    
                if (nrow+1)%4==2:
                    m=segmentImporter.patMatching(line,self.settings.regexDictCompiled[readKey])
                    if m:
                        mdict=m.groupdict()
                        for seg in mdict:
                            parsedSeqDict[seg].append(mdict[seg])
                        component_diff=segment_parsed_set-set(mdict.keys())
                        for seg in component_diff:
                            parsedSeqDict[seg].append("-")
                            parsedQualDict[seg].append("-")
                    else:
                        for seg in segment_parsed:
                            parsedSeqDict[seg].append("-")
                            parsedQualDict[seg].append("-")

                if (nrow+1)%4==0:
                    if m:
                        for component in mdict:
                            extractedQual=line[m.span(component)[0]:m.span(component)[1]]
                            parsedQualDict[component].append(extractedQual)
            

                if (nrow+1)%4000000==0:
                    n_chunk=int((nrow+1)/4000000) #n_chunk>0
                    with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcSeq.pkl"]),mode="wb") as p:
                        pickle.dump(parsedSeqDict,p)
                    with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcQual.pkl"]),mode="wb") as p:
                        pickle.dump(parsedQualDict,p)


                    for i in self.settings.barcodes:
                        if parsedSeqDict.get(i):
                            if i in counterDict:
                                counterDict_tmp=collections.Counter(parsedSeqDict[i])
                                counterDict[i].update(counterDict_tmp)
                            else:
                                counterDict[i]=collections.Counter(parsedSeqDict[i])

                    parsedSeqDict=collections.defaultdict(list)
                    parsedQualDict=collections.defaultdict(list)
                    print(str(int((nrow+1)/4))+" reads were processed for "+readKey,flush=True)

            try:
                n_chunk
                n_chunk+=1
            except NameError:
                n_chunk=0

            if not (nrow+1)%4000000==0:
                for i in self.settings.barcodes:
                    if parsedSeqDict.get(i):
                        counterDict_tmp=collections.Counter(parsedSeqDict[i])
                        if n_chunk==0:
                            counterDict[i]=counterDict_tmp
                        else:
                            counterDict[i].update(counterDict_tmp)
                
                with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcSeq.pkl"]),mode="wb") as p:
                    pickle.dump(parsedSeqDict,p)
                with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcQual.pkl"]),mode="wb") as p:
                    pickle.dump(parsedQualDict,p)
                print(str(int((nrow+1)/4))+" reads were processed for "+readKey,flush=True)
            
            parsedSeqDict=collections.defaultdict(list)
            parsedQualDict=collections.defaultdict(list)
            self.n_chunk=n_chunk
            del n_chunk

        self.counterDict=counterDict
        # self.today_now=today_now
        self.tmpdir=tmpdir+"/"


    def exportExtractedComponents(self):
        iter_num=-1
        for chunk_now in range(self.n_chunk+1):
            filepaths_seq=glob.glob("_".join([self.tmpdir,"*",str(chunk_now),"srcSeq.pkl"]))
            filepaths_qual=glob.glob("_".join([self.tmpdir,"*",str(chunk_now),"srcQual.pkl"]))
            
            if filepaths_seq:
                iter_num+=1
            else:
                continue
                
            if filepaths_seq:
                filepaths_seq.sort()
                filepaths_qual.sort()
                filepaths=[filepaths_seq,filepaths_qual]
                merge_filename="_".join([self.tmpdir,"merge"])
                print("merging file for chunk",str(iter_num),flush=True)
                for n_read,path in enumerate(filepaths):
                    dict_merged=collections.defaultdict(list)
                    for chunk in path:
                        print("chunk:",chunk)
                        if self.settings.flash:
                            try:
                                to_be_processed
                            except NameError:
                                to_be_processed=chunk
                                continue
                        
                        with open(chunk,mode="rb") as pchunk:
                            parsedDict_chunk=pickle.load(pchunk)

                        dict_merged.update(parsedDict_chunk)
                        

                    if self.settings.flash:
                        nrow_uncombined=len(dict_merged["Header"])
                        with open(to_be_processed,mode="rb") as pchunk:
                            parsedDict_chunk=pickle.load(pchunk)
                        for component in self.settings.components:
                            if component not in dict_merged:
                                dict_merged[component]=["-"]*nrow_uncombined
                            dict_merged[component]+=parsedDict_chunk[component]
                    
                
                    if iter_num==0 and n_read==0:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcSeq.tsv.gz",mode="wt",encoding="utf-8")
                    elif iter_num>0 and n_read==0:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcSeq.tsv.gz",mode="at",encoding="utf-8")
                    elif iter_num==0 and n_read==1:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcQual.tsv.gz",mode="wt",encoding="utf-8")
                    elif iter_num>0 and n_read==1:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcQual.tsv.gz",mode="at",encoding="utf-8")

                    d_order = {k:c for c,k in enumerate(["Header"]+self.settings.components)}
                    dict_merged_sorted=sorted(dict_merged.items(), key=lambda x: d_order[x[0]])
                    dict_merged_key=[i[0] for i in dict_merged_sorted]
                    dict_merged_val=[i[1] for i in dict_merged_sorted]
                    if iter_num==0:
                        csvwriter_index=csv.writer(f,delimiter="\t")
                        csvwriter_index.writerow(dict_merged_key)
                    dict_merged_val=list(map(list,zip(*dict_merged_val)))
                    dict_merged_val=["\t".join(i) for i in dict_merged_val]
                    dict_merged_val="\n".join(dict_merged_val)+"\n"
                    f.write(dict_merged_val)
                    f.close()
                
        shutil.rmtree(self.tmpdir)
        print("exporting count dictionary...",flush=True)
        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCount.pkl.gz",mode="wb") as p:
            pickle.dump(self.counterDict,p)
