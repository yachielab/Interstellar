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
import glob
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
        cfg_value_ext,dict_to_terminal = settingImporter.config_extract_value_ext(cfg)
        
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

        regerx_arg_list=[]
        for pat in cfg_value_ext:
            if regex.search(r"^READ1_STRUCTURE|^READ2_STRUCTURE|^INDEX1_STRUCTURE|^INDEX2_STRUCTURE|^READ_FLASH_STRUCTURE",pat):
                regerx_arg_list.append(pat)
        regerx_arg_list = sorted(regerx_arg_list)

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

        self.flash_gzipped_reads,self.src_readPathDict=segmentImporter.merge_reads_flash2(readPathDict,self.flash,self.input_fastq_gzipped,tmpdir,cfg_value_ext,cfg["general"]["NUM_CORES"])
        
        regexDictCompiled={}
        for readKey in regexDict:
            regexDictCompiled[readKey]=[regex.compile(i) for i in regexDict[readKey]]
        self.regexDictCompiled=regexDictCompiled
        func_dict=settingImporter.func_check(cfg_value_ext)
        self.barcodes=func_dict["barcode_list"]
        self.ncore = int(cfg["general"]["NUM_CORES"])

        
class BARISTA_IMPORT(object):
    def __init__(self,settings):
        self.settings=settings

    def importFastq(self):
        fastqDict={}
        for readKey in self.settings.src_readPathDict:
            fq_path=self.settings.src_readPathDict[readKey]
            # if readKey in self.settings.flash_gzipped_reads:
            #     fastqDict[readKey]=segmentImporter.sequenceGenerator(fq_path,self.settings,from_flash=True)
            # else:
            #     fastqDict[readKey]=segmentImporter.sequenceGenerator(fq_path,self.settings)
            if readKey in self.settings.flash_gzipped_reads:
                fastqDict[readKey]=segmentImporter.splitSequenceGenerator(fq_path,self.settings,chunksize=4000000,from_flash=True)
            else:
                fastqDict[readKey]=segmentImporter.splitSequenceGenerator(fq_path,self.settings,chunksize=4000000)
        self.fastqDict=fastqDict


    def extractComponents(self):
        parsedSeqDict=collections.defaultdict(list)
        parsedQualDict=collections.defaultdict(list)
        counterDict={}
        numSeqDict={}
        readKeys = []
        tmpdir=self.settings.tmpdir
        headerSplitRegex=regex.compile(r"[ |\t]+")
        
        # Go through READ1, READ2, ...
        for nread,readKey in enumerate(self.settings.src_readPathDict):
            print("Extracting segments: "+readKey,flush=True)
            segment_parsed=segmentImporter.parseSegmentFromRegex(self.settings.regexDict[readKey])
            segment_parsed_set=set(segment_parsed)
            readKeys.append(readKey)
            
            if self.settings.flash:
                merge_components=segmentImporter.parseSegmentFromRegex(self.settings.regexDict["merge_src"])
                merge_components_set=set(merge_components)

            # Go through chunks of FASTQ
            # Chunks are again split by 4 lines to handle one record with a single job
            for n_chunk,fastq_chunk in enumerate(self.fastqDict[readKey]):
                print("Processing file chunk",n_chunk)
                parsedSeqDict,parsedQualDict = segmentImporter.segmentation_parallel_wrapper(
                    fastq_chunk= fastq_chunk,
                    settings= self.settings,
                    headerSplitRegex= headerSplitRegex,
                    readKey= readKey,
                    segment_parsed= segment_parsed,
                    segment_parsed_set= segment_parsed_set,
                    ncore = self.settings.ncore)
            
            
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
                
                print(1000000*(n_chunk)+len(parsedSeqDict["Header"]),"reads were processed for",readKey,flush=True)

                # if (nrow+1)%4000000==0:
                #     n_chunk=int((nrow+1)/4000000) #n_chunk>0
                #     with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcSeq.pkl"]),mode="wb") as p:
                #         pickle.dump(parsedSeqDict,p)
                #     with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcQual.pkl"]),mode="wb") as p:
                #         pickle.dump(parsedQualDict,p)


                #     for i in self.settings.barcodes:
                #         if parsedSeqDict.get(i):
                #             if i in counterDict:
                #                 counterDict_tmp=collections.Counter(parsedSeqDict[i])
                #                 counterDict[i].update(counterDict_tmp)
                #             else:
                #                 counterDict[i]=collections.Counter(parsedSeqDict[i])

                #     parsedSeqDict=collections.defaultdict(list)
                #     parsedQualDict=collections.defaultdict(list)
                #     print(str(int((nrow+1)/4))+" reads were processed for "+readKey,flush=True)
                
            numSeqDict[readKey] = len(parsedSeqDict["Header"])
            
            # Sequence number check - Interstellar requires the paired end reads and index reads are all sorted and correspond each other.
            if len(parsedSeqDict["Header"]) != numSeqDict[readKeys[0]]:
                if not self.settings.flash:
                    errmsg="Numbers of sequences between input files are inconsistent! Please check all the sequences are sorted in the same order across the input files."
                    raise InputError(errmsg)

            # try:
            #     n_chunk
            #     n_chunk+=1
            # except NameError:
            #     n_chunk=0

            # try:
            #     nrow
            # except NameError:
            #     nrow=0
            
            # if not (nrow+1)%4000000==0 and nrow != 0:
            #     for i in self.settings.barcodes:
            #         if parsedSeqDict.get(i):
            #             counterDict_tmp=collections.Counter(parsedSeqDict[i])
            #             if n_chunk==0:
            #                 counterDict[i]=counterDict_tmp
            #             else:
            #                 counterDict[i].update(counterDict_tmp)
                
            #     with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcSeq.pkl"]),mode="wb") as p:
            #         pickle.dump(parsedSeqDict,p)
            #     with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcQual.pkl"]),mode="wb") as p:
            #         pickle.dump(parsedQualDict,p)
            #     print(str(int((nrow+1)/4))+" reads were processed for "+readKey,flush=True)
            
            # # # print(parsedSeqDict["Header"][:5])
            # # print("_".join([tmpdir+"/",readKey,str(n_chunk),"srcSeq.pkl"]),"\n\n")
            # parsedSeqDict=collections.defaultdict(list)
            # parsedQualDict=collections.defaultdict(list)
            self.n_chunk=n_chunk
            # del n_chunk

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
                        for component in ["Header"]+self.settings.components:
                            if component not in dict_merged:
                                dict_merged[component]=["-"]*nrow_uncombined
                            # else:
                            #     print(parsedDict_chunk[component][:5])
                            #     print("ok\n")
                            dict_merged[component]+=parsedDict_chunk[component]                    
                
                    if iter_num==0 and n_read==0:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcSeq.tsv.gz",mode="wt",encoding="utf-8")
                    elif iter_num>0 and n_read==0:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcSeq.tsv.gz",mode="at",encoding="utf-8")
                    elif iter_num==0 and n_read==1:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcQual.tsv.gz",mode="wt",encoding="utf-8")
                    elif iter_num>0 and n_read==1:
                        f=gzip.open(self.settings.outFilePath_and_Prefix+"_srcQual.tsv.gz",mode="at",encoding="utf-8")

                    # d_order = {k:c for c,k in enumerate(["Header"]+self.settings.components)}
                    # print(d_order)
                    # print(dict_merged.keys())
                    # dict_merged_sorted=sorted(dict_merged.items(), key=lambda x: d_order[x[0]])

                    dict_merged_key=["Header"]+self.settings.components
                    dict_merged_val=[dict_merged[i] for i in ["Header"]+self.settings.components]
                    # if iter_num==0:
                    #     csvwriter_index=csv.writer(f,delimiter="\t")
                    #     csvwriter_index.writerow(dict_merged_key)
                    dict_merged_val=list(map(list,zip(*dict_merged_val)))
                    dict_merged_val=["\t".join(i) for i in dict_merged_val]
                    dict_merged_val="\n".join(dict_merged_val)+"\n"
                    if iter_num==0:
                        dict_merged_val="\t".join(dict_merged_key)+"\n"+dict_merged_val
                    f.write(dict_merged_val)
                    f.close()
                
        shutil.rmtree(self.tmpdir)
        print("Exporting count dictionary...",flush=True)
        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCount.pkl.gz",mode="wb") as p:
            pickle.dump(self.counterDict,p)
