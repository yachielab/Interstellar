import copy
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
import pandas as pd
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
        self.readPathDict = readPathDict
        
        regexDictCompiled={}
        for readKey in regexDict:
            regexDictCompiled[readKey]=[regex.compile(i) for i in regexDict[readKey]]
        self.regexDictCompiled=regexDictCompiled
        func_dict=settingImporter.func_check(cfg_value_ext)
        self.barcodes=func_dict["barcode_list"]
        self.ncore = int(self.opt.ncore)
        self.flash_min=str(cfg_value_ext["FLASH_MIN_OVERLAP"])
        self.flash_max=str(cfg_value_ext["FLASH_MAX_OVERLAP"])
        self.flash_core = str(cfg["general"]["NUM_CORES"])

        
class BARISTA_IMPORT(object):
    def __init__(self,settings):
        self.settings=settings

    def importFastq(self):
        fastqDict={}
        for readKey in self.settings.readPathDict:
            fq_path=self.settings.readPathDict[readKey]
            # if readKey in self.settings.flash_gzipped_reads:
            #     fastqDict[readKey]=segmentImporter.sequenceGenerator(fq_path,self.settings,from_flash=True)
            # else:
            #     fastqDict[readKey]=segmentImporter.sequenceGenerator(fq_path,self.settings)
            if readKey in self.settings.flash_gzipped_reads:
                fastqDict[readKey]=segmentImporter.splitSequenceGenerator(fq_path,self.settings,chunksize=2000000,from_flash=True)
            else:
                fastqDict[readKey]=segmentImporter.splitSequenceGenerator(fq_path,self.settings,chunksize=2000000)
        self.fastqDict=fastqDict


    def extractComponents(self):
        counterDict={}
        numSeqDict=collections.defaultdict(int)
        readKeys = []
        tmpdir=self.settings.tmpdir
        headerSplitRegex=regex.compile(r"[ |\t]+")
        n_chunk_max = 0

        if not self.settings.flash:
            # Go through READ1, READ2, ...
            for nread,readKey in enumerate(self.settings.readPathDict):
                print("Extracting segments: "+readKey,flush=True)
                segment_parsed=segmentImporter.parseSegmentFromRegex(self.settings.regexDict[readKey])
                segment_parsed_set=set(segment_parsed)
                readKeys.append(readKey)
                
                # Go through chunks of FASTQ
                # Chunks will be again split by 4 lines in the wrapper function to handle one record with a single job
                for n_chunk,fastq_chunk in enumerate(self.fastqDict[readKey]):

                    # Making tmp directory for each chunk
                    # This is for the sake of reducing the number of files in a directory
                    outdir_tmp = tmpdir+"/Chunk"+str(n_chunk)
                    os.makedirs(outdir_tmp, exist_ok = True)

                    print("Processing file chunk",n_chunk)
                    n_records,counterDict_tmp = segmentImporter.segmentation_parallel_wrapper(
                        fastq_chunk= fastq_chunk,
                        settings= self.settings,
                        headerSplitRegex= headerSplitRegex,
                        readKey= readKey,
                        segment_parsed= segment_parsed,
                        outdir= outdir_tmp,
                        ncore = self.settings.ncore)
                
                    for i in self.settings.barcodes:
                        if counterDict_tmp.get(i):
                            if i in counterDict:
                                counterDict[i].update(counterDict_tmp[i])
                            else:
                                counterDict[i]=copy.deepcopy(counterDict_tmp[i])
                    
                    numSeqDict[readKey] += n_records
                    print(500000*(n_chunk)+n_records,"reads were processed for",readKey,flush=True)                
                    
                
                # Sequence number check - Interstellar requires the paired end reads and index reads are all sorted and correspond each other.
                if numSeqDict[readKey] != numSeqDict[readKeys[0]]:
                    if not self.settings.flash:
                        errmsg="Numbers of sequences between input files are inconsistent! Please check all the sequences are sorted in the same order across the input files."
                        raise InputError(errmsg)

                if n_chunk >= n_chunk_max:
                    self.n_chunk = n_chunk
                    n_chunk_max = n_chunk


        # FLASH
        if self.settings.flash:
            flash=self.settings.flash.replace("READ","Read")
            flash=flash.replace("INDEX","Index")
            merge_reads=[i+"_src" for i in flash.split("-")]

            # Go through READ1, READ2, ... and store FLASH-merged reads as fastq.
            for readKey in merge_reads:
                # print("Extracting segments: "+readKey,flush=True)
                readKeys.append(readKey)
                
                # Go through chunks of FASTQ
                # Chunks will be again split by 4 lines in the wrapper function to handle one record with a single job
                for n_chunk,fastq_chunk in enumerate(self.fastqDict[readKey]):

                    # Making tmp directory for each chunk
                    # This is for the sake of reducing the number of files in a directory
                    outdir_tmp = tmpdir+"/Chunk"+str(n_chunk)
                    os.makedirs(outdir_tmp, exist_ok = True)
                    os.makedirs(outdir_tmp+"/FLASH", exist_ok = True)

                    if readKey in merge_reads:
                        with open(outdir_tmp+"/FLASH/for_flash_"+readKey+".fastq", mode = "wt") as w:
                            for l in fastq_chunk:
                                w.write(l+"\n")
                        # numSeqDict[readKey] += len(fastq_chunk) / 4
                        continue # On hold for reads that will be FLASHed later
                    
                    # # Non-FLASH reads are processed as they are
                    # print("Processing file chunk",n_chunk)
                    # n_records,counterDict_tmp = segmentImporter.segmentation_parallel_wrapper(
                    #     fastq_chunk= fastq_chunk,
                    #     settings= self.settings,
                    #     headerSplitRegex= headerSplitRegex,
                    #     readKey= readKey,
                    #     segment_parsed= segment_parsed,
                    #     outdir= outdir_tmp,
                    #     ncore = self.settings.ncore)
                
                    # for i in self.settings.barcodes:
                    #     if counterDict_tmp.get(i):
                    #         if i in counterDict:
                    #             counterDict[i].update(counterDict_tmp[i])
                    #         else:
                    #             counterDict[i]=copy.deepcopy(counterDict_tmp[i])
                    
                    # numSeqDict[readKey] += n_records
                    # print(500000*(n_chunk)+n_records,"reads were processed for",readKey,flush=True)                
                
                # Store number of chunks
                if n_chunk >= n_chunk_max:
                    self.n_chunk = copy.deepcopy(n_chunk)
                    n_chunk_max = copy.deepcopy(n_chunk)

                # Sequence number check - Interstellar requires the paired end reads and index reads are all sorted and correspond each other.
                if numSeqDict[readKey] != numSeqDict[readKeys[0]]:
                    if not self.settings.flash:
                        errmsg="Numbers of sequences between input files are inconsistent! Please check all the sequences are sorted in the same order across the input files."
                        raise InputError(errmsg)
            



            # FLASH and segmentation for each chunk
            for chunk_now in range(self.n_chunk+1):
                outdir = tmpdir+"/Chunk"+str(chunk_now)
                FLASHed_filepath = segmentImporter.run_flash2(outdir+"/FLASH/for_flash_"+merge_reads[0]+".fastq", outdir+"/FLASH/for_flash_"+merge_reads[1]+".fastq", self.settings.flash_min, self.settings.flash_max, outdir+"/FLASH", self.settings.flash_core)

                for idx,readKey in enumerate(merge_reads+["merge_src"]):
                    fastq_chunk = []
                    with gzip.open(FLASHed_filepath[idx], mode = "rt") as r:
                        for l in r:
                            l = l.replace("\n", "")
                            fastq_chunk.append(l)
                    
                    segment_parsed=segmentImporter.parseSegmentFromRegex(self.settings.regexDict[readKey])
                    print("Processing file chunk",chunk_now)
                    n_records,counterDict_tmp = segmentImporter.segmentation_parallel_wrapper(
                        fastq_chunk= fastq_chunk,
                        settings= self.settings,
                        headerSplitRegex= headerSplitRegex,
                        readKey= readKey,
                        segment_parsed= segment_parsed,
                        outdir= outdir,
                        ncore = self.settings.ncore)
                
                    for i in self.settings.barcodes:
                        if counterDict_tmp.get(i):
                            if i in counterDict:
                                counterDict[i].update(counterDict_tmp[i])
                            else:
                                counterDict[i]=copy.deepcopy(counterDict_tmp[i])
                    
                    print(n_records)
                    numSeqDict[readKey] += n_records
                    print(numSeqDict[readKey],"reads were processed for",readKey,"\n",flush=True)

                # Extract and store headers to sort the other non-FLASH reads
                header_original = segmentImporter.get_header(outdir+"/FLASH/for_flash_"+merge_reads[0]+".fastq", headerSplitRegex, gz=False)
                orig2idx_dict = {val:idx for idx, val in enumerate(header_original)}

                # Fastq headers need to be unique
                if not len(orig2idx_dict) == len(header_original):
                    raise settingImporter.InputError("FASTQ headers need to be unique")

                header_ordered_uncomb = segmentImporter.get_header(FLASHed_filepath[0], headerSplitRegex, gz=True)
                header_ordered_idx_uncomb = [orig2idx_dict[i] for i in header_ordered_uncomb]
                with open("_".join([outdir+"/FLASH/headers.uncomb.pkl"]),mode="wb") as p:
                    pickle.dump(header_ordered_idx_uncomb, p)

                header_ordered_merged = segmentImporter.get_header(FLASHed_filepath[2], headerSplitRegex, gz=True)
                header_ordered_idx_merged = [orig2idx_dict[i] for i in header_ordered_merged]
                with open("_".join([outdir+"/FLASH/headers.merged.pkl"]),mode="wb") as p:
                    pickle.dump(header_ordered_idx_merged, p)


            

            # Go through READ1, READ2, ... which are not used for read merging.
            for readKey in self.settings.readPathDict:
                if readKey in merge_reads:
                    continue

                print("Extracting segments: "+readKey,flush=True)
                segment_parsed=segmentImporter.parseSegmentFromRegex(self.settings.regexDict[readKey])
                readKeys.append(readKey)
                
                # Go through chunks of FASTQ
                # Chunks will be again split by 4 lines in the wrapper function to handle one record with a single job
                for n_chunk,fastq_chunk in enumerate(self.fastqDict[readKey]):
                    # Non-FLASH reads are processed as they are
                    print("Processing file chunk",n_chunk)
                    outdir_tmp = tmpdir+"/Chunk"+str(n_chunk)
                    header_ordered_idx_uncomb = pd.read_pickle(outdir_tmp + "/FLASH/headers.uncomb.pkl")
                    header_ordered_idx_merged = pd.read_pickle(outdir_tmp + "/FLASH/headers.merged.pkl")
                                        
                    n_records,counterDict_tmp = segmentImporter.segmentation_parallel_wrapper(
                        fastq_chunk= fastq_chunk,
                        settings= self.settings,
                        headerSplitRegex= headerSplitRegex,
                        readKey= readKey,
                        segment_parsed= segment_parsed,
                        outdir= outdir_tmp,
                        ncore = self.settings.ncore,
                        sort_idx_order_uncomb=header_ordered_idx_uncomb,
                        sort_idx_order_merged=header_ordered_idx_merged)
                
                    for i in self.settings.barcodes:
                        if counterDict_tmp.get(i):
                            if i in counterDict:
                                counterDict[i].update(counterDict_tmp[i])
                            else:
                                counterDict[i]=copy.deepcopy(counterDict_tmp[i])

                    numSeqDict[readKey] += n_records
                    print(numSeqDict[readKey],"reads were processed for",readKey,"\n",flush=True)

                # Store number of chunks
                if n_chunk >= n_chunk_max:
                    self.n_chunk = copy.deepcopy(n_chunk)
                    n_chunk_max = copy.deepcopy(n_chunk)

                # Sequence number check - Interstellar requires the paired end reads and index reads are all sorted and correspond each other.
                if numSeqDict[readKey] != numSeqDict[readKeys[0]]:
                    if not self.settings.flash:
                        errmsg="Numbers of sequences between input files are inconsistent! Please check all the sequences are sorted in the same order across the input files."
                        raise InputError(errmsg)
                



        self.counterDict=counterDict
        self.tmpdir=tmpdir+"/"


    def exportExtractedComponents(self):
        # iter_num=-1
        for chunk_now in range(self.n_chunk+1):
            print("Merging file chunk",str(chunk_now),flush=True)

            segmentImporter.merge_parsed_data_parallel_wrapper(
                input_dir= self.tmpdir+"/Chunk"+str(chunk_now),
                settings= self.settings,
                ncore = self.settings.ncore,
                n_chunk= chunk_now
            )

        # shutil.rmtree(self.tmpdir)
        print("Exporting count dictionary...",flush=True)
        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCount.pkl.gz",mode="wb") as p:
            pickle.dump(self.counterDict,p)

