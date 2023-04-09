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
        self.ncore = int(self.opt.ncore)

        
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
                fastqDict[readKey]=segmentImporter.splitSequenceGenerator(fq_path,self.settings,chunksize=2000000,from_flash=True)
            else:
                fastqDict[readKey]=segmentImporter.splitSequenceGenerator(fq_path,self.settings,chunksize=2000000)
        self.fastqDict=fastqDict


    def extractComponents(self):
        parsedSeqDict=collections.defaultdict(list)
        parsedQualDict=collections.defaultdict(list)
        counterDict={}
        numSeqDict=collections.defaultdict(int)
        readKeys = []
        tmpdir=self.settings.tmpdir
        headerSplitRegex=regex.compile(r"[ |\t]+")
        n_chunk_max = 0

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
            # Chunks will be again split by 4 lines in the wrapper function to handle one record with a single job
            for n_chunk,fastq_chunk in enumerate(self.fastqDict[readKey]):

                # Making tmp directory for each chunk
                # This is for the sake of reducing the number of files in a directory
                outdir_tmp = tmpdir+"/Chunk"+str(n_chunk)
                os.makedirs(outdir_tmp,exist_ok = True)

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
            # del n_chunk

        self.counterDict=counterDict
        # self.today_now=today_now
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

