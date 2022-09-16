from . import settingImporter
from . import segmentImporter
from . import settingRequirementCheck
import csv
import pandas as pd
import collections
import gzip
import pickle
import shutil

class settings_qc(object):
    def __init__(self,opt):
        self.opt=opt

    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext,dict_to_terminal = settingImporter.config_extract_value_ext(cfg)
        func_dict=settingImporter.func_check(cfg_value_ext)
        self.qc_targets=[func_dict[val]["QUALITY_FILTER"]["source"] for val in func_dict if "QUALITY_FILTER" in func_dict[val]]
        if self.opt.mode_local:
            self.seq_list = settingImporter.parseInputFileList(self.opt.rawSeq)
            self.qual_list = settingImporter.parseInputFileList(self.opt.rawQual)
            outname = settingImporter.parseInputFileList(self.opt.outname)
        else:
            self.seq_list = [self.opt.rawSeq]
            self.qual_list = [self.opt.rawQual]
            outname = [self.opt.outname]

        self.barcodes=func_dict["barcode_list"]
        self.qscore_dict=settingImporter.getQscoreDict(func_dict)
        self.ncore = int(self.opt.ncore)
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix_list=[outdir+"/"+i for i in outname]

        
class BARISTA_QC(object):
    def __init__(self,settings):
        self.settings=settings
    def qualityCheck(self):
        n_chunk = 0
        for seq,qual,outprefix in zip(self.settings.seq_list,self.settings.qual_list,self.settings.outFilePath_and_Prefix_list):
            seq_chunk=pd.read_pickle(seq)
            qual_chunk=pd.read_pickle(qual)
            counterDict={}

            print("Processing file chunk",n_chunk)
            n_chunk += 1
            
            # Quality filtering
            seq_chunk = segmentImporter.qfilter_parallel_wrapper(seq_chunk, qual_chunk, self.settings.qc_targets, self.settings.qscore_dict, self.settings.ncore)
            
            # Re-creating the count up data
            for col in self.settings.barcodes:
                counterDict[col] = collections.Counter(seq_chunk[col])
            
            # Export the quality filtered segment table to TSV
            seq_chunk.to_pickle(outprefix+"_srcSeq.QC.pkl")
                            
            # Export the count up data
            with gzip.open(outprefix+"_srcCount.QC.pkl.gz",mode="wb") as p:
                pickle.dump(counterDict,p)
                        