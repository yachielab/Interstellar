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
        self.seq=self.opt.rawSeq
        self.qual=self.opt.rawQual
        self.barcodes=func_dict["barcode_list"]
        self.qscore_dict=settingImporter.getQscoreDict(func_dict)
        self.ncore = int(self.opt.ncore)
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=outdir+"/"+outname
        
class BARISTA_QC(object):
    def __init__(self,settings):
        self.settings=settings
    def qualityCheck(self):
        seq_raw=pd.read_csv(self.settings.seq, sep='\t',chunksize=1000000)
        seq_qual=pd.read_csv(self.settings.qual, sep="\t",chunksize=1000000,quoting=csv.QUOTE_NONE)

        counterDict={}
        for n_chunk,seq_qual_zip in enumerate(zip(seq_raw,seq_qual)):
            seq_chunk=seq_qual_zip[0]
            qual_chunk=seq_qual_zip[1]
            
            # Quality filtering
            seq_chunk = segmentImporter.qfilter_parallel_wrapper(seq_chunk, qual_chunk, self.settings.qc_targets, self.settings.qscore_dict, self.settings.ncore)
            
            # Re-creating the count up data
            for col in self.settings.barcodes:
                counterDict_tmp=collections.Counter(seq_chunk[col])
                if n_chunk == 0:
                    counterDict[col] = counterDict_tmp
                else:
                    counterDict[col].update(counterDict_tmp)
            
            # Export the quality filtered segment table to TSV
            if n_chunk==0:
                seq_chunk.to_csv(self.settings.outFilePath_and_Prefix+"_srcSeq.QC.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
            else:
                seq_chunk.to_csv(self.settings.outFilePath_and_Prefix+"_srcSeq.QC.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
        
        # Export the count up data
        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCount.QC.pkl.gz",mode="wb") as p:
            pickle.dump(counterDict,p)
        
        # Copy the quality score table from the import folder to the qc folder to make the file structure consistent
        shutil.copyfile(self.settings.qual,self.settings.outFilePath_and_Prefix+"_srcQual.QC.tsv.gz")
