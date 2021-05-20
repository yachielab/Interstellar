from . import settingImporter
from . import segmentImporter
import regex
import pandas as pd
import collections
import gzip
import pickle
import shutil

class settings_qc(object):
    def __init__(self,opt):
        self.opt=opt
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_qc=settingImporter.readconfig(cfgPath)["import"]
        cfg_qc=settingImporter.configClean(cfg_qc)
        self.qc_targets=cfg_qc["qc_targets"].split(",")
        self.seq=self.opt.rawSeq
        self.qual=self.opt.rawQual
        self.barcodes=cfg_qc["barcodes"].split(",")
        self.min_base_quality=int(cfg_qc["min_base_quality"])
        self.min_avg_quality=int(cfg_qc["min_avg_quality"])
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)
        
class BARISTA_QC(object):
    def __init__(self,settings):
        self.settings=settings
    def qualityCheck(self):
        print("Minimum base quality threshold:",self.settings.min_base_quality)
        print("Average base quality threshold:",self.settings.min_avg_quality)
        seq_raw=pd.read_csv(self.settings.seq, sep='\t',chunksize=1000000)
        seq_qual=pd.read_csv(self.settings.qual, sep="\t",chunksize=1000000)
        counterDict={}
        for n_chunk,seq_qual_zip in enumerate(zip(seq_raw,seq_qual)):
            print("quality filtering for chunk",n_chunk)
            seq_chunk=seq_qual_zip[0]
            qual_chunk=seq_qual_zip[1]
            for component in self.settings.qc_targets:
                bool_filtered=qual_chunk[component].apply(segmentImporter.qualityFiltering,min_base_quality=self.settings.min_base_quality,min_avg_quality=self.settings.min_avg_quality)
                seq_chunk[component][bool_filtered]="-"
                
            for col in self.settings.barcodes:
                counterDict_tmp=collections.Counter(seq_chunk[col])
                if n_chunk==0:
                    counterDict[col]=counterDict_tmp
                else:
                    counterDict[col].update(counterDict_tmp)
            
            if n_chunk==0:
                seq_chunk.to_csv(self.settings.outFilePath_and_Prefix+"_srcSeq.QC.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
            else:
                seq_chunk.to_csv(self.settings.outFilePath_and_Prefix+"_srcSeq.QC.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCount.QC.pkl.gz",mode="wb") as p:
            pickle.dump(counterDict,p)
        shutil.copyfile(self.settings.qual,self.settings.outFilePath_and_Prefix+"_srcQual.QC.tsv.gz")
