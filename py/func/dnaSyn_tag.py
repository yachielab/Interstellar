from . import settingImporter
from . import barcodeConverter
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import time
import os
import datetime

class settings_tag(object):
    def __init__(self,opt):
        self.opt=opt
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_tag=settingImporter.readconfig(cfgPath)["tag"]
        cfg_tag=settingImporter.configClean(cfg_tag)
        self.exportReadStructure={}
        self.tags={}
        for i in cfg_tag:
            if i in ["Read1","Read2","Index1","Index2"] and cfg_tag.get(i):
                self.exportReadStructure[i]=cfg_tag[i].split("+")
            elif i in ["tag_Read1","tag_Read2","tag_Index1","tag_Index2"] and cfg_tag.get(i):
                self.tags[i.split("_")[1]]=cfg_tag[i].split(",")
            
        self.path_to_seq=self.opt.correctedSeq
        self.path_to_avg_qval=self.opt.correctedQual
        self.path_to_rawQual=self.opt.rawQual
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)
class BARISTA_TAG(object):
    def __init__(self,settings):
        self.settings=settings
    def tag(self):
        s_seq=pd.read_csv(self.settings.path_to_seq,sep='\t',dtype=str,chunksize=500000)
        s_avg_qual=pd.read_csv(self.settings.path_to_avg_qval,sep='\t',dtype=str,chunksize=500000)
        s_raw_qual=pd.read_csv(self.settings.path_to_rawQual,sep="\t",dtype=str,chunksize=500000)

        cnt_chunk=0
        for s_seq_chunk,s_avg_qual_chunk,s_raw_qual_chunk in zip(s_seq,s_avg_qual,s_raw_qual):
            print("chunk",cnt_chunk)
            s_avg_qual_chunk=s_avg_qual_chunk.drop("Header",axis=1)
            s_raw_qual_chunk=s_raw_qual_chunk.drop("Header",axis=1)
            s_avg_qual_chunk=s_avg_qual_chunk.astype("int8")
            s_seq_chunk=s_seq_chunk[s_seq_chunk!="-"].dropna()
            s_avg_qual_chunk=s_avg_qual_chunk.loc[s_seq_chunk.index]
            s_raw_qual_chunk=s_raw_qual_chunk.loc[s_seq_chunk.index]

            colnames=list(s_seq_chunk.columns)
            raw_component_names=[i.split(":")[0] for i in colnames if not i=="Header"]
            corrected_component_names=[i.split(":")[1] for i in colnames if not i=="Header"]
            s_seq_chunk.columns=["Header"]+corrected_component_names

            for exportRead in self.settings.exportReadStructure:
                tag_now=self.settings.tags[exportRead]
                structure_now=self.settings.exportReadStructure[exportRead]

                export_pd=pd.DataFrame()
                export_pd["Header"]=s_seq_chunk["Header"].str.cat(s_seq_chunk[tag_now],sep="_")
                if len(structure_now)>1:
                    export_pd["seq"]=s_seq_chunk[structure_now[0]].str.cat(s_seq_chunk[structure_now[1:]],sep="")
                else:
                    export_pd["seq"]=s_seq_chunk[structure_now[0]]
                export_pd["3rd"]=["+"]*export_pd.shape[0]
                export_pd["qual"]=[""]*export_pd.shape[0]
                for component in structure_now:
                    if component in s_avg_qual_chunk:
                        df_seq_qual_tmp=s_seq_chunk[component].str.cat(s_avg_qual_chunk[component].astype(str),sep="_")
                        df_seq_qual_tmp=df_seq_qual_tmp.map(barcodeConverter.getConvQual_ver2)
                        export_pd["qual"]=export_pd["qual"].str.cat(df_seq_qual_tmp,sep="")
                    else:
                        component_raw=raw_component_names[corrected_component_names.index(component)]
                        export_pd["qual"]=export_pd["qual"].str.cat(s_raw_qual_chunk[component_raw],sep="")
                export_pd=export_pd.stack()
                export_pd=export_pd.reset_index()
                export_pd=pd.DataFrame(export_pd[0])
                if cnt_chunk==0:
                    export_pd.to_csv(self.settings.outFilePath_and_Prefix+exportRead+".fastq.gz",mode="w",compression="gzip",sep="\t",index=False,header=False)
                else:
                    export_pd.to_csv(self.settings.outFilePath_and_Prefix+exportRead+".fastq.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
            cnt_chunk+=1