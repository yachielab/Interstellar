from . import settingImporter
from . import barcodeConverter
from . import settingRequirementCheck
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import time
import os
import datetime
import csv

class settings_annotate_header(object):
    def __init__(self,opt):
        self.opt=opt
    
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_annotate_header=settingImporter.readconfig(cfgPath)["annotate_header"]
        cfg_annotate_header=settingImporter.configClean(cfg_annotate_header)

        cfg_raw=settingImporter.readconfig(cfgPath)
        cfg_raw={k:settingImporter.configClean(cfg_raw[k]) for k in cfg_raw}
        cfg_raw=settingRequirementCheck.setDefaultConfig(cfg_raw)
        cfg_value_ext,dict_to_terminal=settingImporter.config_extract_value_ext(cfg_raw)
        
        self.exportReadStructure={}
        self.annotate_headers={}
        for i in cfg_annotate_header:
            if i in ["READ1_STRUCTURE","READ2_STRUCTURE","INDEX1_STRUCTURE","INDEX2_STRUCTURE"] and cfg_annotate_header.get(i):
                self.exportReadStructure[i]=[dict_to_terminal[x] for x in cfg_annotate_header[i].split("+")]
            elif i in ["READ1_TAG","READ2_TAG","INDEX1_TAG","INDEX2_TAG"] and cfg_annotate_header.get(i):
                self.annotate_headers[i.split("_")[0]+"_STRUCTURE"]=cfg_annotate_header[i].split(",")
            
        self.path_to_seq=self.opt.correctedSeq
        self.path_to_avg_qval=self.opt.correctedQual
        self.path_to_rawQual=self.opt.rawQual
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=outdir+"/"+outname

class BARISTA_annotate_header(object):
    def __init__(self,settings):
        self.settings=settings
    def annotate_header(self):
        s_seq=pd.read_csv(self.settings.path_to_seq,sep='\t',dtype=str,chunksize=500000)
        s_avg_qual=pd.read_csv(self.settings.path_to_avg_qval,sep='\t',dtype=str,chunksize=500000)
        s_raw_qual=pd.read_csv(self.settings.path_to_rawQual,sep="\t",dtype=str,chunksize=500000,quoting=csv.QUOTE_NONE)

        cnt_chunk=0
        for s_seq_chunk,s_avg_qual_chunk,s_raw_qual_chunk in zip(s_seq,s_avg_qual,s_raw_qual):
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
                if exportRead in self.settings.annotate_headers:
                    annotate_header_now=self.settings.annotate_headers[exportRead]
                else:
                    annotate_header_now=""
                structure_now=self.settings.exportReadStructure[exportRead]

                export_pd=pd.DataFrame()

                #Deal with non-tag reads
                if not annotate_header_now=="":
                    export_pd["Header"]=s_seq_chunk["Header"].str.cat(s_seq_chunk[annotate_header_now],sep="_")
                else:
                    export_pd["Header"]=s_seq_chunk["Header"]

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

                if exportRead=="READ1_STRUCTURE":
                    read_iden="R1"
                elif exportRead=="READ2_STRUCTURE":
                    read_iden="R2"
                elif exportRead=="INDEX1_STRUCTURE":
                    read_iden="I1"
                elif exportRead=="INDEX2_STRUCTURE":
                    read_iden="I2"

                if cnt_chunk==0:
                    export_pd.to_csv(self.settings.outFilePath_and_Prefix+"_"+read_iden+".fastq.gz",mode="w",compression="gzip",sep="\t",index=False,header=False)
                else:
                    export_pd.to_csv(self.settings.outFilePath_and_Prefix+"_"+read_iden+".fastq.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
            cnt_chunk+=1
            