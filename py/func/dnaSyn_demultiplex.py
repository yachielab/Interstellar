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

class settings_demultiplex(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_demulti=settingImporter.readconfig(cfgPath)["demultiplex"]
        cfg_demulti=settingImporter.configClean(cfg_demulti)
        self.key=cfg_demulti["key"].split(",")
        self.out_components=cfg_demulti["out_components"].split(",")
        self.exportReadStructure={}
        for i in cfg_demulti:
            if i in ["Read1","Read2","Index1","Index2"] and cfg_demulti.get(i):
                self.exportReadStructure[i]=cfg_demulti[i].split("+")
        self.path_to_seq=self.opt.correctedSeq
        self.path_to_avg_qval=self.opt.correctedQual
        self.path_to_rawQual=self.opt.rawQual
        self.export_tsv=self.opt.export_tsv
        outname=self.opt.outname
        outdir=self.opt.outdir

        today_now=str(datetime.datetime.today())
        today_now=regex.sub(r"\W","",today_now)
        new_outdir="_".join([regex.sub("/$","",str(outdir))+"/"+"barista_demulti",outname,today_now])
        flag=os.path.isdir(new_outdir)
        while flag:
            today_now=str(datetime.datetime.today())
            today_now=regex.sub(r"\W","",today_now)
            new_outdir="_".join([regex.sub("/$","",str(outdir))+"/"+"barista_demulti",outname,today_now])
            flag=os.path.isdir(new_outdir)
        os.mkdir(new_outdir)
        self.outFilePath_and_Prefix=new_outdir+"/"+str(outname)
        self.today_now=today_now

class BARISTA_DEMULTIPLEX(object):
    def __init__(self,settings):
        self.settings=settings
    def demultiplex(self):
        s_seq=pd.read_csv(self.settings.path_to_seq,sep='\t',dtype=str,chunksize=500000)
        s_avg_qual=pd.read_csv(self.settings.path_to_avg_qval,sep='\t',dtype=str,chunksize=500000)
        s_raw_qual=pd.read_csv(self.settings.path_to_rawQual,sep="\t",dtype=str,chunksize=500000)
        demulti_key=set()
        cnt_chunk=0
        for s_seq_chunk,s_avg_qual_chunk,s_raw_qual_chunk in zip(s_seq,s_avg_qual,s_raw_qual):
            print("start demultiplexing for chunk",cnt_chunk)
            s_avg_qual_chunk=s_avg_qual_chunk.drop("Header",axis=1)
            s_raw_qual_chunk=s_raw_qual_chunk.drop("Header",axis=1)

            
            colnames=list(s_seq_chunk.columns)
            raw_component_names=[i.split(":")[0] for i in colnames if not i=="Header"]
            corrected_component_names=[i.split(":")[1] for i in colnames if not i=="Header"]
            s_seq_chunk.columns=["Header"]+corrected_component_names            

            if self.settings.export_tsv:
                key_series=s_seq_chunk[self.settings.key].apply("_".join,axis=1)
                demulti_key|=set(key_series)
                s_seq_chunk=s_seq_chunk[["Header"]+self.settings.out_components]
                for eachkey in demulti_key:
                    if not "-" in eachkey:
                        export_pd_tmp=s_seq_chunk[key_series==eachkey]
                        outfilename="_".join([self.settings.outFilePath_and_Prefix,eachkey])+".tsv.gz"
                        if not os.path.isfile(outfilename):
                            export_pd_tmp.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=True)
                        else:
                            export_pd_tmp.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)
            else:
                s_seq_chunk=s_seq_chunk[s_seq_chunk!="-"].dropna()
                s_avg_qual_chunk=s_avg_qual_chunk.loc[s_seq_chunk.index]
                s_raw_qual_chunk=s_raw_qual_chunk.loc[s_seq_chunk.index]
                s_avg_qual_chunk=s_avg_qual_chunk.astype("int8")
                key_series=s_seq_chunk[self.settings.key].apply("_".join,axis=1)
                demulti_key|=set(key_series)
                export_pd=pd.DataFrame()
                export_pd["Header"]=s_seq_chunk["Header"]
                s_seq_chunk=s_seq_chunk.drop("Header",axis=1)
                    
                for exportReadNum in self.settings.exportReadStructure:
                    export_pd["3rd"]=["+"]*export_pd.shape[0]
                    export_pd["qual"]=[""]*export_pd.shape[0]
                    structure_now=self.settings.exportReadStructure[exportReadNum]
                    if len(structure_now)>1:
                        export_pd["seq"]=s_seq_chunk[structure_now[0]].str.cat(s_seq_chunk[structure_now[1:]],sep="")
                    else:
                        export_pd["seq"]=s_seq_chunk[structure_now[0]]

                    for component in structure_now:
                        if component in s_avg_qual_chunk:
                            df_seq_qual_tmp=s_seq_chunk[component].str.cat(s_avg_qual_chunk[component].astype(str),sep="_")
                            df_seq_qual_tmp=df_seq_qual_tmp.map(barcodeConverter.getConvQual_ver2)
                            export_pd["qual"]=export_pd["qual"].str.cat(df_seq_qual_tmp,sep="")
                        else:
                            component_raw=raw_component_names[corrected_component_names.index(component)]
                            export_pd["qual"]=export_pd["qual"].str.cat(s_raw_qual_chunk[component_raw],sep="")

                    export_pd=export_pd[["Header","seq","3rd","qual"]]
                    for eachkey in demulti_key:
                        export_pd_tmp=export_pd[key_series==eachkey]
                        outfilename="_".join([self.settings.outFilePath_and_Prefix,eachkey,exportReadNum])+".fastq.gz"
                        export_pd_tmp=export_pd_tmp.stack()
                        export_pd_tmp=export_pd_tmp.reset_index()
                        export_pd_tmp=pd.DataFrame(export_pd_tmp[0])
                        if not os.path.isfile(outfilename):
                            export_pd_tmp.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=False)
                        else:
                            export_pd_tmp.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)
            cnt_chunk+=1
                        
