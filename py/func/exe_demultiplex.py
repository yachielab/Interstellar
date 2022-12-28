from func.annotate_header import UnknownError
from . import settingImporter
from . import barcodeConverter
from . import settingRequirementCheck
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import csv
import os
import datetime

class settings_demultiplex(object):
    def __init__(self,opt):
        self.opt=opt 

    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_demulti=settingImporter.config_extract_value_demulti(cfg)
        cfg_value_ext,dict_to_terminal=settingImporter.config_extract_value_ext(cfg)

        self.key=cfg_demulti["KEY"].split(",")
        self.run_demulti = True if self.key else False
        self.target=cfg_demulti["TARGET"].split(",")
        if cfg_demulti["TARGET"]=="":
            self.target=cfg_value_ext["value_segment"]
        self.exportReadStructure={}
        self.annotate_headers={}
        for i in cfg_demulti:
            if i in ["READ1_STRUCTURE","READ2_STRUCTURE","INDEX1_STRUCTURE","INDEX2_STRUCTURE"] and cfg_demulti.get(i):
                self.exportReadStructure[i]=[]
                read_structure_list = cfg_demulti[i].split("+")
                for x in read_structure_list:
                    if x in dict_to_terminal:
                        self.exportReadStructure[i].append(dict_to_terminal[x])
                    elif regex.search(r'^\"[^\"]+\"',x):
                        seq=regex.sub('\"',"",x)
                        self.exportReadStructure[i].append(seq)

            elif i in ["READ1_HEADER_ADDITION","READ2_HEADER_ADDITION","INDEX1_HEADER_ADDITION","INDEX2_HEADER_ADDITION"] and cfg_demulti.get(i):
                self.annotate_headers[i.split("_")[0]+"_STRUCTURE"]=cfg_demulti[i].split(",")

        # for i in cfg_demulti:
        #     if i in ["READ1_STRUCTURE","READ2_STRUCTURE","INDEX1_STRUCTURE","INDEX2_STRUCTURE"] and cfg_demulti.get(i):
        #         self.exportReadStructure[i]=cfg_demulti[i].split("+")
        
        if self.opt.mode_local:
            self.path_to_seq = settingImporter.parseInputFileList(self.opt.correctedSeq)
            self.path_to_avg_qval = settingImporter.parseInputFileList(self.opt.correctedQual)
            self.path_to_rawQual = settingImporter.parseInputFileList(self.opt.rawQual)
            outname = settingImporter.parseInputFileList(self.opt.outname)
        else:
            self.path_to_seq = [self.opt.correctedSeq]
            self.path_to_avg_qval = [self.opt.correctedQual]
            self.path_to_rawQual = [self.opt.rawQual]
            outname = [self.opt.outname]

        self.export_tsv=self.opt.export_tsv
        outdir=self.opt.outdir

        today_now=str(datetime.datetime.today())
        today_now=regex.sub(r"\W","",today_now)
        new_outdir_list=["_".join([outdir+"/demulti",i,today_now]) for i in outname]
        for new_outdir in new_outdir_list:
            flag=os.path.isdir(new_outdir)
            while flag:
                today_now=str(datetime.datetime.today())
                today_now=regex.sub(r"\W","",today_now)
                new_outdir="_".join([outdir+"/demulti",outname,today_now])
                flag=os.path.isdir(new_outdir)
            os.mkdir(new_outdir)
        self.outFilePath_and_Prefix_list=[new_outdir+"/"+i for new_outdir,i in zip(new_outdir_list,outname)]
        self.today_now=today_now

        if not self.run_demulti and self.export_tsv:
            raise UnknownError("KEY is required if FORMAT=TSV.")
        
        self.ncore = int(self.opt.ncore)


class BARISTA_DEMULTIPLEX(object):
    def __init__(self,settings):
        self.settings=settings
    def demultiplex(self):
        cnt = 0
        for sseq_path,avg_qual_path,raw_qual_path,prefix in zip(self.settings.path_to_seq,self.settings.path_to_avg_qval,self.settings.path_to_rawQual,self.settings.outFilePath_and_Prefix_list):
            if regex.search(r"\.pkl$",sseq_path):
                s_seq_chunk = pd.read_pickle(sseq_path)
            else:
                s_seq_chunk = pd.read_csv(sseq_path,sep='\t',dtype=str)
            # s_avg_qual=pd.read_csv(avg_qual_path,sep='\t',dtype=str,chunksize=500000)
            # s_raw_qual=pd.read_csv(raw_qual_path,sep="\t",dtype=str,chunksize=500000,quoting=csv.QUOTE_NONE)
            s_avg_qual_chunk = pd.read_pickle(avg_qual_path)
            s_raw_qual_chunk = pd.read_pickle(raw_qual_path)

            demulti_key=set()
            key_iden_list=[]
            print("Start demultiplexing chunk",cnt,flush=True)
            s_avg_qual_chunk=s_avg_qual_chunk.drop("Header",axis=1)
            s_raw_qual_chunk=s_raw_qual_chunk.drop("Header",axis=1)

            colnames=list(s_seq_chunk.columns)
            raw_component_names=[i.split(":")[0] for i in colnames if not i=="Header"]
            corrected_component_names=[i.split(":")[1] for i in colnames if not i=="Header"]
            s_seq_chunk.columns=["Header"]+corrected_component_names       

            if self.settings.export_tsv:
                key_series=s_seq_chunk[self.settings.key].apply("_".join,axis=1)
                demulti_key|=set(key_series)
                s_seq_chunk=s_seq_chunk[["Header"]+self.settings.target]
                key_iden_list_now = barcodeConverter.demultiplex_tsv_parallel_wrapper(s_seq_chunk,key_series,prefix,self.settings.ncore)
                # for eachkey in demulti_key:
                #     if not "-" in eachkey:
                #         export_pd_tmp=s_seq_chunk[key_series==eachkey]
                #         outfilename="_".join([self.settings.outFilePath_and_Prefix,eachkey])+".tsv.gz"
                #         if not os.path.isfile(outfilename):
                #             export_pd_tmp.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=True)
                #         else:
                #             export_pd_tmp.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)
                #         key_iden_list.append("_"+eachkey+".tsv.gz")
            else:
                # Collect segment names to be used
                use_segment_list=set(["Header"])
                for exportReadNum in self.settings.exportReadStructure:
                    if exportReadNum in self.settings.annotate_headers:
                        use_segment_list |= set(self.settings.annotate_headers[exportReadNum])
                    if exportReadNum in self.settings.exportReadStructure:
                        use_segment_list |= set(self.settings.exportReadStructure[exportReadNum])
                use_segment_list |= set(self.settings.key)
                use_segment_list=list(use_segment_list)

                # Constant sequence
                for exportReadNum in self.settings.exportReadStructure:
                    structure_now=self.settings.exportReadStructure[exportReadNum]
                    for component in structure_now:
                        if not component in s_seq_chunk.columns:
                            s_seq_chunk[component]=component
                
                # Cut off unused columns
                s_seq_chunk=s_seq_chunk[use_segment_list]
                
                # Drop missing rows containing segments
                s_seq_chunk=s_seq_chunk[s_seq_chunk!="-"].dropna()
                s_avg_qual_chunk=s_avg_qual_chunk.loc[s_seq_chunk.index]
                s_raw_qual_chunk=s_raw_qual_chunk.loc[s_seq_chunk.index]
                s_avg_qual_chunk=s_avg_qual_chunk.astype("int8")
                if s_seq_chunk.shape[0]==0:
                    raise UnknownError("All rows were dropped out.")
                
                # Generating demultiplex keys
                key_series=s_seq_chunk[self.settings.key].apply("_".join,axis=1)
                demulti_key|=set(key_series)

                export_pd=pd.DataFrame()
                export_pd["Header"]=s_seq_chunk["Header"]
                    
                for exportReadNum in self.settings.exportReadStructure:
                    export_pd["3rd"]=["+"]*export_pd.shape[0]
                    export_pd["qual"]=[""]*export_pd.shape[0]

                    # Annotating header
                    if exportReadNum in self.settings.annotate_headers:
                        annotate_header_now=self.settings.annotate_headers[exportReadNum]
                    else:
                        annotate_header_now=""

                    if not annotate_header_now=="":
                        #tagging
                        export_pd["Header"]=s_seq_chunk["Header"].str.cat(s_seq_chunk[annotate_header_now],sep="_")
                    else:
                        #Deal with non-tag reads
                        export_pd["Header"]=s_seq_chunk["Header"]
                
                    structure_now=self.settings.exportReadStructure[exportReadNum]
                    
                    # Sequence segment concatenation
                    if len(structure_now)>1:
                        export_pd["seq"]=s_seq_chunk[structure_now[0]].str.cat(s_seq_chunk[structure_now[1:]],sep="")
                    else:
                        export_pd["seq"]=s_seq_chunk[structure_now[0]]

                    # Quality segment concatenation
                    for component in structure_now:
                        if component in s_avg_qual_chunk:
                            df_seq_qual_tmp=s_seq_chunk[component].str.cat(s_avg_qual_chunk[component].astype(str),sep="_")
                            df_seq_qual_tmp=df_seq_qual_tmp.map(barcodeConverter.getConvQual_ver2)
                            export_pd["qual"]=export_pd["qual"].str.cat(df_seq_qual_tmp,sep="")
                        elif component in corrected_component_names:
                            component_raw=raw_component_names[corrected_component_names.index(component)]
                            export_pd["qual"]=export_pd["qual"].str.cat(s_raw_qual_chunk[component_raw],sep="")
                        else:
                            baseQuality=40 #hard coded
                            const_quality=[chr(baseQuality+33)*len(component)]*s_raw_qual_chunk.shape[0]
                            export_pd["qual"]=export_pd["qual"].str.cat(const_quality,sep="")

                    export_pd=export_pd[["Header","seq","3rd","qual"]]

                    if exportReadNum == "READ1_STRUCTURE": readIden="R1"
                    if exportReadNum == "READ2_STRUCTURE": readIden="R2"
                    if exportReadNum == "INDEX1_STRUCTURE": readIden="I1"
                    if exportReadNum == "INDEX2_STRUCTURE": readIden="I2"

                    #FASTQ exportation
                    if self.settings.run_demulti:
                        #demultiplex
                        key_iden_list_now = barcodeConverter.demultiplex_fastq_parallel_wrapper(export_pd,key_series,prefix,readIden,self.settings.ncore)
                        # for eachkey in demulti_key:
                        #     export_pd_tmp=export_pd[key_series==eachkey]
                        #     outfilename="_".join([self.settings.outFilePath_and_Prefix,eachkey,readIden])+".fastq.gz"
                        #     export_pd_tmp=export_pd_tmp.stack()
                        #     export_pd_tmp=export_pd_tmp.reset_index()
                        #     export_pd_tmp=pd.DataFrame(export_pd_tmp[0])
                        #     if not os.path.isfile(outfilename):
                        #         export_pd_tmp.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=False)
                        #     else:
                        #         export_pd_tmp.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)
                        #     key_iden_list.append("_"+eachkey+"_"+readIden+".fastq.gz")
                    else:
                        #just tagging
                        outfilename="_".join([prefix,"annotated",readIden])+".fastq.gz"
                        export_pd_tmp=export_pd.stack()
                        export_pd_tmp=export_pd_tmp.reset_index()
                        export_pd_tmp=pd.DataFrame(export_pd_tmp[0])
                        if not os.path.isfile(outfilename):
                            export_pd_tmp.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=False)
                        else:
                            export_pd_tmp.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)
                        key_iden_list_now = ["_annotated_"+readIden+".fastq.gz"]

            cnt += 1
            key_iden_list += key_iden_list_now

            key_iden_list=list(set(key_iden_list))
            with open(prefix+"_demulti_key_list.txt",mode="wt") as w:
                for i in key_iden_list:
                    w.write(i+"\n")
