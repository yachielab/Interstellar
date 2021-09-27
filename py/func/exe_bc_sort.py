from . import settingImporter
from . import barcodeConverter
from . import settingRequirementCheck
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import time

class BARISTA_BC_SORT(object):
    def __init__(self,opt):
        self.opt=opt

    def bc_sort_settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext,dict_to_terminal=settingImporter.config_extract_value_ext(cfg)
        cfg_value_trans = settingImporter.config_extract_value_trans(cfg)
        func_dict=settingImporter.func_check_trans(cfg_value_trans,dict_to_terminal)
        self.dest_segments=cfg_value_trans["dest_segment"]       
        self.value_segment=cfg_value_ext["value_segment"]
        # cvrtOptDict={}
        # self.dest_components=cfg_cvrt["dest_components"].split(",")    
        # for i in self.dest_components:
        #     cvrtOption_now=cfg_cvrt[i]
        #     dict_now=settingImporter.convertOptionParse(cvrtOption_now)
        #     cvrtOptDict[i]=dict_now
        self.func_dict=func_dict

        #extract available src components
        src_components=[]
        for i in self.dest_segments:
            fun=func_dict[i]["func_ordered"][0]
            src_components.append(func_dict[i][fun]["source"])

        conversion_table=pd.read_csv(self.opt.table,sep="\t",header=0)
        source_terminal_pool=[dict_to_terminal[i] for i in conversion_table.columns if i in dict_to_terminal]
        print(list(conversion_table.columns)+source_terminal_pool)
        table_src=[i for i in self.value_segment if i in list(conversion_table.columns)+source_terminal_pool]
        table_dest=[i for i in conversion_table.columns if i not in self.value_segment]
        # print("Warning: On the header line of the table, the barcode name included in the source reference file is regarded as ")
        print("source barcode(s):",table_src)
        print("destination barcode(s):",table_dest)

        self.conversion_table=conversion_table
        self.table_src=table_src
        self.table_dest=table_dest

        # self.component_now=self.opt.destComponent
        self.sseq_to_svalue=self.opt.sseq_to_svalue
        self.tree=self.opt.tree
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=outdir+"/"+outname

    def bc_sort(self):
        #Only for a single global value
        #Cannot be used with samplemerge function
        func_dict=self.func_dict
        with gzip.open(self.sseq_to_svalue,mode="rb") as p:
            sseq_to_svalue=pickle.load(p)
        with gzip.open(self.tree,mode="rb") as p:
            tree=pickle.load(p)

        #Find corresponding pairs of given source components and destination components
        dval_to_sval_relationship=barcodeConverter.dval_to_sval_relationship(func_dict,self.dest_segments)
        for component in dval_to_sval_relationship:
            sval_components=dval_to_sval_relationship[component].split("+")
            if set(sval_components) == set(self.table_src):
                d_component=component
                self.table_src=sval_components
                break
        try:
            d_component
        except:
            raise ValueError("Souce variables: "+self.table_src+" do not exist in conversion process.") #kokomade

        #sseq in conversion table to svalue
        for c,src_component in enumerate(self.table_src):
            if c==0:
                src_value_series=self.conversion_table[src_component].map(lambda x: barcodeConverter.get_svalue(x,component=src_component,ref=sseq_to_svalue)).astype(str)
            else:
                src_value_series_tmp=self.conversion_table[src_component].map(lambda x: barcodeConverter.get_svalue(x,component=src_component,ref=sseq_to_svalue)).astype(str)
                src_value_series=src_value_series.str.cat(src_value_series_tmp,sep="+")
        
        src_value_series.to_csv(self.outFilePath_and_Prefix+"_series.tsv",sep="\t")

        #reinex source sequences in the conversion table into optimized value list
        src_value_prime=src_value_series.apply(barcodeConverter.get_reindex,d_component="+".join(self.table_src),tree=tree)

        src_value_prime.to_csv(self.outFilePath_and_Prefix+"_sval_prime.tsv",sep="\t")
        
        #start dealing with dest whitelist
        dest_df=pd.DataFrame(self.conversion_table[self.table_dest])

        #add optimized svalues to dest table
        dest_df["reindex"]=src_value_prime
        dest_df=dest_df.dropna(how="any")
        dest_df=dest_df.astype({"reindex":int})

        #sort the table by optimized source values
        dest_df=dest_df.sort_values("reindex")
        for i in dest_df.columns:
            if not i=="reindex":
                dest_df[i].to_csv(self.outFilePath_and_Prefix+"_"+i+"_sorted_whitelist.tsv",mode="w",sep="\t",index=False,header=False)

        dest_df.to_csv(self.outFilePath_and_Prefix+"_dest_df.tsv",sep="\t")
