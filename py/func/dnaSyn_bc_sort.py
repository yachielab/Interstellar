from . import settingImporter
from . import barcodeConverter
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
        cfgPath=self.opt.config
        cfg_cvrt=settingImporter.readconfig(cfgPath)["convert"]
        cfg_cvrt=settingImporter.configClean(cfg_cvrt)
        cvrtOptDict={}
        self.dest_components=cfg_cvrt["dest_components"].split(",")    
        for i in self.dest_components:
            cvrtOption_now=cfg_cvrt[i]
            dict_now=settingImporter.convertOptionParse(cvrtOption_now)
            cvrtOptDict[i]=dict_now
        self.cvrtOptDict=cvrtOptDict

        #extract available src components
        src_components=[]
        for i in cvrtOptDict:
            for j in cvrtOptDict[i]:
                if "src_corrected_component" in j:
                    src_components.append(cvrtOptDict[i][j])

        conversion_table=pd.read_csv(self.opt.table,sep="\t",header=0)
        table_src=[i for i in conversion_table.columns if i in src_components]
        table_dest=[i for i in conversion_table.columns if i not in src_components]
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
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)

    def bc_sort(self):
        #Only for a single global value
        #Cannot be used with samplemerge function
        convertOptions=self.cvrtOptDict
        with gzip.open(self.sseq_to_svalue,mode="rb") as p:
            sseq_to_svalue=pickle.load(p)
        with gzip.open(self.tree,mode="rb") as p:
            tree=pickle.load(p)

        #Find corresponding pairs of given source components and destination components
        dval_to_sval_relationship=barcodeConverter.dval_to_sval_relationship(convertOptions)
        for component in dval_to_sval_relationship:
            sval_components=dval_to_sval_relationship[component].split("+")
            if set(sval_components) == set(self.table_src):
                d_component=component
                self.table_src=sval_components
                break
        try:
            d_component
        except:
            raise ValueError("Souce barcodes: "+self.table_src+" do not exist in conversion process.")

        #sseq in conversion table to svalue
        for c,src_component in enumerate(self.table_src):
            if c==0:
                src_value_series=self.conversion_table[src_component].map(lambda x: barcodeConverter.get_svalue(x,component=src_component,ref=sseq_to_svalue)).astype(str)
            else:
                src_value_series_tmp=self.conversion_table[src_component].map(lambda x: barcodeConverter.get_svalue(x,component=src_component,ref=sseq_to_svalue)).astype(str)
                src_value_series=src_value_series.str.cat(src_value_series_tmp,sep="+")
        src_value_prime=src_value_series.apply(barcodeConverter.get_reindex,d_component=d_component,tree=tree)
        
        dest_df=pd.DataFrame(self.conversion_table[self.table_dest])
        dest_df["reindex"]=src_value_prime
        dest_df=dest_df.dropna(how="any")
        dest_df=dest_df.astype({"reindex":int})
        dest_df=dest_df.sort_values("reindex")
        for i in dest_df.columns:
            if not i=="reindex":
                dest_df[i].to_csv(self.outFilePath_and_Prefix+"_"+i+"_sorted_whitelist.tsv",mode="w",sep="\t",index=False,header=False)
