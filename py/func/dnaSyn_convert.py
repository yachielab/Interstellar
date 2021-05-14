from . import settingImporter
from . import barcodeConverter
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import time

class settings_convert(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_cvrt=settingImporter.readconfig(cfgPath)["convert"]
        cfg_cvrt=settingImporter.configClean(cfg_cvrt)
        self.dest_components=cfg_cvrt["dest_components"].split(",")    
        self.tree=self.opt.tree    
        # self.correctedSrc=settingImporter.getCorrectedSrc(self.opt.srcValue,self.opt.srcQuality,self.opt.reference)
        cvrtOptDict={}
        for i in self.dest_components:
            cvrtOption_now=cfg_cvrt[i]
            dict_now=settingImporter.convertOptionParse(cvrtOption_now)
            cvrtOptDict[i]=dict_now
            
        self.cvrtOptDict=cvrtOptDict
        self.samplemerge=self.opt.samplemerge
        self.samplesheet=self.opt.samplesheet
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)
        self.path_to_sval=self.opt.srcValue
        self.path_to_qual=self.opt.srcQuality

class BARISTA_CONVERT(object):
    def __init__(self,settings):
        self.settings=settings
    def convert(self):
        convertOptions=self.settings.cvrtOptDict
        dval_to_sval_relationship = barcodeConverter.dval_to_sval_relationship(convertOptions)
        roots,edge_dict,globalComponents = barcodeConverter.parse_constraint(convertOptions)

        if self.settings.samplemerge:
            samplesheet=pd.read_csv(self.settings.samplesheet,sep="\t",header=None,dtype=str)
            sampledict=dict(zip(samplesheet[0],samplesheet[1]))
            try:
                sample_now=sampledict[self.settings.path_to_sval]
            except:
                raise ValueError("File path",self.settings.path_to_sval,"is not found in the samplesheet!")
        
        with gzip.open(self.settings.tree,mode="rb") as p:
            tree=pickle.load(p)
        s_val =pd.read_csv(self.settings.path_to_sval,sep='\t',dtype=str,chunksize=1000000)
        s_qual=pd.read_csv(self.settings.path_to_qual,sep='\t',chunksize=1000000)

        #s-value conversion
        print("Start converting s-values...")
        for n_chunk,s_val_chunk in enumerate(s_val):
            print("Processing chunk",n_chunk,"...",flush=True)
            s_val_chunk=s_val_chunk.astype(str)
            headers=s_val_chunk["Header"]
            s_val_chunk=pd.concat([headers,barcodeConverter.to_svalue_prime(s_val_chunk,dval_to_sval_relationship)],axis=1)

            compressed_svalue = pd.DataFrame({"Header":s_val_chunk["Header"]})

            if self.settings.samplemerge:
                for component in globalComponents+roots:
                    s_val_chunk[component]=s_val_chunk[component]+":"+sample_now+":"

            if globalComponents:
                for component in globalComponents:
                    compressed_svalue[component] = s_val_chunk[component].map(lambda x: barcodeConverter.compression_global(x,component=component,Tree=tree))
            
            #Local convert, Depth-First Search
            if roots:
                for root_now in roots:
                    tasklist=[root_now]
                    compressed_svalue[root_now] = s_val_chunk[root_now].map(lambda x: barcodeConverter.compression_global(x,component=root_now,Tree=tree))
                    
                    while True:
                        if not tasklist:
                            break

                        parent=tasklist[0]
                        tasklist=tasklist[1:]
                        childs=[edge_dict["child"][idx] for idx,i in enumerate(edge_dict["parent"]) if i==parent]
                        tasklist+=childs

                        if not childs:
                            continue
                            
                        for child in childs:
                            ancestor_list = barcodeConverter.getAncestor(child,edge_dict)
                            parent_series = s_val_chunk[ancestor_list[0]].str.cat(s_val_chunk[ancestor_list[1:]],sep="_")
                            parent_series.name = parent

                            s_val_tmp = parent_series.str.cat(s_val_chunk[child],sep=",")
                            
                            compressed_svalue[child] = s_val_tmp.map(lambda x: barcodeConverter.compression_local(x,component=child,Tree=tree))
            if n_chunk==0:
                compressed_svalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_value.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
            else:
                compressed_svalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_value.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)

        #Quality conversion
        print("Start converting read qualities...")
        for n_chunk,q_val_chunk in enumerate(s_qual):
            print("Processing chunk",n_chunk,"...",flush=True)
            compressed_qvalue=pd.DataFrame({"Header":q_val_chunk["Header"]})

            for component in dval_to_sval_relationship:
                if "+" in dval_to_sval_relationship[component]:
                    svalues=dval_to_sval_relationship[component].split("+")
                    average_qvalue=q_val_chunk[svalues].mean(axis=1).round().astype(int)
                else:
                    average_qvalue=q_val_chunk[dval_to_sval_relationship[component]].astype(int)
                
                compressed_qvalue[component]=average_qvalue
            
            if n_chunk==0:
                compressed_qvalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_qual.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
            else:
                compressed_qvalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_qual.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
       