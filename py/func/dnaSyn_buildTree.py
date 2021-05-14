from . import settingImporter
from . import barcodeConverter
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import time

class settings_buildTree(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_cvrt=settingImporter.readconfig(cfgPath)["convert"]
        cfg_cvrt=settingImporter.configClean(cfg_cvrt)
        self.dest_components=cfg_cvrt["dest_components"].split(",")        
        # self.correctedSrc=settingImporter.getCorrectedSrc(self.opt.srcValue,self.opt.srcQuality,self.opt.reference)
        cvrtOptDict={}
        for i in self.dest_components:
            cvrtOption_now=cfg_cvrt[i]
            dict_now=settingImporter.convertOptionParse(cvrtOption_now)
            cvrtOptDict[i]=dict_now
            # if dict_now["is_combination"] and dict_now.get("combination_group") and dict_now.get("combination_group") not in cvrtOptDict:
            #     cvrtOptDict[dict_now["combination_group"]]=dict_now
        self.cvrtOptDict=cvrtOptDict
        self.samplesheet=self.opt.samplesheet
        self.samplemerge=self.opt.samplemerge
        self.path_to_sval=self.opt.srcValue
        # self.path_to_qual=self.opt.srcQuality
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)
        

class BARISTA_BUILDTREE(object):
    def __init__(self,settings):
        self.settings=settings
    def buildTree(self):
        if self.settings.samplemerge:
            samplesheet=pd.read_csv(self.settings.samplesheet,sep="\t",header=None,dtype=str)
            sampledict=dict(zip(samplesheet[0],samplesheet[1]))
            try:
                sample_now=sampledict[self.settings.path_to_sval]
            except:
                raise ValueError("File path",self.settings.path_to_sval,"is not found in the samplesheet!")
            
        convertOptions=self.settings.cvrtOptDict
        dval_to_sval_relationship = barcodeConverter.dval_to_sval_relationship(convertOptions)
        roots,edge_dict,globalComponents = barcodeConverter.parse_constraint(convertOptions)
        
        count_tree={}
        s_val=pd.read_csv(self.settings.path_to_sval,sep='\t',dtype=str,chunksize=1000000)
        for n_chunk,s_val_chunk in enumerate(s_val):
            print("Building a count tree for chunk",n_chunk,"start...",flush=True)
            s_val_chunk=s_val_chunk.astype(str)
            s_val_chunk=s_val_chunk.replace("-1",np.nan).dropna()
            s_val_chunk=barcodeConverter.to_svalue_prime(s_val_chunk,dval_to_sval_relationship)

            if self.settings.samplemerge:
                for component in globalComponents+roots:
                    s_val_chunk[component]=s_val_chunk[component]+":"+sample_now+":"

            if globalComponents:
                count_tree=barcodeConverter.buildTree_global(s_val_chunk,globalComponents,Tree=count_tree)
            count_tree = barcodeConverter.buildTree(s_val_chunk,roots,edge_dict,Tree=count_tree)

        with gzip.open(self.settings.outFilePath_and_Prefix+"_Tree.pkl.gz",mode="wb") as p:
                pickle.dump(count_tree,p)
        