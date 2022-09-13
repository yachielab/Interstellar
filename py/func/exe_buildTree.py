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

class settings_buildTree(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext,dict_to_terminal=settingImporter.config_extract_value_ext(cfg)
        func_dict_ext=settingImporter.func_check(cfg_value_ext)
        cfg_value_trans = settingImporter.config_extract_value_trans(cfg)
        func_dict=settingImporter.func_check_trans(cfg_value_trans,dict_to_terminal)
        self.child2parent_val=settingImporter.getAllocation(func_dict_ext,cfg_value_trans,cfg_value_ext)
        self.dest_segments=cfg_value_trans["dest_segment"]       
        self.value_segment=cfg_value_ext["value_segment"]
        self.func_dict=func_dict
        self.samplesheet=os.path.expanduser(self.opt.samplesheet)
        self.samplemerge=self.opt.samplemerge
        if self.opt.mode_local:
            self.path_to_sval = settingImporter.parseInputFileList(self.opt.srcValue)
            outname = settingImporter.parseInputFileList(self.opt.outname)
        else:
            self.path_to_sval = [self.opt.srcValue]
            outname = [self.opt.outname]
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix_list=[outdir+"/"+i for i in outname]

        value_variables=[]
        for i in self.value_segment:
            if "SEQ2VALUE" in func_dict_ext[i]["func_ordered"]:
                value_variables.append(i)
        self.value_variables=value_variables
        self.ncore = int(self.opt.ncore)
        

class BARISTA_BUILDTREE(object):
    def __init__(self,settings):
        self.settings=settings
    def buildTree(self):
        for val_file_path,prefix in zip(self.settings.path_to_sval,self.settings.outFilePath_and_Prefix_list):
            if self.settings.samplemerge:
                samplesheet=pd.read_csv(self.settings.samplesheet,sep="\t",header=None,dtype=str)
                sampledict=dict(zip(samplesheet[0],samplesheet[1]))
                try:
                    sample_now=sampledict[val_file_path]
                except:
                    raise ValueError("File path",val_file_path,"is not found in the samplesheet!")
                
            func_dict=self.settings.func_dict
            d2s_dict = barcodeConverter.dval_to_sval_relationship(func_dict,self.settings.dest_segments)
            values_in_destarg=list(d2s_dict.values())
            
            # Get value2seq segments
            value2seq_sources=[]
            for dest_segment in func_dict:
                if "RANDSEQ_ASSIGNMENT" in func_dict[dest_segment]["func_ordered"] or "WHITELIST_ASSIGNMENT" in func_dict[dest_segment]["func_ordered"]:
                    value2seq_sources.append(d2s_dict[dest_segment]) 
            
            roots,edge_dict,globalComponents = barcodeConverter.parse_constraint(self.settings.value_segment,values_in_destarg,self.settings.child2parent_val,self.settings.value_variables,value2seq_sources)

            count_tree={}

            # s_val=pd.read_csv(val_file_path,sep='\t',dtype=str,chunksize=1000000)
            s_val=[pd.read_pickle(val_file_path)]
            for n_chunk,s_val_chunk in enumerate(s_val):
                print("Sample:",sample_now,"| Building a count tree for chunk",n_chunk,"start...",flush=True)
                count_tree = barcodeConverter.build_count_tree_parallel_wrapper(s_val_chunk,roots,edge_dict,d2s_dict,self.settings,globalComponents,sample_now,count_tree,self.settings.ncore)
                # s_val_chunk=s_val_chunk.astype(str)
                # s_val_chunk=s_val_chunk.replace("-1",np.nan).dropna()
                # s_val_chunk=barcodeConverter.to_svalue_prime(s_val_chunk,dval_to_sval_relationship,roots,edge_dict)

                # if self.settings.samplemerge:
                #     for component in globalComponents+roots:
                #         s_val_chunk[component]=s_val_chunk[component]+":"+sample_now+":"

                # if globalComponents:
                #     count_tree=barcodeConverter.buildTree_global(s_val_chunk,globalComponents,Tree=count_tree)
                # count_tree = barcodeConverter.buildTree(s_val_chunk,roots,edge_dict,Tree=count_tree)

            with gzip.open(prefix+"_Tree.pkl.gz",mode="wb") as p:
                    pickle.dump(count_tree,p)
        