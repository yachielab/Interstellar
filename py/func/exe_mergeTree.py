from . import settingImporter
from . import barcodeConverter
from . import settingRequirementCheck
import regex
import collections
import pandas as pd
import numpy as np
import gzip
import pickle
import time

class settings_mergetree(object):
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
        self.localTreePkl=self.opt.localTreePkl
        self.samplesheet=self.opt.samplesheet
        self.samplemerge=self.opt.samplemerge
        self.func_dict=func_dict
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=outdir+"/"+outname
        value_variables=[]
        for i in self.value_segment:
            if "SEQ2VALUE" in func_dict_ext[i]["func_ordered"]:
                value_variables.append(i)
        self.value_variables=value_variables

class BARISTA_MERGETREE(object):
    def __init__(self,settings):
        self.settings=settings

    def mergeTree_v2(self):
        pkl=self.settings.localTreePkl
        func_dict=self.settings.func_dict
        dval_to_sval_relationship = barcodeConverter.dval_to_sval_relationship(func_dict,self.settings.dest_segments)
        values_in_destarg=list(dval_to_sval_relationship.values())
        value2seq_sources=[]
        for dest_segment in func_dict:
            if "RANDSEQ_ASSIGNMENT" in func_dict[dest_segment]["func_ordered"] or "WHITELIST_ASSIGNMENT" in func_dict[dest_segment]["func_ordered"]:
                value2seq_sources.append(dval_to_sval_relationship[dest_segment]) 

        roots,edge_dict,globalComponents = barcodeConverter.parse_constraint(self.settings.value_segment,values_in_destarg,self.settings.child2parent_val,self.settings.value_variables,value2seq_sources)
        

        #Pile up counttree
        mergedCountTree={}
        for cnt,p in enumerate(pkl):
            print("Piling up count trees...",cnt+1,flush=True)
            #load counttree
            with gzip.open(p,mode="rb") as rp:
                count_tree=pickle.load(rp)

            #go through the component in the counttree
            for component in count_tree:
                if not component in mergedCountTree: #1st counttree
                    mergedCountTree[component]=count_tree[component]
                else: #2nd~
                    if type(count_tree[component])==collections.Counter: #For global value
                        mergedCountTree[component].update(count_tree[component])
                    else: #For local value
                        for parent in count_tree[component]:
                            child_freq_now=count_tree[component][parent]
                            if parent in mergedCountTree[component]:
                                mergedCountTree[component][parent].update(child_freq_now)
                            else:
                                mergedCountTree[component][parent]=child_freq_now
        
        # barcodeConverter.printDictTest(mergedCountTree)
        #Sort piled up tree by count, and add indices
        print("Start sorting...",flush=True)
        Tree={}
        iterate_keys=list(mergedCountTree.keys())
        for component in iterate_keys:
            
            Tree[component]={}
            if type(mergedCountTree[component])==collections.Counter: #For global value
                Tree[component]=barcodeConverter.sortCounter(mergedCountTree[component])
            else: #For local value
                for parent in mergedCountTree[component]:
                    child_freq_now=mergedCountTree[component][parent]
                    Tree[component][parent]=barcodeConverter.sortCounter(child_freq_now)
            
            #to minimize memory consumption
            del mergedCountTree[component]
        
        # for i in Tree:
        #     for n,k in enumerate(Tree[i]):
        #         if n<10:
        #             print(i,k,Tree[i][k])
        
        
        #Generate a merged tree 
        if self.settings.samplemerge:
            """
            If multiple samples were merged, 
            - The same values from different samples are converted into different d-values 
            - Since the tree containing multiple samples is too large, the tree is split by each sample
            """
            samplesheet=pd.read_csv(self.settings.samplesheet,sep="\t",header=None,dtype=str)
            sample_names=list(set(samplesheet[1]))
            iterate_keys=list(Tree.keys())
            for sample_name in sample_names:
                #Tree is split by each sample
                subTree={}
                for component in iterate_keys:
                    sample_matching_pattern=":"+sample_name+":"
                    
                    if type(Tree[component])==collections.Counter: #For global value
                        counttree_tmp={k:Tree[component][k] for k in Tree[component] if sample_matching_pattern in k}
                        subTree[component]=counttree_tmp
                    else: #For local value
                        counttree_tmp={k:Tree[component][k] for k in Tree[component] if sample_matching_pattern in k}
                        subTree[component]={}
                        for parent in counttree_tmp:
                            child_freq_now=counttree_tmp[parent]
                            subTree[component][parent]=child_freq_now
                
                print("Sorting done for "+sample_name+". Compressing and exporting tree...",flush=True)
                with gzip.open(self.settings.outFilePath_and_Prefix+"_"+sample_name+"_mergeTree.pkl.gz",mode="wb") as p:
                        pickle.dump(subTree,p)

                #Maximum information space
                info_dic={}
                for component in subTree:
                    if component in roots+globalComponents: #For global value
                        info_dic[component]=len(Tree[component]) #maximum size of the merged tree
                    else:
                        len_list=[len(subTree[component][i]) for i in subTree[component]]
                        info_dic[component]=max(len_list)

                with gzip.open(self.settings.outFilePath_and_Prefix+"_"+sample_name+"_size_info.pkl.gz",mode="wb") as p:
                        pickle.dump(info_dic,p)

        else:
            print("Sorting done. Compressing and exporting tree...",flush=True)
            with gzip.open(self.settings.outFilePath_and_Prefix+"_mergeTree.pkl.gz",mode="wb") as p:
                    pickle.dump(Tree,p)

            #Maximum information space
            info_dic={}
            for component in Tree:
                if component in roots+globalComponents: #For global value
                    info_dic[component]=len(Tree[component])
                else:
                    len_list=[len(Tree[component][i]) for i in Tree[component]]
                    info_dic[component]=max(len_list)

            with gzip.open(self.settings.outFilePath_and_Prefix+"_size_info.pkl.gz",mode="wb") as p:
                    pickle.dump(info_dic,p)