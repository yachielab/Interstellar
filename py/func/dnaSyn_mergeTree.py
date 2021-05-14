from . import settingImporter
from . import barcodeConverter
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
        cfgPath=self.opt.config
        cfg_cvrt=settingImporter.readconfig(cfgPath)["convert"]
        cfg_cvrt=settingImporter.configClean(cfg_cvrt)
        self.dest_components=cfg_cvrt["dest_components"].split(",")
        self.localTreePkl=self.opt.localTreePkl
        self.samplesheet=self.opt.samplesheet
        self.samplemerge=self.opt.samplemerge
        cvrtOptDict={}
        for i in self.dest_components:
            cvrtOption_now=cfg_cvrt[i]
            dict_now=settingImporter.convertOptionParse(cvrtOption_now)
            cvrtOptDict[i]=dict_now
            # if dict_now["is_combination"] and dict_now.get("combination_group") and dict_now.get("combination_group") not in cvrtOptDict:
            #     cvrtOptDict[dict_now["combination_group"]]=dict_now
        self.cvrtOptDict=cvrtOptDict
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)

class BARISTA_MERGETREE(object):
    def __init__(self,settings):
        self.settings=settings

    def mergeTree_v2(self):
        pkl=self.settings.localTreePkl
        num_pkl=len(pkl)

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
                
                print("Sorting done for",sample_name,". Compressing and exporting tree...",flush=True)
                with gzip.open(self.settings.outFilePath_and_Prefix+"_"+sample_name+"_mergeTree.pkl.gz",mode="wb") as p:
                        pickle.dump(subTree,p)

                #Maximum information space
                info_dic={}
                for component in subTree:
                    if self.settings.cvrtOptDict[component]["type"]=="globalConvert":
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
                if self.settings.cvrtOptDict[component]["type"]=="globalConvert":
                    info_dic[component]=len(Tree[component])
                else:
                    len_list=[len(Tree[component][i]) for i in Tree[component]]
                    info_dic[component]=max(len_list)

            with gzip.open(self.settings.outFilePath_and_Prefix+"_size_info.pkl.gz",mode="wb") as p:
                    pickle.dump(info_dic,p)