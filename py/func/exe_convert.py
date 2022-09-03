from . import settingImporter
from . import barcodeConverter
from . import settingRequirementCheck
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import os

class settings_convert(object):
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
        self.tree=self.opt.tree
        self.func_dict=func_dict
        self.samplemerge=self.opt.samplemerge
        self.samplesheet=os.path.expanduser(self.opt.samplesheet)
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=outdir+"/"+outname
        self.path_to_sval=self.opt.srcValue
        self.path_to_qual=self.opt.srcQuality
        value_variables=[]
        for i in self.value_segment:
            if "SEQ2VALUE" in func_dict_ext[i]["func_ordered"]:
                value_variables.append(i)
        self.value_variables=value_variables
        self.ncore = int(self.opt.ncore)


class BARISTA_CONVERT(object):
    def __init__(self,settings):
        self.settings=settings
    def convert(self):
        if self.settings.samplemerge:
            samplesheet=pd.read_csv(self.settings.samplesheet,sep="\t",header=None,dtype=str)
            sampledict=dict(zip(samplesheet[0],samplesheet[1]))
            try:
                sample_now=sampledict[self.settings.path_to_sval]
            except:
                raise ValueError("File path",self.settings.path_to_sval,"is not found in the samplesheet!")
        
        func_dict=self.settings.func_dict
        d2s_dict = barcodeConverter.dval_to_sval_relationship(func_dict,self.settings.dest_segments)
        values_in_destarg=list(d2s_dict.values())
        value2seq_sources=[]
        for dest_segment in func_dict:
            if "RANDSEQ_ASSIGNMENT" in func_dict[dest_segment]["func_ordered"] or "WHITELIST_ASSIGNMENT" in func_dict[dest_segment]["func_ordered"]:
                value2seq_sources.append(d2s_dict[dest_segment]) 
        roots,edge_dict,globalComponents = barcodeConverter.parse_constraint(self.settings.value_segment,values_in_destarg,self.settings.child2parent_val,self.settings.value_variables,value2seq_sources)
        
        global_used_in_dest=list(set(globalComponents) & set(d2s_dict.values()))
        destname_dict={d2s_dict[k]:k for k in d2s_dict}

        with gzip.open(self.settings.tree,mode="rb") as p:
            tree=pickle.load(p)
        s_val =pd.read_csv(self.settings.path_to_sval,sep='\t',dtype=str,chunksize=1000000)
        s_qual=pd.read_csv(self.settings.path_to_qual,sep='\t',chunksize=1000000)

        #s-value conversion
        print("Start converting s-values...")
        for n_chunk,s_val_chunk in enumerate(s_val):
            print("Processing chunk",n_chunk,"...",flush=True)
            # s_val_chunk=s_val_chunk.astype(str)
            # s_val_chunk=barcodeConverter.to_svalue_prime(s_val_chunk,dval_to_sval_relationship,roots,edge_dict)
            # compressed_svalue = pd.DataFrame({"Header":s_val_chunk["Header"]})

            # if self.settings.samplemerge:
            #     for component in globalComponents+roots:
            #         s_val_chunk[component]=s_val_chunk[component]+":"+sample_now+":"

            # if global_used_in_dest:
            #     for component in global_used_in_dest:
            #         compressed_svalue[destname_dict[component]] = s_val_chunk[component].map(lambda x: barcodeConverter.compression_global(x,component=component,Tree=tree))
            
            # #Local convert, Depth-First Search
            # if roots:
            #     for root_now in roots:
            #         tasklist=[root_now]
            #         compressed_svalue[destname_dict[root_now]] = s_val_chunk[root_now].map(lambda x: barcodeConverter.compression_global(x,component=root_now,Tree=tree))    
            #         while True:
            #             if not tasklist:
            #                 break
            #             parent=tasklist[0]
            #             tasklist=tasklist[1:]
            #             childs=[edge_dict["child"][idx] for idx,i in enumerate(edge_dict["parent"]) if i==parent]
            #             tasklist+=childs

            #             if not childs:
            #                 continue
                            
            #             for child in childs:
            #                 if not child in destname_dict:
            #                     continue
            #                 ancestor_list = barcodeConverter.getAncestor(child,edge_dict)
            #                 # ancestor_list = [i.replace(",","+") for i in ancestor_list]
            #                 # parent=parent.replace(",","+")
            #                 # child=child.replace(",","+")
            #                 parent_series = s_val_chunk[ancestor_list[0]].str.cat(s_val_chunk[ancestor_list[1:]],sep="_")
            #                 parent_series.name = parent
            #                 s_val_tmp = parent_series.str.cat(s_val_chunk[child],sep=",")  
            #                 compressed_svalue[destname_dict[child]] = s_val_tmp.map(lambda x: barcodeConverter.compression_local(x,component=child,Tree=tree))
            compressed_svalue = barcodeConverter.optimize_value_table_parallel_wrapper(s_val_chunk,d2s_dict,roots,edge_dict,self.settings,globalComponents,sample_now,global_used_in_dest,destname_dict,tree,self.settings.ncore)
            if n_chunk==0:
                compressed_svalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_value.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
            else:
                compressed_svalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_value.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)

        #Quality conversion
        print("Start converting read qualities...")
        dest_value_segments=list(set(roots+globalComponents+list(edge_dict["parent"])+list(edge_dict["child"])))
        components_used=list(set(dest_value_segments) & set(d2s_dict.values()))
        for n_chunk,q_val_chunk in enumerate(s_qual):
            print("Processing chunk",n_chunk,"...",flush=True)
            # compressed_qvalue=pd.DataFrame({"Header":q_val_chunk["Header"]})
            # for component in components_used:
            #     if "+" in component:
            #         svalues=component.split("+")
            #         average_qvalue=q_val_chunk[svalues].mean(axis=1).round().astype(int)
            #     else:
            #         average_qvalue=q_val_chunk[component].astype(int)
                
            #     compressed_qvalue[destname_dict[component]]=average_qvalue
            compressed_qvalue = barcodeConverter.quality_conversion_parallel_wrapper(q_val_chunk,components_used,destname_dict,self.settings.ncore)

            if n_chunk==0:
                compressed_qvalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_qual.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
            else:
                compressed_qvalue.to_csv(self.settings.outFilePath_and_Prefix+"_converted_qual.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
       