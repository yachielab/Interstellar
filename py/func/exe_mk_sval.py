from . import settingImporter
from . import barcodeCorrecter
from . import settingRequirementCheck
import csv
import pandas as pd
import gzip
import pickle
import sys

class settings_make_s(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext,dict_to_terminal = settingImporter.config_extract_value_ext(cfg)
        func_dict=settingImporter.func_check(cfg_value_ext)

        self.corrected_components=cfg_value_ext["value_segment"]
        self.correctionDictPkl=self.opt.correctedPickle
        if self.opt.mode_local:
            self.rawFastqPath = dict(seq=settingImporter.parseInputFileList(self.opt.rawSeq), qual=settingImporter.parseInputFileList(self.opt.rawQual))
            outname = settingImporter.parseInputFileList(self.opt.outname)
        else:
            self.rawFastqPath = dict(seq=[self.opt.rawSeq], qual=[self.opt.rawQual])
            outname = [self.opt.outname]
        self.correctOptDict=func_dict
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix_list=[outdir+"/"+i for i in outname]
        self.ncore = int(self.opt.ncore)
        

class BARISTA_MAKE_S(object):
    def __init__(self,settings):
        self.settings=settings

    def make_s_seq(self):
        with gzip.open(self.settings.correctionDictPkl) as p:
            correctionDictionaries=pickle.load(p)
        
        cnt = 0
        for seq,outprefix in zip(self.settings.rawFastqPath["seq"], self.settings.outFilePath_and_Prefix_list):
            # parsedSeq_raw_chunk=pd.read_csv(seq, sep='\t', chunksize=1000000)
            parsedSeq_raw_df = pd.read_pickle(seq)
            components_raw=[]
            for k in self.settings.corrected_components:
                func_tmp=self.settings.correctOptDict[k]["func_ordered"][0]
                components_raw.append(self.settings.correctOptDict[k][func_tmp]["source"])
            self.components_raw=components_raw

            print("Making a sequence table for chunk",cnt,"...",flush=True)
            cnt += 1
            parsedSeq_raw_df = barcodeCorrecter.seqCleanUp_parallel_wrapper(parsedSeq_raw_df,components_raw,correctionDictionaries,self.settings,self.settings.ncore)

                # cols_ordered=["Header"]+[i for i in sorted(parsedSeq_raw_df.columns) if not i=="Header"]
                # cols_ordered_final=[]
                # for ncol,component_raw_now in enumerate(cols_ordered):
                #     if ncol>0:
                #         if not component_raw_now in components_raw:
                #             parsedSeq_raw_df=parsedSeq_raw_df.drop(component_raw_now,axis=1)
                #             continue
                #         col_index=components_raw.index(component_raw_now)
                #         component_corrected_now=self.settings.corrected_components[col_index]
                #         opt_now=self.settings.correctOptDict[component_corrected_now]
                        
                #         if (component_corrected_now in correctionDictionaries) and ("I2M_CORRECTION" in opt_now["func_ordered"] or "M2A_CORRECTION" in opt_now["func_ordered"] or "BARTENDER_CORRECTION" in opt_now["func_ordered"] or "CUSTOM_CORRECTION" in opt_now["func_ordered"]):
                #             parsedSeq_raw_df[component_raw_now]=parsedSeq_raw_df[component_raw_now].map(lambda x: barcodeCorrecter.seq_correct_and_write(x,reference=correctionDictionaries[component_corrected_now]["correctionDict"]))
                #         parsedSeq_raw_df.rename(columns={component_raw_now:component_raw_now+":"+component_corrected_now},inplace=True)
                #         cols_ordered_final.append(component_raw_now+":"+component_corrected_now)
                #     else:
                #         cols_ordered_final.append("Header")
                # parsedSeq_raw_df = parsedSeq_raw_df[cols_ordered_final]
                
            parsedSeq_raw_df.to_pickle(outprefix+"_correct_result.pkl")
        
    
    def make_s_value(self):
        with gzip.open(self.settings.correctionDictPkl) as p:
            correctionDictionaries=pickle.load(p)
        
        ref_dic={}
        for component_corrected_now in correctionDictionaries:
            ref_tup=tuple(correctionDictionaries[component_corrected_now]["reference"])
            ref_dic[component_corrected_now]={k:v for k,v in zip(ref_tup,range(len(ref_tup)))}

        cnt = 0
        for qual,outprefix in zip(self.settings.rawFastqPath["qual"], self.settings.outFilePath_and_Prefix_list):
            parsedSeq_raw_chunk=pd.read_pickle(outprefix+"_correct_result.pkl")
            # parsedQual_raw_chunk=pd.read_csv(qual, sep='\t',chunksize=1000000,quoting=csv.QUOTE_NONE)
            parsedQual_raw_chunk=pd.read_pickle(qual)
            components_raw=self.components_raw

            # with gzip.open(outprefix+"_sseq_to_svalue.pkl.gz",mode="wb") as p:
            #         pickle.dump(ref_dic,p)
                                    
            for cat in ["seq","qual"]:  
                if cat=="seq":
                    parsedSeq_df=parsedSeq_raw_chunk
                else:
                    parsedSeq_df=parsedQual_raw_chunk

                print("Category:",cat,"| Making a value table for chunk",cnt)
                res = barcodeCorrecter.gen_value_table_parallel_wrapper(parsedSeq_df,cat,components_raw,self.settings,correctionDictionaries,ref_dic,self.settings.ncore)
                # res=pd.DataFrame()
                # for ncol,component in enumerate(parsedSeq_df.columns):                    
                #     if ncol==0:
                #         res["Header"]=parsedSeq_df[component]
                #     if ncol>0: 
                #         if cat=="seq":
                #             component_raw_now=component.split(":")[0]
                #             component_corrected_now=component.split(":")[1]
                #         else:
                #             component_raw_now=component
                #             if not component_raw_now in components_raw:
                #                 continue
                #             col_index=components_raw.index(component_raw_now)
                #             component_corrected_now=self.settings.corrected_components[col_index]

                #         if component_corrected_now in correctionDictionaries:
                #             # print(component_corrected_now,flush=True)
                #             if cat=="seq":
                #                 res[component_corrected_now]=parsedSeq_df[component].map(lambda x: barcodeCorrecter.seq_to_val_ver2(x,dic=ref_dic[component_corrected_now]))
                #             else:
                #                 res[component_corrected_now]=parsedSeq_df[component_raw_now].map(barcodeCorrecter.qualityProcessing)
                if cat=="seq":
                    # res.to_csv(outprefix+"_correct_srcValue.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
                    res.to_pickle(outprefix+"_correct_srcValue.pkl")
                else:
                    # res.to_csv(outprefix+"_correct_srcQual.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
                    res.to_pickle(outprefix+"_correct_srcQual.pkl")
                # elif cnt > 0:
                #     if cat=="seq":
                #         res.to_csv(outprefix+"_correct_srcValue.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
                #     else:
                #         res.to_csv(outprefix+"_correct_srcQual.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
            cnt += 1