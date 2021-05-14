from . import settingImporter
from . import barcodeCorrecter
import regex
import pandas as pd
import gzip
import pickle
import time

class settings_make_s(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_correct=settingImporter.readconfig(cfgPath)["correct"]
        cfg_correct=settingImporter.configClean(cfg_correct)
        self.corrected_components=cfg_correct["src_corrected_components"].split(",")
        self.correctionDictPkl=self.opt.correctedPickle
        self.rawFastqPath=dict(seq=self.opt.rawSeq,qual=self.opt.rawQual)
        
        correctOptDict={}
        corrected_components_sorted=[]
        for i in self.corrected_components:
            correctOption_now=cfg_correct[i]
            if not regex.search("^combination",correctOption_now):
                correctOptDict[i]=settingImporter.correctOptionParse(correctOption_now)
                corrected_components_sorted.append(i)
        for i in self.corrected_components:
            correctOption_now=cfg_correct[i]
            if regex.search("^combination",correctOption_now):
                correctOptDict[i]=settingImporter.correctOptionParse(correctOption_now)
                corrected_components_sorted.append(i)
        self.correctOptDict=correctOptDict
        self.corrected_components_sorted=corrected_components_sorted
        self.resultonly=self.opt.resultonly
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)
        # self.mk_s_value_components=cfg_correct["make_s_value"].split(",")

class BARISTA_MAKE_S(object):
    def __init__(self,settings):
        self.settings=settings

    def make_s_seq(self):
        with gzip.open(self.settings.correctionDictPkl) as p:
            correctionDictionaries=pickle.load(p)
        parsedSeq_raw_chunk=pd.read_csv(self.settings.rawFastqPath["seq"], sep='\t',chunksize=1000000)
        options=[self.settings.correctOptDict[i] for i in self.settings.corrected_components_sorted]
        components_raw=[dic["src_raw_components"] for dic in options]
        for cnt,parsedSeq_raw_df in enumerate(parsedSeq_raw_chunk):
            print("Making a sequence table for chunk",cnt,"...",flush=True)
            for ncol,component_raw_now in enumerate(parsedSeq_raw_df.columns):
                if ncol>0:
                    if not component_raw_now in components_raw:
                        parsedSeq_raw_df=parsedSeq_raw_df.drop(component_raw_now,axis=1)
                        continue
                    col_index=components_raw.index(component_raw_now)
                    component_corrected_now=self.settings.corrected_components_sorted[col_index]
                    opt_now=self.settings.correctOptDict[component_corrected_now]
                    
                    if (component_corrected_now in correctionDictionaries) and (opt_now["method"]=="correct" or opt_now["method"]=="from_bt"):
                        parsedSeq_raw_df[component_raw_now]=parsedSeq_raw_df[component_raw_now].map(lambda x: barcodeCorrecter.seq_correct_and_write(x,reference=correctionDictionaries[component_corrected_now]["correctionDict"]))
                    parsedSeq_raw_df.rename(columns={component_raw_now:component_raw_now+":"+component_corrected_now},inplace=True)

            if cnt==0:
                parsedSeq_raw_df.to_csv(self.settings.outFilePath_and_Prefix+"_correct_result.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
            else:
                parsedSeq_raw_df.to_csv(self.settings.outFilePath_and_Prefix+"_correct_result.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
        
    
    def make_s_value(self):
        t0=time.time()
        with gzip.open(self.settings.correctionDictPkl) as p:
            correctionDictionaries=pickle.load(p)

        parsedSeq_raw_chunk=pd.read_csv(self.settings.outFilePath_and_Prefix+"_correct_result.tsv.gz", sep='\t',chunksize=1000000)
        parsedQual_raw_chunk=pd.read_csv(self.settings.rawFastqPath["qual"], sep='\t',chunksize=1000000)
        options=[self.settings.correctOptDict[i] for i in self.settings.corrected_components_sorted]
        components_raw=[dic["src_raw_components"] for dic in options]
        # print(options,components_raw)

        
        ref_dic={}
        for component_corrected_now in correctionDictionaries:
            ref_tup=tuple(correctionDictionaries[component_corrected_now]["reference"])
            ref_dic[component_corrected_now]={k:v for k,v in zip(ref_tup,range(len(ref_tup)))}
        with gzip.open(self.settings.outFilePath_and_Prefix+"_sseq_to_svalue.pkl.gz",mode="wb") as p:
                pickle.dump(ref_dic,p)
                                
        for cat in ["seq","qual"]:  
            t0=time.time()
            if cat=="seq":
                df_chunk=parsedSeq_raw_chunk
            else:
                df_chunk=parsedQual_raw_chunk

            for cnt,parsedSeq_df in enumerate(df_chunk):
                print(cat,"Making a value table for chunk",cnt)
                res=pd.DataFrame()
                for ncol,component in enumerate(parsedSeq_df.columns):                    
                    if ncol==0:
                        res["Header"]=parsedSeq_df[component]
                    if ncol>0: 
                        if cat=="seq":
                            component_raw_now=component.split(":")[0]
                            component_corrected_now=component.split(":")[1]
                        else:
                            component_raw_now=component
                            if not component_raw_now in components_raw:
                                continue
                            col_index=components_raw.index(component_raw_now)
                            component_corrected_now=self.settings.corrected_components_sorted[col_index]

                        if component_corrected_now in correctionDictionaries:
                            # print(component_corrected_now,flush=True)
                            if cat=="seq":
                                res[component_corrected_now]=parsedSeq_df[component].map(lambda x: barcodeCorrecter.seq_to_val_ver2(x,dic=ref_dic[component_corrected_now]))
                            else:
                                res[component_corrected_now]=parsedSeq_df[component_raw_now].map(barcodeCorrecter.qualityProcessing)
                print(cat,time.time()-t0)

                if cnt==0:
                    if cat=="seq":
                        res.to_csv(self.settings.outFilePath_and_Prefix+"_correct_srcValue.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
                    else:
                        res.to_csv(self.settings.outFilePath_and_Prefix+"_correct_srcQual.tsv.gz",mode="w",compression="gzip",sep="\t",index=False)
                else:
                    if cat=="seq":
                        res.to_csv(self.settings.outFilePath_and_Prefix+"_correct_srcValue.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
                    else:
                        res.to_csv(self.settings.outFilePath_and_Prefix+"_correct_srcQual.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
                print(cat,time.time()-t0)
