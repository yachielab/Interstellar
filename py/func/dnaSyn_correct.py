from . import settingImporter
from . import barcodeCorrecter
import regex
import collections
import gzip
import pickle
import csv
import pandas as pd

class settings_correct(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_correct=settingImporter.readconfig(cfgPath)["correct"]
        cfg_correct=settingImporter.configClean(cfg_correct)
        self.corrected_components=cfg_correct["src_corrected_components"].split(",")
        self.importPkl=self.opt.pickle
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
        self.yaxis_scale=self.opt.yaxis_scale
        self.show_summary=not self.opt.no_show_summary
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)
        # self.mk_s_value_components=cfg_correct["make_s_value"].split(",")
        
class BARISTA_CORRECT(object):
    def __init__(self,settings):
        self.settings=settings

    def importExtractedComponents(self):
        counterDict_pileup=collections.defaultdict(lambda: collections.Counter())
        for i in self.settings.importPkl:
            print("start unzip",flush=True)
            with gzip.open(i,mode="rb") as p:
                counter_tmp=pickle.load(p)
            print("start pileup",flush=True)
            for component in counter_tmp:
                counterDict_pileup[component].update(counter_tmp[component])
        counterDict=dict(counterDict_pileup)   
        # print("start export",flush=True)
        # with gzip.open(self.settings.outFilePath_and_Prefix+"_srcPileup.pkl.gz",mode="wb") as p:
        #     pickle.dump(counterDict,p)     

        self.counterDict=counterDict


    def correct_component(self):
        correspondicgDict={k:self.settings.correctOptDict[k]["src_raw_components"] for k in self.settings.correctOptDict}
        corresponding_key=list(correspondicgDict.keys())
        corresponding_val=list(correspondicgDict.values())
        correctionDictionaries={}
        # correctedQualDict={}
        for rawSegment in self.counterDict:
            if rawSegment not in corresponding_val:
                continue
            correctedComponent=corresponding_key[corresponding_val.index(rawSegment)]
            print("correction started for "+correctedComponent,flush=True)
            correctOpt=self.settings.correctOptDict[correctedComponent]
            #component_raw_corresponding=correctOpt["src_raw_components"]
            if correctOpt["method"]=="correct":
                correctedTables=barcodeCorrecter.bcCorrect(correctOpt,self.counterDict,self.settings.yaxis_scale,self.settings.show_summary,self.settings.outFilePath_and_Prefix)
                correctionDictionaries[correctedComponent]=correctedTables
            elif correctOpt["method"]=="equal":
                correctionDictionaries[correctedComponent]={}
                correctionDictionaries[correctedComponent]["reference"]=list(self.counterDict[correctOpt["src_raw_components"]].keys())
            elif correctOpt["method"]=="from_bt":
                correctionDictionaries[correctedComponent]={}
                bc_file=correctOpt["barcode_csv"]
                clstr_file=correctOpt["cluster_csv"]
                df_clstr=pd.read_csv(clstr_file)
                correctionDictionaries[correctedComponent]["reference"]=list(df_clstr["Center"])
                correctionDictionaries[correctedComponent]["correctionDict"]=barcodeCorrecter.gen_bt_dict(bc_file,clstr_file)
            elif correctOpt["method"]=="from_starcode":
                correctionDictionaries[correctedComponent]={}
                #need to parse starcode output to fill correctionDictionaries[correctedComponent]["reference"] and correctionDictionaries[correctedComponent]["correctionDict"]

        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCorrect.pkl.gz",mode="wb") as p:
            pickle.dump(correctionDictionaries,p)

        ref={i:correctionDictionaries[i]["reference"] for i in correctionDictionaries}
        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcReference.tsv.gz",mode="wt",encoding="utf-8") as wref:
            csvwriter_ref=csv.writer(wref,delimiter="\t")
            for component in ref:
                col1=component
                col2=",".join(ref[component])
                csvwriter_ref.writerow([col1,col2])

