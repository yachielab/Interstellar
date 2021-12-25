from . import settingImporter
from . import barcodeCorrecter
from . import settingRequirementCheck
import regex
import collections
import gzip
import pickle
import csv
import pandas as pd
import os
import re


class settings_correct(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext,dict_to_terminal = settingImporter.config_extract_value_ext(cfg)
        func_dict=settingImporter.func_check(cfg_value_ext)
        self.correctOptDict=func_dict
        self.corrected_components=cfg_value_ext["value_segment"]
        self.importPkl=self.opt.pickle
        # correctOptDict={}
        # for i in self.corrected_components:
        #    correctOption_now=cfg_correct[i]
        #    if not regex.search("^combination",correctOption_now):
        #        correctOptDict[i]=settingImporter.correctOptionParse(correctOption_now)
        
        # self.correctOptDict=correctOptDict
        self.yaxis_scale=self.opt.yaxis_scale
        self.show_summary=not self.opt.no_show_summary
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outdir=outdir
        self.outFilePath_and_Prefix=outdir+"/"+outname
        # self.mk_s_value_components=cfg_correct["make_s_value"].split(",")

        
class BARISTA_CORRECT(object):
    def __init__(self,settings):
        self.settings=settings


    def importExtractedComponents(self):
        counterDict_pileup=collections.defaultdict(lambda: collections.Counter())
        for i in self.settings.importPkl:
            print("Merging file: "+i,flush=True)
            with gzip.open(i,mode="rb") as p:
                counter_tmp=pickle.load(p)
            for component in counter_tmp:
                counterDict_pileup[component].update(counter_tmp[component])
        counterDict=dict(counterDict_pileup)   

        self.counterDict=counterDict


    def correct_component(self):
        correspondicgDict={}
        for k in self.settings.corrected_components:
            func_tmp=self.settings.correctOptDict[k]["func_ordered"][0]
            correspondicgDict[k]=self.settings.correctOptDict[k][func_tmp]["source"]
        corresponding_key=list(correspondicgDict.keys())
        corresponding_val=list(correspondicgDict.values())
        correctionDictionaries={}
        # correctedQualDict={}
        for rawSegment in self.counterDict:
            if rawSegment not in corresponding_val:
                continue
            correctedComponent=corresponding_key[corresponding_val.index(rawSegment)]
            correctOpt=self.settings.correctOptDict[correctedComponent]
            #component_raw_corresponding=correctOpt["src_raw_components"]
            if "BARTENDER" in correctOpt["func_ordered"]:
                bartender_path=os.path.dirname(re.sub(r"\/$","",self.settings.outdir))+"/to_bt/to_bt"+"_"+rawSegment+"_bartender"
                correctionDictionaries[correctedComponent]={}
                bc_file=bartender_path+"_barcode.csv"
                clstr_file=bartender_path+"_cluster.csv"
                df_clstr=pd.read_csv(clstr_file)
                correctionDictionaries[correctedComponent]["reference"]=list(df_clstr["Center"])
                correctionDictionaries[correctedComponent]["correctionDict"]=barcodeCorrecter.gen_bt_dict(bc_file,clstr_file)
            elif "KNEE_CORRECT" in correctOpt["func_ordered"] or "WHITELIST_CORRECT" in correctOpt["func_ordered"]:
                correctedTables=barcodeCorrecter.bcCorrect(correctOpt,self.counterDict,self.settings.yaxis_scale,self.settings.show_summary,self.settings.outFilePath_and_Prefix)
                correctionDictionaries[correctedComponent]=correctedTables
            else:
                correctionDictionaries[correctedComponent]={}
                correctionDictionaries[correctedComponent]["reference"]=list(self.counterDict[rawSegment].keys())
                        # elif correctOpt["method"]=="from_starcode":
            #     correctionDictionaries[correctedComponent]={}
            #     #need to parse starcode output to fill correctionDictionaries[correctedComponent]["reference"] and correctionDictionaries[correctedComponent]["correctionDict"]

        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcCorrect.pkl.gz",mode="wb") as p:
            pickle.dump(correctionDictionaries,p)

        ref={i:correctionDictionaries[i]["reference"] for i in correctionDictionaries}
        with gzip.open(self.settings.outFilePath_and_Prefix+"_srcReference.tsv.gz",mode="wt",encoding="utf-8") as wref:
            csvwriter_ref=csv.writer(wref,delimiter="\t")
            for component in ref:
                col1=component
                col2=",".join(ref[component])
                csvwriter_ref.writerow([col1,col2])

