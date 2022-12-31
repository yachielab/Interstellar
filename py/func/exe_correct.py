from . import settingImporter
from . import barcodeCorrecter
from . import settingRequirementCheck
from . import interstellar_setup
import regex
import collections
import gzip
import pickle
import csv
import pandas as pd
import os
import re
import glob
import subprocess
import sys


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
        self.importPkl = settingImporter.parseInputFileList(self.opt.pickle)

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
        self.ncore = int(self.opt.ncore)
        # self.mk_s_value_components=cfg_correct["make_s_value"].split(",")

        
class BARISTA_CORRECT(object):
    def __init__(self,settings):
        self.settings=settings


    # Pile up the count information
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


    # Sequence error correction
    def correct_component(self):

        # Get correspondence between source segment name and corrected source segment name
        correspondingDict={}
        for k in self.settings.corrected_components:
            func_tmp=self.settings.correctOptDict[k]["func_ordered"][0]
            correspondingDict[k]=self.settings.correctOptDict[k][func_tmp]["source"]
        corresponding_key=list(correspondingDict.keys()) # corrected source segment name
        corresponding_val=list(correspondingDict.values()) # raw source segment name
        correctionDictionaries={}

        # Setup custom corection
        is_qc=interstellar_setup.checkRequiredFile("_srcSeq.QC.pkl",glob.glob(os.path.dirname(re.sub(r"\/$","",self.settings.outdir))+"/qc/*"))
        custom_correction_segments=list()
        for rawSegment in self.counterDict:
            if rawSegment not in corresponding_val:
                continue

            correctedComponent=corresponding_key[corresponding_val.index(rawSegment)]
            correctOpt=self.settings.correctOptDict[correctedComponent]
            if "CUSTOM_CORRECTION" in correctOpt["func_ordered"]:
                custom_correction_segments.append(rawSegment)
        
        if is_qc:
            src_seq_paths=glob.glob(os.path.dirname(re.sub(r"\/$","",self.settings.outdir))+"/qc/*_srcSeq.QC.pkl")
        else:
            src_seq_paths=glob.glob(os.path.dirname(re.sub(r"\/$","",self.settings.outdir))+"/import/*_srcSeq.pkl")
        barcodeCorrecter.custom_correction_setup(custom_correction_segments,self.settings.outdir,src_seq_paths)

        #Correction        
        for rawSegment in self.counterDict:
            if rawSegment not in corresponding_val:
                continue
            correctedComponent=corresponding_key[corresponding_val.index(rawSegment)]
            correctOpt=self.settings.correctOptDict[correctedComponent]

            if "CUSTOM_CORRECTION" in correctOpt["func_ordered"]:         
                input_filename=self.settings.outdir+"/custom"+"_"+rawSegment+".input.csv"
                output_filename=self.settings.outdir+"/custom"+"_"+rawSegment+".output.csv"
                cmd=["bash",correctOpt["CUSTOM_CORRECTION"]["shell_script"],input_filename,output_filename]
                cmdline=" ".join(cmd)
                s=subprocess.run(cmdline,shell=True)
                if s.returncode != 0:
                    print("Custom correction failed', file=sys.stderr")
                    sys.exit(1)

                correctionDictionaries[correctedComponent]={}
                df_corrected = pd.read_csv(output_filename,header=None,sep="\t")
                counter_dict=collections.Counter(df_corrected[1])
                correctionDictionaries[correctedComponent]["reference"]=sorted(counter_dict.keys(),key=counter_dict.__getitem__,reverse=True)
                correctionDictionaries[correctedComponent]["correctionDict"]=barcodeCorrecter.gen_custom_dict(df_corrected)

            elif "BARTENDER_CORRECTION" in correctOpt["func_ordered"]:
                bartender_path=os.path.dirname(re.sub(r"\/$","",self.settings.outdir))+"/to_bt/to_bt"+"_"+rawSegment+"_bartender"
                correctionDictionaries[correctedComponent]={}
                bc_file=bartender_path+"_barcode.csv"
                clstr_file=bartender_path+"_cluster.csv"
                df_clstr=pd.read_csv(clstr_file)
                correctionDictionaries[correctedComponent]["reference"]=list(df_clstr["Center"])
                correctionDictionaries[correctedComponent]["correctionDict"]=barcodeCorrecter.gen_bt_dict(bc_file,clstr_file)

            elif "I2M_CORRECTION" in correctOpt["func_ordered"] or "M2A_CORRECTION" in correctOpt["func_ordered"]:
                min_num_reads = correctOpt["I2M_CORRECTION"]["min_num_reads"] if "I2M_CORRECTION" in correctOpt["func_ordered"] else correctOpt["M2A_CORRECTION"]["min_num_reads"]
                correctedTables=barcodeCorrecter.bcCorrect(correctOpt,self.counterDict,self.settings.yaxis_scale,self.settings.show_summary,self.settings.outFilePath_and_Prefix,self.settings.ncore,min_num_reads)
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

