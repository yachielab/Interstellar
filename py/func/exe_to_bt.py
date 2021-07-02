from . import settingImporter
from . import settingRequirementCheck
import regex
import pandas as pd
import numpy as np
import subprocess
import sys


class settings_to_bt(object):
    def __init__(self,opt):
        self.opt=opt
    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext = settingImporter.config_extract_value_ext(cfg)
        func_dict=settingImporter.func_check(cfg_value_ext)
        self.func_dict=func_dict
        self.bt_targets=[func_dict[val]["BARTENDER"]["source"] for val in func_dict if "BARTENDER" in func_dict[val]]
        # if cfg_to_bt["need_match"]:
        #     self.need_match=cfg_to_bt["need_match"].split(",")
        # else:
        self.need_match=None
        self.seq=self.opt.rawSeq
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)

class BARISTA_TO_BT(object):
    def __init__(self,settings):
        self.settings=settings
    def to_bt(self):
        dic_rowcount={k:0 for k in self.settings.bt_targets}
        for processCount,input_filepath in enumerate(self.settings.seq):
            df_tmp=pd.read_csv(input_filepath,sep="\t")

            if self.settings.need_match:
                df_selected=df_tmp[self.settings.need_match].replace("-",np.nan)
                df_selected=df_selected.dropna()
                row_index_survive=list(df_selected.index)
                print("Survived index",row_index_survive)
                df_tmp=df_tmp.loc[row_index_survive]

            for bc in self.settings.bt_targets:
                if processCount==0:
                    with open(self.settings.outFilePath_and_Prefix+"_"+bc+".csv",mode="wt") as w:
                        for entry in list(df_tmp[bc]):
                            if not entry=="-":
                                dic_rowcount[bc]+=1
                                w.write(entry+","+str(dic_rowcount[bc])+"\n")
                else:
                    with open(self.settings.outFilePath_and_Prefix+"_"+bc+".csv",mode="at") as w:
                        for entry in list(df_tmp[bc]):
                            if not entry=="-":
                                dic_rowcount[bc]+=1
                                w.write(entry+","+str(dic_rowcount[bc])+"\n")
    def bartender(self):
        func_dict_bartender=dict()
        for bc in self.settings.bt_target:
            for val in self.settings.func_dict:
                if "BARTENDER" in self.settings.func_dict[val]:
                    if self.settings.func_dict[val]["BARTENDER"]["source"]==bc:
                        func_dict_bartender[bc]=self.settings.func_dict[val]["BARTENDER"]

        for bc in self.settings.bt_target:
            input_filename=self.settings.outFilePath_and_Prefix+"_"+bc+".csv"
            cmd=["bartender_single_com","-f",input_filename,"-o",self.settings.outFilePath_and_Prefix+"_"+bc+"_bartender"]
            for key in func_dict_bartender:
                if not key=="source":
                    cmd+=[key,func_dict_bartender[key]]
        
            cmdline=" ".join(cmd)
            s=subprocess.run(cmdline,shell=True)
            if s.returncode != 0:
                print("qsub failed: import', file=sys.stderr")
                sys.exit(1)
