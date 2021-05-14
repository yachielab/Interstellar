from . import settingImporter
import regex
import pandas as pd
import numpy as np

class settings_to_bt(object):
    def __init__(self,opt):
        self.opt=opt
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_to_bt=settingImporter.readconfig(cfgPath)["to_bt"]
        cfg_to_bt=settingImporter.configClean(cfg_to_bt)
        self.bt_targets=cfg_to_bt["bt_targets"].split(",")
        if cfg_to_bt["need_match"]:
            self.need_match=cfg_to_bt["need_match"].split(",")
        else:
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