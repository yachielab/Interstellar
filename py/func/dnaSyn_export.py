from . import settingImporter
from . import barcodeConverter
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import time
import os
import collections
import itertools as it

class settings_export(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfgPath=self.opt.config
        cfg_export=settingImporter.readconfig(cfgPath)["export"]
        cfg_export=settingImporter.configClean(cfg_export)
        self.export_components=cfg_export["export_components"].split(",")
        self.destValue=self.opt.destValue
        self.destQual=self.opt.destQual
        self.rawSeq=self.opt.rawSeq
        self.rawQual=self.opt.rawQual
        self.is_barcodelist=self.opt.export_bclist

        size_info=self.opt.size_info
        with gzip.open(size_info,mode="rb") as p:
            self.size_info=pickle.load(p)
        
        exportOptDict={}
        for i in self.export_components:
            exportOption_now=cfg_export[i]
            dict_now=settingImporter.exportOptionParse(exportOption_now)
            exportOptDict[i]=dict_now
        
        exportReadStructure={}
        for i in cfg_export:
            if i in ["Read1_dest","Read2_dest","Index1_dest","Index2_dest"] and cfg_export.get(i):
                exportReadStructure[i]=cfg_export[i].split("+")
        self.exportOptDict=exportOptDict
        self.exportReadStructure=exportReadStructure
        
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=regex.sub("/$","",str(outdir))+"/"+str(outname)
class BARISTA_EXPORT(object):
    def __init__(self,settings):
        self.settings=settings

    def generateReferences(self):
        print("generating reference...",flush=True)
        exportOptDict=self.settings.exportOptDict
        referenceDict=collections.defaultdict(list)
        #typeDict={}
        export_components_list=self.settings.export_components
        processed_components_list=[]
        for export_component in export_components_list:
            if export_component in processed_components_list:
                continue
            print(export_component)
            processed_components_list.append(export_component)
            opt_now=exportOptDict[export_component]

            #typeDict[export_component]=opt_now["type"]
            if opt_now["type"]=="whitelistExport" or opt_now["type"]=="randomExport":
                # if opt_now.get("whitelist_correspondence"):
                #     continue

                if opt_now.get("is_sorted"):
                    with open(opt_now["whitelist"],mode="rt",encoding="utf-8") as f:
                        referenceDict[export_component]=[regex.sub("\n","",i) for i in f]
                    continue
                
                ref_now=barcodeConverter.buildReference(opt_now,self.settings.size_info)

                if opt_now.get("combination_group"):
                    d_val_now=opt_now["d_val"]
                    reference_group=[ref_now]
                    component_group=[export_component]
                    group_now=opt_now["combination_group"]
                    for otherComponent in export_components_list:
                        opt_other=exportOptDict[otherComponent]
                        if opt_other.get("combination_group") and opt_other.get("combination_group")==group_now and not otherComponent==export_component:
                            processed_components_list.append(otherComponent)
                            component_group.append(otherComponent)
                            ref_other=barcodeConverter.buildReference(opt_other,self.settings.size_info)
                            reference_group.append(ref_other)
                    
                    reference_product=[]
                    maxsize=self.settings.size_info[d_val_now]
                    if maxsize<200000:
                        maxsize=200000
                    for cnt,i in enumerate(it.product(*reference_group)):
                        if cnt >= maxsize:
                            break
                        reference_product.append("_".join(i))

                    for combinedReference in reference_product:
                        splittedReference=combinedReference.split("_")
                        for idx,eachComponent in enumerate(component_group):
                            referenceDict[eachComponent].append(splittedReference[idx])

                else:
                    referenceDict[export_component]=ref_now
                    
                print("Destination barcode library size:",len(referenceDict[export_component]),flush=True)
        #self.typeDict=typeDict
        # for i in referenceDict:
        #     L="\n".join(referenceDict[i])+"\n"
        #     with open(self.settings.outFilePath_and_Prefix+"_"+i+".txt",mode="wt") as w:
        #         w.write(L)
        self.referenceDict=referenceDict

    def exportSequence(self):
        print("exporting sequence...",flush=True)
        
        exportOptDict=self.settings.exportOptDict
        exportReadStructure=self.settings.exportReadStructure
        
        referenceDict=self.referenceDict
        export_read_exist=list(exportReadStructure.keys())
        # with gzip.open(self.settings.reindex,mode="rb") as p:
        #     reindex_dict=pickle.load(p)

        d_val=pd.read_csv(self.settings.destValue,sep='\t',dtype=str,chunksize=500000)
        d_qual=pd.read_csv(self.settings.destQual,sep='\t',dtype=str,chunksize=500000)
        s_seq=pd.read_csv(self.settings.rawSeq,sep='\t',dtype=str,chunksize=500000)
        s_qual=pd.read_csv(self.settings.rawQual,sep='\t',dtype=str,chunksize=500000)

        chunkCount=0
        for d_val_chunk,d_qual_chunk,s_seq_chunk,s_qual_chunk in zip(d_val,d_qual,s_seq,s_qual):
            survived_idx_dict={}
            barcode_correspondence=pd.DataFrame()
            barcode_correspondence["Header"]=s_seq_chunk["Header"]
            for read_now in export_read_exist:
                print("start processing for",read_now,"and chunk",chunkCount,"...",flush=True)
                fastq_parse=pd.DataFrame()
                fastq_parse["Header"]=s_seq_chunk["Header"]
                fastq_parse["Third"]=["+"]*fastq_parse.shape[0]
                
                readStructure_now=exportReadStructure[read_now]
                for cnt_comp,component in enumerate(readStructure_now):
                    t0=time.time()
                    # print("\n"+component,flush=True)
                    opt_now=exportOptDict[component]

                    if opt_now["type"]=="equalExport":
                        s_seq_component=opt_now["src_raw_components"]
                        seq_export_tmp=s_seq_chunk[s_seq_component].apply(barcodeConverter.genEqSeq,length=opt_now.get("length"),datatype="seq",add_nuc=opt_now.get("add_nucleotide"))
                        qual_export_tmp=s_qual_chunk[s_seq_component].apply(barcodeConverter.genEqSeq,length=opt_now.get("length"),datatype="qual",baseQuality=None)
                    
                    elif opt_now["type"]=="whitelistExport" or opt_now["type"]=="randomExport":
                        d_val_component=opt_now["d_val"]
                        # if d_val_component in reindex_dict:
                        #     reindex_now=reindex_dict[d_val_component]
                        # else:
                        #     reindex_now=None

                        reference_now=referenceDict[component]
                        
                        t_map=time.time()
                        d_val_chunk[d_val_component]=d_val_chunk[d_val_component].map(int)
                        d_qual_chunk[d_val_component]=d_qual_chunk[d_val_component].map(int)
                        # print("mappint time:",time.time()-t_map)
                        t_convseq=time.time()
                        # seq_export_tmp=d_val_chunk[d_val_component].apply(barcodeConverter.getConvSeq,reference=reference_now,reindex=reindex_now)
                        seq_export_tmp=d_val_chunk[d_val_component].apply(barcodeConverter.getConvSeq,reference=reference_now)
                        # print("convert time:",time.time()-t_convseq)
                        df_tmp_cat=seq_export_tmp.str.cat(d_qual_chunk[d_val_component].astype(str),sep="_")
                        t_getqual=time.time()
                        qual_export_tmp=df_tmp_cat.map(barcodeConverter.getConvQual_ver2)
                        # print("get qual seq:",time.time()-t_getqual)

                        if self.settings.is_barcodelist:
                            barcode_correspondence[component]=seq_export_tmp

                    elif opt_now["type"]=="constantExport":
                        seq_export_tmp=[opt_now["sequence"]]*d_val_chunk.shape[0]
                        length_now=len(opt_now["sequence"])
                        if opt_now.get("baseQuality"):
                            baseQuality=opt_now["baseQuality"]
                        else:
                            baseQuality=30
                        const_quality=chr(baseQuality+33)*length_now
                        qual_export_tmp=[const_quality]*d_val_chunk.shape[0]

                    if cnt_comp==0:
                        fastq_parse["seq"]=seq_export_tmp
                        fastq_parse["qual"]=qual_export_tmp
                    else:
                        fastq_parse["seq"]=fastq_parse["seq"].str.cat(seq_export_tmp,sep="")
                        fastq_parse["qual"]=fastq_parse["qual"].str.cat(qual_export_tmp,sep="")
                    print("processing for",component,"end:",time.time()-t0,flush=True)

                # print("start NA trimming...\n",flush=True)
                fastq_parse=fastq_parse.dropna(how="any")
                survived_idx=fastq_parse.index

                survived_idx_dict[read_now]=set(survived_idx)

                # for eachline in ["Header","seq","Third","qual"]:
                #     if eachline=="Header":
                #         fastq=fastq_parse[eachline]
                #     else:
                #         fastq=fastq.str.cat(fastq_parse[eachline],sep="\n")
                fastq_parse=fastq_parse[["Header","seq","Third","qual"]]
                with open(self.settings.outFilePath_and_Prefix+read_now+"tmp.pkl",mode="wb") as p:
                    pickle.dump(fastq_parse,p)

            if self.settings.is_barcodelist:
                barcode_correspondence=barcode_correspondence.fillna("-")
                if chunkCount==0:
                    barcode_correspondence.to_csv(self.settings.outFilePath_and_Prefix+"_bclist.tsv.gz",mode="w",compression="gzip",sep="\t",index=False,header=True)
                else:
                    barcode_correspondence.to_csv(self.settings.outFilePath_and_Prefix+"_bclist.tsv.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
              
            print("merging fastq...",flush=True)
            merged_idx=set()
            for c,key in enumerate(survived_idx_dict):
                if c==0:
                    merged_idx=survived_idx_dict[key]
                else:
                    merged_idx&=survived_idx_dict[key]

            for read_now in export_read_exist:
                if read_now=="Read1_dest":
                    read_now_out="R1"
                elif read_now=="Read2_dest":
                    read_now_out="R2"
                elif read_now=="Index1_dest":
                    read_now_out="I1"
                elif read_now=="Index2_dest":
                    read_now_out="I2"

                with open(self.settings.outFilePath_and_Prefix+read_now+"tmp.pkl",mode="rb") as p:
                    fastq_parse=pickle.load(p)
                fastq_parse=fastq_parse.loc[merged_idx]
                fastq_parse=fastq_parse.stack()
                fastq_parse=fastq_parse.reset_index()
                fastq_parse=pd.DataFrame(fastq_parse[0])
                print(read_now_out,"exporting...")

                if chunkCount==0:
                    fastq_parse.to_csv(self.settings.outFilePath_and_Prefix+"_"+read_now_out+".fastq.gz",mode="w",compression="gzip",sep="\t",index=False,header=False)
                else:
                    fastq_parse.to_csv(self.settings.outFilePath_and_Prefix+"_"+read_now_out+".fastq.gz",mode="a",compression="gzip",sep="\t",index=False,header=False)
            chunkCount+=1
            print("\n")

        for read_now in export_read_exist:
            os.remove(self.settings.outFilePath_and_Prefix+read_now+"tmp.pkl")