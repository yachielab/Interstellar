from shutil import copy
from . import settingImporter
from . import barcodeConverter
from . import settingRequirementCheck
import regex
import pandas as pd
import numpy as np
import gzip
import pickle
import time
import os
import collections
import itertools as it
import copy

class settings_export(object):
    def __init__(self,opt):
        self.opt=opt 
    def settingGetter(self):
        cfg=settingImporter.readconfig(self.opt.config)
        cfg={k:settingImporter.configClean(cfg[k]) for k in cfg}
        cfg=settingRequirementCheck.setDefaultConfig(cfg)
        cfg_value_ext=settingImporter.config_extract_value_ext(cfg)
        func_dict_ext=settingImporter.func_check(cfg_value_ext)
        cfg_value_trans = settingImporter.config_extract_value_trans(cfg)
        func_dict=settingImporter.func_check_trans(cfg_value_trans)
        self.dest_segments=cfg_value_trans["dest_segment"]

        self.destValue=self.opt.destValue
        self.destQual=self.opt.destQual
        self.rawSeq=self.opt.rawSeq
        self.rawQual=self.opt.rawQual
        self.is_barcodelist=self.opt.export_bclist
        size_info=self.opt.size_info
        with gzip.open(size_info,mode="rb") as p:
            self.size_info=pickle.load(p)
        
        # exportOptDict={}
        # for i in self.export_components:
        #     exportOption_now=cfg_export[i]
        #     dict_now=settingImporter.exportOptionParse(exportOption_now)
        #     exportOptDict[i]=dict_now
        
        exportReadStructure={}
        for i in cfg_value_trans:
            if i in ["READ1_STRUCTURE","READ2_STRUCTURE","INDEX1_STRUCTURE","INDEX2_STRUCTURE"] and cfg_value_trans.get(i):
                exportReadStructure[i]=cfg_value_trans[i].split("+")
        self.func_dict=func_dict
        self.func_dict_ext=func_dict_ext
        self.exportReadStructure=exportReadStructure

        self.dest_barcode_segment=[i for i in self.dest_segments if func_dict[i]["func_ordered"][0]=="WHITELIST_ASSIGNMENT" or func_dict[i]["func_ordered"][0]=="RANDSEQ_ASSIGNMENT"]
        outname=self.opt.outname
        outdir=self.opt.outdir
        self.outFilePath_and_Prefix=outdir+"/"+outname


class BARISTA_EXPORT(object):
    def __init__(self,settings):
        self.settings=settings

    def generateReferences(self):
        print("generating reference...",flush=True)
        func_dict=self.settings.func_dict
        func_dict_cp=copy.deepcopy(func_dict)
        referenceDict=collections.defaultdict(list)
        dest_segments=self.settings.dest_segments
        processed_components_list=[]
        for seg in self.settings.dest_barcode_segment:
            if seg in processed_components_list:
                continue
            processed_components_list.append(seg)
            opt_now=func_dict_cp[seg]
            func_now=opt_now["func_ordered"][0]
            opt_now[func_now]["source"]="+".join(opt_now[func_now]["source"]) #concatenate source value (stored in a list as a default)

            seg_split=seg.split(",")

            """
            seg = destseg1,destseg2,destseg3
            seg_split = [destseg1,destseg2,destseg3]
            seg_each = destseg1
            """
            reference_group=[]
            for n_each,seg_each in enumerate(seg_split):
                if func_now=="WHITELIST_ASSIGNMENT" or func_now=="RANDSEQ_ASSIGNMENT":
                    #bc sort
                    if func_now=="WHITELIST_ASSIGNMENT" and opt_now["WHITELIST_ASSIGNMENT"].get("correspondence_table"):
                        wl_now=opt_now["WHITELIST_ASSIGNMENT"]["path"][n_each]
                        with open(wl_now,mode="rt",encoding="utf-8") as f:
                            referenceDict[seg_each]=[regex.sub("\n","",i) for i in f]
                        continue
                    
                    #normal
                    ref_now=barcodeConverter.buildReference(opt_now,n_each,self.settings.size_info)
                    reference_group.append(ref_now)

                # d_val_now=opt_now["d_val"]
                # component_group=[export_component]
                # group_now=opt_now["combination_group"]
                # for otherComponent in export_components_list:
                #     opt_other=exportOptDict[otherComponent]
                #     if opt_other.get("combination_group") and opt_other.get("combination_group")==group_now and not otherComponent==export_component:
                #         processed_components_list.append(otherComponent)
                #         component_group.append(otherComponent)
                #         ref_other=barcodeConverter.buildReference(opt_other,self.settings.size_info)
                #         reference_group.append(ref_other)
            reference_product=[]
            maxsize=self.settings.size_info[opt_now[func_now]["source"]]
            if maxsize<200000:
                maxsize=200000
            for cnt,i in enumerate(it.product(*reference_group)):
                if cnt >= maxsize:
                    break
                reference_product.append("_".join(i))

            for combinedReference in reference_product:
                splittedReference=combinedReference.split("_")
                for idx,eachComponent in enumerate(seg_split):
                    referenceDict[eachComponent].append(splittedReference[idx])
            
            # else:
            #     referenceDict[export_component]=ref_now
            for eachComponent in seg_split:
                print("Destination barcode library size: ",eachComponent,len(referenceDict[eachComponent]),flush=True)
 
        self.referenceDict=referenceDict


    def exportSequence(self):
        print("exporting sequence...",flush=True)
        
        func_dict=self.settings.func_dict
        funcdict_key_list=list(func_dict.keys())
        exportReadStructure=self.settings.exportReadStructure
        
        referenceDict=self.referenceDict
        export_read_exist=list(exportReadStructure.keys())
        
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
                    for k in funcdict_key_list:
                        if component+"," in k or ","+component in k or component==k:
                            funcdict_key=copy.copy(k)
                    try:
                        funcdict_key
                    except:
                        raise(component+" was not found in func_dict keys.")
                    

                    opt_now=copy.deepcopy(func_dict[funcdict_key])
                    func_now=opt_now["func_ordered"][0]
                    opt_now[func_now]["source"]="+".join(opt_now[func_now]["source"])

                    if func_now=="PASS":
                        s_seq_component_clean=opt_now[func_now]["source"]
                        func_tmp=self.settings.func_dict_ext[s_seq_component_clean]["func_ordered"][0]
                        s_seq_component_raw=self.settings.func_dict_ext[s_seq_component_clean][func_tmp]["source"]
                        s_seq_component=s_seq_component_raw+":"+s_seq_component_clean
                        seq_export_tmp=s_seq_chunk[s_seq_component].apply(barcodeConverter.genEqSeq,length=opt_now[func_now].get("length"),datatype="seq",add_nuc=opt_now.get("add_nucleotide"))
                        qual_export_tmp=s_qual_chunk[s_seq_component_raw].apply(barcodeConverter.genEqSeq,length=opt_now[func_now].get("length"),datatype="qual",baseQuality=None)
                    
                    elif func_now=="WHITELIST_ASSIGNMENT" or func_now=="RANDSEQ_ASSIGNMENT":
                        d_val_component=opt_now[func_now]["source"]

                        reference_now=referenceDict[component]
                        d_val_chunk[funcdict_key]=d_val_chunk[funcdict_key].map(int)
                        d_qual_chunk[funcdict_key]=d_qual_chunk[funcdict_key].map(int)
                        seq_export_tmp=d_val_chunk[funcdict_key].apply(barcodeConverter.getConvSeq,reference=reference_now)
                        df_tmp_cat=seq_export_tmp.str.cat(d_qual_chunk[funcdict_key].astype(str),sep="_")
                        qual_export_tmp=df_tmp_cat.map(barcodeConverter.getConvQual_ver2)

                        if self.settings.is_barcodelist:
                            barcode_correspondence[component]=seq_export_tmp

                    elif func_now=="CONSTANT":
                        seq_export_tmp=[opt_now[func_now]["sequence"]]*d_val_chunk.shape[0]
                        length_now=len(opt_now[func_now]["sequence"])
                        baseQuality=30 #hard coded
                        const_quality=chr(baseQuality+33)*length_now
                        qual_export_tmp=[const_quality]*d_val_chunk.shape[0]
                        
                    
                    if cnt_comp==0:
                        fastq_parse["seq"]=seq_export_tmp
                        fastq_parse["qual"]=qual_export_tmp
                    else:
                        fastq_parse["seq"]=fastq_parse["seq"].str.cat(seq_export_tmp,sep="")
                        fastq_parse["qual"]=fastq_parse["qual"].str.cat(qual_export_tmp,sep="")
                    print("processing for",component,"end:",time.time()-t0,flush=True)

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
            merged_idx=sorted(list(merged_idx))

            for read_now in export_read_exist:
                print(read_now)
                if read_now=="READ1_STRUCTURE":
                    read_now_out="R1"
                elif read_now=="READ2_STRUCTURE":
                    read_now_out="R2"
                elif read_now=="INDEX1_STRUCTURE":
                    read_now_out="I1"
                elif read_now=="INDEX2_STRUCTURE":
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