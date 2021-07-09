import configparser
import json
from argparse import ArgumentError, ArgumentParser
from types import FunctionType
from pandas._config.config import option_context
import regex
import gzip
import collections
import os
from distutils.util import strtobool
import pandas as pd
import re
import copy


def readconfig(input_cfg):
    cfg=configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(input_cfg)

    return cfg


def configClean(cfgDict):
    new_cfgDict={}
    for key in cfgDict:
        new_cfgDict[key]=regex.sub("\\t|\\n| ","",cfgDict[key])
    return new_cfgDict


def genSampleDir(proj_dir,samplesheet):
    os.makedirs(proj_dir+"/_flg",exist_ok=True)
    if samplesheet=="":
        os.makedirs(proj_dir+"/sample1",exist_ok=True)
        sampledir_list=[proj_dir+"/sample1"]
    else:
        samplesheet=pd.read_csv(samplesheet,sep="\t",header=None)
        samples_uniq=set(samplesheet[1])
        for i in samples_uniq:
            os.makedirs(proj_dir+"/"+i,exist_ok=True)
        sampledir_list=[proj_dir+"/"+i for i in samples_uniq]
        os.makedirs(proj_dir+"/_multisample",exist_ok=True)
        os.makedirs(proj_dir+"/_multisample/mergeTree",exist_ok=True)
        os.makedirs(proj_dir+"/_multisample/samplesheet",exist_ok=True)

    return sampledir_list


def config_extract_value_ext(cfg_raw):
    cfg=copy.deepcopy(cfg_raw)
    cfg_value_ext=cfg["value_extraction"]
    #sequence segment label extractioon
    segments=[]
    for key in cfg_value_ext:
        if "READ1_STRUCTURE" in key or "READ2_STRUCTURE" in key or "INDEX1_STRUCTURE" in key or "INDEX2_STRUCTURE" in key or "READ_FLASH" in key:
            segments+=re.findall(r"\?P\<([^\>]+)\>",cfg_value_ext[key])
    cfg_value_ext["segments"]=list(set(segments))

    #parental-local allocation
    cfg_value_ext["parent"]={}
    for key in cfg_value_ext:
        if not type(cfg_value_ext[key])==str:
            continue
        if re.search(r"CHILD_OF\([^\)]+\)",cfg_value_ext[key]):
            parent=cfg_value_ext[key].replace("CHILD_OF(","").replace(")","")
            cfg_value_ext["parent"][key]=parent
    # for seg in cfg_value_ext["segments"]:
    #     if cfg_value_ext.get(seg):
    #         parents_list=cfg_value_ext[seg].replace("CHILD_OF(","").replace(")","").split(",")
    #         cfg_value_ext["parent"][seg]=parents_list
    
    #value segment label extraction
    all_keys =set(cfg_value_ext.keys())
    keys_tmp =set([i for i in cfg_value_ext if "READ1_STRUCTURE" in i or "READ2_STRUCTURE" in i or "INDEX1_STRUCTURE" in i or "INDEX2_STRUCTURE" in i or "READ_FLASH" in i])
    keys_tmp|=set(["READ1_DIR","READ2_DIR","INDEX1_DIR","INDEX2_DIR"])
    keys_tmp|=set(["FLASH","FLASH_MIN_OVERLAP","FLASH_MAX_OVERLAP"])
    keys_tmp|=set(cfg_value_ext["segments"]+["segments","parent"]+list(cfg_value_ext["parent"].keys()))
    cfg_value_ext["value_segment"]=sorted(list(all_keys-keys_tmp))

    #detecting functions
    functions_used=[]
    for key in cfg_value_ext["value_segment"]:
        # print(cfg_value_ext[key])
        commandline=cfg_value_ext[key].split(">>")
        functions_used+=[re.sub(r"\(.*\)","",i) for i in commandline]
    cfg_value_ext["functions_used"]=set(functions_used)

    return cfg_value_ext



def funcIllegalOptionCheck(val,func_dict,funcname,available_option):
    if funcname in func_dict[val]["func_ordered"]:
        d_now=func_dict[val][funcname]
        for key in d_now:
            if not key in available_option:
                raise KeyError("Option "+key+" is not available for "+funcname)



def funcRequiredOptionCheck(val,func_dict,funcname,required_list):
    if funcname in func_dict[val]["func_ordered"]:
        d_now=func_dict[val][funcname]
        for key in required_list:
            if not key in d_now:
                raise KeyError("Option "+key+" is required for "+funcname)


def funcSetDefault(val,func_dict,funcname,option_list,default_list,opt_int=[]):
    if funcname in func_dict[val]["func_ordered"]:
        d_now=func_dict[val][funcname]
        for n,opt in enumerate(option_list):
            try:
                d_now[opt]
                if d_now[opt]=="":
                    raise KeyError()
            except:
                func_dict[val][funcname][opt]=default_list[n]

        if funcname=="WHITELIST_CORRECT":
            func_dict[val][funcname]["path"]=os.path.expanduser(func_dict[val][funcname]["path"])
        elif funcname=="WHITELIST_ASSIGNMENT":
            if "path" in func_dict[val][funcname]:
                func_dict[val][funcname]["path"]=list(map(os.path.expanduser,func_dict[val][funcname]["path"]))
            func_dict[val][funcname]["correspondence_table"]=os.path.expanduser(func_dict[val][funcname]["correspondence_table"])

        for i in opt_int:
            if type(func_dict[val][funcname][i])==list:
                func_dict[val][funcname][i]=list(map(int,func_dict[val][funcname][i]))
            else:
                func_dict[val][funcname][i]=int(func_dict[val][funcname][i])
    return func_dict


def func_parse(func_line,func_collection,dest=False):
    func_line=func_line.split(">>")
    d=dict()
    funcname_list=[]
    for func in func_line:
        m=re.search(r"([^\(]+)\(([^\)]*)\)",func)
        funcname=m.group(1)
        if not funcname in func_collection:
            raise KeyError("Function '"+funcname+"' doesn't exist.")
        contents=m.group(2).split(",")

        if not dest:
            contents={i.split(":")[0]:i.split(":")[1] for i in contents if not i==""}
            d[funcname]=contents
        else:
            opt=""
            dict_out=dict()
            for i in contents:
                if opt=="source" and len(i.split(":"))==1:
                    dict_out["source"].append(i)
                elif opt=="path" and len(i.split(":"))==1:
                    dict_out["path"].append(i)
                elif funcname=="RANDSEQ_ASSIGNMENT" and opt=="length":
                    dict_out["length"].append(i)
                else:
                    opt=i.split(":")[0]
                    val=i.split(":")[1]
                    if opt=="source" or opt=="path" or (funcname=="RANDSEQ_ASSIGNMENT" and opt=="length"):
                        dict_out[opt]=[val]
                    else:
                        dict_out[opt]=val
            d[funcname]=dict_out
        funcname_list.append(funcname)
    d["func_ordered"]=funcname_list
    return d



def parse_function_value_ext(cfg):
    func_collection=["QUALITY_FILTER","PASS","KNEE_CORRECT","WHITELIST_CORRECT","BARTENDER","VALUE"]
    value_names=cfg["value_segment"]
    func_dict=dict()
    for val in value_names:
        func_line=cfg[val]
        func_dict[val]=func_parse(func_line,func_collection)
    return func_dict



def parse_function_value_trans(cfg):
    func_collection=["WHITELIST_ASSIGNMENT","RANDSEQ_ASSIGNMENT","PASS","CONSTANT"]
    destseg_names=cfg["dest_segment"]
    func_dict=dict()
    for seg in destseg_names:
        func_line=cfg[seg]
        func_dict[seg]=func_parse(func_line,func_collection,dest=True)
    return func_dict



def func_check(cfg):
    func_dict=parse_function_value_ext(cfg)
    barcode_list=[]
    for val in func_dict:
        segment_now=[]
        func_oredered=func_dict[val]["func_ordered"]
        for func in func_oredered:
            content_dict=func_dict[val][func]
            if "source" in content_dict and not content_dict["source"]=="":
                segment_now.append(content_dict["source"])
        
        #Inconsistent source segment name
        if len(set(segment_now))>1:
            raise KeyError("Source segment names are inconsistent in "+val)
        if len(segment_now)==0:
            raise KeyError("Source segment should be specified in "+val)

        #Order check
        if "PASS" in func_oredered and len(func_oredered)>1:
            raise ArgumentError("Functions PASS cannot be performed with other functions.")      
        if ("BARTENDER" in func_oredered and "KNEE_CORRECT" in func_oredered) or ("BARTENDER" in func_oredered and "WHITELIST_CORRECT" in func_oredered):
            raise ArgumentError("Functions BARTENDER cannot be performed with KNEE_CORRECT or WHITELIST_CORRECT.")

        #QUALITY_FILTER Correct option check
        funcIllegalOptionCheck(val,func_dict,"QUALITY_FILTER",["source","min_base_score","min_avg_score"])
        func_dict=funcSetDefault(val,func_dict,"QUALITY_FILTER",["min_base_score","min_avg_score"],[5,20],opt_int=["min_base_score","min_avg_score"])

        #KNEE Correct option check
        funcIllegalOptionCheck(val,func_dict,"KNEE_CORRECT",["source","rank","dist"])
        func_dict=funcSetDefault(val,func_dict,"KNEE_CORRECT",["rank","dist"],["auto",1],opt_int=["dist"])

        #WHITELIST Correct option check
        funcIllegalOptionCheck(val,func_dict,"WHITELIST_CORRECT",["source","path","dist"])
        funcRequiredOptionCheck(val,func_dict,"WHITELIST_CORRECT",["path"])
        func_dict=funcSetDefault(val,func_dict,"WHITELIST_CORRECT",["dist"],[1],opt_int=["dist"])

        #Bartender option check
        funcIllegalOptionCheck(val,func_dict,"BARTENDER",["source","-c","-d","-l","-z","-s"])

        #VALUE option check
        funcIllegalOptionCheck(val,func_dict,"VALUE",["source"])

        #PASS Correct option check
        funcIllegalOptionCheck(val,func_dict,"PASS",["source","length"])
        funcRequiredOptionCheck(val,func_dict,"PASS",["source"])
        func_dict=funcSetDefault(val,func_dict,"PASS",["length"],[0],opt_int=["length"])

        #Values
        if "VALUE" in func_oredered:
            barcode_list.append(segment_now[0])
        for func in func_dict[val]:
            if not func=="func_ordered":
                func_dict[val][func]["source"]=segment_now[0]
    func_dict["barcode_list"]=barcode_list
    return func_dict



def config_extract_value_trans(cfg_raw):
    cfg=copy.deepcopy(cfg_raw)
    cfg_value_ext=config_extract_value_ext(cfg)
    cfg_value_trans=cfg["value_translation"]
    cfg_value_trans["parent"]=cfg_value_ext["parent"]
    
    #value segment label extraction
    all_keys =set(cfg_value_trans.keys())
    keys_tmp =set([i for i in cfg_value_trans if "READ1_STRUCTURE" in i or "READ2_STRUCTURE" in i or "INDEX1_STRUCTURE" in i or "INDEX2_STRUCTURE" in i])
    keys_tmp|=set(["parent"])
    cfg_value_trans["dest_segment"]=list(all_keys-keys_tmp)

    #detecting functions
    functions_used=[]
    for key in cfg_value_trans["dest_segment"]:
        commandline=cfg_value_trans[key]
        functions_used.append(re.sub(r"\(.*\)","",commandline))
    cfg_value_trans["functions_used"]=set(functions_used)

    #sort bc?
    cfg_value_trans["bc_sort"]=False
    cfg_value_trans["val2table"]={}
    func_dict=parse_function_value_trans(cfg_value_trans)
    for dest_seg in func_dict:
        for fun in func_dict[dest_seg]["func_ordered"]:
            if "correspondence_table" in func_dict[dest_seg][fun]:
                if not func_dict[dest_seg][fun]["correspondence_table"]=="":
                    cfg_value_trans["bc_sort"]=True
                    cfg_value_trans["val2table"][dest_seg]=func_dict[dest_seg][fun]["correspondence_table"]
    
    #available segments
    cfg_value_trans["available_seg"]=[]
    for L in cfg_value_trans["dest_segment"]:
        L=L.split(",")
        cfg_value_trans["available_seg"]+=L
   
    return cfg_value_trans


def func_check_trans(cfg):
    func_dict=parse_function_value_trans(cfg)
    for val in func_dict:
        # segment_now=[]
        # func_oredered=func_dict[val]["func_ordered"]
        # for func in func_oredered:
        #     content_dict=func_dict[val][func]
        #     if "source" in content_dict and not content_dict["source"]=="":
        #         segment_now+=content_dict["source"]
        
        #WHITELIST_ASSIGNMENT Correct option check
        funcIllegalOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source","path","correspondence_table"])
        func_dict=funcSetDefault(val,func_dict,"WHITELIST_ASSIGNMENT",["correspondence_table"],[""])
        if "WHITELIST_ASSIGNMENT" in func_dict[val] and func_dict[val]["WHITELIST_ASSIGNMENT"]["correspondence_table"]=="":
            funcRequiredOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source","path"])
        else:
            funcRequiredOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source"])

        #RANDSEQ_ASSIGNMENT Correct option check
        funcIllegalOptionCheck(val,func_dict,"RANDSEQ_ASSIGNMENT",["source","length"])
        funcRequiredOptionCheck(val,func_dict,"RANDSEQ_ASSIGNMENT",["source","length"])
        func_dict=funcSetDefault(val,func_dict,"RANDSEQ_ASSIGNMENT",[],[],opt_int=["length"])

        #CONSTANT Correct option check
        funcIllegalOptionCheck(val,func_dict,"CONSTANT",["sequence"])
        funcRequiredOptionCheck(val,func_dict,"CONSTANT",["sequence"])
        func_dict=funcSetDefault(val,func_dict,"CONSTANT",["source"],[["constant"]])

        #PASS Correct option check
        funcIllegalOptionCheck(val,func_dict,"PASS",["source","length","add"])
        funcRequiredOptionCheck(val,func_dict,"PASS",["source"])
        func_dict=funcSetDefault(val,func_dict,"PASS",["length","add"],[0,"A"],opt_int=["length"])

    return func_dict



def config_extract_value_demulti(cfg_raw):
    cfg=copy.deepcopy(cfg_raw)
    cfg_value_ext=config_extract_value_ext(cfg)
    cfg_value_demulti=cfg["demultiplex"]
    
    #value segment label extraction
    # all_keys =set(cfg_value_trans.keys())
    # keys_tmp =set([i for i in cfg_value_trans if "READ1_STRUCTURE" in i or "READ2_STRUCTURE" in i or "INDEX1_STRUCTURE" in i or "INDEX2_STRUCTURE" in i])
    # keys_tmp|=set(["parent"])
    # cfg_value_trans["dest_segment"]=list(all_keys-keys_tmp)

    #detecting functions
    # functions_used=[]
    # for key in cfg_value_trans["dest_segment"]:
    #     commandline=cfg_value_trans[key]
    #     functions_used.append(re.sub(r"\(.*\)","",commandline))
    # cfg_value_trans["functions_used"]=set(functions_used)

    # #sort bc?
    # cfg_value_trans["bc_sort"]=False
    # cfg_value_trans["val2table"]={}
    # func_dict=parse_function_value_trans(cfg_value_trans)
    # for dest_seg in func_dict:
    #     for fun in func_dict[dest_seg]["func_ordered"]:
    #         if "correspondence_table" in func_dict[dest_seg][fun]:
    #             if not func_dict[dest_seg][fun]["correspondence_table"]=="":
    #                 cfg_value_trans["bc_sort"]=True
    #                 cfg_value_trans["val2table"][dest_seg]=func_dict[dest_seg][fun]["correspondence_table"]
    
    #available segments
    cfg_value_demulti["available_seg"]=[]
    for i in cfg_value_ext["value_segment"]:
        cfg_value_demulti["available_seg"].append(i)
   
    return cfg_value_demulti



def getQscoreDict(func_dict):
    qscore_dict=dict()
    for val in func_dict:
        if "QUALITY_FILTER" in func_dict[val]:
            segment=func_dict[val]["QUALITY_FILTER"]["source"]
            # if "min_base_score" in func_dict[val]["QUALITY_FILTER"]:
            #     min_base=int(func_dict[val]["QUALITY_FILTER"]["min_base_score"])
            # else:
            #     min_base=5
            
            # if "min_avg_score" in func_dict[val]["QUALITY_FILTER"]:
            #     min_avg=int(func_dict[val]["QUALITY_FILTER"]["min_avg_score"])
            # else:
            #     min_avg=20
            qscore_dict[segment]={"min_base":func_dict[val]["QUALITY_FILTER"]["min_base_score"],"min_avg":func_dict[val]["QUALITY_FILTER"]["min_avg_score"]}
    return qscore_dict


def getValFromRawsegment(raw_list,cfg,func_dict,segment_name):
    outlist=[]
    for rawseg in raw_list:
        flg=False
        for val in cfg[segment_name]:
            tmp_funcname=func_dict[val]["func_ordered"][0]
            corresponding_raw_segname=func_dict[val][tmp_funcname]["source"]
            if rawseg==corresponding_raw_segname:
                outlist.append(val)
                flg=True
        if not flg:
            raise KeyError("The segment ",rawseg," is not found in the value configuration. See the config file.")
    outlist=sorted(outlist)
    outname=",".join(outlist)
    return outname


def getDestFromValsegment(raw_list,cfg,func_dict,segment_name):
    outlist=[]
    for rawseg in raw_list:
        flg=False
        for val in cfg[segment_name]:
            tmp_funcname=func_dict[val]["func_ordered"][0]
            corresponding_raw_segname=func_dict[val][tmp_funcname]["source"]
            if rawseg==corresponding_raw_segname:
                outlist.append(val)
                flg=True
        if not flg:
            raise KeyError("The segment ",rawseg," is not found in the value configuration. See the config file.")
    outlist=sorted(outlist)
    outname=",".join(outlist)
    return outname


def getAllocation(func_dict_ext,cfg,cfg_ext):
    raw_child2parent_dict=cfg["parent"]

    #Parental-local allocations in value segment
    val_child2parent_dict=dict()
    for child in raw_child2parent_dict:
        parent=raw_child2parent_dict[child]
        child_split=sorted(child.split(","))
        parent_split=sorted(parent.split(","))
        val_child=getValFromRawsegment(child_split,cfg_ext,func_dict_ext,segment_name="value_segment")
        val_parent=getValFromRawsegment(parent_split,cfg_ext,func_dict_ext,segment_name="value_segment")
        val_child2parent_dict[val_child]=val_parent
    
    return val_child2parent_dict

    # #Parental-local alloications in dest segment
    # dest_child2parent_dict=dict()
    # for child in val_child2parent_dict:
    #     parent=val_child2parent_dict[child]
    #     val_child=",".join(sorted(child.split(",")))
    #     val_parent=",".join(sorted(parent.split(",")))
    #     dest_child=getDestFromValsegment(val_child,cfg,func_dict,segment_name="dest_segment")
    #     dest_parent=getDestFromValsegment(val_parent,cfg,func_dict,segment_name="dest_segment")
    #     val_child2parent_dict[dest_child]=dest_parent

    



def getImportSeqFromFiles(file_seq_gz,file_qual_gz):
    seq={}
    qual={}
    with gzip.open(file_seq_gz,mode="rt") as seqf:
        for rownum,l_seq in enumerate(seqf):
            l_seq=regex.sub("\n","",l_seq)
            l_seq_split=l_seq.split("\t")
            if rownum==0:
                fileheader=l_seq_split
            else:
                for colnum,component in enumerate(l_seq_split):
                    if colnum==0:
                        seq[component]={}
                        fastqHeader_now=component
                    seq[fastqHeader_now][fileheader[colnum]]=component
    with gzip.open(file_qual_gz,mode="rt") as qualf:
        for rownum,l_qual in enumerate(qualf):
            l_qual=regex.sub("\n","",l_qual)
            l_qual_split=l_qual.split("\t")
            if rownum==0:
                fileheader=l_qual_split
            else:
                for colnum,component in enumerate(l_qual_split):
                    if colnum==0:
                        qual[component]={}
                        fastqHeader_now=component
                    qual[fastqHeader_now][fileheader[colnum]]=component
    importSeq=dict(seq=seq,qual=qual)
    return importSeq

def seqYielder(file_seq_gz):
    with gzip.open(file_seq_gz,mode="rt") as seqf:
        for i in seqf:
            yield i


def getCorrectedSrc(index,quality,reference=""):
    with gzip.open(index,mode="rt") as f:
        correctedIndexDict=collections.defaultdict(list)
        for rownum,l in enumerate(f):
            l=regex.sub("\n","",l)
            l_split=l.split("\t")
            if rownum==0:
                header=l_split
            else:
                for colnum,component in enumerate(l_split):
                    correctedIndexDict[header[colnum]].append(component)
    with gzip.open(quality,mode="rt") as f:
        qualityDict=collections.defaultdict(list)
        for rownum,l in enumerate(f):
            l=regex.sub("\n","",l)
            l_split=l.split("\t")
            if rownum==0:
                header=l_split
            else:
                for colnum,component in enumerate(l_split):
                    if colnum==0:
                        qualityDict[header[colnum]].append(component)
                    else:
                        qualityDict[header[colnum]].append(int(component))
    if reference:
        with gzip.open(reference,mode="rt") as f:
            references=collections.defaultdict(list)
            for l in f:
                l=regex.sub("\n","",l)
                l_split_tab=l.split("\t")
                l_split_comma=l_split_tab[1].split(",")
                references[l_split_tab[0]]=l_split_comma
    else:
        references=""
    correctedSrc=dict(correctedIndexDict=correctedIndexDict,qualityDict=qualityDict,references=references)
    return correctedSrc


def correctOptionParse(correctOption):
    correctOptDict=dict()
    m=regex.search(r"^[^\(]+\((?P<option>[^\)]*)\)",correctOption)
    correctOpt=m.groupdict()["option"]
    correctOptSplit=correctOpt.split(",")
    if regex.search(r"^correct",correctOption):
        correctOptDict["method"]="correct" 
    elif regex.search(r"^equal",correctOption):
        correctOptDict["method"]="equal"
    elif regex.search(r"^from_bt",correctOption):
        correctOptDict["method"]="from_bt"

    if regex.search(r"^correct|^equal|^from_bt",correctOption):
        for eachOpt in correctOptSplit:
            if "maxDist" in eachOpt.split(":")[0] or "seed_min" in eachOpt.split(":")[0]:
                correctOptDict[eachOpt.split(":")[0]]=int(eachOpt.split(":")[1])    
            elif "analyzed_portion" in eachOpt.split(":")[0]:
                correctOptDict[eachOpt.split(":")[0]]=float(eachOpt.split(":")[1])
            else:
                correctOptDict[eachOpt.split(":")[0]]=eachOpt.split(":")[1]  
        if not correctOptDict.get("analyzed_portion"):
            correctOptDict["analyzed_portion"]=1
        if not correctOptDict.get("rank_threshold"):
            correctOptDict["rank_threshold"]=""
    elif regex.search(r"^combination",correctOption):
        correctOptDict["method"]="combination"
        correctOptDict["combinatorialComponents"]=correctOptSplit    
    return correctOptDict

def convertOptionParse(convertOption):
    convertOptDict=dict()
    m=regex.search(r"^[^\(]+\((?P<option>[^\)]*)\)",convertOption)
    convertOpt=m.groupdict()["option"]
    convertOptSplit=convertOpt.split(",")
    convertOptDict={}
    if regex.search(r"^globalConvert",convertOption):
        convertOptDict["type"]="globalConvert" 
    elif regex.search(r"^localConvert",convertOption):
        convertOptDict["type"]="localConvert"

    for eachOpt in convertOptSplit:
        eachOptSplit=eachOpt.split(":")
        convertOptDict[eachOptSplit[0]]=eachOptSplit[1] 
        # if eachOptSplit[0]=="constraint":
        #     convertOptDict[eachOptSplit[0]]=eachOptSplit[1].split("+")
        if eachOptSplit[0]=="is_combination":
            convertOptDict[eachOptSplit[0]]=strtobool(eachOptSplit[1])
    
    return convertOptDict

def exportOptionParse(exportOption):
    exportOptDict=dict()
    m=regex.search(r"^[^\(]+\((?P<option>[^\)]*)\)",exportOption)
    exportOpt=m.groupdict()["option"]
    exportOptSplit=exportOpt.split(",")
    exportOptDict={}
    if regex.search(r"^equalExport",exportOption):
        exportOptDict["type"]="equalExport" 
    elif regex.search(r"^whitelistExport",exportOption):
        exportOptDict["type"]="whitelistExport"
    elif regex.search(r"^randomExport",exportOption):
        exportOptDict["type"]="randomExport"
    elif regex.search(r"^constantExport",exportOption):
        exportOptDict["type"]="constantExport"

    for eachOpt in exportOptSplit:
        if ("length" in eachOpt.split(":")[0] or "baseQuality" in eachOpt.split(":")[0]) and eachOpt.split(":")[1]:
            exportOptDict[eachOpt.split(":")[0]]=int(eachOpt.split(":")[1])
        elif eachOpt.split(":")[0]=="is_sorted":
            exportOptDict[eachOpt.split(":")[0]]=strtobool(eachOpt.split(":")[1])
        else:
            exportOptDict[eachOpt.split(":")[0]]=eachOpt.split(":")[1]    
    return exportOptDict

