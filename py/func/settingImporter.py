import configparser
from dataclasses import replace
import json
from argparse import ArgumentError, ArgumentParser
from random import sample
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

from . import settingRequirementCheck

class InputError(Exception):
    pass



def readconfig(input_cfg):
    cfg=configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(input_cfg)

    return cfg


def configClean(cfgDict,qconf=True,query=""):
    if qconf:
        new_cfgDict={}
        for key in cfgDict:
            if key=="QOPTION":
                opt_now=cfgDict[key]
                opt_now=regex.sub(r"^\"|\"$","",opt_now)
                settingRequirementCheck.check_qoption(opt_now)
                new_cfgDict[key]=opt_now
            else:
                new_cfgDict[key]=regex.sub("\\t|\\n| ","",cfgDict[key])
        return new_cfgDict
    else:
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
        samplesheet=samplesheet.dropna(how="all")
        samples_uniq=set(samplesheet[1])
        fileprefix_uniq=set(samplesheet[0])

        # Rule 1: A same prefix can't be used for different samples
        sample_to_fileprefix=collections.defaultdict(list)
        for s,pfx in zip(samplesheet[1],samplesheet[0]):
            sample_to_fileprefix[s].append(pfx)
        for s1 in sample_to_fileprefix:
            for s2 in sample_to_fileprefix:
                if s1 != s2:
                    p1 = set(sample_to_fileprefix[s1])
                    p2 = set(sample_to_fileprefix[s2])
                    if len(p1 & p2) > 0:
                        errmsg = "A same prefix can't be used for different samples! This error was caused by the samples "+s1+" and "+s2+".  Check the sample sheet."
                        raise InputError(errmsg)

        # Rule 2: File prefixes should be completely unique each other
        for pfx1 in fileprefix_uniq:
            for pfx2 in fileprefix_uniq:
                if pfx1 != pfx2:
                    m1 = re.search(pfx1,pfx2)
                    m2 = re.search(pfx2,pfx1)
                    if m1 or m2:
                        errmsg = "File prefixes should be completely unique each other! This error was caused by the file prefixes "+pfx1+" and "+pfx2+" Check the sample sheet."
                        raise InputError(errmsg)

        for i in samples_uniq:
            os.makedirs(proj_dir+"/"+i,exist_ok=True)
        sampledir_list=[proj_dir+"/"+i for i in samples_uniq]
        os.makedirs(proj_dir+"/_multisample",exist_ok=True)
        os.makedirs(proj_dir+"/_multisample/mergeTree",exist_ok=True)
        os.makedirs(proj_dir+"/_multisample/samplesheet",exist_ok=True)

    return sampledir_list


def concat_extended_process(cfg_value_ext,untreated_segment_pool):
    additional_value_dict={}
    parents=set()
    process_inherit_dict=collections.defaultdict(dict)
    for seg in untreated_segment_pool:
        line_now=cfg_value_ext[seg]
        func_blocks=line_now.split(">>")
        first_block=func_blocks[0]
        fisrt_block_contents=re.search(r"[^\(]+\(([^\)]*)\)",first_block)
        fisrt_block_contents=fisrt_block_contents.group(1)
        for i in fisrt_block_contents.split(","):
            key=i.split(":")[0]
            val=i.split(":")[1]
            if key=="source":
                process_inherit_dict[seg]["parent_segment"]=val
                process_inherit_dict[seg]["func_line"]=re.sub(r"source:"+val+",|,source:"+val+"|source:"+val,"",line_now)
                parents.add(val)
    
    # ancestors=parents-set(untreated_segment_pool)
    ancestors=cfg_value_ext["segments"]
    dict_to_terminal={}
    # print(ancestors)
    for seg_anc in ancestors:
        function_list=[]
        target_anc=copy.copy(seg_anc)
        segment_lineage=[target_anc]
        while True:
            terminal=True
            for seg in process_inherit_dict:
                if target_anc == process_inherit_dict[seg]["parent_segment"]:
                    terminal=False
                    function_list.append(process_inherit_dict[seg]["func_line"])
                    target_anc=copy.copy(seg)
                    segment_lineage.append(seg)
                    break
            if not function_list: #In case no processes happened after regex segmentation for the segment
                idx=0
                while True:
                    if not seg_anc+".clean"+str(idx) in (parents | set(untreated_segment_pool)):
                        new_key=seg_anc+".clean"+str(idx)
                        break
                    idx+=1

                segment_lineage.append(new_key)
                target_anc=new_key
                if seg_anc in set(list(cfg_value_ext["parent"].keys())+list(cfg_value_ext["parent"].values())):
                    additional_value_dict[new_key]="SEQ2VALUE(source:"+seg_anc+")"
                else:
                    additional_value_dict[new_key]="PASS(source:"+seg_anc+")"
                break
            else: #In case some processes done for the segment
                if terminal:
                    first_block_split=function_list[0].split(">>")
                    if re.search("\(\)",first_block_split[0]):
                        fisrt_block_contents=re.sub("\(","(source:"+seg_anc,first_block_split[0])
                    else:
                        fisrt_block_contents=re.sub("\(","(source:"+seg_anc+",",first_block_split[0])
                    first_block_split[0]=fisrt_block_contents
                    function_list[0]=">>".join(first_block_split)
                    process_streamline=">>".join(function_list)
                    if seg_anc in set(list(cfg_value_ext["parent"].keys())+list(cfg_value_ext["parent"].values())):
                        if not "SEQ2VALUE" in process_streamline:
                            process_streamline+=">>SEQ2VALUE()"
                    additional_value_dict[target_anc]=process_streamline
                    break
        for i in segment_lineage:
            dict_to_terminal[i]=target_anc
    # print(additional_value_dict)
    
    return additional_value_dict,dict_to_terminal


def config_extract_value_ext(cfg_raw):
    cfg=copy.deepcopy(cfg_raw)
    cfg_value_ext=cfg["value_extraction"]
    #sequence segment label extractioon
    segments=[]
    for key in cfg_value_ext:
        if "READ1_STRUCTURE" in key or "READ2_STRUCTURE" in key or "INDEX1_STRUCTURE" in key or "INDEX2_STRUCTURE" in key or "READ_FLASH_STRUCTURE" in key:
            segments+=re.findall(r"\?P\<([^\>]+)\>",cfg_value_ext[key])
    cfg_value_ext["segments"]=sorted(list(set(segments)))

    #parental-local allocation
    cfg_value_ext["parent"]={}
    for child in cfg_value_ext:
        if not type(cfg_value_ext[child])==str:
            continue
        if re.search(r"CHILD_OF\([^\)]+\)",cfg_value_ext[child]):
            parent=cfg_value_ext[child].replace("CHILD_OF(","").replace(")","")
            cfg_value_ext["parent"][child]=parent
    # for seg in cfg_value_ext["segments"]:
    #     if cfg_value_ext.get(seg):
    #         parents_list=cfg_value_ext[seg].replace("CHILD_OF(","").replace(")","").split(",")
    #         cfg_value_ext["parent"][seg]=parents_list
    
    #value segment label extraction
    all_keys =set(cfg_value_ext.keys())
    keys_tmp =set([i for i in cfg_value_ext if "READ1_STRUCTURE" in i or "READ2_STRUCTURE" in i or "INDEX1_STRUCTURE" in i or "INDEX2_STRUCTURE" in i or "READ_FLASH" in i])
    keys_tmp|=set(["READ1_PATH","READ2_PATH","INDEX1_PATH","INDEX2_PATH"])
    keys_tmp|=set(["FLASH","FLASH_MIN_OVERLAP","FLASH_MAX_OVERLAP"])
    keys_tmp|=set(cfg_value_ext["segments"]+["segments","parent"]+list(cfg_value_ext["parent"].keys()))

    untreated_segment_pool=list(all_keys-keys_tmp)

    process_dict,dict_to_terminal=concat_extended_process(cfg_value_ext,untreated_segment_pool)

    cfg_value_ext["value_segment"]=sorted(list(process_dict.keys()))
    for k in process_dict:
        cfg_value_ext[k]=process_dict[k]

    #detecting functions
    functions_used=[]
    for key in cfg_value_ext["value_segment"]:
        commandline=cfg_value_ext[key].split(">>")
        functions_used+=[re.sub(r"\(.*\)","",i) for i in commandline]
    cfg_value_ext["functions_used"]=set(functions_used)

    return cfg_value_ext,dict_to_terminal



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

        if funcname=="M2A_CORRECTION":
            func_dict[val][funcname]["path"]=os.path.expanduser(func_dict[val][funcname]["path"])
        elif funcname=="SEQ2SEQ":
            func_dict[val][funcname]["conversion_table"]=os.path.expanduser(func_dict[val][funcname]["conversion_table"])
        elif funcname=="WHITELIST_ASSIGNMENT":
            if "allowlist_path" in func_dict[val][funcname]:
                func_dict[val][funcname]["allowlist_path"]=list(map(os.path.expanduser,func_dict[val][funcname]["allowlist_path"]))
            # func_dict[val][funcname]["conversion_table"]=os.path.expanduser(func_dict[val][funcname]["conversion_table"])

        for i in opt_int:
            if type(func_dict[val][funcname][i])==list:
                func_dict[val][funcname][i]=list(map(int,func_dict[val][funcname][i]))
            else:
                func_dict[val][funcname][i]=int(func_dict[val][funcname][i])
    return func_dict


def func_parse(func_line,func_collection,dest=False,dict_to_terminal=""):
    func_line=func_line.split(">>")
    d=dict()
    funcname_list=[]
    for func in func_line:
        m=re.search(r"([^\(]+)\(([^\)]*)\)",func)
        if m:
            funcname=m.group(1)
            if not funcname in func_collection:
                raise KeyError("Function '"+funcname+"' doesn't exist.")
            contents=m.group(2).split(",")
        elif re.search(r"^\".+\"|\'.+\'$", func):
            m2=re.search(r"^\"(.+)\"|\'(.+)\'$", func)
            seq=m2.group(1)
            funcname="CONSTANT"
            contents=["sequence:"+seq]
        else:
            funcname="PASS"
            contents=["source:"+func]

        #Packing func options into a dict
        if not dest:
            contents={i.split(":")[0]:i.split(":")[1] for i in contents if not i==""}
            d[funcname]=contents
        else:                
            opt=""
            dict_out=dict()
            for i in contents:
                if opt=="source" and len(i.split(":"))==1:
                    dict_out["source"].append(dict_to_terminal[i]) #translate the given source segment name into terminal source name

                elif (opt=="allowlist_path" or opt=="randseq_pattern") and len(i.split(":"))==1:
                    dict_out[opt].append(i)

                else:
                    opt=i.split(":")[0]
                    val=i.split(":")[1] 
                    if opt=="source":
                        dict_out[opt]=[dict_to_terminal[val]] #translate the given source segment name into terminal source name
                    elif opt=="allowlist_path" or opt=="randseq_pattern":
                        dict_out[opt]=[val]
                    else:
                        dict_out[opt]=val

                if "source" in dict_out and str(dict_out["source"])==list:
                    dict_out["source"]=sorted(dict_out["source"])

            # Convert VALUE2SEQ into WHITELIST or RANDSEQ
            if funcname=="VALUE2SEQ":
                # if dict_out.get("allowlist_path") or dict_out.get("conversion_table"):
                if dict_out.get("allowlist_path"):
                    funcname="WHITELIST_ASSIGNMENT"
                elif dict_out.get("randseq_pattern"):
                    funcname="RANDSEQ_ASSIGNMENT"
                else:
                    # raise KeyError("Whitelist path, sequence conversion table or random sequence length are required.")
                    raise KeyError("Allowlist path or random sequence length are required.")
            d[funcname]=dict_out
        funcname_list.append(funcname)
    d["func_ordered"]=funcname_list

    return d



def parse_function_value_ext(cfg):
    func_collection=["QUALITY_FILTER","PASS","I2M_CORRECTION","M2A_CORRECTION","BARTENDER_CORRECTION","CUSTOM_CORRECTION","SEQ2VALUE"]
    value_names=cfg["value_segment"]
    func_dict=dict()
    for val in value_names:
        func_line=cfg[val]
        func_dict[val]=func_parse(func_line,func_collection)
    return func_dict


#Around here, treatment of adding "PASS" func / VALUE2SEQ classification into WHITE / RAND
def parse_function_value_trans(cfg,dict_to_terminal):
    func_collection=["SEQ2SEQ","VALUE2SEQ","PASS","CONSTANT"]
    destseg_names=cfg["dest_segment"]
    func_dict=dict()
    for seg in destseg_names:
        func_line=cfg[seg]
        func_dict[seg]=func_parse(func_line,func_collection,dest=True,dict_to_terminal=dict_to_terminal)
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
            raise ArgumentError("Function PASS cannot be performed with other functions.")      
        if ("BARTENDER_CORRECTION" in func_oredered) and ("I2M_CORRECTION" in func_oredered or "M2A_CORRECTION" in func_oredered or "CUSTOM_CORRECTION" in func_oredered):
            raise ArgumentError("Function BARTENDER_CORRECTION cannot be performed with other correction functions.")
        if ("CUSTOM_CORRECTION" in func_oredered) and ("I2M_CORRECTION" in func_oredered or "M2A_CORRECTION" in func_oredered or "BARTENDER_CORRECTION" in func_oredered):
            raise ArgumentError("Function CUSTOM_CORRECTION cannot be performed with other correction functions.")

        #QUALITY_FILTER Correct option check
        funcIllegalOptionCheck(val,func_dict,"QUALITY_FILTER",["source","min_nucleotide_Q-score","min_avg_Q-score"])
        func_dict=funcSetDefault(val,func_dict,"QUALITY_FILTER",["min_nucleotide_Q-score","min_avg_Q-score"],[5,20],opt_int=["min_nucleotide_Q-score","min_avg_Q-score"])

        #I2M CORRECTION option check
        funcIllegalOptionCheck(val,func_dict,"I2M_CORRECTION",["source","rank","levenshtein_distance","min_num_reads"])
        func_dict=funcSetDefault(val,func_dict,"I2M_CORRECTION",["rank","levenshtein_distance","min_num_reads"],["auto",1,1],opt_int=["levenshtein_distance","min_num_reads"])

        #M2A CORRECTION option check
        funcIllegalOptionCheck(val,func_dict,"M2A_CORRECTION",["source","path","levenshtein_distance","min_num_reads"])
        funcRequiredOptionCheck(val,func_dict,"M2A_CORRECTION",["path"])
        func_dict=funcSetDefault(val,func_dict,"M2A_CORRECTION",["levenshtein_distance","min_num_reads"],[1,1],opt_int=["levenshtein_distance","min_num_reads"])

        #Bartender option check
        funcIllegalOptionCheck(val,func_dict,"BARTENDER_CORRECTION",["source","-c","-d","-l","-z","-s"])

        #CUSTOM_CORRECTION option check
        funcIllegalOptionCheck(val,func_dict,"CUSTOM_CORRECTION",["source","shell_script"])
        funcRequiredOptionCheck(val,func_dict,"CUSTOM_CORRECTION",["shell_script"])

        #SEQ2VALUE option check
        funcIllegalOptionCheck(val,func_dict,"SEQ2VALUE",["source"])

        #PASS Correct option check
        funcIllegalOptionCheck(val,func_dict,"PASS",["source"])
        funcRequiredOptionCheck(val,func_dict,"PASS",["source"])
        func_dict=funcSetDefault(val,func_dict,"PASS",["length"],[0],opt_int=["length"])

        #Values
        if ("SEQ2VALUE" in func_oredered) or ("I2M_CORRECTION" in func_oredered) or ("M2A_CORRECTION" in func_oredered) or ("BARTENDER_CORRECTION" in func_oredered) or ("CUSTOM_CORRECTION" in func_oredered):
            barcode_list.append(segment_now[0])
        for func in func_dict[val]:
            if not func=="func_ordered":
                func_dict[val][func]["source"]=segment_now[0]
    func_dict["barcode_list"]=barcode_list
    return func_dict



def config_extract_value_trans(cfg_raw):
    cfg=copy.deepcopy(cfg_raw)
    cfg_value_ext,dict_to_terminal=config_extract_value_ext(cfg)
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
    # cfg_value_trans["bc_sort"]=False
    # cfg_value_trans["val2table"]={}
    # func_dict=parse_function_value_trans(cfg_value_trans,dict_to_terminal)
    # for dest_seg in func_dict:
    #     for fun in func_dict[dest_seg]["func_ordered"]:
    #         if "conversion_table" in func_dict[dest_seg][fun]:
    #             if not func_dict[dest_seg][fun]["conversion_table"]=="":
    #                 cfg_value_trans["bc_sort"]=True
    #                 cfg_value_trans["val2table"][dest_seg]=func_dict[dest_seg][fun]["conversion_table"]
    
    #available segments
    cfg_value_trans["available_seg"]=[]
    for L in cfg_value_trans["dest_segment"]:
        L=L.split(",")
        cfg_value_trans["available_seg"]+=L

    cfg_value_trans["available_seg"] += cfg_value_ext["value_segment"]
    cfg_value_trans["available_seg"] += cfg_value_ext["segments"]
   
   
    return cfg_value_trans


def func_check_trans(cfg,dict_to_terminal):
    func_dict=parse_function_value_trans(cfg,dict_to_terminal)
    for val in func_dict:
        # segment_now=[]
        # func_oredered=func_dict[val]["func_ordered"]
        # for func in func_oredered:
        #     content_dict=func_dict[val][func]
        #     if "source" in content_dict and not content_dict["source"]=="":
        #         segment_now+=content_dict["source"]
        
        #WHITELIST_ASSIGNMENT Correct option check
        # funcIllegalOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source","allowlist_path","conversion_table"])
        funcIllegalOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source","allowlist_path"])
        func_dict=funcSetDefault(val,func_dict,"WHITELIST_ASSIGNMENT",[""],[""])
        # if "WHITELIST_ASSIGNMENT" in func_dict[val] and func_dict[val]["WHITELIST_ASSIGNMENT"]["conversion_table"]=="":
        #     funcRequiredOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source","allowlist_path"])
        # else:
        #     funcRequiredOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source"])
        funcRequiredOptionCheck(val,func_dict,"WHITELIST_ASSIGNMENT",["source","allowlist_path"])
        
        #SEQ2SEQ Correct option check
        funcIllegalOptionCheck(val,func_dict,"SEQ2SEQ",["source","conversion_table"])
        funcRequiredOptionCheck(val,func_dict,"SEQ2SEQ",["source","conversion_table"])

        #RANDSEQ_ASSIGNMENT Correct option check
        funcIllegalOptionCheck(val,func_dict,"RANDSEQ_ASSIGNMENT",["source","randseq_pattern"])
        funcRequiredOptionCheck(val,func_dict,"RANDSEQ_ASSIGNMENT",["source","randseq_pattern"])
        func_dict=funcSetDefault(val,func_dict,"RANDSEQ_ASSIGNMENT",[],[])

        #CONSTANT Correct option check
        funcIllegalOptionCheck(val,func_dict,"CONSTANT",["sequence"])
        funcRequiredOptionCheck(val,func_dict,"CONSTANT",["sequence"])
        func_dict=funcSetDefault(val,func_dict,"CONSTANT",["source"],[["constant"]])

        #PASS Correct option check
        funcIllegalOptionCheck(val,func_dict,"PASS",["source"])
        funcRequiredOptionCheck(val,func_dict,"PASS",["source"])
        func_dict=funcSetDefault(val,func_dict,"PASS",["length","add"],[0,"A"],opt_int=["length"])

    return func_dict



def config_extract_value_demulti(cfg_raw):
    cfg=copy.deepcopy(cfg_raw)
    cfg_value_ext,dict_to_terminal=config_extract_value_ext(cfg)
    cfg_value_demulti=cfg["demultiplex"]
      
    #available segments
    cfg_value_demulti["available_seg"]=[]
    for i in cfg_value_ext["value_segment"]:
        cfg_value_demulti["available_seg"].append(i)
    cfg_value_demulti["available_seg"] += cfg_value_ext["segments"]
    
    return cfg_value_demulti


def config_extract_value_annotate_header(cfg_raw):
    cfg=copy.deepcopy(cfg_raw)
    cfg_value_ext,dict_to_terminal=config_extract_value_ext(cfg)
    cfg_value_annotate_header=cfg["annotate_header"]
       
    #available segments
    cfg_value_annotate_header["available_seg"]=[]
    for i in cfg_value_ext["value_segment"]:
        cfg_value_annotate_header["available_seg"].append(i)
    cfg_value_annotate_header["available_seg"] += cfg_value_ext["segments"]
   
    return cfg_value_annotate_header


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
            qscore_dict[segment]={"min_base":func_dict[val]["QUALITY_FILTER"]["min_nucleotide_Q-score"],"min_avg":func_dict[val]["QUALITY_FILTER"]["min_avg_Q-score"]}
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


def parseInputFileList(in_path):
    L = []
    with open(in_path,mode="rt") as f:
        for x in f:
            L.append(x.replace("\n","")) 
    return L