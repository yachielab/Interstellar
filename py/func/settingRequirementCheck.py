import glob
import os
import re

class EmptyError(Exception):
    pass

def pathExistCheck(check_list):
    for path in check_list:
        actual_path=glob.glob(path)
        try:
            actual_path[0]
        except IndexError:
            raise EmptyError(path+" is not found.")


# QOPTION should include <mem> and <num_cores> to specify the resource information when submitting jobs
def check_qoption(qoption):
    if (not "<mem>" in qoption) or (not "<num_cores>" in qoption):
        err_msg = "'<mem>' and '<num_cores>' are both required to specify the resource information."
        raise ValueError(err_msg) 


# Filling default values
def tryAndFill(cfg,section,key,defaultvalue,required=False,choice="",parse_home=False):
    try:
        cfg[section][key]
        if cfg[section][key]=="":
            raise EmptyError
    except (KeyError,EmptyError):
        if required:
            err_msg="Parameter '"+key+"' is required."
            raise KeyError(err_msg)
        else:
            return defaultvalue

    if choice:
        choice_list=choice.split(",")
        if not cfg[section][key] in choice_list:
            raise KeyError("Parameter '"+key+"' must be selected from "+choice)
    
    if parse_home:
        cfg[section][key]=os.path.expanduser(cfg[section][key])

    if section=="demultiplex":
        if key=="FORMAT":
            if cfg[section][key]=="FASTQ":
                cfg[section][key]="fastq"
            elif cfg[section][key]=="TSV":
                cfg[section][key]="tsv"
    
    # if section=="qsub":
    #     if key=="QOPTION":
    #         cfg[section][key]=

    return cfg[section][key]

def setMemory(maxmem,minmem,scale_factor):
    maxmem=float(maxmem)
    minmem=float(minmem)
    scale_factor=float(scale_factor)
    mem_now=round(minmem*scale_factor)

    if mem_now > maxmem:
        return str(maxmem)
    else:
        return str(mem_now)

def setDefaultConfig(cfg):
    #Section=general
    cfg["general"]["PROJECT_NAME"]=tryAndFill(cfg,"general","PROJECT_NAME","Interstellar")
    cfg["general"]["PROJECT_DIR"] =tryAndFill(cfg,"general","PROJECT_DIR","",required=True,parse_home=True)
    cfg["general"]["SAMPLESHEET"] =tryAndFill(cfg,"general","SAMPLESHEET","",required=True,parse_home=True)
    cfg["general"]["SET_SHELL_ENV"]=tryAndFill(cfg,"general","SET_SHELL_ENV","",required=True,parse_home=True)
    cfg["general"]["NUM_CORES"]   =tryAndFill(cfg,"general","NUM_CORES","1")
    
    #Section=value_extraction
    if "value_extraction" in cfg:
        cfg["value_extraction"]["READ1_PATH"]=tryAndFill(cfg,"value_extraction","READ1_PATH","",parse_home=True)
        cfg["value_extraction"]["READ2_PATH"]=tryAndFill(cfg,"value_extraction","READ2_PATH","",parse_home=True)
        cfg["value_extraction"]["INDEX1_PATH"]=tryAndFill(cfg,"value_extraction","INDEX1_PATH","",parse_home=True)
        cfg["value_extraction"]["INDEX2_PATH"]=tryAndFill(cfg,"value_extraction","INDEX2_PATH","",parse_home=True)
        cfg["value_extraction"]["FLASH"]=tryAndFill(cfg,"value_extraction","FLASH","",choice="READ1-READ2,READ2-READ1,READ1-INDEX1,INDEX1-READ1,READ1-INDEX2,INDEX2-READ1,READ2-INDEX1,INDEX1-READ2,READ2-INDEX2,INDEX2-READ2,INDEX1-INDEX2,INDEX2-INDEX1")
        cfg["value_extraction"]["FLASH_MIN_OVERLAP"]=tryAndFill(cfg,"value_extraction","FLASH_MIN_OVERLAP","20")
        cfg["value_extraction"]["FLASH_MAX_OVERLAP"]=tryAndFill(cfg,"value_extraction","FLASH_MAX_OVERLAP","30")
        cfg["value_extraction"]["READ1_STRUCTURE"]=tryAndFill(cfg,"value_extraction","READ1_STRUCTURE","")
        cfg["value_extraction"]["READ2_STRUCTURE"]=tryAndFill(cfg,"value_extraction","READ2_STRUCTURE","")
        cfg["value_extraction"]["INDEX1_STRUCTURE"]=tryAndFill(cfg,"value_extraction","INDEX1_STRUCTURE","")
        cfg["value_extraction"]["INDEX2_STRUCTURE"]=tryAndFill(cfg,"value_extraction","INDEX2_STRUCTURE","")

    #Section=value_translation
    if "value_translation" in cfg:
        cfg["value_translation"]["READ1_STRUCTURE"]=tryAndFill(cfg,"value_translation","READ1_STRUCTURE","")
        cfg["value_translation"]["READ2_STRUCTURE"]=tryAndFill(cfg,"value_translation","READ2_STRUCTURE","")
        cfg["value_translation"]["INDEX1_STRUCTURE"]=tryAndFill(cfg,"value_translation","INDEX1_STRUCTURE","")
        cfg["value_translation"]["INDEX2_STRUCTURE"]=tryAndFill(cfg,"value_translation","INDEX2_STRUCTURE","")
    
    #Section=demultiplex
    if "demultiplex" in cfg:
        cfg["demultiplex"]["KEY"]=tryAndFill(cfg,"demultiplex","KEY","")
        cfg["demultiplex"]["FORMAT"]=tryAndFill(cfg,"demultiplex","FORMAT","",required=True,choice="FASTQ,TSV")
        cfg["demultiplex"]["TARGET"]=tryAndFill(cfg,"demultiplex","TARGET","")
        cfg["demultiplex"]["READ1_STRUCTURE"] =tryAndFill(cfg,"demultiplex","READ1_STRUCTURE","")
        cfg["demultiplex"]["READ2_STRUCTURE"] =tryAndFill(cfg,"demultiplex","READ2_STRUCTURE","")
        cfg["demultiplex"]["INDEX1_STRUCTURE"]=tryAndFill(cfg,"demultiplex","INDEX1_STRUCTURE","")
        cfg["demultiplex"]["INDEX2_STRUCTURE"]=tryAndFill(cfg,"demultiplex","INDEX2_STRUCTURE","")
        cfg["demultiplex"]["READ1_HEADER_ADDITION"]=tryAndFill(cfg,"demultiplex","READ1_HEADER_ADDITION","")
        cfg["demultiplex"]["READ2_HEADER_ADDITION"]=tryAndFill(cfg,"demultiplex","READ2_HEADER_ADDITION","")
        cfg["demultiplex"]["INDEX1_HEADER_ADDITION"]=tryAndFill(cfg,"demultiplex","INDEX1_HEADER_ADDITION","")
        cfg["demultiplex"]["INDEX2_HEADER_ADDITION"]=tryAndFill(cfg,"demultiplex","INDEX2_HEADER_ADDITION","")
    
    #Section=annotate_header
    if "annotate_header" in cfg:
        cfg["annotate_header"]["READ1_HEADER_ADDITION"]=tryAndFill(cfg,"annotate_header","READ1_HEADER_ADDITION","")
        cfg["annotate_header"]["READ2_HEADER_ADDITION"]=tryAndFill(cfg,"annotate_header","READ2_HEADER_ADDITION","")
        cfg["annotate_header"]["INDEX1_HEADER_ADDITION"]=tryAndFill(cfg,"annotate_header","INDEX1_HEADER_ADDITION","")
        cfg["annotate_header"]["INDEX2_HEADER_ADDITION"]=tryAndFill(cfg,"annotate_header","INDEX2_HEADER_ADDITION","")
        cfg["annotate_header"]["READ1_STRUCTURE"] =tryAndFill(cfg,"annotate_header","READ1_STRUCTURE","")
        cfg["annotate_header"]["READ2_STRUCTURE"] =tryAndFill(cfg,"annotate_header","READ2_STRUCTURE","")
        cfg["annotate_header"]["INDEX1_STRUCTURE"]=tryAndFill(cfg,"annotate_header","INDEX1_STRUCTURE","")
        cfg["annotate_header"]["INDEX2_STRUCTURE"]=tryAndFill(cfg,"annotate_header","INDEX2_STRUCTURE","")
    
    return cfg


def setDefaultQConfig(cfg):
    cfg["qsub"]["MEM_MAX"]=tryAndFill(cfg,"qsub","MEM_MAX","","128")
    cfg["qsub"]["MEM_MIN"] =tryAndFill(cfg,"qsub","MEM_MIN","","6")
    cfg["qsub"]["NUM_READS"]=tryAndFill(cfg,"qsub","NUM_READS","2000000")
    cfg["qsub"]["QOPTION"]=tryAndFill(cfg,"qsub","QOPTION","",required=True)
    cfg["qsub"]["mem_import"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],1)
    cfg["qsub"]["mem_qc"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],1)
    cfg["qsub"]["mem_to_bt"]=cfg["qsub"]["MEM_MAX"]
    cfg["qsub"]["mem_correct"]=cfg["qsub"]["MEM_MAX"]
    cfg["qsub"]["mem_mk_sval"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],2)
    cfg["qsub"]["mem_buildTree"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],1)
    cfg["qsub"]["mem_mergeTree"]=cfg["qsub"]["MEM_MAX"]
    cfg["qsub"]["mem_convert"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],3)
    cfg["qsub"]["mem_bc_sort"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],1)
    cfg["qsub"]["mem_export"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],3)
    cfg["qsub"]["mem_demultiplex"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],2)
    cfg["qsub"]["mem_annotate_header"]=setMemory(cfg["qsub"]["MEM_MAX"],cfg["qsub"]["MEM_MIN"],2)
    return cfg["qsub"]


# def setDefaultFunc_ext(func_dict):
#     for val in func_dict:
#         for func in func_dict[val]["func_oredered"]:
#             func_options=func_dict[val][func]
#             if func=="QUALITY_FILTER":


def setDefaultValueInConfig(cmd,config,distribute=False):
    if cmd=="import":
        required_parameters=["src_raw_components"]
        for i in required_parameters:
            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                err_msg="The parameter "+i+" is not found in the config file."
                raise KeyError(err_msg)
        flash_options=["flash_min_overlap","flash_max_overlap"]
        for i in flash_options:
            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                if i=="flash_min_overlap":
                    config[i]=10
                elif i=="flash_max_overlap":
                    config[i]=65
        return config
    
    elif cmd=="streamline":
        required_parameters=["pipeline","qoption","template_shellscript"]
        for i in required_parameters:
            if i=="qoption" and not distribute:
                continue
            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                err_msg="The parameter "+i+" is not found in the config file."
                raise KeyError(err_msg)
        config["template_shellscript"]=os.path.expanduser(config["template_shellscript"])
        
        for i in ["mem_max","chunksize"]:
            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                if i=="mem_max":
                    config["mem_max"]='128'
                if i=="chunksize":
                    config["chunksize"]='2000000'
        
        for i in config["pipeline"].split(","):
            try:
                    config["mem_"+i]
                    if config["mem_"+i]=="":
                        raise EmptyError
            except (KeyError,EmptyError):
                if i in ["import","qc","to_bt","bc_sort"]:
                    config["mem_"+i]='8'
                elif i in ["mk_sval","buildTree","convert","export","annotate_header","demultiplex"]:
                    config["mem_"+i]='16'
                elif i in ["mergeTree","correct"]:
                    config["mem_"+i]=config["mem_max"]

        return config