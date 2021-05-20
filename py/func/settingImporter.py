import configparser
import json
from argparse import ArgumentParser
import regex
import gzip
import collections
import csv
from distutils.util import strtobool

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

