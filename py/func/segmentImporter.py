import regex
import gzip
import random
import string
import subprocess
import collections
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from itertools import zip_longest



def merge_reads_flash2(readPathDict,flash,gzipped,tmpdir,config,n_core):
    if not flash:
        return [],readPathDict
    else:
        flash=flash.replace("READ","Read")
        flash=flash.replace("INDEX","Index")
        merge_reads=[i+"_src" for i in flash.split("-")]
        if gzipped:
            print("gzipped read decompressing...",flush=True)
            decompressed_fastq=[]
            for i in merge_reads:
                cmdlist_gunzip=["gunzip","-c",readPathDict[i],">",tmpdir+"/"+i+".fastq"]
                cmd=" ".join(cmdlist_gunzip)
                subprocess.run(cmd,shell=True)
                decompressed_fastq.append(tmpdir+"/"+i+".fastq")
        else:
            decompressed_fastq=[readPathDict[i] for i in merge_reads]

        flash_min=str(config["FLASH_MIN_OVERLAP"])
        flash_max=str(config["FLASH_MAX_OVERLAP"])
        print("Merging reads using Flash2...",flush=True)
        cmdlist_flash=["flash2","-z","-t",str(n_core),"-m",flash_min,"-M",flash_max,"-d",tmpdir,"-o","merge",decompressed_fastq[0],decompressed_fastq[1]]
        cmd=" ".join(cmdlist_flash)
        subprocess.run(cmd,shell=True)
        print("Merging done!",flush=True)

        for cnt,i in enumerate(merge_reads):
            if cnt==0:
                readPathDict[i]=tmpdir+"/"+"merge.notCombined_1.fastq.gz"
            else:
                readPathDict[i]=tmpdir+"/"+"merge.notCombined_2.fastq.gz"
        readPathDict["merge_src"]=tmpdir+"/"+"merge.extendedFrags.fastq.gz"
        gzipped=merge_reads+["merge_src"]
        return gzipped,readPathDict


def sequenceGenerator(fq_path,settings,from_flash=False):
    if settings.input_fastq_gzipped or from_flash:
        with gzip.open(fq_path,mode="rt") as input_file:
            for i in input_file:
                i=i.replace("\n","")
                yield i
    else:
        with open(fq_path,mode="rt") as input_file:
            for i in input_file:
                i=i.replace("\n","")
                yield i


def splitSequenceGenerator(fq_path,settings,chunksize,from_flash=False):  # Generator function.
    FILLER = object()  # Unique object
    
    if settings.input_fastq_gzipped or from_flash:
        with gzip.open(fq_path,mode="rt") as input_file:
            lines = (line.rstrip() for line in input_file)
            for group in zip_longest(*([iter(lines)]*chunksize), fillvalue=FILLER):
                limit = group.index(FILLER) if group[-1] is FILLER else len(group)
                yield group[:limit]  # Sliced to remove any filler.
    else:
        with open(fq_path,mode="rt") as input_file:
            lines = (line.rstrip() for line in input_file)
            for group in zip_longest(*([iter(lines)]*chunksize), fillvalue=FILLER):
                limit = group.index(FILLER) if group[-1] is FILLER else len(group)
                yield group[:limit]  # Sliced to remove any filler.



def randomname(n):
   return ''.join(random.choices(string.ascii_letters + string.digits, k=n))


def parseSegmentFromRegex(regexList):
    parse_res=set()
    for pat_now in regexList:
        pat_parsed_raw=regex.findall(r"\?P\<[^\>]+\>",pat_now)
        pat_parsed=[regex.sub(r"\?P\<|\>","",x) for x in pat_parsed_raw]
        parse_res|=set(pat_parsed)
    parse_res=list(parse_res)
    return parse_res


def patMatching(readNow,patternList):
    for patNow in patternList:
        mtch=patNow.search(readNow)
        if mtch:
            return mtch
    return None


def qualityFiltering(qualNow,min_base_quality,min_avg_quality):
    qual_list=[ord(i)-33 for i in qualNow]
    qual_min_now=min(qual_list)
    qual_avg_now=sum(qual_list)/len(qual_list)

    if qual_min_now < min_base_quality or qual_avg_now < min_avg_quality:
        return True
    else:
        return False


def qualityFilteringForDataFrame(df_set,qc_targets,qscore_dict):
    seq_chunk = df_set[0]
    qual_chunk = df_set[1]

    for seg in qc_targets:
        bool_filtered=qual_chunk[seg].map(lambda x: qualityFiltering(x,min_base_quality=qscore_dict[seg]["min_base"],min_avg_quality=qscore_dict[seg]["min_avg"]))
        seq_chunk[seg][bool_filtered]="-"
    
    return seq_chunk
    

def split_fastq_per_lines(lines, n=4):
    for idx in range(0, len(lines), n):
        yield lines[idx:idx + n]
    

def fastq_segmentation(subchunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set):
    parsedSeqDict_tmp=collections.defaultdict(list)
    parsedQualDict_tmp=collections.defaultdict(list)
    for lines in split_fastq_per_lines(subchunk):
        line_header = lines[0]
        line_sequence = lines[1]
        line_quality = lines[3]
        
        # Header treatment
        header=headerSplitRegex.split(line_header)[0]
        parsedSeqDict_tmp["Header"].append(header)
        parsedQualDict_tmp["Header"].append(header)
        
        # Sequence segmentation
        m=patMatching(line_sequence,settings.regexDictCompiled[readKey])
        if m:
            mdict=m.groupdict()
            for seg in mdict:
                parsedSeqDict_tmp[seg].append(mdict[seg])

            component_diff=segment_parsed_set-set(mdict.keys())
            for seg in component_diff:
                parsedSeqDict_tmp[seg].append("-")
                parsedQualDict_tmp[seg].append("-")
        else: # If the regex doesn't match, quality segments should also be all missing
            for seg in segment_parsed:
                parsedSeqDict_tmp[seg].append("-")
                parsedQualDict_tmp[seg].append("-")

        # Quality segmentation
        if m:
            for component in mdict:
                extractedQual=line_quality[m.span(component)[0]:m.span(component)[1]]
                parsedQualDict_tmp[component].append(extractedQual)
    
    return parsedSeqDict_tmp,parsedQualDict_tmp


# Multithreading implementation for read segmentation
def segmentation_parallel_wrapper(fastq_chunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set,ncore):
    out = Parallel(n_jobs=ncore,verbose=3)(
        delayed(fastq_segmentation)(subchunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set) 
            for subchunk in split_fastq_per_lines(fastq_chunk,n=160000)) # Each CPU processes 40000 FASTQ records
    for n,dict_tuple in enumerate(out):
        if n == 0:
            parsedSeqDict = dict_tuple[0]
            parsedQualDict = dict_tuple[1]
        else:
            for seg in parsedSeqDict:
                parsedSeqDict[seg] += dict_tuple[0][seg]
                parsedQualDict[seg] += dict_tuple[1][seg]
    return parsedSeqDict,parsedQualDict


# Multithreading implementation for quality filtering
def qfilter_parallel_wrapper(seq_chunk, qual_chunk, qc_targets, qscore_dict, ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    seq_subchunks = np.array_split(seq_chunk,ncore)
    qual_subchunks = np.array_split(qual_chunk,ncore)

    retLst = Parallel(n_jobs=ncore,verbose=3)(delayed(qualityFilteringForDataFrame)(df_set, qc_targets, qscore_dict) for df_set in zip(seq_subchunks,qual_subchunks))
    return pd.concat(retLst)