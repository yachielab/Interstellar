import sys
import regex
import gzip
import random
import string
import subprocess
import collections
import math
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from itertools import zip_longest
import time
import pickle
import glob

# New function of FLASH
def run_flash2(path_read_1, path_read_2, flash_min, flash_max, outdir, n_core):
    print("Merging reads using Flash2...",flush=True)
    cmdlist_flash=["flash2", "-z", "-t", str(n_core), "-m", flash_min, "-M", flash_max, "-d", outdir, "-o", "merge", path_read_1, path_read_2]
    cmd=" ".join(cmdlist_flash)
    subprocess.run(cmd,shell=True)
    print("Merging done!",flush=True)

    filepath_out = [outdir+"/merge.notCombined_1.fastq.gz", outdir+"/merge.notCombined_2.fastq.gz", outdir+"/merge.extendedFrags.fastq.gz"]
    return filepath_out


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
        print(seg, flush=True)
        print(qual_chunk[seg].head(), flush=True)
        bool_filtered=qual_chunk[seg].map(lambda x: qualityFiltering(x,min_base_quality=qscore_dict[seg]["min_base"],min_avg_quality=qscore_dict[seg]["min_avg"]))
        seq_chunk[seg][bool_filtered]="-"
    
    return seq_chunk
    

def split_yield_fastq_per_lines(lines, n=4):
    for idx in range(0, len(lines), n):
        yield lines[idx:idx + n]


# def split_fastq_per_lines(lines, n=4):
#     L = [lines[idx:idx + n] for idx in range(0, len(lines), n)]
#     return L
    

def fastq_segmentation(cpu_idx,subchunk,headerSplitRegex,segment_parsed,regex_pattern_now,outdir,readKey,barcodes):
    t0 = time.time()
    parsedSeqDict_tmp={k:["-"]*int(len(subchunk)/4) for k in (["Header"]+segment_parsed)}
    parsedQualDict_tmp={k:["-"]*int(len(subchunk)/4)for k in (["Header"]+segment_parsed)}
    
    # print("Formatting",time.time()-t0,"sec", flush=True)

    t_getheader=0
    t_regex = 0
    t_packdic = 0

    for idx,lines in enumerate(split_yield_fastq_per_lines(subchunk,n=4)):
        # line_header = lines[0]
        # line_sequence = lines[1]
        # line_quality = lines[3]
        
        # Header treatment
        T=time.time()
        header=headerSplitRegex.split(lines[0])[0]
        parsedSeqDict_tmp["Header"][idx] = header
        parsedQualDict_tmp["Header"][idx] = header
        t_getheader += time.time() - T
        
        # Sequence segmentation
        T = time.time()
        m=patMatching(lines[1],regex_pattern_now)
        t_regex += time.time() - T

        T = time.time()
        if m:
            mdict=m.groupdict()
            for seg in mdict:
                parsedSeqDict_tmp[seg][idx] = mdict[seg]
                # extractedQual=lines[3][m.span(component)[0]:m.span(component)[1]]
                parsedQualDict_tmp[seg][idx] = lines[3][m.span(seg)[0] : m.span(seg)[1]]
        t_packdic += time.time() - T

            # component_diff=segment_parsed_set-set(mdict.keys())
            # for seg in component_diff:
            #     parsedSeqDict_tmp[seg].append("-")
            #     parsedQualDict_tmp[seg].append("-")
        # else: # If the regex doesn't match, quality segments should also be all missing
        #     for seg in segment_parsed:
        #         parsedSeqDict_tmp[seg].append("-")
        #         parsedQualDict_tmp[seg].append("-")

        # Quality segmentation
        # if m:
        #     for component in mdict:
        #         extractedQual=lines[3][m.span(component)[0]:m.span(component)[1]]
        #         parsedQualDict_tmp[component].append(extractedQual)
    
    with open("_".join([outdir+"/CPU"+str(cpu_idx),readKey,"srcSeq.pkl"]),mode="wb") as p:
        pickle.dump(parsedSeqDict_tmp,p)
    with open("_".join([outdir+"/CPU"+str(cpu_idx),readKey,"srcQual.pkl"]),mode="wb") as p:
        pickle.dump(parsedQualDict_tmp,p)
    
    counterDict_tmp = {}
    for bc in barcodes:
        if bc in parsedSeqDict_tmp:
            counterDict_tmp[bc] = collections.Counter(parsedSeqDict_tmp[bc])

    # print("Get header",t_getheader,"sec", flush=True)
    # print("Regex search",t_regex,"sec", flush=True)
    # print("Packing dict",t_packdic,"sec", flush=True)
    # print("Job time",time.time()-t0,"sec",flush=True)
    return counterDict_tmp


# def fastq_segmentation2(record,settings,headerSplitRegex,readKey,segment_parsed):
#     # parsedSeqDict_tmp=collections.defaultdict(list)
#     # parsedQualDict_tmp=collections.defaultdict(list)
    
#     # Header treatment
#     out_seq = [headerSplitRegex.split(record[0])[0]]
#     out_qual = [headerSplitRegex.split(record[0])[0]]
    
#     # Sequence segmentation
#     m=patMatching(record[1],settings.regexDictCompiled[readKey])
#     if m:
#         mdict=m.groupdict()
#         for seg in segment_parsed:
#             if seg in mdict:
#                 out_seq.append(mdict[seg])
#                 out_qual.append(record[3][m.span(seg)[0] : m.span(seg)[1]])
#             else:
#                 out_seq.append("-")
#                 out_qual.append("-")
#     else:
#         out_seq += ["-"] * len(segment_parsed)
#         out_qual += ["-"] * len(segment_parsed)

#     return out_seq,out_qual


# Multithreading implementation for read segmentation
def segmentation_parallel_wrapper(fastq_chunk,settings,headerSplitRegex,readKey,segment_parsed,outdir,ncore):
    n_lines = len(fastq_chunk)
    n_records = n_lines / 4
    n_records_per_CPU = math.ceil(n_records / ncore)
    n_lines_per_CPU = 4 * n_records_per_CPU
    regex_pattern_now = settings.regexDictCompiled[readKey]

    out = Parallel(n_jobs=ncore,verbose=2,backend="multiprocessing")(
        delayed(fastq_segmentation)(cpu_idx,subchunk,headerSplitRegex,segment_parsed,regex_pattern_now,outdir,readKey,settings.barcodes) 
            for cpu_idx,subchunk in enumerate(split_yield_fastq_per_lines(fastq_chunk,n = n_lines_per_CPU))) 

    # Pile up the count dictionaries from all CPUs
    counterDict = {}
    for n,counterDict_tmp in enumerate(out):
        if n == 0:
            for seg in counterDict_tmp:
                counterDict[seg] = counterDict_tmp[seg]
        else:
            for seg in counterDict_tmp:
                counterDict[seg].update(counterDict_tmp[seg])

    return n_records,counterDict


# Merging parsed seq dicts
def merge_parsed_data_process(input_dir,cpu_idx,settings):
    
    filepaths_seq = glob.glob("_".join([input_dir+"/CPU"+str(cpu_idx),"*","srcSeq.pkl"]))
    filepaths_qual = glob.glob("_".join([input_dir+"/CPU"+str(cpu_idx),"*","srcQual.pkl"]))

    # if filepaths_seq:
    #     iter_num+=1
    # else:
    #     continue
        
    # if filepaths_seq:
    filepaths_seq.sort()
    filepaths_qual.sort()
    filepaths=[filepaths_seq,filepaths_qual]

    seq_qual_tsv_list = []

    for n_read,pathlist in enumerate(filepaths):
        dict_merged=collections.defaultdict(list)
        to_be_processed = ""
        
        for path in pathlist:
            if settings.flash:
                if "_merge_src_srcSeq.pkl" in path or "_merge_src_srcQual.pkl" in path:
                    to_be_processed=path
                    continue
            
            # Load segment seq/qual dictionaries from paired end reads and merge them into a single dictionary
            with open(path,mode="rb") as pchunk:
                parsedDict_chunk=pickle.load(pchunk)
            dict_merged.update(parsedDict_chunk)
        
        # FLASH data merging
        if settings.flash and to_be_processed != "":
            nrow_uncombined=len(dict_merged["Header"])

            # Load FLASHed segment sequences
            with open(to_be_processed,mode="rb") as pchunk:
                parsedDict_chunk = pickle.load(pchunk)
            
            # Puck the FLASHed segments but also put segments from uncombined segments
            for component in ["Header"]+settings.components:
                if component in parsedDict_chunk:
                    if component not in dict_merged:
                        dict_merged[component] = ["-"]*nrow_uncombined
                    dict_merged[component] += parsedDict_chunk[component]

                if component not in parsedDict_chunk:
                    continue
        
        # # FLASH data merging; if the chunk only contains non-flash reads
        # if settings.flash and to_be_processed == "":
        #     nrow_uncombined=len(dict_merged["Header"])

        #     # Puck the non-FLASH segments
        #     for component in ["Header"]+settings.components:
        #         if component not in dict_merged:

                    

                                 

        dict_merged_key=["Header"]+settings.components
        dict_merged_val=[dict_merged[i] for i in ["Header"]+settings.components]
        dict_merged_val=list(map(list,zip(*dict_merged_val)))
        dict_merged_val = pd.DataFrame(dict_merged_val,columns=dict_merged_key)
        
        seq_qual_tsv_list.append(dict_merged_val)
        
    return seq_qual_tsv_list
        # f.write(dict_merged_val)
        # f.close()


# Multithreading implementation for file merging
def merge_parsed_data_parallel_wrapper(input_dir, settings, ncore, n_chunk):
    retLst = Parallel(n_jobs=ncore,verbose=2)(
        delayed(merge_parsed_data_process)(input_dir,cpu_idx,settings) for cpu_idx in range(ncore))

    pd.concat([i[0] for i in retLst]).reset_index(drop=True).to_pickle(settings.outFilePath_and_Prefix+"_Chunk"+str(n_chunk)+"_srcSeq.pkl")
    pd.concat([i[1] for i in retLst]).reset_index(drop=True).to_pickle(settings.outFilePath_and_Prefix+"_Chunk"+str(n_chunk)+"_srcQual.pkl")
    
    # with gzip.open(settings.outFilePath_and_Prefix+"_Chunk"+str(n_chunk)+"_srcSeq.tsv.gz",mode="wt",encoding="utf-8") as w:
    #     w.write("\n".join([i[0] for i in retLst]))
    
    # with gzip.open(settings.outFilePath_and_Prefix+"_Chunk"+str(n_chunk)+"_srcQual.tsv.gz",mode="wt",encoding="utf-8") as w:
    #     w.write("\n".join([i[1] for i in retLst]))



# Multithreading implementation for quality filtering
def qfilter_parallel_wrapper(seq_chunk, qual_chunk, qc_targets, qscore_dict, ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    seq_subchunks = np.array_split(seq_chunk,ncore)
    qual_subchunks = np.array_split(qual_chunk,ncore)

    retLst = Parallel(n_jobs=ncore,verbose=2)(delayed(qualityFilteringForDataFrame)(df_set, qc_targets, qscore_dict) for df_set in zip(seq_subchunks,qual_subchunks))
    return pd.concat(retLst)


# Multithreading implementation for merging segment files
# def merge_count_tree_parallel_wrapper(pkl_path_list,ncore):
#     # group the path list by pairs
#     pkl_path_pairs = [pkl_path_list[idx:idx + 2] for idx in range(0, len(pkl_path_list), 2)]

#     # First round, input = path
#     count_tree_list = Parallel(n_jobs=ncore,require='sharedmem',verbose=8)(
#         delayed(merge_count_tree)(pkl_path_pair,path=True) for pkl_path_pair in pkl_path_pairs)

#     # 2nd~ round, iterate the pair merging processes until everything is merged
#     while True:
#         if len(count_tree_list) == 1:
#             # Exit condition: everything is merged
#             return count_tree_list[0]
#         else:
#             # Group count tree list by pairs if there are more than one trees in the list
#             count_tree_list = [count_tree_list[idx:idx + 2] for idx in range(0, len(count_tree_list), 2)]
#             count_tree_list = Parallel(n_jobs=ncore,require='sharedmem',verbose=8)(
#                 delayed(merge_count_tree)(tree_pair, path=False) for tree_pair in count_tree_list)
