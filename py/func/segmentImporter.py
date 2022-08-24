import regex
import gzip
import random
import string
import subprocess
import collections
import pickle
from joblib import Parallel, delayed
import multiprocessing
from itertools import zip_longest



def merge_reads_flash2(readPathDict,flash,gzipped,tmpdir,config):
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
        cmdlist_flash=["flash2","-z","-m",flash_min,"-M",flash_max,"-d",tmpdir,"-o","merge",decompressed_fastq[0],decompressed_fastq[1]]
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


# def helper_map(df,conv_dict,seg):
#     time.sleep(.5)
#     out = df[seg].map(lambda x: test_conversion(x,conv_dict))
#     return out

# def applyParallel(df, func, ncore, conv_dict, seg):
#     retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(func)(df_chunk,conv_dict,seg) for df_chunk in df)
#     return pd.concat(retLst)

# for ncore in range(1,11):
#     df_split = pd.read_csv("/Users/yusukekijima/work/yachie/droParser/tool_test/Interteslaar/test220812_debug_segname/process_ex1/E10.5_10x/value_extraction/_work/import/E105_10X_R1.fastq.gz_srcSeq.tsv.gz",
#                     sep="\t",
#                     header=0,
#                     chunksize=10000)
#     start = time.time()
#     result = applyParallel(df_split, helper_map, ncore=ncore, conv_dict=conversion_dict, seg="src_segment1")
#     end = time.time()
#     print("Ncores =",ncore,"Elapsed time:",end - start,"secs")


def split_fastq_per_lines(lines, n=4):
    for idx in range(0, len(lines), n):
        yield lines[idx:idx + n]
    

def segmentation_parallel_wrapper(fastq_chunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set,ncore):
    parsedSeqDict=collections.defaultdict(list)
    parsedQualDict=collections.defaultdict(list)
    Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(
        delayed(fastq_segmentation)(parsedSeqDict,parsedQualDict,subchunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set) 
            for subchunk in split_fastq_per_lines(fastq_chunk,n=40000))
    return parsedSeqDict,parsedQualDict


def fastq_segmentation(parsedSeqDict,parsedQualDict,subchunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set):
    for lines in split_fastq_per_lines(subchunk):
        line_header = lines[0]
        line_sequence = lines[1]
        line_quality = lines[3]
        
        # Header treatment
        header=headerSplitRegex.split(line_header)[0]
        parsedSeqDict["Header"].append(header)
        parsedQualDict["Header"].append(header)
        
        # Sequence segmentation
        m=patMatching(line_sequence,settings.regexDictCompiled[readKey])
        if m:
            mdict=m.groupdict()
            for seg in mdict:
                parsedSeqDict[seg].append(mdict[seg])

            component_diff=segment_parsed_set-set(mdict.keys())
            for seg in component_diff:
                parsedSeqDict[seg].append("-")
                parsedQualDict[seg].append("-")
        else: # If the regex doesn't match, quality segments should also be all missing
            for seg in segment_parsed:
                parsedSeqDict[seg].append("-")
                parsedQualDict[seg].append("-")

        # Quality segmentation
        if m:
            for component in mdict:
                extractedQual=line_quality[m.span(component)[0]:m.span(component)[1]]
                parsedQualDict[component].append(extractedQual)


    # if (nrow+1)%4000000==0:
    #     n_chunk=int((nrow+1)/4000000) #n_chunk>0
    #     with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcSeq.pkl"]),mode="wb") as p:
    #         pickle.dump(parsedSeqDict,p)
    #     with open("_".join([tmpdir+"/",readKey,str(n_chunk),"srcQual.pkl"]),mode="wb") as p:
    #         pickle.dump(parsedQualDict,p)


    #     for i in settings.barcodes:
    #         if parsedSeqDict.get(i):
    #             if i in counterDict:
    #                 counterDict_tmp=collections.Counter(parsedSeqDict[i])
    #                 counterDict[i].update(counterDict_tmp)
    #             else:
    #                 counterDict[i]=collections.Counter(parsedSeqDict[i])

    #     parsedSeqDict=collections.defaultdict(list)
    #     parsedQualDict=collections.defaultdict(list)
    #     print(str(int((nrow+1)/4))+" reads were processed for "+readKey,flush=True)
    
    # numSeqDict[readKey] = nrow