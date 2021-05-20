import regex
import gzip
import random
import string
import subprocess

def merge_reads_flash2(readPathDict,flash,gzipped,tmpdir,config):
    if not flash:
        return [],readPathDict
    else:
        flash=flash.replace("R","Read")
        flash=flash.replace("I","Index")
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

        flash_min=str(config["flash_min_overlap"])
        flash_max=str(config["flash_max_overlap"])
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