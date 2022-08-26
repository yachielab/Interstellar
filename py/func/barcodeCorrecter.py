from kneed import KneeLocator
import regex
import math
import matplotlib.pyplot as plt
from symspellpy.symspellpy import SymSpell, Verbosity
import pandas as pd
import time
import numpy as np
import os
import datetime
import random
import string
import pandas as pd
from joblib import Parallel, delayed



def getRawIndex(seq,rawReference):
    return rawReference.index(seq)+1


# Helper function for error correction
def sequenceCorrection(seq_chunk,suggestion_verbosity,correctOpt,correction_method,symspelldb,wlset):
    seq_chunk = seq_chunk.map(lambda x:findClosestCandidate(x,suggestion_verbosity,correctOpt[correction_method]["levenshtein_distance"],symspelldb,wlset))
    return seq_chunk

# Return closest hit
def findClosestCandidate(seq,suggestion_verbosity,dist_thresh,symspelldb,wlset):
    if wlset:
        if seq in wlset:
            return seq
    suggestions = symspelldb.lookup(seq,suggestion_verbosity,dist_thresh,transfer_casing=True)
    if len(suggestions)>=1:
        top_candidate=suggestions[0].term.upper()
        return top_candidate
    else:
        return "-"

# Multithreading implementation for error correction
def correct_parallel_wrapper(seq_series,suggestion_verbosity,correctOpt,correction_method,symspelldb,ncore,wlset=None):
    # Split the sequence series into chunks by the number of CPUs
    seq_chunks = np.array_split(seq_series,ncore)
    
    retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(sequenceCorrection)(seq_chunk,suggestion_verbosity,correctOpt,correction_method,symspelldb,wlset) for seq_chunk in seq_chunks)
    return pd.concat(retLst)


def bcCorrect(correctOpt,counterDict,yaxis_scale,show_summary,outname,ncore):
    suggestion_verbosity = Verbosity.TOP
    func_tmp=correctOpt["func_ordered"][0]
    srcComponent=correctOpt[func_tmp]["source"]
    srcCounter=counterDict[srcComponent]

    seqCount_sort=sorted(srcCounter.items(),key=lambda x:x[1],reverse=True)
    seqCountSummary=dict(seq=list(),rank=list(),logrank=list(),count=list())
    for idx,countSet in enumerate(seqCount_sort):
        if countSet[0]=="-":
            continue
        seqCountSummary["seq"].append(countSet[0])
        seqCountSummary["count"].append(countSet[1])
        seqCountSummary["rank"].append(idx+1)
        seqCountSummary["logrank"].append(math.log10(idx+1))

    analyzedPosition=len(seqCountSummary["count"])
    kneepoint_idx=len(seqCountSummary["count"])
    # if kneepoint_idx ==0:
    #     kneepoint_idx = 1
    seqlen_min=min([len(i) for i in seqCountSummary["seq"]])
    seedlen=math.floor(seqlen_min/2)

    seq_minority=[]
    seq_discarded=[]

    if "I2M_CORRECTION" in correctOpt:
        rank_threshold=correctOpt["I2M_CORRECTION"]["rank"]
    else:
        rank_threshold=len(seqCountSummary["count"])

    # Imputation-to-majority
    if "I2M_CORRECTION" in correctOpt["func_ordered"]:

        if rank_threshold=="auto":
            knee=KneeLocator(seqCountSummary["rank"],seqCountSummary["count"],S=10,curve="convex",direction="decreasing",interp_method='interp1d')
            print("Knee point of",srcComponent,"was identified.",flush=True)
            kneepoint_idx=seqCountSummary["rank"].index(knee.knee)
        else:
            kneepoint_idx=int(rank_threshold)
        seq_majority=seqCountSummary["seq"][:kneepoint_idx]
        seq_minority=seqCountSummary["seq"][kneepoint_idx:analyzedPosition]
        seq_discarded=seqCountSummary["seq"][analyzedPosition:]
        print("Majority length",str(len(seq_majority)),flush=True)

        #store majority dictionary
        correctionDict_maj={}
        seqlen_min=min([len(i) for i in seq_majority if not i == "-"])
        seedlen=math.floor(seqlen_min/2)
        
        today_now=str(datetime.datetime.today())
        today_now=regex.sub(r"\.|:| ","-",today_now)
        random_string=''.join(random.choices(string.ascii_letters + string.digits, k=15))
        fname="_".join([outname,srcComponent,today_now,random_string,"seq_majority.tmp"])
        if seq_minority:
            with open(fname,mode="wt") as w:
                for seqm in seq_majority:
                    w.write(seqm+"\n")

            symspelldb=SymSpell(correctOpt["I2M_CORRECTION"]["levenshtein_distance"],seedlen)
            symspelldb.create_dictionary(fname)

            t0=time.time()
            print("Imputation-to-majority correction has been started...",flush=True)
            seq_minority_pd=pd.Series(seq_minority)
            seq_minority_pd_corrected = correct_parallel_wrapper(seq_minority_pd,suggestion_verbosity,correctOpt,"I2M_CORRECTION",symspelldb,ncore)
            # seq_minority_pd_corrected=seq_minority_pd.map(lambda x:findMostFeasibleCandidate(x,suggestion_verbosity,correctOpt["I2M_CORRECTION"]["levenshtein_distance"],symspelldb))

            correctionDict_maj={k:v for k,v in zip(list(seq_minority_pd),list(seq_minority_pd_corrected))}
            print("Done.")
            os.remove(fname)

        for seq in seq_majority:
            correctionDict_maj[seq]=seq
        for seq in seq_discarded:
            correctionDict_maj[seq]="-"
        correctionDict=dict(correctionDict=correctionDict_maj,reference=seq_majority)

    # Mapping-to-allowlist
    if "M2A_CORRECTION" in correctOpt["func_ordered"]:
        if not "I2M_CORRECTION" in correctOpt["func_ordered"]:
            seq_majority=seqCountSummary["seq"][:kneepoint_idx]
            correctionDict_maj={}
            for seq in seq_majority:
                correctionDict_maj[seq]=seq
            correctionDict=dict(correctionDict=correctionDict_maj,reference=seq_majority)

        print("Map-to-allowlist correction has been started...",flush=True)
        print("Building reference...",flush=True)
        with open(correctOpt["M2A_CORRECTION"]["path"],mode="rt",encoding="utf-8") as f:
            wl=[regex.sub("\n","",i) for i in f]

        wlset=set(wl)
        print("Seed length: ",seedlen,flush=True)
        symspelldb=SymSpell(correctOpt["M2A_CORRECTION"]["levenshtein_distance"],seedlen)
        symspelldb.create_dictionary(correctOpt["M2A_CORRECTION"]["path"])
        seq_majority_pd=pd.Series(seq_majority)
        print("Correct...",flush=True)
        seq_majority_pd_corrected = correct_parallel_wrapper(seq_majority_pd,suggestion_verbosity,correctOpt,"M2A_CORRECTION",symspelldb,ncore,wlset)
        # seq_majority_pd_corrected=seq_majority_pd.map(lambda x:findMostFeasibleCandidate(x,suggestion_verbosity,correctOpt["M2A_CORRECTION"]["levenshtein_distance"],symspelldb,wlset))
        print("Done.",flush=True)
        correctionDict_wl={k:v for k,v in zip(list(seq_majority_pd),list(seq_majority_pd_corrected))}
        # seq_success=seq_majority_pd[seq_majority_pd_corrected!="-"]
        # seq_fail=seq_majority_pd[seq_majority_pd_corrected=="-"]

        for seq in seq_minority+seq_discarded:
            correctionDict_wl[seq]="-"
            if not correctionDict_maj[seq]=="-":
                corrected_with_majority=correctionDict_maj[seq]
                if not correctionDict_wl[corrected_with_majority]=="-":
                    correctionDict_wl[seq]=correctionDict_wl[corrected_with_majority]
        t2=time.time()

        correctionDict["correctionDict"]=correctionDict_wl
        correctionDict["reference"]=wl

    total_count=0
    t3=time.time()
    # print("make plotting data...",flush=True)
    rawSeqIndexDict={k:v for v,k in enumerate(seqCountSummary["seq"])}
    for i in correctionDict["correctionDict"]:
        if correctionDict["correctionDict"][i]!="-":
            total_count+=seqCountSummary["count"][rawSeqIndexDict[i]]
    # print("plotting data was made",round(time.time()-t3),flush=True)
    # print("plotting...\n",flush=True)
    fig=plt.figure()
    p1=fig.add_subplot(1,1,1)
    if "-" in seqCountSummary["seq"]:
        idx_nohit=seqCountSummary["seq"].index("-")
        del seqCountSummary["count"][idx_nohit]
        del seqCountSummary["rank"][len(seqCountSummary["rank"])-1]
    p1.plot(seqCountSummary["rank"],seqCountSummary["count"],c="pink")
    if show_summary:
        p1.text(1,2,"detected barcode: "+str(sum(seqCountSummary["count"])))
        p1.text(1,1,"corrected barcode: "+str(total_count))
    if not rank_threshold=="all":
        plt.vlines(kneepoint_idx,min(seqCountSummary["count"]),max(seqCountSummary["count"]),linestyles="dashed")
        p1.text(kneepoint_idx+1,max(seqCountSummary["count"]),kneepoint_idx)
    plt.xscale("log")
    if yaxis_scale=="log":
        plt.yscale("log")

    plt.savefig(outname+"_barcodeSummary_"+srcComponent+".png")

    if not "-" in correctionDict["correctionDict"]:
        correctionDict["correctionDict"]["-"]="-"

    if not "-" in correctionDict["reference"]:
        correctionDict["reference"].append("-")

    return correctionDict



def qualityProcessing(qualitySequence):
    qualSum=0
    for c in qualitySequence:
        qualSum+=ord(c)
    qualAvg=round(qualSum/len(qualitySequence))
    return qualAvg


def seq_correct_and_write(seq,reference):
    if seq in reference:
        return reference[seq]
    else:
        return "-"


# def seq_to_val(seq,d_ref):
#     if seq in d_ref["ref_set"]:
#         # print(seq,d_ref["ref_np"].shape)
#         #return np.where(d_ref["ref_np"]==seq)
#         return d_ref["ref_tup"].index(seq)
#     else:
#         return -1


def seq_to_val_ver2(seq,dic):
    if seq in dic:
        if seq=="-":
            return -1
        return(dic[seq])
    else:
        return -1


def gen_bt_dict(bt_barcode,bt_cluster):
    df_barcode=pd.read_csv(bt_barcode)
    df_cluster=pd.read_csv(bt_cluster)
    df_barcode=df_barcode[["Unique.reads","Cluster.ID"]]
    df_cluster=df_cluster[["Cluster.ID","Center"]]
    df_barcode=df_barcode.merge(df_cluster,on="Cluster.ID")

    correctionDict={k:df_barcode["Center"][i] for i,k in enumerate(df_barcode["Unique.reads"])}
    return correctionDict


def custom_correction_setup(target_segment_list,outdir,src_seq_paths):
    dic_rowcount={k:0 for k in target_segment_list}

    # Export the input files for a custom error correction code per barcode
    for processCount,input_filepath in enumerate(src_seq_paths):
        df_tmp=pd.read_csv(input_filepath,sep="\t")

        for bc in target_segment_list:
            seg_seqs=[i for i in list(df_tmp[bc]) if i != "-"] # Segment sequence list without NA characters
            row_list=list(map(str,list(dic_rowcount[bc]+np.array(range(len(seg_seqs)))+1))) # Making row numbers
            dic_rowcount[bc]+=len(seg_seqs)

            # Export
            if processCount==0:
                with open(outdir+"/custom"+"_"+bc+".input.csv",mode="wt") as w:
                    for seq,row in zip(seg_seqs,row_list):
                        w.write(seq+"\t"+row+"\n")

            else:
                with open(outdir+"/custom"+"_"+bc+".input.csv",mode="at") as w:
                    for seq,row in zip(seg_seqs,row_list):
                        w.write(seq+"\t"+row+"\n")


def gen_custom_dict(df_corrected):
    seq_raw=df_corrected[0]
    seq_corrected=df_corrected[1]
    correctionDict={k:seq_corrected[i] for i,k in enumerate(seq_raw)}
    return correctionDict


# Create a clean segment sequence table
def gen_clean_segment_seq(df_subchunk,segments_raw,correctionDictionaries,settings):
    cols_ordered=["Header"]+[i for i in sorted(df_subchunk.columns) if not i=="Header"]
    cols_ordered_final=[]

    for ncol,segment_raw_now in enumerate(cols_ordered):
        if ncol>0:
            if not segment_raw_now in segments_raw:
                df_subchunk=df_subchunk.drop(segment_raw_now,axis=1)
                continue
            col_index=segments_raw.index(segment_raw_now)
            component_corrected_now=settings.corrected_components[col_index]
            opt_now=settings.correctOptDict[component_corrected_now]
            
            if (component_corrected_now in correctionDictionaries) and ("I2M_CORRECTION" in opt_now["func_ordered"] or "M2A_CORRECTION" in opt_now["func_ordered"] or "BARTENDER_CORRECTION" in opt_now["func_ordered"] or "CUSTOM_CORRECTION" in opt_now["func_ordered"]):
                df_subchunk[segment_raw_now]=df_subchunk[segment_raw_now].map(lambda x: seq_correct_and_write(x,reference=correctionDictionaries[component_corrected_now]["correctionDict"]))
            df_subchunk.rename(columns={segment_raw_now:segment_raw_now+":"+component_corrected_now},inplace=True)
            cols_ordered_final.append(segment_raw_now+":"+component_corrected_now)
        else:
            cols_ordered_final.append("Header")

    df_subchunk = df_subchunk[cols_ordered_final]
    return df_subchunk


# Multithreading implementation for cleas segment sequence table generation
def seqCleanUp_parallel_wrapper(df_chunk,segments_raw,correctionDictionaries,settings,ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    df_subchunks = np.array_split(df_chunk,ncore)

    retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(gen_clean_segment_seq)(df_subchunk,segments_raw,correctionDictionaries,settings) for df_subchunk in df_subchunks)
    return pd.concat(retLst)


# Sequence to value, sequence or quality segment
def seq_to_value_table(df_subchunk,cat,segments_raw,settings,correctionDictionaries,ref_dic):
    res=pd.DataFrame()
    for ncol,component in enumerate(df_subchunk.columns):                    
        if ncol==0:
            res["Header"]=df_subchunk[component]
        if ncol>0: 
            if cat=="seq":
                component_raw_now=component.split(":")[0]
                component_corrected_now=component.split(":")[1]
            else:
                component_raw_now=component
                if not component_raw_now in segments_raw:
                    continue
                col_index=segments_raw.index(component_raw_now)
                component_corrected_now=settings.corrected_components[col_index]

            if component_corrected_now in correctionDictionaries:
                if cat=="seq":
                    res[component_corrected_now]=df_subchunk[component].map(lambda x: seq_to_val_ver2(x,dic=ref_dic[component_corrected_now]))
                else:
                    res[component_corrected_now]=df_subchunk[component_raw_now].map(qualityProcessing)
    return res


# Multithreading implementation for value table generation
def gen_value_table_parallel_wrapper(df_chunk,cat,segments_raw,settings,correctionDictionaries,ref_dic,ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    df_subchunks = np.array_split(df_chunk,ncore)

    retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(seq_to_value_table)(df_subchunk,cat,segments_raw,settings,correctionDictionaries,ref_dic) for df_subchunk in df_subchunks)
    return pd.concat(retLst)