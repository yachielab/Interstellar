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

def findMostFeasibleCandidate(seq,suggestion_verbosity,dist_thresh,symspelldb,wlset=None):
    if wlset:
        if seq in wlset:
            return seq
    suggestions = symspelldb.lookup(seq,suggestion_verbosity,dist_thresh,transfer_casing=True)
    if len(suggestions)>=1:
        top_candidate=suggestions[0].term.upper()
        return top_candidate
    else:
        return "-"

    
def getRawIndex(seq,rawReference):
    return rawReference.index(seq)+1

def bcCorrect(correctOpt,counterDict,yaxis_scale,show_summary,outname):
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

    analyzedPosition=len(seqCountSummary["count"])-1        
    
    if "KNEE_CORRECT" in correctOpt["func_ordered"]:
        rank_threshold=correctOpt["KNEE_CORRECT"]["rank"]
        # seed_min=correctOpt.get("seed_min")
            
        if rank_threshold=="auto":
            knee=KneeLocator(seqCountSummary["rank"],seqCountSummary["count"],S=10,curve="convex",direction="decreasing",interp_method='interp1d')
            #knee=KneeLocator(seqCountSummary["logrank"],seqCountSummary["count"],S=1,curve="convex",direction="decreasing",interp_method='polynomial')
            print("Finding knee for",srcComponent,"was done!",flush=True)
            kneepoint_idx=seqCountSummary["rank"].index(knee.knee)
            #kneepoint_idx=seqCountSummary["logrank"].index(knee.knee)
        # elif rank_threshold=="all":
        #     kneepoint_idx=analyzedPosition
        else:
            kneepoint_idx=int(rank_threshold)
        seq_majority=seqCountSummary["seq"][:kneepoint_idx]
        seq_minority=seqCountSummary["seq"][kneepoint_idx:analyzedPosition]
        seq_discarded=seqCountSummary["seq"][analyzedPosition:]
        print("majority length",str(len(seq_majority)),flush=True)

        #store majority dictionary
        correctionDict_maj={}
        seqlen_min=min([len(i) for i in seq_majority if not i == "-"])
        seedlen=math.floor(seqlen_min/2)
        # if seed_min and seedlen<seed_min:
        #     seedlen=seed_min
        suggestion_verbosity = Verbosity.TOP
        today_now=str(datetime.datetime.today())
        today_now=regex.sub(r"\.|:| ","-",today_now)
        random_string=''.join(random.choices(string.ascii_letters + string.digits, k=15))
        fname="_".join([outname,srcComponent,today_now,random_string,"seq_majority.tmp"])
        if seq_minority:
            with open(fname,mode="wt") as w:
                for seqm in seq_majority:
                    w.write(seqm+"\n")

            symspelldb=SymSpell(correctOpt["KNEE_CORRECT"]["dist"],seedlen)
            symspelldb.create_dictionary(fname)

            t0=time.time()
            print("Minority correction has been started...",flush=True)
            seq_minority_pd=pd.Series(seq_minority)
            seq_minority_pd_corrected=seq_minority_pd.map(lambda x:findMostFeasibleCandidate(x,suggestion_verbosity,correctOpt["KNEE_CORRECT"]["dist"],symspelldb))

            correctionDict_maj={k:v for k,v in zip(list(seq_minority_pd),list(seq_minority_pd_corrected))}
            t1=time.time()
            print("Minority correction done",round(t1-t0),"sec")
            seq_success=seq_minority_pd[seq_minority_pd_corrected!="-"]
            seq_fail=seq_minority_pd[seq_minority_pd_corrected=="-"]
            os.remove(fname)
            
        for seq in seq_majority:
            correctionDict_maj[seq]=seq
        for seq in seq_discarded:
            correctionDict_maj[seq]="-"
        correctionDict=dict(correctionDict=correctionDict_maj,reference=seq_majority)

    if "WHITELIST_CORRECT" in correctOpt["func_ordered"]:
        if not "KNEE_CORRECT" in correctOpt["func_ordered"]:
            seq_majority=seqCountSummary["seq"][:kneepoint_idx]
            correctionDict_maj={}
            for seq in seq_majority:
                correctionDict_maj[seq]=seq
            correctionDict=dict(correctionDict=correctionDict_maj,reference=seq_majority)

        print("Whitelist correction has been started...",flush=True)
        print("Reference building...",flush=True)
        with open(correctOpt["WHITELIST_CORRECT"]["path"],mode="rt",encoding="utf-8") as f:
            wl=[regex.sub("\n","",i) for i in f]

        wlset=set(wl)
        print("Seed length: ",seedlen,flush=True)
        symspelldb=SymSpell(correctOpt["WHITELIST_CORRECT"]["dist"],seedlen)
        symspelldb.create_dictionary(correctOpt["WHITELIST_CORRECT"]["path"])
        print("Reference build done.",flush=True)
        seq_majority_pd=pd.Series(seq_majority)
        print("correct...",flush=True)
        t0=time.time()
        seq_majority_pd_corrected=seq_majority_pd.map(lambda x:findMostFeasibleCandidate(x,suggestion_verbosity,correctOpt["WHITELIST_CORRECT"]["dist"],symspelldb,wlset))
        t1=time.time()
        print("Whitelist correction done",t1-t0,"sec",flush=True)
        correctionDict_wl={k:v for k,v in zip(list(seq_majority_pd),list(seq_majority_pd_corrected))}
        seq_success=seq_majority_pd[seq_majority_pd_corrected!="-"]
        seq_fail=seq_majority_pd[seq_majority_pd_corrected=="-"]

        for seq in seq_minority+seq_discarded:
            correctionDict_wl[seq]="-"
            if not correctionDict_maj[seq]=="-":
                corrected_with_majority=correctionDict_maj[seq]
                if not correctionDict_wl[corrected_with_majority]=="-":
                    correctionDict_wl[seq]=correctionDict_wl[corrected_with_majority]
        t2=time.time()
        print("Correction dictionary building done",round(t2-t1),"sec",flush=True)
                
        correctionDict["correctionDict"]=correctionDict_wl
        correctionDict["reference"]=wl

    total_count=0
    t3=time.time()
    print("make plotting data...",flush=True)
    rawSeqIndexDict={k:v for v,k in enumerate(seqCountSummary["seq"])}
    for i in correctionDict["correctionDict"]:
        if correctionDict["correctionDict"][i]!="-":
            total_count+=seqCountSummary["count"][rawSeqIndexDict[i]]
    print("plotting data was made",round(time.time()-t3),flush=True)
    print("plotting...\n",flush=True)
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

