import collections
import regex
import itertools as it
import random
import pandas as pd
import numpy as np
import gzip,pickle
import copy
from joblib import Parallel, delayed
import os
import sys
import subprocess


class UnknownError(Exception):
    pass
# random.seed(0)


def dval_to_sval_relationship(func_dict,dest_segment):
    dval_to_sval_relationship={}
    for val in dest_segment:
        opt_now=func_dict[val][func_dict[val]["func_ordered"][0]]
        s_values=opt_now["source"]
        
        s_value="+".join(s_values)
        dval_to_sval_relationship[val]=s_value

    return dval_to_sval_relationship


def parse_constraint(value_segment,values_in_destarg,child2parent_val,value_variables,value2seq_sources):
    noParentSegments=[]
    edge_dict=dict(parent=[],child=[])
    for component in child2parent_val:
        edge_dict["child"].append(component)
        edge_dict["parent"].append(child2parent_val[component])
    for val in value_segment:
        if not val in edge_dict["child"]:
            noParentSegments.append(val)
    
    roots=list(set(edge_dict["parent"])-set(edge_dict["child"]))
    # tips =list(set(edge_dict["child"])-set(edge_dict["parent"]))
    globalComponents=list(set(noParentSegments+values_in_destarg)-set(roots)-set(edge_dict["child"]))
    globalComponents=[i for i in globalComponents if i in value_variables+value2seq_sources]

    roots=[i.replace(",","+") for i in roots]
    globalComponents=[i.replace(",","+") for i in globalComponents]
    edge_dict["child"]=[k.replace(",","+") for k in edge_dict["child"]]
    edge_dict["parent"]=[k.replace(",","+") for k in edge_dict["parent"]]
    
    return roots,edge_dict,globalComponents


def to_svalue_prime(s_val_chunk,dval_to_sval_relationship,roots,edge_dict):
    #update s-value by interpreting the combination
    combi_idx=0
    for component in dval_to_sval_relationship:
        component_svalue=dval_to_sval_relationship[component]
        # Generate series of s-values 
        if "+" in component_svalue:
            combi_idx+=1
            component_svalue_split = component_svalue.split("+")
            series_svalue = s_val_chunk[component_svalue_split[0]].str.cat(s_val_chunk[component_svalue_split[1:]],sep="+")
            series_svalue.name=component_svalue
            if combi_idx==1:
                s_val_chunk_prime = pd.DataFrame(series_svalue)
            else:
                s_val_chunk_prime = pd.concat([s_val_chunk_prime,series_svalue],axis=1)
    if combi_idx>0:
        s_val_chunk=pd.concat([s_val_chunk,s_val_chunk_prime],axis=1)

    combi_idx=0
    input_list=list(set(roots+edge_dict["parent"]+edge_dict["child"]))
    for val in input_list:
        # Generate series of s-values 
        if "+" in val:
            # val=val.replace(",","+")
            if not val in s_val_chunk.columns:
                combi_idx+=1
                values = val.split("+")
                series_svalue = s_val_chunk[values[0]].str.cat(s_val_chunk[values[1:]],sep="+")
                series_svalue.name=val
                if combi_idx==1:
                    s_val_chunk_prime = pd.DataFrame(series_svalue)
                else:
                    s_val_chunk_prime = pd.concat([s_val_chunk_prime,series_svalue],axis=1)
    if combi_idx>0:
        s_val_chunk=pd.concat([s_val_chunk,s_val_chunk_prime],axis=1)
    
    return s_val_chunk


# def to_svalue_prime_shrink(s_val_chunk,dval_to_sval_relationship):
#     #convert s-value into new s-value interpreting the combination
#     combi_idx=0
#     for component in dval_to_sval_relationship:
#         component_svalue=dval_to_sval_relationship[component]

#         # Generate series of s-values 
#         if "+" in component_svalue:
#             combi_idx+=1
#             component_svalue_split = component_svalue.split("+")
#             series_svalue = s_val_chunk[component_svalue_split[0]].str.cat(s_val_chunk[component_svalue_split[1:]],sep="+")
#             series_svalue.name=component_svalue
#             if combi_idx==1:
#                 s_val_chunk_prime = pd.DataFrame(series_svalue)
#             else:
#                 s_val_chunk_prime = pd.concat([s_val_chunk_prime,series_svalue],axis=1)
#     s_val_chunk=pd.concat([s_val_chunk["Header"],s_val_chunk_prime],axis=1)
#     return s_val_chunk


def getAncestor(target,edge_dict):
    ancestor_list=[]
    while True:
        if not target in edge_dict["child"]:
            break
        parent_now=edge_dict["parent"][edge_dict["child"].index(target)]
        ancestor_list.append(parent_now)
        target=parent_now
    return ancestor_list[::-1] #Sort from the root to leaf


def update_tree(component,Tree_tmp,Tree,root=False):
    if root:
        # for n,i in enumerate(Tree_tmp):
        #     if n<3:
        #         print(component,i,Tree_tmp[i])
        if not component in Tree:
            Tree[component]=Tree_tmp
        else:
            Tree[component].update(Tree_tmp)
    else:
        if not component in Tree:
            Tree[component] = Tree_tmp
        else:
            for k in Tree_tmp:
                if k in Tree[component]:
                    Tree[component][k].update(Tree_tmp[k])
                else:
                    Tree[component][k]=Tree_tmp[k]
    # print(component)
    # for i in Tree:
    #     for n,k in enumerate(Tree[i]):
    #         if n<3:
    #             print(i,k,Tree[i][k])
    # print("\n")
    return Tree


def buildTree_global(s_val_chunk,globalComponents,Tree):
    for component in globalComponents:
        component=component.replace(",","+")
        Tree=update_tree(component,collections.Counter(s_val_chunk[component]),Tree,root=True)
    return Tree


def buildTree(s_val_chunk,roots,edge_dict,Tree):
    # Let the read structure like below
    # ROOT──C1──C3
    #   └───C2
    #
    # TREE structure
    # TREE[ROOT]={10:0, 2:1, 3:2, ...}
    # TREE[C1] = {0:{5:0, 2:1, 10:2, ...},    1:{100:0, 23:1, 310:2, ...},  2:{15:0, 112:1, 130:2, ...}, ...}
    # TREE[C2] = {0_5:{2:0, 1:1, 100:2, ...}, 0_2:{120:0, 27:1, 10:2, ...}, 0_10:{125:0, 1:1, 7:2, ...}, ...}
    # TREE[C3] = {0_5_2:{1:0, 8:1, 4:2, ...}, 0_5_1:{14:0, 12:1, 3:2, ...}, 0_5_100:{6:0, 7:1, 1:2, ...}, ...}

    for root_now in roots:
        #Add / Update the root component of Tree 
        Tree=update_tree(root_now,collections.Counter(s_val_chunk[root_now]),Tree,root=True)
        tasklist=[root_now] #starting from the root
        while True:
            if not tasklist:
                break

            parent=tasklist[0] #current parent (1)
            tasklist=tasklist[1:] #others

            childs=[edge_dict["child"][idx] for idx,i in enumerate(edge_dict["parent"]) if i==parent] #current childs for the parent (1): val1 / val2,val3 / ...
            tasklist+=childs #the current childs are gonna be the next parents

            if not childs:
                continue

                
            for child in childs:
                ancestor_list = getAncestor(child,edge_dict) #from the root to leaf
                # ancestor_list = [i.replace(",","+") for i in ancestor_list]
                # parent=parent.replace(",","+")
                # child=child.replace(",","+")
                parent_series = s_val_chunk[ancestor_list[0]].str.cat(s_val_chunk[ancestor_list[1:]],sep="_")
                parent_series.name = parent
                s_val_tmp = pd.concat([parent_series,s_val_chunk[child]],axis=1)
                s_val_tmp = s_val_tmp.groupby(parent)
                s_val_tmp = s_val_tmp[child].apply(list).to_dict()
                Tree_tmp = {k:collections.Counter(s_val_tmp[k]) for k in s_val_tmp}

                
                Tree=update_tree(child,Tree_tmp,Tree)
    return Tree


# Data processing function for buildTree
def buildTree_process(df_subchunk,roots,edge_dict,d2s_dict,settings,globalComponents,sample_now):
    df_subchunk=df_subchunk.astype(str)
    df_subchunk=df_subchunk.replace("-1",np.nan).dropna()
    df_subchunk=to_svalue_prime(df_subchunk,d2s_dict,roots,edge_dict)

    if settings.samplemerge:
        for component in globalComponents+roots:
            df_subchunk[component]=df_subchunk[component]+":"+sample_now+":"

    count_tree_subchunk = {}
    if globalComponents:
        count_tree_subchunk = buildTree_global(df_subchunk,globalComponents,Tree=count_tree_subchunk)
    count_tree_subchunk = buildTree(df_subchunk,roots,edge_dict,Tree=count_tree_subchunk)

    return count_tree_subchunk


# Multithreading implementation for building a value tree
def build_count_tree_parallel_wrapper(df_chunk,roots,edge_dict,d2s_dict,settings,globalComponents,sample_now,count_tree,ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    df_subchunks = np.array_split(df_chunk,ncore)

    count_tree_list = Parallel(n_jobs=ncore,verbose=8,backend="multiprocessing")(
        delayed(buildTree_process)(df_subchunk,roots,edge_dict,d2s_dict,settings,globalComponents,sample_now) for df_subchunk in df_subchunks)
    
    # Merge each subchunk tree generated by multiprocessing into the main tree
    for count_tree_subchunk in count_tree_list:
        for segment in count_tree_subchunk:
            if not segment in count_tree: #If the count tree was formatted
                count_tree[segment]=count_tree_subchunk[segment]
            else: #2nd~
                if type(count_tree_subchunk[segment])==collections.Counter: #For global value
                    count_tree[segment].update(count_tree_subchunk[segment])
                else: #For local value
                    for parent in count_tree_subchunk[segment]:
                        child_freq_now=count_tree_subchunk[segment][parent]
                        if parent in count_tree[segment]:
                            count_tree[segment][parent].update(child_freq_now)
                        else:
                            count_tree[segment][parent]=child_freq_now
    return count_tree
                

def compression_global(s_value,component,Tree):
    if s_value in Tree[component]:
        return Tree[component][s_value]
    else:
        return "-1"


def compression_local(s_value_combined,component,Tree):
    constraint = s_value_combined.split(",")[0]
    s_value=s_value_combined.split(",")[1]

    if constraint in Tree[component] and s_value in Tree[component][constraint]:
        return Tree[component][constraint][s_value]
    else:
        return "-1"


# Count tree merging process
def merge_count_tree(tree_pair,path=False):
    if path:
        # The input data is a pair of paths and needs to be loaded
        tree_pair_new = []
        for p in tree_pair:
            with gzip.open(p,mode="rb") as rp:
                count_tree_tmp=pickle.load(rp)
                tree_pair_new.append(count_tree_tmp)
        del tree_pair
        tree_pair = tree_pair_new

    # Merge tree pairs into a single count tree
    if len(tree_pair) == 1:
        return tree_pair[0]
    else:
        tree_query = tree_pair[0]
        tree_reference = tree_pair[1]
        for component in tree_query:
            if not component in tree_reference:
                tree_reference[component]=tree_query[component]
            else: 
                if type(tree_query[component])==collections.Counter: #For global value
                    tree_reference[component].update(tree_query[component])
                else: #For local value
                    for parent in tree_query[component]:
                        child_freq_now=tree_query[component][parent]
                        if parent in tree_reference[component]:
                            tree_reference[component][parent].update(child_freq_now)
                        else:
                            tree_reference[component][parent]=child_freq_now
        return tree_reference


# Multithreading implementation for merging count trees
def merge_count_tree_parallel_wrapper(pkl_path_list,ncore):
    # group the path list by pairs
    pkl_path_pairs = [pkl_path_list[idx:idx + 2] for idx in range(0, len(pkl_path_list), 2)]

    # First round, input = path
    count_tree_list = Parallel(n_jobs=ncore,backend='multiprocessing',verbose=8)(
        delayed(merge_count_tree)(pkl_path_pair,path=True) for pkl_path_pair in pkl_path_pairs)

    # 2nd~ round, iterate the pair merging processes until everything is merged
    while True:
        if len(count_tree_list) == 1:
            # Exit condition: everything is merged
            return count_tree_list[0]
        else:
            # Group count tree list by pairs if there are more than one trees in the list
            count_tree_list = [count_tree_list[idx:idx + 2] for idx in range(0, len(count_tree_list), 2)]
            count_tree_list = Parallel(n_jobs=ncore,backend='multiprocessing',verbose=8)(
                delayed(merge_count_tree)(tree_pair, path=False) for tree_pair in count_tree_list)


def printDictTest(d):
   if d:
    for i in d:
        for n,k in enumerate(d[i]):
            if n<10:
                print(i,k,d[i][k])


def localJudge(func_dict,constraint_tmp,constraintGroup):
    if type(constraint_tmp)==list:
        constraintGroup+=constraint_tmp
        constraint_type="non_local"
        for i in constraint_tmp:
            if func_dict[i]["type"]=="localConvert":
                constraint_type="localConvert"
    else:
        constraintGroup.append(constraint_tmp)
        constraint_type=func_dict[constraint_tmp]["type"]
    return constraintGroup,constraint_type


def sortCounter(counterDict):
    dict_tmp={}
    local_count_sort=sorted(counterDict.items(),key=lambda x:x[1],reverse=True)
    local_k_sorted=[c[0] for c in local_count_sort]
    if "-1" in local_k_sorted:
        del local_k_sorted[local_k_sorted.index("-1")]
        dict_tmp["-1"]=-1
    for idx,k in enumerate(local_k_sorted):
        dict_tmp[k]=idx
    return dict_tmp


# Value optimization process
def optimize_value_process(df_subchunk,d2s_dict,roots,edge_dict,settings,globalComponents,sample_now,global_used_in_dest,destname_dict,tree):
    df_subchunk=df_subchunk.astype(str)
    df_subchunk=to_svalue_prime(df_subchunk,d2s_dict,roots,edge_dict)
    compressed_svalue = pd.DataFrame({"Header":df_subchunk["Header"]})

    if settings.samplemerge:
        for component in globalComponents+roots:
            df_subchunk[component]=df_subchunk[component]+":"+sample_now+":"

    if global_used_in_dest:
        for component in global_used_in_dest:
            compressed_svalue[destname_dict[component]] = df_subchunk[component].map(lambda x: compression_global(x,component=component,Tree=tree))
    
    #Local convert, Depth-First Search
    if roots:
        for root_now in roots:
            tasklist=[root_now]
            compressed_svalue[destname_dict[root_now]] = df_subchunk[root_now].map(lambda x: compression_global(x,component=root_now,Tree=tree))    
            while True:
                if not tasklist:
                    break
                parent=tasklist[0]
                tasklist=tasklist[1:]
                childs=[edge_dict["child"][idx] for idx,i in enumerate(edge_dict["parent"]) if i==parent]
                tasklist+=childs

                if not childs:
                    continue
                    
                for child in childs:
                    if not child in destname_dict:
                        continue
                    ancestor_list = getAncestor(child,edge_dict)
                    parent_series = df_subchunk[ancestor_list[0]].str.cat(df_subchunk[ancestor_list[1:]],sep="_")
                    parent_series.name = parent
                    s_val_tmp = parent_series.str.cat(df_subchunk[child],sep=",")  
                    compressed_svalue[destname_dict[child]] = s_val_tmp.map(lambda x: compression_local(x,component=child,Tree=tree))
    return compressed_svalue


# Multithreading implementation for optimizing value table
def optimize_value_table_parallel_wrapper(df_chunk,d2s_dict,roots,edge_dict,settings,globalComponents,sample_now,global_used_in_dest,destname_dict,tree,ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    df_subchunks = np.array_split(df_chunk,ncore)

    retLst = Parallel(n_jobs=ncore,backend='multiprocessing',verbose=8)(
        delayed(optimize_value_process)(df_subchunk,d2s_dict,roots,edge_dict,settings,globalComponents,sample_now,global_used_in_dest,destname_dict,tree) for df_subchunk in df_subchunks)
    return pd.concat(retLst)


# Quality data handling process in conversion process
def quality_conversion_process(q_val_subchunk,segments_used,destname_dict):
    compressed_qvalue=pd.DataFrame({"Header":q_val_subchunk["Header"]})
    for segment in segments_used:
        if "+" in segment:
            svalues=segment.split("+")
            average_qvalue=q_val_subchunk[svalues].mean(axis=1).round().astype(int)
        else:
            average_qvalue=q_val_subchunk[segment].astype(int)
        
        compressed_qvalue[destname_dict[segment]] = average_qvalue
    return compressed_qvalue


# Multithreading implementation for optimizing value table
def quality_conversion_parallel_wrapper(df_chunk,segments_used,destname_dict,ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    df_subchunks = np.array_split(df_chunk,ncore)

    retLst = Parallel(n_jobs=ncore,verbose=8)(
        delayed(quality_conversion_process)(df_subchunk,segments_used,destname_dict) for df_subchunk in df_subchunks)
    return pd.concat(retLst)


def IUPAC(x,seed):
    random.seed(seed)
    if x in ["A","T","G","C"]:
        return x
    elif x=="R":
        y = "AG"
    elif x=="Y":
        y =  "CT"
    elif x=="S":
        y =  "GC"
    elif x=="W":
        y =  "AT"
    elif x=="K":
        y =  "GT"
    elif x=="M":
        y =  "AC"
    elif x=="B":
        y =  "CGT"
    elif x=="D":
        y =  "AGT"
    elif x=="H":
        y =  "ACT"
    elif x=="V":
        y =  "ACG"
    elif x=="N":
        y =  "ATGC"
    else:
        raise UnknownError("The chracter",x,"is not in the IUPAC code list!")

    random.seed(seed)
    return ''.join(random.sample(y ,len(y)))


def buildReference(opt_now,seg_idx,info_dic,max_default=7000000): #n_each?
    func_now=opt_now["func_ordered"][0]
    if func_now=="WHITELIST_ASSIGNMENT":
        with open(opt_now[func_now]["allowlist_path"][seg_idx],mode="rt",encoding="utf-8") as f:
            reference_now=[regex.sub("\n","",i) for i in f]
            random.seed(0)
            random.shuffle(reference_now)
    elif func_now=="RANDSEQ_ASSIGNMENT":
        val_segment=opt_now[func_now]["source"]
        maxsize=info_dic[val_segment] 
        if maxsize < max_default:
            maxsize=max_default

        # seqlength=int(opt_now[func_now]["randseq_length"][seg_idx])
        randseq=opt_now[func_now]["randseq_pattern"][seg_idx]
        reference_iter=it.product(*([IUPAC(i,seed) for seed,i in enumerate(randseq)]))
        reference_now=[]
        for cnt,i in enumerate(reference_iter):
            if cnt>=maxsize:
                break
            reference_now.append(i)
        reference_now=["".join(i) for i in reference_now]
        random.seed(0)
        random.shuffle(reference_now)
    return reference_now


def get_index(idx,ref):
    if "-1" in idx:
        return -1
    else:
        return ref[idx]


def mean(x):
    return round(sum(x)/len(x))


def getLocalIndex(row,constraint,localComponent,tree):
    convMotherIdx=row[0]
    srcLocalIdx=row[1]
    if "-1" in convMotherIdx or "-1" in srcLocalIdx:
        return -1
    else:
        return tree["constraintTree"][constraint][localComponent][convMotherIdx][srcLocalIdx]


def genEqSeq(query,datatype,length=None,add_nuc=None,baseQuality=None):
    if not length or length==0:
        length=len(query)
    if not add_nuc:
        add_nuc="A"   
    additionalSeqLength=0
    raw_length=len(query)
    if length > raw_length:
        additionalSeqLength=length-raw_length
    
    if query=="-":
        return np.nan

    if datatype=="seq":
        return query[:length]+add_nuc*additionalSeqLength
    elif datatype=="qual":
        baseQuality=round(sum(map(ord,query[:length]))/int(len(query[:length])))
        return query[:length]+chr(baseQuality)*additionalSeqLength 


def getConvSeq(query,reference):
    if query >= len(reference) or query < 0:
        return None
    else:
        return reference[query]


def getConvQual_ver2(query):
    if not pd.isna(query):
        L=query.split("_")
        length=len(L[0])
        avg_qual=int(L[1])
        return chr(avg_qual)*length
    else:
        return None


def get_svalue(query,component,ref):
    if query in ref[component]:
        return ref[component][query]
    else:
        return "-1"


def get_reindex(query,d_component,tree):
    reindex_dict=tree[d_component]
    if query not in reindex_dict:
        return None
    else:
        return reindex_dict[query]


def build_SEQ2SEQ(conversion_table,dict_to_terminal,value_segment,ref_dict):
    source_terminal_pool=[dict_to_terminal[i] for i in conversion_table.columns if i in dict_to_terminal]
    table_src =[i for i in value_segment if i in list(conversion_table.columns)+source_terminal_pool]
    table_dest=[i for i in conversion_table.columns if i not in value_segment]

    print("Source segment(s):",table_src)
    print("Destination segment(s):",table_dest)

    conversion_table_src = conversion_table[sorted(table_src)]
    conversion_table_dest= conversion_table[table_dest]

    if len(table_src) > 1:
        conversion_table_src = conversion_table_src[conversion_table_src.columns[0]].str.cat(conversion_table_src[conversion_table_src.columns[1:]],sep="_")
    else:
        conversion_table_src = conversion_table_src[conversion_table_src.columns[0]]

    for dest in conversion_table_dest.columns:
        ser = conversion_table_dest[dest]
        for idx,srcseq in enumerate(conversion_table_src):
            ref_dict[dest][srcseq]=ser[idx]
        ref_dict[dest]["-"]="-"
    
    return ref_dict


def SEQ2SEQ(seq,dest_seg,dic):
    if "-" in seq:
        return "-"
    else:
        return dic[dest_seg][seq]


def toQscore(query):
    baseQuality=round(sum(map(ord,query))/int(len(query)))
    return baseQuality


# Destination segment sequence generation process
def gen_dest_seq_process(subchunk_zip,exportReadStructure,read_now,settings,funcdict_key_list,func_dict,seq2seq_refs,referenceDict):
    d_val_chunk = subchunk_zip[0]
    d_qual_chunk = subchunk_zip[1]
    s_seq_chunk = subchunk_zip[2]
    s_qual_chunk = subchunk_zip[3]
    
    fastq_parse=pd.DataFrame()
    fastq_parse["Header"]=s_seq_chunk["Header"]
    fastq_parse["Third"]=["+"]*fastq_parse.shape[0]

    barcode_correspondence = pd.DataFrame()
    
    readStructure_now=exportReadStructure[read_now]
    for cnt_comp,component in enumerate(readStructure_now):
        if regex.search(r'^\"[^\"]+\"',component):
            # Implementation for constant sequences in the read structure configuration
            seq=regex.sub('\"',"",component)
            seq_export_tmp=[seq]*d_val_chunk.shape[0]
            length_now=len(seq)
            baseQuality=40 #hard coded
            const_quality=chr(baseQuality+33)*length_now
            qual_export_tmp=[const_quality]*d_val_chunk.shape[0]
            
        else:
            if component in settings.src_segments_available:
                # If source segment names are directly provided to the read structures
                s_seq_component_clean = settings.dict_to_terminal[component]
                func_tmp=settings.func_dict_ext[s_seq_component_clean]["func_ordered"][0]
                s_seq_component_raw=settings.func_dict_ext[s_seq_component_clean][func_tmp]["source"]
                s_seq_component=s_seq_component_raw+":"+s_seq_component_clean
                seq_export_tmp=s_seq_chunk[s_seq_component].apply(genEqSeq,datatype="seq")
                qual_export_tmp=s_qual_chunk[s_seq_component_raw].apply(genEqSeq,datatype="qual",baseQuality=None)
            
            else:
                # If the destination segment names are specified
                for k in funcdict_key_list:
                    if (component+"," in k) or (","+component in k) or (component==k):
                        funcdict_key=copy.copy(k)
                try:
                    funcdict_key
                except:
                    raise(component+" was not found in func_dict keys.")

                opt_now=copy.deepcopy(func_dict[funcdict_key])
                func_now=opt_now["func_ordered"][0]
                opt_now[func_now]["source"]="+".join(opt_now[func_now]["source"])

                if func_now=="PASS":
                    s_seq_component_clean=opt_now[func_now]["source"]
                    func_tmp=settings.func_dict_ext[s_seq_component_clean]["func_ordered"][0]
                    s_seq_component_raw=settings.func_dict_ext[s_seq_component_clean][func_tmp]["source"]
                    s_seq_component=s_seq_component_raw+":"+s_seq_component_clean
                    seq_export_tmp=s_seq_chunk[s_seq_component].apply(genEqSeq,length=opt_now[func_now].get("length"),datatype="seq",add_nuc=opt_now.get("add_nucleotide"))
                    qual_export_tmp=s_qual_chunk[s_seq_component_raw].apply(genEqSeq,length=opt_now[func_now].get("length"),datatype="qual",baseQuality=None)
                
                elif func_now=="SEQ2SEQ": #Fix this part next to apply a new dictionary format
                    s_seq_component_clean=opt_now[func_now]["source"]
                    s_seq_component_clean_list=sorted(s_seq_component_clean.split("+"))
                    func_tmp_list=[settings.func_dict_ext[i]["func_ordered"][0] for i in s_seq_component_clean_list]
                    s_seq_component_raw_list=[settings.func_dict_ext[i][func_tmp_list[idx]]["source"] for idx,i in enumerate(s_seq_component_clean_list)]
                    s_seq_component_list=[i+":"+s_seq_component_clean_list[idx] for idx,i in enumerate(s_seq_component_raw_list)]
                    seq_export_tmp=s_seq_chunk[s_seq_component_list]
                    qual_export_tmp=s_qual_chunk[s_seq_component_raw_list]

                    if len(seq_export_tmp.columns) > 1:
                        seq_export_tmp=seq_export_tmp[seq_export_tmp.columns[0]].str.cat(seq_export_tmp[seq_export_tmp.columns[1:]],sep="_")
                        qual_export_tmp=qual_export_tmp[qual_export_tmp.columns[0]].str.cat(qual_export_tmp[qual_export_tmp.columns[1:]])
                    else:
                        seq_export_tmp=seq_export_tmp[seq_export_tmp.columns[0]]
                        qual_export_tmp=qual_export_tmp[qual_export_tmp.columns[0]]

                    seq_export_tmp=seq_export_tmp.apply(SEQ2SEQ,dest_seg=component,dic=seq2seq_refs)
                    qual_export_tmp=qual_export_tmp.apply(toQscore).astype(str)
                    qual_export_tmp=seq_export_tmp.str.cat(qual_export_tmp,sep="_")
                    qual_export_tmp=qual_export_tmp.map(getConvQual_ver2)

                
                elif func_now=="WHITELIST_ASSIGNMENT" or func_now=="RANDSEQ_ASSIGNMENT":
                    d_val_component=opt_now[func_now]["source"]

                    reference_now=referenceDict[component]
                    d_val_chunk[funcdict_key]=d_val_chunk[funcdict_key].map(int)
                    d_qual_chunk[funcdict_key]=d_qual_chunk[funcdict_key].map(int)
                    seq_export_tmp=d_val_chunk[funcdict_key].apply(getConvSeq,reference=reference_now)
                    df_tmp_cat=seq_export_tmp.str.cat(d_qual_chunk[funcdict_key].astype(str),sep="_")
                    qual_export_tmp=df_tmp_cat.map(getConvQual_ver2)

                    if settings.is_barcodelist:
                        barcode_correspondence[component]=seq_export_tmp

                elif func_now=="CONSTANT":
                    seq_export_tmp=[opt_now[func_now]["sequence"]]*d_val_chunk.shape[0]
                    length_now=len(opt_now[func_now]["sequence"])
                    baseQuality=40 #hard coded
                    const_quality=chr(baseQuality+33)*length_now
                    qual_export_tmp=[const_quality]*d_val_chunk.shape[0]

        if cnt_comp==0:
            fastq_parse["seq"]=seq_export_tmp
            fastq_parse["qual"]=qual_export_tmp
        else:
            fastq_parse["seq"]=fastq_parse["seq"].str.cat(seq_export_tmp,sep="")
            fastq_parse["qual"]=fastq_parse["qual"].str.cat(qual_export_tmp,sep="")
        
        return fastq_parse,barcode_correspondence
    


# Multithreading implementation for generating a destination sequence table
def gen_dest_seq_parallel_wrapper(df_zip,exportReadStructure,read_now,settings,funcdict_key_list,func_dict,seq2seq_refs,referenceDict,barcode_correspondence,ncore):
    # Further split the data chunks into subchunks by the number of CPUs
    d_val_subchunks = np.array_split(df_zip[0],ncore)
    d_qual_subchunks = np.array_split(df_zip[1],ncore)
    s_seq_subchunks = np.array_split(df_zip[2],ncore)
    s_qual_subchunks = np.array_split(df_zip[3],ncore)

    retLst = Parallel(n_jobs=ncore,backend='multiprocessing',verbose=8)(
        delayed(gen_dest_seq_process)(subchunk_zip,exportReadStructure,read_now,settings,funcdict_key_list,func_dict,seq2seq_refs,referenceDict) for subchunk_zip in zip(d_val_subchunks,d_qual_subchunks,s_seq_subchunks,s_qual_subchunks))
    
    fastq_parse = [i[0] for i in retLst]
    barcode_correspondence_tmp = [i[1] for i in retLst]
    barcode_correspondence = pd.concat([barcode_correspondence,pd.concat(barcode_correspondence_tmp)], axis=1)

    return pd.concat(fastq_parse),barcode_correspondence


# Demultiplex to tsv file, process
# def demultiplex_tsv_process(s_seq_chunk,eachkey,key_series,settings):
    # if not "-" in eachkey:
    #     export_pd_tmp=s_seq_chunk[key_series==eachkey]
    #     outfilename="_".join([settings.outFilePath_and_Prefix,eachkey])+".tsv.gz"
    #     if not os.path.isfile(outfilename):
    #         export_pd_tmp.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=True)
    #     else:
    #         export_pd_tmp.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)

    #     return "_"+eachkey+".tsv.gz"
    # else:
    #     return ""
def demultiplex_tsv_process(sub_df,eachkey,settings):
    outfilename="_".join([settings.outFilePath_and_Prefix,eachkey])+".tsv.gz"
    if not os.path.isfile(outfilename):
        sub_df.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=True)
    else:
        sub_df.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)

    return "_"+eachkey+".tsv.gz"


# Demultiplex to FASTQ file, process
def demultiplex_fastq_process(export_pd,eachkey,settings,readIden):
    export_pd_tmp=export_pd[["Header","seq","3rd","qual"]]
    outfilename="_".join([settings.outFilePath_and_Prefix,eachkey,readIden])+".fastq.gz"
    export_pd_tmp=export_pd_tmp.stack()
    export_pd_tmp=export_pd_tmp.reset_index()
    export_pd_tmp=pd.DataFrame(export_pd_tmp[0])
    if not os.path.isfile(outfilename):
        export_pd_tmp.to_csv(outfilename,mode="w",compression="gzip",sep="\t",index=False,header=False)
    else:
        export_pd_tmp.to_csv(outfilename,mode="a",compression="gzip",sep="\t",index=False,header=False)
    return "_"+eachkey+"_"+readIden+".fastq.gz"


# Multithreading implementation for demultiplexing tsv file
def demultiplex_tsv_parallel_wrapper(s_seq_chunk,key_series,settings,ncore):
    s_seq_chunk = s_seq_chunk.copy()
    s_seq_chunk["demulti_keys"] = key_series
    
    # Group by demultiplex keys
    s_seq_chunk = s_seq_chunk.groupby("demulti_keys")
    key_iden_list_now = Parallel(n_jobs=ncore,require=None,verbose=3)(
        delayed(demultiplex_tsv_process)(sub_df,eachkey,settings) for eachkey,sub_df in s_seq_chunk if not "-" in eachkey)
    
    key_iden_list_now = [i for i in key_iden_list_now if not i == ""]
    return key_iden_list_now
    

# Multithreading implementation for demultiplexing tsv file
def demultiplex_fastq_parallel_wrapper(export_pd,key_series,settings,readIden,ncore):
    export_pd = export_pd.copy()
    export_pd["demulti_keys"] = key_series

    # Group by demultiplex keys
    export_pd = export_pd.groupby("demulti_keys")

    # Parallelize the df subsetting and exporting process
    key_iden_list_now = Parallel(n_jobs=ncore,require=None,verbose=3,backend="multiprocessing")(
        delayed(demultiplex_fastq_process)(sub_df,eachkey,settings,readIden) for eachkey,sub_df in export_pd)
    
    key_iden_list_now = [i for i in key_iden_list_now if not i == ""]
    return key_iden_list_now


# Generate shell coommands to concatenate the file chunks with the same key
def gen_cmd_tsv_concat(key,sampledir,target_files):
    key_gunzip=key.replace(".gz","")
    cmd1=["gunzip -c"]+[target_files[0]]+["| head -n1 >",sampledir+"/demultiplex/out/demultiplex.header"+key_gunzip]
    cmd1=" ".join(cmd1)
    cmd2=["echo"]+target_files+["| xargs cat | zgrep -v Header >",sampledir+"/demultiplex/out/demultiplex.content"+key_gunzip]
    cmd2=" ".join(cmd2)
    cmd3=["cat",sampledir+"/demultiplex/out/demultiplex.header"+key_gunzip,sampledir+"/demultiplex/out/demultiplex.content"+key_gunzip+" | gzip -c > ",sampledir+"/demultiplex/out/demultiplex"+key]
    cmd3=" ".join(cmd3)
    cmd4=["rm",sampledir+"/demultiplex/out/demultiplex.header"+key_gunzip,sampledir+"/demultiplex/out/demultiplex.content"+key_gunzip]
    cmd4=" ".join(cmd4)
    return [cmd1,cmd2,cmd3,cmd4]


# TSV file concatenation process
def tsv_concat_process(key,sampledir,demulti_filename_list):
    target_files=[t for t in demulti_filename_list if key in os.path.basename(t)]
    concat_commands = gen_cmd_tsv_concat(key,sampledir,target_files)
    cmd = ";".join(concat_commands)
    s=subprocess.run(cmd,shell=True)
    if s.returncode != 0:
        print("Job failed: Demultiuplexed file merge', file=sys.stderr")
        sys.exit(1)


# FASTQ file concatenation process
def fastq_concat_process(key,sampledir,demulti_filename_list):
    target_files=[t for t in demulti_filename_list if key in os.path.basename(t)]
    cmd=["echo"]+target_files+["| xargs cat >",sampledir+"/demultiplex/out/demultiplex"+key]
    cmd=" ".join(cmd)
    s=subprocess.run(cmd,shell=True)
    if s.returncode != 0:
        print("Job failed: Demultiuplexed file merge', file=sys.stderr")
        sys.exit(1)


# Multithreading implementation for concatenating demultiplexed tsv files
def concat_demultiplex_files_parallel_wrapper(key_set,sampledir,demulti_filename_list,ncore,format_type):
    if format_type=="tsv":
        Parallel(n_jobs=ncore,require=None,verbose=3)(delayed(tsv_concat_process)(key,sampledir,demulti_filename_list) for key in key_set)
    else:
        Parallel(n_jobs=ncore,require=None,verbose=3)(delayed(fastq_concat_process)(key,sampledir,demulti_filename_list) for key in key_set)
    