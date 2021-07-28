import collections
import regex
import itertools as it
import random
import pandas as pd
import numpy as np
import gzip,pickle
# random.seed(0)


def dval_to_sval_relationship(func_dict,dest_segment):
    dval_to_sval_relationship={}
    for val in dest_segment:
        opt_now=func_dict[val][func_dict[val]["func_ordered"][0]]
        s_values=opt_now["source"]
        
        s_value="+".join(s_values)
        dval_to_sval_relationship[val]=s_value

    return dval_to_sval_relationship


def parse_constraint(value_segment,values_in_destarg,child2parent_val,value_variables):
    noParentSegments=[]
    edge_dict=dict(parent=[],child=[])
    for component in child2parent_val:
        edge_dict["child"].append(component)
        edge_dict["parent"].append(child2parent_val[component])
    for val in value_segment:
        if not val in edge_dict["child"]:
            noParentSegments.append(component)

    roots=list(set(edge_dict["parent"])-set(edge_dict["child"]))
    # tips =list(set(edge_dict["child"])-set(edge_dict["parent"]))
    globalComponents=list(set(noParentSegments+values_in_destarg)-set(roots)-set(edge_dict["child"]))
    globalComponents=[i for i in globalComponents if i in value_variables]

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

# def compression_tree(s_val_chunk,component,edge_dict,roots,constraint,Tree):
#     for root_now in roots:
#         tasklist=[root_now]


#         while True:
#             if not tasklist:
#                 break

#             parent=tasklist[0]
#             tasklist=tasklist[1:]

#             childs=[edge_dict["child"][idx] for idx,i in enumerate(edge_dict["parent"]) if i==parent]
#             tasklist+=childs

#             if not childs:
#                 continue
                
#             for child in childs:
#                 ancestor_list = getAncestor(child,edge_dict)
#                 parent_series = s_val_chunk[ancestor_list[0]].str.cat(s_val_chunk[ancestor_list[1:]],sep="_")
#                 parent_series.name = parent            
            


# def buildConstraintTree_singleScan(constraintList,convertedDataframe,func_dict,s_val_chunk,constraintTree):
#     convertedDataframe=convertedDataframe.applymap(str)
#     for cnt,i in enumerate(constraintList):
#         localComponent=i[0]
#         constraintGroup=i[1:]
#         constraint="+".join(constraintGroup)
#         if set(constraintGroup) <= set(convertedDataframe.columns):
#             opt_now=func_dict[localComponent]
            
#             for c,s in enumerate(constraintGroup):
#                 if c==0:
#                     constraint_col_concat=convertedDataframe[s]
#                 else:
#                     constraint_col_concat=constraint_col_concat.str.cat(convertedDataframe[s],sep="+")
#             constraint_col_concat=tuple(constraint_col_concat)

#             if not constraintTree[constraint].get(localComponent):
#                 constraintTree[constraint][localComponent]=collections.defaultdict(collections.Counter)
                
#             if opt_now["is_combination"]:
#                 sources=[opt_now[i] for i in opt_now if "src_corrected_components" in i]
#                 local_src_val=s_val_chunk[sources[0]].str.cat(s_val_chunk[sources[1:]],sep="_")
                
#             else:
#                 sources=opt_now["src_corrected_components"]
#                 local_src_val=s_val_chunk[sources]
            
#             local_src_val=tuple(local_src_val)
#             for cnt,mother in enumerate(constraint_col_concat):
#                 constraintTree[constraint][localComponent][mother].update([local_src_val[cnt]])
            
#     return constraintTree

# def buildConstraintTree_v2(constraintList,convertedDataframe,func_dict,s_val_chunk,constraintTree):
#     processed=[]
#     convertedDataframe=convertedDataframe.applymap(str)
#     for cnt,i in enumerate(constraintList):
#         localComponent=i[0]
#         constraintGroup=i[1:]
#         constraint="+".join(constraintGroup)
#         if set(constraintGroup) <= set(convertedDataframe.columns):
#             processed.append(cnt)
#             opt_now=func_dict[localComponent]
            
#             for c,s in enumerate(constraintGroup):
#                 if c==0:
#                     constraint_col_concat=convertedDataframe[s]
#                 else:
#                     constraint_col_concat=constraint_col_concat.str.cat(convertedDataframe[s],sep="+")
#             constraint_col_concat=tuple(constraint_col_concat)

#             if not constraintTree[constraint].get(localComponent):
#                 constraintTree[constraint][localComponent]=collections.defaultdict(collections.Counter)
#                 # constraintTree[constraint][localComponent]=collections.defaultdict(list)

#             # constraintTree_tmp=collections.defaultdict(dict)
#             # constraintTree_tmp[constraint][localComponent]=collections.defaultdict(list)

#             if opt_now["is_combination"]:
#                 sources=[opt_now[i] for i in opt_now if "src_corrected_components" in i]
#                 local_src_val=s_val_chunk[sources[0]].str.cat(s_val_chunk[sources[1:]],sep="_")
#                 # for c,s in enumerate(sources):
#                 #     if c==0:
#                 #         local_src_val=s_val_chunk[s]
#                 #     else:
#                 #         local_src_val=local_src_val.str.cat(s_val_chunk[s],sep="_")
#             else:
#                 sources=opt_now["src_corrected_components"]
#                 local_src_val=s_val_chunk[sources]
            
#             # for mother in constraint_set:
#             #     pd_tmp=concat_pd.query("mother==@mother")
#             #     constraintTree[constraint][localComponent][mother]|=set(pd_tmp["daughter"])
#             local_src_val=tuple(local_src_val)
#             for cnt,mother in enumerate(constraint_col_concat):
#                 # constraintTree_tmp[constraint][localComponent][mother].append(local_src_val[cnt])
#                 constraintTree[constraint][localComponent][mother].update([local_src_val[cnt]])

#             # for mother in constraint_col_concat:
#             #     if not constraintTree[constraint][localComponent].get(mother):
#             #         constraintTree[constraint][localComponent][mother]=collections.Counter(constraintTree_tmp[constraint][localComponent][mother])
#             #     else:
#             #         constraintTree[constraint][localComponent][mother].update(constraintTree_tmp[constraint][localComponent][mother])

#     return processed,constraintTree
            



# def buildConstraintTree(constraintList,convertedDestination,func_dict,correctedIndexDict):
#     constraintTree=collections.defaultdict(dict)
#     processed=[]
#     for cnt,i in enumerate(constraintList):
#         localComponent=i[0]
#         constraintGroup=i[1:]
#         constraint="+".join(constraintGroup)
#         if set(constraintGroup) <= set(convertedDestination):
#             processed.append(cnt)
#             opt_now=func_dict[localComponent]
#             convertedIdxGroup=[convertedDestination[i] for i in constraintGroup]
#             convertedIdxGroup=list(map(list,zip(*convertedIdxGroup)))
#             convertedIdxGroup=["+".join(map(str,i)) for i in convertedIdxGroup]
#             constraintTree[constraint][localComponent]=collections.defaultdict(list)
#             if opt_now["is_combination"]:
#                 sources=[opt_now[i] for i in opt_now if "src_corrected_components" in i]
#                 sources_idx=[correctedIndexDict[i] for i in sources]
#                 sources_idx=list(map(list,zip(*sources_idx)))
#                 sources_idx=["_".join(map(str,i)) for i in sources_idx]
#                 for num,i in enumerate(convertedIdxGroup):
#                     constraintTree[constraint][localComponent][i].append(sources_idx[num])
#             else:
#                 sources=opt_now["src_corrected_components"]
#                 sources_idx=correctedIndexDict[sources]
#                 for num,i in enumerate(convertedIdxGroup):
#                     constraintTree[constraint][localComponent][i].append(sources_idx[num])
#             for i in convertedIdxGroup:
#                 constraintTree[constraint][localComponent][i]=collections.Counter([constraint][localComponent][i])

#     for rmv in processed:
#         del constraintList[rmv]  
#     return constraintList,constraintTree

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

def buildReference(opt_now,seg_idx,info_dic,max_default=7000000): #n_each?
    func_now=opt_now["func_ordered"][0]
    if func_now=="WHITELIST_ASSIGNMENT":
        with open(opt_now[func_now]["path"][seg_idx],mode="rt",encoding="utf-8") as f:
            reference_now=[regex.sub("\n","",i) for i in f]
            random.seed(0)
            random.shuffle(reference_now)
    elif func_now=="RANDSEQ_ASSIGNMENT":
        val_segment=opt_now[func_now]["source"]
        maxsize=info_dic[val_segment] 
        if maxsize < max_default:
            maxsize=max_default

        seqlength=int(opt_now[func_now]["length"][seg_idx])
        reference_iter=it.product(*(["ATGC"]*seqlength))
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
    # print("start")
    convMotherIdx=row[0]
    srcLocalIdx=row[1]
    # print(constraint,localComponent,convMotherIdx,srcLocalIdx)
    # print(row[1])
    if "-1" in convMotherIdx or "-1" in srcLocalIdx:
        return -1
    else:
        #print(constraint,localComponent,convMotherIdx,srcLocalIdx)
        #print(type(constraint),type(localComponent),type(convMotherIdx),type(srcLocalIdx))
        #print(tree)
        return tree["constraintTree"][constraint][localComponent][convMotherIdx][srcLocalIdx]

def genEqSeq(query,length,datatype,add_nuc=None,baseQuality=None):
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

# def getConvSeq(query,reference,reindex):
#     if reindex:
#         query_now=reindex[query]
#     else:
#         query_now=query
def getConvSeq(query,reference):
    if query >= len(reference) or query < 0:
        return None
    else:
        return reference[query]

# def getConvQual(query):
#     if query[0]:
#         length=len(query[0])
#         avg_qual=query[1]
#         return chr(avg_qual)*length
#     else:
#         return None

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

# def get_reindex(query,src_index_dict,reindex_dict):
#     if query not in src_index_dict:
#         return None
#     else:
#         if src_index_dict[query] not in reindex_dict:
#             #print(query,src_index_dict[query],"something wrong")
#             return None
#         else:
#             return reindex_dict[src_index_dict[query]]

def get_reindex(query,d_component,tree):
    reindex_dict=tree[d_component]
    if query not in reindex_dict:
        #print(query,src_index_dict[query],"something wrong")
        return None
    else:
        return reindex_dict[query]