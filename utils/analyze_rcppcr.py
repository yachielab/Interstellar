from argparse import ArgumentParser
import pandas as pd
import os
import re
import glob
import string
import collections
import pyper

class parseError(Exception):
    pass

def get_option():
    argparser = ArgumentParser()
    #argparser.add_argument('-libtype', '--libtype', type=str,required=True,choices=["AD","DB"],help="library type, AD or DB")
    argparser.add_argument('-TList', '--TList', type=str,default="",help="tsv files containing the position info")
    argparser.add_argument('-tag', '--tag_assign', type=str,default="",help="csv files containing the position info")
    argparser.add_argument('-bc', '--bc', type=str,required=True,help="path to a directory containing demultiplexed bc files")
    argparser.add_argument('-lox', '--lox', type=str,required=True,help="path to a directory containing demultiplexed lox files")
    argparser.add_argument('-d', '--outdir', type=str,default=".",help="path to demultiplexed lox files")
    argparser.add_argument('-o', '--outname', type=str,default="rcppcr",help="path to demultiplexed lox files")
    argparser.add_argument('-no_plot', '--no_plot', action="store_true",help="no plot")
    argparser.add_argument('-Dec2018', '--Dec2018', action="store_true",help="format of Dec2018")
    argparser.add_argument('-ADDB', '--ADDB', default="",help="AD or DB")
    argparser.add_argument('-thresh_bc1', '--thresh_bc1', type=float,default=0.7,help="path to demultiplexed lox files")
    argparser.add_argument('-thresh_bc2', '--thresh_bc2', type=float,default=0.7,help="path to demultiplexed lox files")
    argparser.add_argument('-thresh_loxp', '--thresh_loxp', type=float,default=0.7,help="path to demultiplexed lox files")
    argparser.add_argument('-thresh_lox2272', '--thresh_lox2272', type=float,default=0.7,help="path to demultiplexed lox files")
    argparser.add_argument('-thresh_lox_bc1', '--thresh_lox_bc1', type=float,default=0.5,help="path to demultiplexed lox files")
    argparser.add_argument('-thresh_lox_bc2', '--thresh_lox_bc2', type=float,default=0.5,help="path to demultiplexed lox files")
    return argparser.parse_args()

def searchDirectory(bc_dir,lox_dir):
    bc_files=glob.glob(bc_dir)
    lox_files=glob.glob(lox_dir)

    if not bc_files or not lox_dir:
        raise parseError("Direcotry does not exist or empty")

def get_PRC(series_now):
    sample_name=series_now["Type"]
    plate_name=series_now["GS_plate_384"]
    plate_index=series_now["Plate_PCR_index"]
    row_column_index_set=series_now["RC_PCR_index"]
    position=series_now["Position"]
    
    plate=re.search(r"^P([^\-]+)\-.+",plate_index)
    plate=int(plate.group(1))
    if row_column_index_set=="BL":
        row_set=[2*(i+1) for i in range(8)]
    elif row_column_index_set=="TL":
        row_set=[2*i+1 for i in range(8)]
    col_set=[2*i+1 for i in range(12)]

    PRC_set=[]
    for row in row_set:
        for col in col_set:
            PRC_set.append("plate"+str(plate)+"_row"+str(row)+"_col"+str(col)+".tsv")

    return sample_name,plate_name,position,row_column_index_set,PRC_set


def get_PRC_2018(series_now,addb,row_column_index_set):
    sample_name=series_now["sample"]
    plate_name=series_now["sample"].split("-")[0]
    plate_index=series_now["plate"]
    # if addb =="AD":
    #     row_column_index_set="TL" #Row=Odd / Col=Odd in AD
    # else:
    #     row_column_index_set="BR" #Row=Odd / Col=Odd in DB
    position=series_now["sample"].split("-")[1]
    
    print(row_column_index_set)
    plate=re.search(r"^P([^\-]+)\-.+",plate_index)
    plate=int(plate.group(1))
    if row_column_index_set=="BR":
        row_set=[2*(i+1) for i in range(8)]
        col_set=[2*(i+1) for i in range(12)]
    elif row_column_index_set=="TL":
        row_set=[2*i+1 for i in range(8)]
        col_set=[2*i+1 for i in range(12)]

    PRC_set=[]
    for row in row_set:
        for col in col_set:
            PRC_set.append("plate"+str(plate)+"_row"+str(row)+"_col"+str(col)+".tsv")

    return sample_name,plate_name,position,row_column_index_set,PRC_set


def get_384_position(position,row_idx,col_idx,rowset_384,colset_384,row_column_index_set):
    if row_column_index_set=="BL":
        row_384_idx=row_idx/2-1
        col_384_idx=(col_idx-1)/2
    elif row_column_index_set=="BR":
        row_384_idx=row_idx/2-1
        col_384_idx=col_idx/2-1
    elif row_column_index_set=="TL":
        row_384_idx=(row_idx-1)/2
        col_384_idx=(col_idx-1)/2
    elif row_column_index_set=="TR":
        row_384_idx=(row_idx-1)/2
        col_384_idx=col_idx/2-1

    print(row_column_index_set,row_idx,col_idx,position)

    if position=="TL":
        row_384_idx=2*row_384_idx
        col_384_idx=2*col_384_idx
    elif position=="TR":
        row_384_idx=2*row_384_idx
        col_384_idx=2*col_384_idx+1
    elif position=="BL":
        row_384_idx=2*row_384_idx+1
        col_384_idx=2*col_384_idx
    elif position=="BR":
        row_384_idx=2*row_384_idx+1
        col_384_idx=2*col_384_idx+1

    row_384=rowset_384[int(row_384_idx)]
    col_384=colset_384[int(col_384_idx)]

    return row_384,col_384

def counter_parse(dict_counter):
    counter_sorted=sorted(dict_counter.items(),key=lambda x:x[1],reverse=True)
    rank,bc,cnt=list(),list(),list()
    for n,i in enumerate(counter_sorted):
        rank.append(n+1)
        bc.append(i[0])
        cnt.append(i[1])
    counter_df=pd.DataFrame({"Rank":rank,"BC":bc,"Count":cnt})
    return counter_df

def get_bc_profile(bc_df,out_figure_prefix,no_plot,thresh_bc1,thresh_bc2):
    bc_qcdict={}
    bc1_counter=collections.Counter(bc_df["s_bc1"])
    bc1_counter_df=counter_parse(bc1_counter)
    bc_qcdict["BC1"]=bc1_counter_df["BC"][0]
    bc_qcdict["Num_dominant_bc1"]=bc1_counter_df["Count"][0]
    bc_qcdict["Ratio_dominant_bc1"]=bc1_counter_df["Count"][0]/sum(bc1_counter_df["Count"])
    if bc_qcdict["Ratio_dominant_bc1"]>thresh_bc1:
        bc_qcdict["Quality_filter_bc1"]="Pass"
    else:
        bc_qcdict["Quality_filter_bc1"]="Drop"
    
    bc2_counter=collections.Counter(bc_df["s_bc2"])
    bc2_counter_df=counter_parse(bc2_counter)
    bc_qcdict["BC2"]=bc2_counter_df["BC"][0]
    bc_qcdict["Num_dominant_bc2"]=bc2_counter_df["Count"][0]
    bc_qcdict["Ratio_dominant_bc2"]=bc2_counter_df["Count"][0]/sum(bc2_counter_df["Count"])
    if bc_qcdict["Ratio_dominant_bc2"]>thresh_bc2:
        bc_qcdict["Quality_filter_bc2"]="Pass"
    else:
        bc_qcdict["Quality_filter_bc2"]="Drop"

    if not no_plot:
        r = pyper.R(use_pandas='True')
        r.assign("bc1_df", bc1_counter_df)
        r.assign("bc2_df", bc2_counter_df)
        r.assign("out_figure_prefix",out_figure_prefix)
        r("source(file='"+os.path.dirname(__file__)+"/plot_rcppcr.R')")
        
    return bc_qcdict

def get_lox_profile(lox_df,out_figure_prefix,thresh_bc1,thresh_bc2,thresh_loxp,thresh_lox2272):
    lox_qcdict={}
    bc_qcdict={}
    lox_qcdict["Num_lox_reads"]=lox_df.shape[0]

    loxP_counter=collections.Counter(lox_df["s_loxp"])
    loxP_counter_df=counter_parse(loxP_counter)
    lox_qcdict["LoxP"]=loxP_counter_df["BC"][0]
    lox_qcdict["Num_dominant_LoxP"]=loxP_counter_df["Count"][0]
    lox_qcdict["Ratio_dominant_LoxP"]=loxP_counter_df["Count"][0]/sum(loxP_counter_df["Count"])
    if lox_qcdict["Ratio_dominant_LoxP"]>thresh_loxp and lox_qcdict["LoxP"]=="ATAACTTCGTATAATGTATGCTATACGAAGTTAT":
        lox_qcdict["Quality_filter_LoxP"]="Pass"
    else:
        lox_qcdict["Quality_filter_LoxP"]="Drop"
    
    lox2272_counter=collections.Counter(lox_df["s_lox2272"])
    lox2272_counter_df=counter_parse(lox2272_counter)
    print(lox2272_counter_df)
    lox_qcdict["Lox2272"]=lox2272_counter_df["BC"][0]
    lox_qcdict["Num_dominant_Lox2272"]=lox2272_counter_df["Count"][0]
    lox_qcdict["Ratio_dominant_Lox2272"]=lox2272_counter_df["Count"][0]/sum(lox2272_counter_df["Count"])
    if lox_qcdict["Ratio_dominant_Lox2272"]>thresh_lox2272 and lox_qcdict["Lox2272"]=="ATAACTTCGTATAGGATACTTTATACGAAGTTAT":
        lox_qcdict["Quality_filter_Lox2272"]="Pass"
    else:
        lox_qcdict["Quality_filter_Lox2272"]="Drop"

    bc1_counter=collections.Counter(lox_df["s_bc1"])
    bc1_counter_df=counter_parse(bc1_counter)
    bc_qcdict["Ratio_dominant_bc1"]=bc1_counter_df["Count"][0]/sum(bc1_counter_df["Count"])
    if bc_qcdict["Ratio_dominant_bc1"]>thresh_bc1:
        bc_qcdict["BC1"]=bc1_counter_df["BC"][0]
    else:
        bc_qcdict["BC1"]="-"
    
    bc2_counter=collections.Counter(lox_df["s_bc2"])
    bc2_counter_df=counter_parse(bc2_counter)
    bc_qcdict["Num_dominant_bc2"]=bc2_counter_df["Count"][0]
    bc_qcdict["Ratio_dominant_bc2"]=bc2_counter_df["Count"][0]/sum(bc2_counter_df["Count"])
    if bc_qcdict["Ratio_dominant_bc2"]>thresh_bc2:
        bc_qcdict["BC2"]=bc2_counter_df["BC"][0]
    else:
        bc_qcdict["BC2"]="-"

    return bc_qcdict,lox_qcdict

def process_bc(TList,res_df,no_plot,bc_dir,out_prefix,thresh_bc1,thresh_bc2):
    TList_bc = TList.query("RCP_PCR_type=='BC'")
    pat_parse_filename=re.compile(r"^.+plate([^_]+)_row([^_]+)_col([^\.]+)\.tsv")
    rowset_384=[i for i in string.ascii_uppercase[:16]]
    colset_384=[i+1 for i in range(24)]
    for ind in TList_bc.index:
        series_now=TList_bc.iloc[ind]
        sample_name,plate_name,position,row_column_index_set,PRC_set=get_PRC(series_now)
        print("Processing plate index="+str(ind+1)+" for BC data...")
        for well in PRC_set:
            res_dict={}

            bc_file=glob.glob(bc_dir+"/*"+well)
            if not bc_file:
                continue
            m=pat_parse_filename.search(bc_file[0])
            original_idx="plate"+m.group(1)+"_row"+m.group(2)+"_col"+m.group(3)
            row_idx=int(m.group(2))
            col_idx=int(m.group(3))
            bc_df=pd.read_csv(bc_file[0],sep="\t")

            row_384,col_384=get_384_position(position,row_idx,col_idx,rowset_384,colset_384,row_column_index_set)
            pos_384=str(row_384)+str(col_384)

            res_dict["Original_index"]=original_idx
            res_dict["Sample_name"]=sample_name
            res_dict["Plate_name"]=plate_name
            res_dict["Position_384"]=pos_384
            res_dict["Row"]=row_384
            res_dict["Column"]=col_384
            res_dict["Num_bc_reads"]=bc_df.shape[0]

            out_figure_prefix=out_prefix+sample_name+"_"+plate_name+"_"+pos_384
            bc_qcdict=get_bc_profile(bc_df,out_figure_prefix,no_plot,thresh_bc1,thresh_bc2)
            res_dict.update(bc_qcdict)

            res_df.loc[sample_name+"_"+plate_name+"_"+pos_384]="-"
            for key in res_dict:
                res_df.loc[sample_name+"_"+plate_name+"_"+pos_384][key]=res_dict[key]

    return res_df

def process_bc_Dec2018(tag_assignment_df,addb,res_df,no_plot,bc_dir,out_prefix,thresh_bc1,thresh_bc2):
    tag_assignment_df = tag_assignment_df.query("ADDB==@addb")
    tag_assignment_df = tag_assignment_df.query("type=='BC'")
    tag_assignment_df = tag_assignment_df.reset_index()

    if addb=="AD": #hard coding
        row_column_index_set="TL"
    else:
        row_column_index_set="BR"
    
    pat_parse_filename=re.compile(r"^.+plate([^_]+)_row([^_]+)_col([^\.]+)\.tsv")
    rowset_384=[i for i in string.ascii_uppercase[:16]]
    colset_384=[i+1 for i in range(24)]
    for ind in tag_assignment_df.index:
        series_now=tag_assignment_df.iloc[ind]
        sample_name,plate_name,position,row_column_index_set,PRC_set=get_PRC_2018(series_now,addb,row_column_index_set)

        #hard coding
        print("Processing plate index="+str(ind+1)+" for BC data...")
        for well in PRC_set:
            res_dict={}

            bc_file=glob.glob(bc_dir+"/*"+well)
            if not bc_file:
                continue


            m=pat_parse_filename.search(bc_file[0])
            original_idx="plate"+m.group(1)+"_row"+m.group(2)+"_col"+m.group(3)
            row_idx=int(m.group(2))
            col_idx=int(m.group(3))
            bc_df=pd.read_csv(bc_file[0],sep="\t")

            if bc_df.shape[0]==0:
                continue


            row_384,col_384=get_384_position(position,row_idx,col_idx,rowset_384,colset_384,row_column_index_set)
            pos_384=str(row_384)+str(col_384)

            res_dict["Original_index"]=original_idx
            res_dict["Sample_name"]=sample_name
            res_dict["Plate_name"]=plate_name
            res_dict["Position_384"]=pos_384
            res_dict["Row"]=row_384
            res_dict["Column"]=col_384
            res_dict["Num_bc_reads"]=bc_df.shape[0]

            out_figure_prefix=out_prefix+sample_name+"_"+plate_name+"_"+pos_384
            bc_qcdict=get_bc_profile(bc_df,out_figure_prefix,no_plot,thresh_bc1,thresh_bc2)
            res_dict.update(bc_qcdict)

            res_df.loc[sample_name+"_"+plate_name+"_"+pos_384]="-"
            for key in res_dict:
                res_df.loc[sample_name+"_"+plate_name+"_"+pos_384][key]=res_dict[key]

    return res_df

def process_lox(TList,res_df,lox_dir,out_prefix,thresh_lox_bc1,thresh_lox_bc2,thresh_loxp,thresh_lox2272):
    TList_lox = TList.query("RCP_PCR_type=='Lox'")
    TList_lox.index=range(TList_lox.shape[0])
    pat_parse_filename=re.compile(r"^.+plate([^_]+)_row([^_]+)_col([^\.]+)\.tsv")
    rowset_384=[i for i in string.ascii_uppercase[:16]]
    colset_384=[i+1 for i in range(24)]
    for ind in TList_lox.index:
        series_now=TList_lox.iloc[ind]
        sample_name,plate_name,position,row_column_index_set,PRC_set=get_PRC(series_now)
        print("Processing plate index="+str(ind+1)+" for Lox data...")
        
        for well in PRC_set:
            lox_file=glob.glob(lox_dir+"/*"+well)
            if not lox_file:
                continue
            m=pat_parse_filename.search(lox_file[0])
            plate_idx=m.group(1)
            row_idx=int(m.group(2))
            col_idx=int(m.group(3))
            lox_df=pd.read_csv(lox_file[0],sep="\t")

            row_384,col_384=get_384_position(position,row_idx,col_idx,rowset_384,colset_384,row_column_index_set)
            pos_384=str(row_384)+str(col_384)

            out_figure_prefix=out_prefix+sample_name+"_"+plate_name+"_"+pos_384
            bc_qcdict,lox_qcdict=get_lox_profile(lox_df,out_figure_prefix,thresh_bc1=thresh_lox_bc1,thresh_bc2=thresh_lox_bc2,thresh_loxp=thresh_loxp,thresh_lox2272=thresh_lox2272)
            
            #QC barcode correspondence
            if bc_qcdict["BC1"]==bc_qcdict["BC2"]:
                lox_qcdict["Correspondence_BC_Lox"]="Drop"
            else:
                barcode_name="BC1" if not bc_qcdict["BC1"]=="-" else "BC2"
                if bc_qcdict[barcode_name]==res_df.loc[sample_name+"_"+plate_name+"_"+pos_384][barcode_name]:
                    lox_qcdict["Correspondence_BC_Lox"]="Pass"
                else:
                    lox_qcdict["Correspondence_BC_Lox"]="Drop"
            
            #QC total
            if lox_qcdict["Correspondence_BC_Lox"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif lox_qcdict["Quality_filter_LoxP"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif lox_qcdict["Quality_filter_Lox2272"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif res_df.loc[sample_name+"_"+plate_name+"_"+pos_384]["Quality_filter_bc1"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif res_df.loc[sample_name+"_"+plate_name+"_"+pos_384]["Quality_filter_bc2"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            else:
                lox_qcdict["Quality_filter_total"]="Pass"
            
            #Store in the data frame
            for key in lox_qcdict:
                res_df.loc[sample_name+"_"+plate_name+"_"+pos_384][key]=lox_qcdict[key]

    return res_df

def process_lox_Dec2018(tag_assignment_df,addb,res_df,lox_dir,out_prefix,thresh_lox_bc1,thresh_lox_bc2,thresh_loxp,thresh_lox2272):
    tag_assignment_df = tag_assignment_df.query("ADDB==@addb")
    tag_assignment_df = tag_assignment_df.query("type=='lox'")
    tag_assignment_df = tag_assignment_df.reset_index()
    tag_assignment_df.index=range(tag_assignment_df.shape[0])

    if addb=="AD": #hard coding
        row_column_index_set="TL"
    else:
        row_column_index_set="TL"

    
    pat_parse_filename=re.compile(r"^.+plate([^_]+)_row([^_]+)_col([^\.]+)\.tsv")
    rowset_384=[i for i in string.ascii_uppercase[:16]]
    colset_384=[i+1 for i in range(24)]
    for ind in tag_assignment_df.index:
        series_now=tag_assignment_df.iloc[ind]
        sample_name,plate_name,position,row_column_index_set,PRC_set=get_PRC_2018(series_now,addb,row_column_index_set)
        print("Processing plate index="+str(ind+1)+" for Lox data...")
        
        for well in PRC_set:
            lox_file=glob.glob(lox_dir+"/*"+well)
            if not lox_file:
                continue
            m=pat_parse_filename.search(lox_file[0])
            plate_idx=m.group(1)
            row_idx=int(m.group(2))
            col_idx=int(m.group(3))
            lox_df=pd.read_csv(lox_file[0],sep="\t")

            if lox_df.shape[0]==0:
                continue

            row_384,col_384=get_384_position(position,row_idx,col_idx,rowset_384,colset_384,row_column_index_set)
            pos_384=str(row_384)+str(col_384)

            if not sample_name+"_"+plate_name+"_"+pos_384 in res_df.index:
                continue

            out_figure_prefix=out_prefix+sample_name+"_"+plate_name+"_"+pos_384
            bc_qcdict,lox_qcdict=get_lox_profile(lox_df,out_figure_prefix,thresh_bc1=thresh_lox_bc1,thresh_bc2=thresh_lox_bc2,thresh_loxp=thresh_loxp,thresh_lox2272=thresh_lox2272)
            

            #QC barcode correspondence
            if bc_qcdict["BC1"]==bc_qcdict["BC2"]:
                lox_qcdict["Correspondence_BC_Lox"]="Drop"
            else:
                barcode_name="BC1" if not bc_qcdict["BC1"]=="-" else "BC2"

                if bc_qcdict[barcode_name]==res_df.loc[sample_name+"_"+plate_name+"_"+pos_384][barcode_name]:
                    lox_qcdict["Correspondence_BC_Lox"]="Pass"
                else:
                    lox_qcdict["Correspondence_BC_Lox"]="Drop"
            
            #QC total
            if lox_qcdict["Correspondence_BC_Lox"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif lox_qcdict["Quality_filter_LoxP"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif lox_qcdict["Quality_filter_Lox2272"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif res_df.loc[sample_name+"_"+plate_name+"_"+pos_384]["Quality_filter_bc1"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            elif res_df.loc[sample_name+"_"+plate_name+"_"+pos_384]["Quality_filter_bc2"]=="Drop":
                lox_qcdict["Quality_filter_total"]="Drop"
            else:
                lox_qcdict["Quality_filter_total"]="Pass"
            
            #Store in the data frame
            for key in lox_qcdict:
                res_df.loc[sample_name+"_"+plate_name+"_"+pos_384][key]=lox_qcdict[key]

    return res_df

if __name__ == "__main__":
    opt=get_option()
    TList=opt.TList
    tag_assign=opt.tag_assign
    bc_dir=opt.bc
    lox_dir=opt.lox
    outdir=opt.outdir
    no_plot=opt.no_plot
    outname=opt.outname
    out_prefix=re.sub("/$","",outdir)+"/"+outname
    thresh_bc1=opt.thresh_bc1
    thresh_bc2=opt.thresh_bc2
    thresh_loxp=opt.thresh_loxp
    thresh_lox2272=opt.thresh_lox2272
    thresh_lox_bc1=opt.thresh_lox_bc1
    thresh_lox_bc2=opt.thresh_lox_bc2
    

    
    #start proessing for bc
    res_df=pd.DataFrame({"Sample_name":[], #ok
                         "Original_index":[],
                         "Plate_name":[], #ok
                         "Position_384":[], #ok
                         "Row":[], #ok
                         "Column":[], #ok
                         "Quality_filter_total":[], 
                         "Num_bc_reads":[], #ok
                         "BC1":[], #ok
                         "Num_dominant_bc1":[], #ok
                         "Ratio_dominant_bc1":[], #ok
                         "Quality_filter_bc1":[], #ok
                         "BC2":[], #ok
                         "Num_dominant_bc2":[], #ok
                         "Ratio_dominant_bc2":[], #ok
                         "Quality_filter_bc2":[], #ok
                         "Correspondence_BC_Lox":[],
                         "Num_lox_reads":[], #ok
                         "LoxP":[],
                         "Num_dominant_LoxP":[],
                         "Ratio_dominant_LoxP":[],
                         "Quality_filter_LoxP":[],
                         "Lox2272":[],
                         "Num_dominant_Lox2272":[],
                         "Ratio_dominant_Lox2272":[],
                         "Quality_filter_Lox2272":[]})
    
    if not opt.Dec2018:
        TList=pd.read_csv(TList,sep='\t')

        print("Processing BC-RCP-PCR...")
        res_df=process_bc(TList,res_df,no_plot,bc_dir,out_prefix,thresh_bc1,thresh_bc2)

        print("\nProcessing Lox-RCP-PCR...")
        res_df=process_lox(TList,res_df,lox_dir,out_prefix,thresh_lox_bc1,thresh_lox_bc2,thresh_loxp,thresh_lox2272)
    else:
        addb=opt.ADDB
        if not addb in ["AD","DB"]:
            raise parseError()
        tag_assignment_df = pd.read_csv(tag_assign,header=None)
        tag_assignment_df.columns = ["sample","ADDB","type","plate"]
        
        print("Processing BC-RCP-PCR...")
        res_df=process_bc_Dec2018(tag_assignment_df,addb,res_df,no_plot,bc_dir,out_prefix,thresh_bc1,thresh_bc2)

        print("\nProcessing Lox-RCP-PCR...")
        res_df=process_lox_Dec2018(tag_assignment_df,addb,res_df,lox_dir,out_prefix,thresh_lox_bc1,thresh_lox_bc2,thresh_loxp,thresh_lox2272)

    res_df=res_df.sort_values(["Column","Row"])
    res_df.to_csv(out_prefix+"_summary.tsv",sep="\t",index=False)