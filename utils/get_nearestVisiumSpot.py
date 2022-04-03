#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import numpy as np
import pandas as pd
import math
import itertools as it

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-s_loc', '--srcLocation', type=str,required=True,help="source coordinates")
    argparser.add_argument('-exp_rate', '--expansionRate', type=int,default=1,help="source coordinates")
    argparser.add_argument('-pigeonhole', '--pigeonhole', action="store_true",help="all puck are asigned to different spot")
    argparser.add_argument('-exp_s_loc', '--expandedSrcLocation', default="", type=str,help="expanded source coordinates")
    argparser.add_argument('-d_loc', '--destLocation', type=str,required=True,help='destination coordinates')
    argparser.add_argument('-o', '--outname', type=str,default="barista",help='output file name')
    argparser.add_argument('-d', '--outdir', type=str, default=".", help='output directory')
    return argparser.parse_args()

def coord_expand(s_loc,d_loc,outdir,expansionRate,pigeonhole,outname):
    #s_loc, d_loc: pandas dataframe
    if pigeonhole:
        d_diffs = np.expand_dims(d_loc, axis=1) - np.expand_dims(d_loc, axis=0)
        d_loc_dist = np.sqrt(np.sum(d_diffs ** 2, axis=-1))
        d_loc_dist=list(set(d_loc_dist.flatten()))
        del d_loc_dist[d_loc_dist.index(0)]
        d_loc_dist=sorted(d_loc_dist)
        L=[]
        for i in d_loc_dist:
            if i<min(d_loc_dist)*1.5:
                L.append(i)
            else:
                break
        threshold=math.ceil(math.ceil(max(L))*2/np.sqrt(3))

        s_diffs=np.expand_dims(s_loc, axis=1) - np.expand_dims(s_loc, axis=0)
        s_loc_dist=np.sqrt(np.sum(s_diffs ** 2, axis=-1))
        s_loc_dist_min=None
        for c,ar in enumerate(s_loc_dist):
            ar=np.delete(ar,np.where(ar==0))
            if c==0:
                s_loc_dist_min=min(ar)
            elif min(ar)<s_loc_dist_min:
                s_loc_dist_min=min(ar)
        
        expansion_rate=threshold/s_loc_dist_min
    else:
        expansion_rate=expansionRate

    s_loc=s_loc*expansion_rate
    
    #expanded coordinate of source data
    s_loc.to_csv(outdir+"/"+outname+"_expanded.tsv",sep="\t")
    return s_loc

if __name__ == "__main__":
    opt=get_option()
    s_loc=pd.read_csv(opt.srcLocation,sep=",|\t",engine="python")
    d_loc=pd.read_csv(opt.destLocation,sep=",|\t",engine="python",header=None)
    
    s_loc.columns=["barcode","x","y"]
    d_loc.columns=["barcode","x","y"]

    s_loc.index=s_loc["barcode"]
    s_loc=s_loc[["x","y"]]
    d_loc.index=d_loc["barcode"]
    d_loc=d_loc[["x","y"]]
    
    outdir=opt.outdir
    outname=opt.outname
    if not opt.expandedSrcLocation:
        print("start expanding the source corrdinates...")
        new_coord=coord_expand(s_loc,d_loc,outdir,opt.expansionRate,opt.pigeonhole,outname)
        new_coord["barcode"]=new_coord.index
    else:
        new_coord=pd.read_csv(opt.expandedSrcLocation,sep="\t")
    
    print("start barcode alignment...")
    delta_x=min(d_loc["x"])-min(new_coord["x"])
    delta_y=min(d_loc["y"])-min(new_coord["y"])
    
    #align left bottom
    d_loc["x"]=d_loc["x"]-delta_x
    d_loc["y"]=d_loc["y"]-delta_y
    d_loc.to_csv(outdir+"/"+outname+"_moved_destLocation.tsv",sep="\t")
    
    #calculate section index for each spot in source
    x_d_index=((new_coord["x"]-min(new_coord["x"]))//(max(d_loc["x"])-min(d_loc["x"]))).astype(int).astype(str)
    y_d_index=((new_coord["y"]-min(new_coord["y"]))//(max(d_loc["y"])-min(d_loc["y"]))).astype(int).astype(str)
    d_index=x_d_index.str.cat(y_d_index,sep="_")

    #unit size for x and y
    x_unit=max(d_loc["x"])-min(d_loc["x"])
    y_unit=max(d_loc["y"])-min(d_loc["y"])
    print("unit size x:",x_unit)
    print("unit size y:",y_unit)

    d_loc_base=pd.DataFrame()
    dloc_ymin_pls5=min(d_loc["y"])+5
    df_tmp=d_loc.query("y<@dloc_ymin_pls5")
    blank=(sorted(df_tmp["x"])[1]-sorted(df_tmp["x"])[0])/2
    d_loc_base["x"]=d_loc["x"]-min(d_loc["x"])+blank
    d_loc_base["y"]=d_loc["y"]-min(d_loc["y"])+blank

    d_loc_new=pd.DataFrame()
    new_loc=pd.DataFrame()
    for i in range(max(list(x_d_index.astype(int)))+1):
        for j in range(max(list(y_d_index.astype(int)))+1):
            tmp=pd.DataFrame()
            if i == 0:
                tmp["x"]=d_loc["x"]
            else:
                tmp["x"]= min(d_loc["x"])+x_unit*i+d_loc_base["x"]

            if j == 0:
                tmp["y"]=d_loc["y"]
            else:
                tmp["y"]= min(d_loc["y"])+y_unit*j+d_loc_base["y"]

            tmp["ind"]=str(i)+"_"+str(j)
            
            new_loc=pd.concat([new_loc,tmp],axis=0)
    new_loc.to_csv(outdir+"/"+outname+"_tiled_visium.tsv",sep="\t",header=False,index=False)

    def find_nearest(coord):
        x=float(coord.split("_")[0])
        y=float(coord.split("_")[1])
        index_x=int(coord.split("_")[2])
        index_y=int(coord.split("_")[3])

        # d_loc_new["x"]=min(d_loc["x"])+x_unit*index_x+d_loc_base["x"]
        # d_loc_new["y"]=min(d_loc["y"])+y_unit*index_y+d_loc_base["y"]
        if index_x==0:
            d_loc_new["x"]=d_loc["x"]
        else:
            # d_loc_new["x"]=d_loc["x"]+x_unit*index_x+d_loc_base["x"] #orifinal coord(left bottom aligned) + 
            d_loc_new["x"]=min(d_loc["x"])+x_unit*index_x+d_loc_base["x"]

        if index_y==0:
            d_loc_new["y"]=d_loc["y"]
        else:
            # d_loc_new["y"]=d_loc["y"]+y_unit*index_y+d_loc_base["y"]
            d_loc_new["y"]=min(d_loc["y"])+y_unit*index_y+d_loc_base["y"]
        # d_loc_new["x"]=d_loc["x"]+(max(d_loc["x"])*index_x)
        # d_loc_new["y"]=d_loc["y"]+(max(d_loc["y"])*index_y)

        np_d_loc=np.array(d_loc_new)
        idx=(((np_d_loc[:, 0] - x)**2 + (np_d_loc[:, 1] - y)**2)**0.5).argmin()
        # print([list(d_loc_new.index)[idx]]+list(np_d_loc[idx].astype(str)))
        # print(d_loc_new.iloc[idx])
        return "_".join([list(d_loc_new.index)[idx]]+list(np_d_loc[idx].astype(str))+[str(x),str(y)])

    new_coord_str=new_coord.astype(str)
    new_coord_series=new_coord_str["x"].str.cat(new_coord_str["y"],sep="_")
    new_coord_series=new_coord_series.str.cat(d_index,sep="_")
    
    print("start finding nearest spot...")
    aligned_coord=new_coord_series.map(find_nearest)
    aligned_coord=pd.DataFrame(aligned_coord.str.split("_",expand=True))
    aligned_coord["d_index"]=d_index

    aligned_coord.columns=["dest_barcode","x","y","orig_x","orig_y","d_index"]
    aligned_coord.index=new_coord["barcode"]
    aligned_coord["d_ind_bc"]=aligned_coord["d_index"].str.cat(aligned_coord["dest_barcode"],sep="_")

    print("exporting...")
    aligned_coord.to_csv(outdir+"/"+outname+"_aligned.tsv",sep="\t")