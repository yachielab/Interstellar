#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from argparse import ArgumentParser
import gzip
import re
import numpy as np
import pandas as pd

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-i', '--input', type=str,required=True,help="input fastq")
    argparser.add_argument('-q', '--baseQuality', type=int,default=33,help="base quality of generated fastq")
    argparser.add_argument('-o', '--outname', type=str,default="barista",help='output file name')
    argparser.add_argument('-d', '--outdir', type=str, default=".", help='output directory')
    return argparser.parse_args()

def fastq_generator(input):
    with gzip.open(input,mode="rt") as seqf:
        for c,i in enumerate(seqf):
            if c%4==0:
                yield i.replace("\n","")

def header_processor(fastq_iter,quality,outdir,outname):
    I1=[]
    for c,headline in enumerate(fastq_iter):
        headline=headline.split(" ")
        header=headline[0]
        barcodes=re.sub(":.+$","",headline[1])

        #header
        I1.append(header)
        #seq
        I1.append(barcodes)
        #3rd
        I1.append("+")
        #qual
        I1.append(chr(quality)*len(barcodes))

        if c%4000000==0:
            print(c,"reads have been processed.")
            I1="\n".join(I1)+"\n"
            if c==0:
                with gzip.open(outdir+"/"+outname+"_I1.fastq.gz",mode="wt") as wI1:
                    wI1.write(I1)
            else:
                with gzip.open(outdir+"/"+outname+"_I1.fastq.gz",mode="at") as wI1:
                    wI1.write(I1)
            I1=[]

    print(c,"reads have been processed.")
    I1="\n".join(I1)+"\n"
    with gzip.open(outdir+"/"+outname+"_I1.fastq.gz",mode="at") as wI1:
        wI1.write(I1)


if __name__ == "__main__":
    opt=get_option()
    input_fastq=opt.input
    quality=opt.baseQuality+30

    fastq_iter=fastq_generator(input_fastq)
    header_processor(fastq_iter,quality,opt.outdir,opt.outname)

