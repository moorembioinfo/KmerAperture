#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import argparse
import os, sys

def add_args(a):
    parser = argparse.ArgumentParser(description="Precluster")
    parser.add_argument(
        "--fastas",
        help="Path to fasta dir ",
        required=True,
    )
    parser.add_argument(
        "--threshold",
        help="Cluster based on this level of truncation ",
        required=False,
    )
    args = parser.parse_args(a)
    return args

def precluster(filelist):
    print('Running sourmash sketch...\n')
    scmd = f'sourmash sketch dna -p k=31 {" ".join(filelist)}'
    os.system(scmd)
    #print(scmd)
    print('Running sourmash compare...\n')
    ccmd = 'sourmash compare --csv precluster_sm.csv *sig'
    os.system(ccmd)
    #os.remove('*sig')

def dendro(threshold, ngenomes):

    if threshold:
        #fh = open('fulllinkagearray.txt')
        #Z = fh.readlines()[0].rstrip()
        
        df=pd.read_csv('precluster_sm.csv')
        df.index = df.columns
        Z = linkage(df, 'average')

        #fig, ax = plt.subplots(figsize=(20, 10))
        #plt.title(f'Full dendrogram')
        #plt.ylabel('Distance (UPGMA)')
        #plt = dendrogram(Z, leaf_rotation=90)
        #plt.savefig('fulldendrogram.png', dpi=400)
        #plt.close()

        fig, ax = plt.subplots(figsize=(20, 10))
        plt.title(f'Truncated dendrogram')
        plt.ylabel('Distance (UPGMA)')
        d2 = dendrogram(Z, leaf_rotation=90)
        plt.savefig(f'truncated_level_{threshold}.png', dpi=400)
        plt.close()
    else:
        
        df=pd.read_csv('precluster_sm.csv')
        df.index = df.columns
        Z = linkage(df, 'average')
        outname='fulllinkagearray.txt'
        output=open(outname, "w")
        #output.write(Z)
        output.close()

        w = int(ngenomes /10)
        
        fig, ax = plt.subplots(figsize=(w, 10))                                                                       
        plt.title(f'Full dendrogram')                                                                             
        plt.ylabel('Distance (UPGMA)')                                                                                  
        d1 = dendrogram(Z, leaf_rotation=90)                                                                          
        plt.savefig(f'Full_dendrogram.png', dpi=400)
        plt.close()


if __name__=='__main__':

    args = add_args(sys.argv[1:])
    #refine =args.refine
    threshold = args.threshold
    genomedir = args.fastas
    gList =[]
    for filename in os.listdir(genomedir):
        if filename.endswith('.fna') or filename.endswith('.fasta') or filename.endswith('.fas'):
            gList.append(genomedir + filename)
    print(f"Found {len(gList)} genomes for comparison")

    if threshold:
        dendro(threshold, refine)
    else:
        precluster(gList)
        dendro(threshold, len(gList))
