#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import numpy as np
import argparse
import os, sys


def add_args(a):
    parser = argparse.ArgumentParser(description="Precluster")
    parser.add_argument(
        "--fastas",
        help="Path to fasta dir ",
        required=False,
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
    for f in filelist:
        fn = f.split('/')[-1]
        scmd = f'sourmash sketch dna -p k=31 -o {fn}.sig {f}'
        os.system(scmd)
        #print(scmd)
    pathex = filelist[0]
    path='/'.join(pathex.split('/')[0:-1])
    print('Running sourmash compare...\n')
    ccmd = f'sourmash compare --csv precluster_sm.csv *sig'
    os.system(ccmd)

def dendro(ngenomes):

    df=pd.read_csv('precluster_sm.csv')
    df.index = df.columns
    print(df.head)
    Z = linkage(df, 'average')
    print(Z)
    w = int(ngenomes /20)
    if w<20:
        w=20
    fig, ax = plt.subplots(figsize=(w, 10))                                                                       
    plt.title(f'Full dendrogram (UPGMA)')                                                                             
    plt.ylabel('Distance')                                                                                  
    d1 = dendrogram(Z, leaf_rotation=90, labels=df.index, leaf_font_size=2)  
    plt.savefig(f'Full_dendrogram.png', dpi=400)
    plt.close()

def cluster_genomes(threshold):
    df=pd.read_csv('precluster_sm.csv')
    df.index = df.columns
    Z = linkage(df, 'average')
    labels = fcluster(Z, t=threshold, criterion='distance') 
    print(labels)
    df['clusters'] = labels
    df.to_csv('sourmashclusters.csv')

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
        cluster_genomes(threshold)
    else:
        precluster(gList)
        dendro(len(gList))
