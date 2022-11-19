#!/usr/bin/env python3
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import argparse
import os

def add_args(a):
    parser = argparse.ArgumentParser(description="Precluster")
    parser.add_argument(
        "--fastas",
        help="Path to fasta dir ",
        required=True,
    )
    parser.add_argument(
        "--refine",
        help="Cluster based on this level of truncation ",
        required=False,
    )
    args = parser.parse_args(a)
    return args

def precluster(filelist):
    print('Running sourmash sketch...\n')
    scmd = 'sourmash sketch dna -p k=31'
    os.system(scmd)
    print('Running sourmash compare...\n')
    ccmd = 'sourmash compare --csv precluster_sm.csv *sig'
    os.system(ccmd)
    os.remove('*sig')

def dendro(threshold, refine):

    if refine:
        fh = open('fulllinkagearray.txt')
        Z = fh.readlines()[0].rstrip()
        fig, ax = plt.subplots(figsize=(20, 10))
        plt.title(f'Truncated dendrogram')
        plt.ylabel('Distance (Ward)')
        plt = dendrogram(Z, leaf_rotation=90)
        plt.savefig(f'truncated_level_{refine}.png', dpi=400)

    else:
        df=pd.read_csv('precluster_sm.csv')
        df.index = df.columns
        Z = linkage(df, 'ward')

        fig, ax = plt.subplots(figsize=(20, 10))
        plt.title(f'Full dendrogram')
        plt.ylabel('Distance (Ward)')
        plt = dendrogram(Z, leaf_rotation=90)
        plt.savefig('fulldendrogram.png', dpi=400)
        plt.close()

        fig, ax = plt.subplots(figsize=(20, 10))
        plt.title(f'Truncated dendrogram')
        plt.ylabel('Distance (Ward)')
        plt = dendrogram(Z, leaf_rotation=90)
        plt.savefig(f'truncated_level_{refine}.png', dpi=400)

        outname='fulllinkagearray.txt'
        output=open(outname, "w")
        output.write(Z)
        output.close()

if __name__=='__main__':

    args = add_args(sys.argv[1:])
    refine =args.refine
    threshold = args.threshold
    genomedir = args.fastas
    gList =[]
    for filename in os.listdir(genomedir):
        if filename.endswith('.fna') or filename.endswith('.fasta') or filename.endswith('.fas'):
            gList.append(genomedir + filename)
    print(f"Found {len(gList)} genomes for comparison")

    if refine:
        dendro(threshold, refine)
    else:
        precluster(filelist)
        dendro(threshold, False)
