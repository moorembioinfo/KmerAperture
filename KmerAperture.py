#!/usr/bin/env python3
import pandas as pd
import numpy as np
import screed
import os
import random
from itertools import groupby
import time
import argparse
import sys
import itertools


def add_args(a):
    parser = argparse.ArgumentParser(description="KmerAperture")
    parser.add_argument(
        "--fastas",
        "-f",
        help="Path to directory containing fastas for analysis",
        required=True,
    )
    parser.add_argument(
        "--reference",
        "-r",
        help="Reference genome for comparison",
        required=True,
    )
    parser.add_argument(
        "--kmersize",
        "-k",
        help="Kmer k size. Must be an odd number",
        type=int,
        default=21,
        required=False,
    )

    args = parser.parse_args(a)
    return args

def canonicalise(kmer):
    canonical_kmer=''
    try:
        rc_kmer = screed.rc(kmer)
        if kmer < rc_kmer:
            canonical_kmer = kmer
        else:
            canonical_kmer = rc_kmer
    except:
        pass
    return canonical_kmer

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        c_kmer = canonicalise(kmer)
        if c_kmer:
            if not 'N' in c_kmer:
                kmers.append(c_kmer)
    return kmers

def read_kmers_from_file(filename, ksize):
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence
        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers
    return all_kmers

def get_uniques(kmer1set, kmer2set):
    kmeruniq = kmer1set - kmer2set
    return kmeruniq

def get_ranges(lst):
    pos = (j - i for i, j in enumerate(lst))
    t = 0
    for i, els in groupby(pos):
        l = len(list(els))
        el = lst[t]
        t += l
        yield (el, el+l)

def get_indices(kmeruniq, kmers):
    indexlist = [i for i, e in enumerate(kmers) if e in kmeruniq]
    return indexlist

def get_accessory(kmer1ranges, ksize):
    acclength = 0
    allSNPranges = []
    accranges=[]
    for pair in kmer1ranges:
        rangediff = pair[1] - pair[0]
        if rangediff == ksize:
            allSNPranges.append(pair)
        if (rangediff>=(ksize*2)):
            acclength+=rangediff
            acclength +=(ksize-1)
            accranges.append(pair)
    return allSNPranges, accranges, acclength


def assert_kmer(kmerranges, k, kmers2):
    klist = []
    for pair in kmerranges:
        kp=[]
        startpos = pair[0]
        endpos = pair[1]
        kgap = int((k-1)/2)
        km = kmers2[startpos+kgap]
        km_rc = screed.rc(km)
        mkmer0 = km[:kgap] + km[kgap+1:]
        mkmer1 = km_rc[:kgap] + km_rc[kgap+1:]
        klist.extend([mkmer0, mkmer1])#
    return(klist)

def find_dense_SNP(kmer2ranges, kmer1ranges, k, kmers2, kmers1):

    SNPs=0
    kend = (2*k)
    seriessize = range(k+2,kend)

    for L in seriessize:
        middlekmers1 = []
        middlekmers2 = []
        k2_L_ranges = []
        k1_L_ranges = []
        for pair in kmer1ranges:
            rangediff = pair[1] - pair[0]
            if rangediff == L:
                k1_L_ranges.append(pair)
        for pair in kmer2ranges:
            rangediff = pair[1] - pair[0]
            if rangediff == L:
                k2_L_ranges.append(pair)

        #Space between SNPs at L=k+2  is L-k-1. But for python, +1
        a =[]
        b = []
        spacer = (L-k)
        for pair in k1_L_ranges:
            startpos = pair[0]
            mkmer1 = kmers1[startpos + (k-1)]
            mkmer3 = mkmer1[1:spacer] + mkmer1[spacer+1:]
            km1_rc=screed.rc(mkmer1)
            mkmer2 = km1_rc[1:spacer] + km1_rc[spacer+1:]
            middlekmers1.extend([mkmer3, mkmer2])
            a.extend([mkmer1, km1_rc])
        for pair in k2_L_ranges:
            startpos = pair[0]
            mkmer1 = kmers2[startpos + (k-1)]
            mkmer3 = mkmer3[1:spacer] + mkmer3[spacer+1:]
            km1_rc=screed.rc(mkmer1)
            mkmer2 = km1_rc[1:spacer] + km1_rc[spacer+1:]
            middlekmers2.extend([mkmer3, mkmer2])
            b.extend([mkmer1, km1_rc])

        denseSNPs = len(set(middlekmers1).intersection(set(middlekmers2)))
        pairs_kmers = list(itertools.product(a, b))
        dSNPs = 0
        for pair in pairs_kmers:
            counter =0
            for p, g in zip(pair[0], pair[1]):
                if p==g:
                    counter+=1
            if (counter>=(k-3)) & (counter<k):
                dSNPs+=counter
                SNPs+=counter

    return(SNPs)


def run_KmerAperture(gList, reference, ksize):

    print('Reading in file 1')

    kmers1 = read_kmers_from_file(reference, ksize)
    kmer1set=set(kmers1)

    outname = f'./{reference}_{ksize}.csv'
    output=open(outname, "w")
    output.write('gID,matchedSNP,denseSNPs,acc1,acc2\n')

    outname2 = f'./{reference}_{ksize}_timings.csv'
    output2=open(outname2, "w")
    output2.write('Timetoread,timeforset,timeforSNP\n')


    for genome2 in gList:
        print(f'Reading in query genome {genome2}')

        time0 = time.time()
        kmers2 = read_kmers_from_file(genome2, ksize)
        kmer2set=set(kmers2)
        readtime= (time.time())-time0

        analysistime0 =time.time()
        kmer2uniq = get_uniques(kmer2set, kmer1set)
        kmer2indices = get_indices(kmer2uniq, kmers2)
        kmer2indices.sort()
        kmer2ranges = get_ranges(kmer2indices)
        kmer2ranges_=list(kmer2ranges)
        SNPranges2, accranges2, acclength2 = get_accessory(kmer2ranges_, ksize)

        klist2 = assert_kmer(SNPranges2, ksize, kmers2)

        kmer1uniq = get_uniques(kmer1set, kmer2set)
        kmer1indices = get_indices(kmer1uniq, kmers1)
        kmer1indices.sort()
        kmer1ranges = get_ranges(kmer1indices)
        kmer1ranges_=list(kmer1ranges)
        SNPranges1, accranges1, acclength1 = get_accessory(kmer1ranges_, ksize)

        klist1 = assert_kmer(SNPranges1, ksize, kmers1)

        int(matchedSNPs = len(set(klist1).intersection(set(klist2)))/2)
        denseSNPs = find_dense_SNP(kmer2ranges_, kmer1ranges_, ksize, kmers2, kmers1)
        analysistime = (time.time())-analysistime0

        result =f"{genome2},{matchedSNPs},{denseSNPs},{acclength1},{acclength2}\n"
        output.write(result)
        print(result)

        timeresult =f"{readtime},{analysistime}\n"
        output2.write(timeresult)
        print(timeresult)


if __name__=='__main__':

    args = add_args(sys.argv[1:])
    reference =args.reference
    genomedir = args.fastas
    gList =[]
    for filename in os.listdir(genomedir):
        if filename.endswith('.fna') or filename.endswith('.fasta') or filename.endswith('.fas'):
            gList.append(genomedir + filename)

    print(f"Found {len(gList)} genomes for comparison")

    run_KmerAperture(
        gList,
        reference,
        args.kmersize)
