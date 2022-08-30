#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
import screed
import os
import random
from itertools import groupby
import time
import argparse
import sys


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
    parser.add_argument(
    "--sensitive",
    "-s",
    help="Run in sensitive mode",
    default=False,
    action="store_true"
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
    SNPs =0
    allSNPranges = []
    for pair in kmer1ranges:
        rangediff = pair[1] - pair[0]
        if rangediff == ksize:
            SNPs+=1
            allSNPranges.append(pair)
    return allSNPranges, SNPs

def assert_kmer(kmerranges, k, kmers2):
    SNPs = 0
    for pair in kmerranges:
        kp=[]
        startpos = pair[0]
        endpos = pair[1]

        kf=kmers2[startpos]
        ke=kmers2[endpos-1]
        kf_rc = screed.rc(kf)
        ke_rc = screed.rc(ke)

        kgap = int((k-1)/2)
        km = kmers2[startpos+kgap]
        kt1=(kf[kgap:] + ke[1:kgap+1])
        kt2=(kf_rc[kgap:] + ke[1:kgap+1])
        kt3=(kf[kgap:] + ke_rc[1:kgap+1])
        kt4=(kf_rc[kgap:] + ke_rc[1:kgap+1])
        kp=[kt1,kt2,kt3,kt4]
        if str(km) in kp:
            SNPs+=1
    return(SNPs)

def run_KmerAperture(gList, reference, ksize, sensitive):

    print('Reading in file 1')

    kmers1 = read_kmers_from_file(reference, ksize)
    kmer1set=set(kmers1)

    outname = f'./{reference}_{ksize}.csv'
    output=open(outname, "w")
    output.write('gID,Jaccard,Union,Intersection,SNP1,SNP2\n')

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
        SNPranges2, kmer2SNPs = get_accessory(kmer2ranges, ksize)
        if sensitive:
            kmer2SNPs = assert_kmer(SNPranges2, ksize, kmers2)
        analysistime = (time.time())-analysistime0

        kmer1uniq = get_uniques(kmer1set, kmer2set)
        kmer1indices = get_indices(kmer1uniq, kmers1)
        kmer1indices.sort()
        kmer1ranges = get_ranges(kmer1indices)
        SNPranges1, kmer1SNPs = get_accessory(kmer1ranges, ksize)

        intersection = len(kmer1set.intersection(kmer2set))
        union = len(kmer1set.union(kmer2set))
        jaccard = intersection/union
        setanalysistime = (time.time())-time0

        result =f"{genome2},{jaccard},{union},{intersection},{kmer1SNPs},{kmer2SNPs}\n"
        output.write(result)
        print(result)

        timeresult =f"{readtime},{setanalysistime},{analysistime}\n"
        output2.write(timeresult)
        print(timeresult)


if __name__=='__main__':

    args = add_args(sys.argv[1:])
    reference =args.reference
    gList = list(Path(args.fastas).glob("*.fasta"))
    print(gList)
    run_KmerAperture(
        gList,
        reference,
        args.kmersize,
        args.sensitive)
