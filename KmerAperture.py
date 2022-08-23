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
        help="Kmer k size. Preferably provide an odd number",
        type=int,
        default=21,
        required=False,
    )
    args = parser.parse_args(a)
    return args

def canonicalise(kmer):
    canonical_kmer=''
    rc_kmer = screed.rc(kmer)
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
    return canonical_kmer

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        c_kmer = canonicalise(kmer)
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
    kmer1uniq = kmer1set - kmer2set
    kmer2uniq = kmer2set - kmer1set

    return kmer1uniq, kmer2uniq

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
    unioncorrect =0
    SNPs =0
    for pair in kmer1ranges:
        rangediff = pair[1] - pair[0]
        if rangediff > ksize:
            unioncorrect+=rangediff
        if rangediff == ksize:
            SNPs+=1

    return unioncorrect, SNPs

def run_KmerAperture(gList, reference, ksize):

    print('Reading in file 1')

    kmers1 = read_kmers_from_file(reference, ksize)
    kmer1set=set(kmers1)

    outname = f'{reference}_{ksize}.csv'
    output=open(outname, "w")
    output.write('gID,Jaccard,Union,Intersection,SNP1,SNP2\n')

    outname2 = f'{reference}_{ksize}_timings.csv'
    output2=open(outname2, "w")
    output2.write('Timetoread,timeforset,timeforSNP\n')


    for genome2 in gList:
        print(f'Reading in query genome {genome2}')

        time0 = time.time()
        kmers2 = read_kmers_from_file(genome2, ksize)
        readtime= (time.time())-time0

        analysistime0 =time.time()
        kmer2set=set(kmers2)
        kmer1uniq, kmer2uniq = get_uniques(kmer1set, kmer2set)
        kmer2indices = get_indices(kmer2uniq, kmers2)
        kmer2indices.sort()
        kmer2ranges = get_ranges(kmer2indices)
        kmer2diff, kmer2SNPs = get_accessory(kmer2ranges, ksize)
        analysistime = (time.time())-analysistime0

        kmer1indices = get_indices(kmer1uniq, kmers1)
        kmer1indices.sort()
        kmer1ranges = get_ranges(kmer1indices)
        kmer1diff, kmer1SNPs = get_accessory(kmer1ranges, ksize)

        intersection = len(kmer1set.intersection(kmer2set))
        union = len(kmer1set.union(kmer2set))
        jaccard = intersection/union
        setanalysistime = (time.time())-time0

        result =f"{genome2},{jaccard},{union},{intersection},{kmer1SNPs},{kmer2SNPs}\n"
        output.write(result)
        print(result)

        timeresult =f"{readtime},{setanalysistime},{analysistime}"
        output2.write(timeresult)
        print(timeresult)

if __name__=='__main__':

    args = add_args(sys.argv[1:])
    reference =args.reference
    gList = list(Path(args.fastas).glob("*.[fa][fas][fasta][fna]"))
    run_KmerAperture(
        gList,
        reference,
        args.kmersize)
