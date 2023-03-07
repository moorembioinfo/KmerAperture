from itertools import groupby
import screed
import numpy as np
from Bio.Seq import reverse_complement
import os
import pandas as pd

def compl(base):
    basedict = {'A':'T', 'T':'A', 'G': 'C', 'C':'G'}
    return(basedict.get(base))

def rangedict(accranges):
    rng_dict = {}
    for r in accranges:
        start, end = r
        for i in range(start, end+1):
            rng_dict[i] = '-'
    return rng_dict

def canon(naivekmers):

    allkmers = []
    for kmer in naivekmers:
        canonical_kmer=kmer
        rckmer = str(reverse_complement(kmer))
        if kmer> rckmer:
            canonical_kmer=rckmer
        allkmers.append(canonical_kmer)

    return allkmers

def build_kmers(sequence, ksize):
    kmers = []
    naivekmers = [sequence[x:x+ksize].upper() for x in range(len(sequence) - ksize + 1)]
    kmers =canon(naivekmers)
    return kmers

def read_kmers_from_file(filename, ksize):
    all_kmers = []
    sequence=''
    for record in screed.open(filename):
        sequence += record.sequence
    all_kmers=build_kmers(sequence, ksize)
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
