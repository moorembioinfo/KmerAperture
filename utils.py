from itertools import groupby
import screed
import numpy as np
from Bio.Seq import reverse_complement
import os
import pandas as pd

def canon(naivekmers):
    '''
    Take the first base of the kmer and determine
    the canonical kmer based on a higher or lower
    value than its complement
    '''
    rcdict = {'A': 1, 'C':1, 'T':2, 'G':2}
    allkmers=[]
    count =0
    countrc = 0
    for kmer in naivekmers:
        base = kmer[0]
        val = rcdict.get(base)
        if val == 1:
            allkmers.append(kmer)
        elif val == 2:
            rckmer = str(reverse_complement(kmer))
            allkmers.append(rckmer)
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
