from itertools import groupby
import screed
import numpy as np


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
        kmers.append(c_kmer)
    return kmers

def read_kmers_from_file(filename, ksize):
    all_kmers = []
    sequence=''
    for record in screed.open(filename):
        sequence += record.sequence#+contigspace
        #kmers = build_kmers(sequence, ksize)
        #all_kmers += kmers
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
