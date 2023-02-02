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
    #complement = {'A': 'T', 'C':'G', 'T':'A', 'G':'C'}
    rcdict = {'AA': 1, 'AC':1, 'AT':0, 'AG':1,
    'TT': 2, 'TA':0, 'TC':2, 'TG':2,
    'CA': 1, 'CC':1, 'CG':0, 'CT':2,
    'GC':0, 'GA':1, 'GT':2, 'GG':2}
    allkmers=[]
    count =0
    countrc =0
    for kmer in naivekmers:
        if not 'N' in kmer:
            for pos in kmer:
                firstpos=0
                lastpos = -1
                base = kmer[0] + kmer[-1]
                #if 'N' in base:
                #    continue
                val = rcdict.get(base)
                if val>0:
                    if val == 1:
                        allkmers.append(kmer)
                        count+=1
                    elif val == 2:
                        rckmer = str(reverse_complement(kmer))
                        allkmers.append(rckmer)
                        countrc+=1
                    break
                else:
                    firstpos+=1
                    lastpos-=1
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
