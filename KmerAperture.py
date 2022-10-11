#!/usr/bin/env python3
import pandas as pd
import numpy as np
import screed
import os
import random
import time
import argparse
import sys
import itertools
import math
from utils import canonicalise, build_kmers, read_kmers_from_file
from utils import get_uniques, get_ranges, get_indices


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
        "--polySNP",
        "-p",
        help="Output matrix of polymorphic sites ",
        default=False,
        required=False,
    )
    args = parser.parse_args(a)
    return args



def get_accessory(kmer1ranges, ksize):
    '''
    Find contiguous kmer series indicative of SNPs
    or accessory
    '''
    acclength = 0
    allSNPranges = []
    accranges=[]
    for pair in kmer1ranges:
        rangediff = pair[1] - pair[0]
        if rangediff == ksize:
            allSNPranges.append(pair)
        if (rangediff>=(ksize*2)):
            acclength+=rangediff
            acclength -=(ksize-1)
            accranges.append(pair)
    return allSNPranges, accranges, acclength

def assert_kmer(kmerranges, k, kmers2):
    '''
    Match middle k-mers with middle base removed
    for both possible strands
    '''
    klist = []
    positiondict ={}
    for pair in kmerranges:
        kp=[]
        startpos = pair[0]
        endpos = pair[1]
        kgap = int((k-1)/2)
        km = kmers2[startpos+kgap]
        km_rc = screed.rc(km)
        mkmer0 = km[:kgap] + km[kgap+1:]
        mkmer1 = km_rc[:kgap] + km_rc[kgap+1:]
        klist.extend([mkmer0, mkmer1])

        positiondict[mkmer0] = (startpos+k) #Should be -1, but SNP positions are 1 indexed
        positiondict[mkmer1] = (startpos+k) #-1
    return(klist, positiondict)

def get_SNPs(middlekinter, klist1pos, klist2pos, sequence, qfilename, refdict):
    '''
    Find and extract SNPs estimated by KmerAperture
    '''
    refdict = {}
    querydict = {}
    sequence=''
    for record in screed.open(reference):
        sequence += record.sequence
    qsequence =''
    for record in screed.open(qfilename):
        qsequence += record.sequence

    for mkmer in list(middlekinter):
        refpos = klist1pos.get(mkmer)
        refseq = sequence[refpos-5:refpos+4]
        querypos = klist2pos.get(mkmer)
        queryseq = qsequence[querypos-5:querypos+4]
        rcseq = screed.rc(queryseq)
        km = queryseq[:4]+queryseq[5:]
        refkm = refseq[:4]+refseq[5:]
        rckm=rcseq[:4]+rcseq[5:]
        SNP=''
        if km == refkm:
            SNP=queryseq[4]
        elif rckm == refkm:
            SNP=rcseq[4]
        refbase = refseq[4]

        if not refpos in refdict.keys():
            refdict[refpos] = refbase
        querydict[refpos] = SNP

    return(querydict, refdict)


def run_KmerAperture(gList, reference, ksize, polySNPmat):

    print(f'Reading in reference genome {reference}')
    kmers1 = read_kmers_from_file(reference, ksize)
    kmers1_=[]
    for kmer in kmers1:
        if not 'N' in kmer:
            kmers1_.append(kmer)
    kmer1set=set(kmers1_)

    outname = f'./{reference}_{ksize}.csv'
    output=open(outname, "w")
    output.write('gID,SNP,acc1,acc2\n')

    querynamedict = {}
    refdict = {}

    for genome2 in gList:
        print(f'Reading in query genome {genome2}')
        kmers2 = read_kmers_from_file(genome2, ksize)
        kmers2_=[]
        for kmer in kmers2:
            if not 'N' in kmer:
                kmers2_.append(kmer)
        kmer2set=set(kmers2_)

        kmer2uniq = get_uniques(kmer2set, kmer1set)
        kmer2indices = get_indices(kmer2uniq, kmers2)
        kmer2indices.sort()
        kmer2ranges = get_ranges(kmer2indices)
        kmer2ranges_=list(kmer2ranges)
        SNPranges2, accranges2, acclength2 = get_accessory(kmer2ranges_, ksize)
        klist2, klist2pos = assert_kmer(SNPranges2, ksize, kmers2)

        kmer1uniq = get_uniques(kmer1set, kmer2set)
        kmer1indices = get_indices(kmer1uniq, kmers1)
        kmer1indices.sort()
        kmer1ranges = get_ranges(kmer1indices)
        kmer1ranges_=list(kmer1ranges)
        SNPranges1, accranges1, acclength1 = get_accessory(kmer1ranges_, ksize)
        klist1, klist1pos = assert_kmer(SNPranges1, ksize, kmers1)

        middlekinter = set(klist1).intersection(set(klist2))
        matchedSNPs = int(len(middlekinter)/2)

        querydict, refdict = get_SNPs(middlekinter, klist1pos, klist2pos, reference, genome2, refdict)

        querynamedict[genome2] = querydict

        result =f"{genome2},{matchedSNPs},{acclength1},{acclength2}\n"
        output.write(result)

    #df = pd.DataFrame.from_dict(refdict)

    if polySNPmat:
        print('Generating SNP output...\n')
        df = pd.DataFrame(list(refdict.items()), columns = ['refpos','refbase'])
        df.columns = ['refpos','refbase']

        for key in querynamedict:
            querydict_=querynamedict.get(key)
            df[key] = df['refpos'].map(querydict_)

            df.loc[df[key].isna(),key] = df['refbase']


        df.to_csv('SNPmatrix.polymorphic.csv')
    print('Finished!')

if __name__=='__main__':

    args = add_args(sys.argv[1:])
    reference =args.reference
    genomedir = args.fastas
    kmersize = args.kmersize
    polySNPmat = args.polySNP
    if (kmersize % 2) != 1:

        print('\nPlease enter an odd numbered integer for k\n\nExiting...')
        exit()

    gList =[]
    for filename in os.listdir(genomedir):
        if filename.endswith('.fna') or filename.endswith('.fasta') or filename.endswith('.fas'):
            gList.append(genomedir + filename)

    print(f"Found {len(gList)} genomes for comparison")

    run_KmerAperture(
        gList,
        reference,
        kmersize,
        polySNPmat)
