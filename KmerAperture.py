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
from Bio.Seq import reverse_complement

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

def canon(naivekmers):
    allkmers=[]
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


def find_dense_SNP2(kmer2ranges, kmer1ranges, k, kmers2, kmers1, filename, qfilename):

    '''
    Find contiguous kmer series from size k+2 to n(k-1)+1
    Extract full sequence and pattern match SNPs for count and extract
    '''
    SNPs = 0

    sequence=''
    for record in screed.open(reference):
        sequence += record.sequence
    qsequence =''
    for record in screed.open(qfilename):
        qsequence += record.sequence

    upperboundSNP = 10
    ks = k+2
    ke = (upperboundSNP *(k-1)) + 1
    for L in range(ks, ke):
        k2_L_ranges = []
        k1_L_ranges = []
        for pair, pair2 in zip(kmer1ranges, kmer2ranges):
            rangediff1 = pair[1] - pair[0]
            if rangediff1 == L:
                k1_L_ranges.append(pair[0])
            rangediff2 = pair2[1] - pair2[0]
            if rangediff2 == L:
                k2_L_ranges.append(pair2[0])

        Lseqlen = L+(k-1)
        aseqs = []
        bseqs = []
        for pos in k1_L_ranges:
            exseq = sequence[pos:pos+Lseqlen]
            rcexseq = str(reverse_complement(exseq))
            aseqs.extend([exseq, rcexseq])
        for pos in k2_L_ranges:
            exseq = qsequence[pos:pos+Lseqlen]
            rcexseq = str(reverse_complement(exseq))
            bseqs.extend([exseq, rcexseq])

        if aseqs and bseqs:
            pairs_seqs = list(itertools.product(aseqs, bseqs))
            for pairseq in pairs_seqs:
                #print(pairseq)
                if (len(pairseq[0]) > 1) and (len(pairseq[1]) > 1): #This should always be the case and the condition not be needed
                    pmindex = [index for index, elem in enumerate(pairseq[0]) if elem != pairseq[1][index]]
                    #print(len(pmindex))
                    cutoff = round(Lseqlen*0.15)
                    #print(f'{L}, {Lseqlen}, {cutoff}')
                    if len(pmindex) < cutoff:
                        outr = list(get_ranges(pmindex))

                        noindel =True
                        indellen=0
                        for r in outr:
                            if r[1]-r[0] >1:
                                print('Indel!')
                                noindel=False
                                indellen+=(r[1]-r[0])
                        if noindel:
                            SNPs+=len(pmindex)
                            break
                        else:
                            if len(pmindex) > indellen:
                                SNPs+=(len(pmindex)-indellen)
                            break


    return(SNPs)



def find_dense_SNP(kmer2ranges, kmer1ranges, k, kmers2, kmers1):

    SNPs2=0
    SNPs3=0
    SNPs4=0
    duplicates = []
    for SNPc in [2, 3, 4]:
        ksend = (SNPc *(k-1)) + 1
        ksstart = SNPc+k
        seriesrange = range(ksstart,ksend)
        for L in seriesrange:
            k2_L_ranges = []
            k1_L_ranges = []
            for pair in kmer1ranges:
                rangediff = pair[1] - pair[0]
                if rangediff == L:
                    k1_L_ranges.append(pair)
            for pair2 in kmer2ranges:
                rangediff2 = pair2[1] - pair2[0]
                if rangediff2 == L:
                    k2_L_ranges.append(pair2)
            a =[]
            b = []
            spacer = (L-k)
            for pair in k1_L_ranges:
                startpos = pair[0]
                mkmer1 = kmers1[startpos+(k-1)]
                km1_rc=screed.rc(mkmer1)

                rmkmer1 = kmers1[startpos+(L-k)]
                rkm1_rc = screed.rc(rmkmer1)
                a.extend([mkmer1, km1_rc, rmkmer1, rkm1_rc])

            for pair2 in k2_L_ranges:
                startpos2 = pair2[0]
                mkmer2 = kmers2[startpos2+(k-1)]
                km2_rc=screed.rc(mkmer2)

                rmkmer2=kmers2[startpos2+(L-k)]
                rkm2_rc=screed.rc(rmkmer2)
                b.extend([mkmer2, km2_rc, rmkmer2, rkm2_rc])

            pairs_kmers = list(itertools.product(a, b))
            dSNPs = 0
            for pair in pairs_kmers:
                counter =0
                for p, g in zip(pair[0], pair[1]):
                    if p==g:
                        counter+=1
                snps =k-counter
                if snps==SNPc:
                    if SNPc==2:
                        SNPs2+=snps
                        break
                    if SNPc==3:
                        SNPs3+=snps
                        break
                    if SNPc==4:
                        SNPs4+=SNPc
                        break
    return(SNPs2, SNPs3, SNPs4)


def run_KmerAperture(gList, reference, ksize):

    print('Reading in file 1')

    kmers1 = read_kmers_from_file(reference, ksize)
    kmer1set=set(kmers1)

    outname = f'./{reference}_{ksize}.csv'
    output=open(outname, "w")
    output.write('gID,Jaccard,matchedSNP,denseSNPs2,denseSNPs3,denseSNPs4,newdense,acc1,acc2\n')

    outname2 = f'./{reference}_{ksize}_timings.csv'
    output2=open(outname2, "w")
    output2.write('Timetoread,timeforset,timeforSNP,timefordenseSNP\n')


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

        matchedSNPs = int(len(set(klist1).intersection(set(klist2)))/2)
        analysistime = (time.time())-analysistime0
        denseSNPs2, denseSNPs3, denseSNPs4 = find_dense_SNP(kmer2ranges_, kmer1ranges_, ksize, kmers2, kmers1)
        analysistime2 = (time.time())-analysistime0

        newdense = find_dense_SNP2(kmer2ranges_, kmer1ranges_, ksize, kmers2, kmers1, reference, genome2)


        Jtime0 = time.time()
        J = len(kmer1set.intersection(kmer2set))/len(kmer1set.union(kmer2set))
        jtime = time.time()-Jtime0

        result =f"{genome2},{J},{matchedSNPs},{denseSNPs2},{denseSNPs3},{denseSNPs4},{newdense},{acclength1},{acclength2}\n"
        output.write(result)
        print(result)

        timeresult =f"{readtime},{jtime},{analysistime},{analysistime2}\n"
        output2.write(timeresult)

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
