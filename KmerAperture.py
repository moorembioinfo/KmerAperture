#!/usr/bin/env python3
import pandas as pd
import numpy as np
import screed
import random
import time
import argparse
import sys
import itertools
import math
from utils import canon, build_kmers, read_kmers_from_file
from utils import get_uniques, get_ranges, get_indices
from Bio.Seq import reverse_complement
import os
import subprocess
dir_path = str(os.path.dirname(os.path.realpath(__file__)))



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
        "--pyonly",
        "-py",
        help="Run KmerAperture in python only. Slower than using the ocaml parser",
        default=False,
        action='store_true',
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
            acclength +=(ksize)
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

def find_dense_SNP2(kmer2ranges, kmer1ranges, k, kmers2, kmers1, filename, qfilename):

    '''
    Find contiguous kmer series from size k+2 to n(k-1)+1
    Extract full sequence and pattern match SNPs for count and extract
    '''
    SNPs = 0
    Ns = 0

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
                                #print('Indel!')
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

def get_indels(kmer2ranges, k, kmers2, kmers1):
    '''
    Match kmers flanking series in alternate genome(s)
    if they're contiguous in alternate
    '''
    for pair2 in kmer2ranges:
        rangediff2 = pair2[1] - pair2[0]
        if (rangediff2 >= (k+2)) and (rangediff2 <= (k+50)):
            kmer1 = kmers2[pair2[0]-1]
            kmer2 = kmers2[pair2[1]+1]

            print(rangediff2)
            print([kmer1, kmer2])

            refpos = get_indices([kmer1, kmer2], kmers1)
            print(refpos)
            if (refpos[-1] - refpos[0]) == 4:
                print(f'Indel of size {rangediff} (-k+1)')





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

def outputs(refdict, querynamedict):
    '''
    Generate kmeraperture outputs
    '''
    print('Generating SNP output...\n')
    df = pd.DataFrame(list(refdict.items()), columns = ['refpos','refbase'])
    df.columns = ['refpos','refbase']
    for key in querynamedict:
        querydict_=querynamedict.get(key)
        df[key] = df['refpos'].map(querydict_)
        df.loc[df[key].isna(),key] = df['refbase']
    df.to_csv('SNPmatrix.polymorphic.csv')


    print('\n\n\nFinished, thanks for using KmerAperture!\n\n\n')



def run_KmerAperture(gList, reference, ksize, pyonly):

    print(f'Reading in reference genome {reference}')
    if pyonly:
        kmers1 = read_kmers_from_file(reference, ksize)
    else:

        ocamlparser = dir_path+'/parser/KmerApertureParser'
        proc = subprocess.Popen([ocamlparser, reference, str(ksize)], stdout=subprocess.PIPE, encoding='utf8')
        kmers1 = proc.stdout.read().split()
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
        #Generate canonical kmers for query genome
        print(f'Reading in query genome {genome2}')
        if pyonly:
            kmers2 = read_kmers_from_file(genome2, ksize)
        else:
            proc = subprocess.Popen([ocamlparser, genome2, str(ksize)], stdout=subprocess.PIPE, encoding='utf8')
            kmers2 = proc.stdout.read().split()
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

        newdense = find_dense_SNP2(kmer2ranges_, kmer1ranges_, ksize, kmers2, kmers1, reference, genome2)

        #get_indels(kmer2ranges_, ksize, kmers2, kmers1)


        SNPs = matchedSNPs+newdense

        result =f"{genome2},{SNPs},{acclength1},{acclength2}\n"
        output.write(result)

    #df = pd.DataFrame.from_dict(refdict)

    outputs(
    refdict,
    querynamedict
    )

if __name__=='__main__':

    args = add_args(sys.argv[1:])
    reference =args.reference
    genomedir = args.fastas
    kmersize = args.kmersize
    pyonly = args.pyonly

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
        pyonly)
