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
import json
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
        if (rangediff>(ksize)):
            acclength+=((rangediff-ksize)+1)
            accranges.append(pair)
    return allSNPranges, accranges, acclength

def filter_accessory_positions(accranges, indelranges, denseranges, ksize):
    '''
    For all k-mer series >k assume they're query-accessory
    Filter out those that instead represent dense SNPs or indels
    '''
    indelranges.extend(denseranges)
    acconly =[]
    for ranges in accranges:
        if ranges not in indelranges:
            acconly.append(ranges)

    acclength =0
    for ranges in acconly:
        rangediff = ranges[1]-ranges[0]
        acclength+=((rangediff-ksize)+1)
    return(acconly, acclength)

def get_core(refacclists):
    '''
    Take the positions that are accessory to the reference
    (absent in a query genome)
    Find if they're absent in >5% of query genomes
    '''
    positions = set()
    for r in refacclists:
        for start, end in r:
            positions.update(range(start, end+1))
    positions = np.array(sorted(positions))
    max_len = max(len(r) for r in refacclists)
    matrix = np.zeros((len(refacclists), len(positions)))
    for i, r in enumerate(refacclists):
        for j, p in enumerate(positions):
            if any(start <= p <= end for start, end in r):
                matrix[i, j] = 1
    freq = matrix.sum(axis=0) / len(refacclists)

    notcoresites = []
    for p in positions[freq > 0.05]:
        notcoresites.append(p)
    return(notcoresites)


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

    rdict ={}
    qdict ={}

    SNPs = 0
    denseSNP_L = 0

    denseranges = []
    denseranges_ref = []

    sequence=''
    for record in screed.open(reference):
        sequence += record.sequence
    qsequence =''
    for record in screed.open(qfilename):
        qsequence += record.sequence

    upperboundSNP = 100
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
                if (len(pairseq[0]) > 1) and (len(pairseq[1]) > 1): #This should always be the case and the condition not be needed
                    pmindex = [index for index, elem in enumerate(pairseq[0]) if elem != pairseq[1][index]]
                    outr = list(get_ranges(pmindex))

                    noindel =True
                    indellen=0
                    for r in outr:
                        if r[1]-r[0] >1:
                            noindel=False
                    if noindel:
                        SNPs+=len(pmindex)
                        ri = sequence.index(pairseq[0]) +k
                        qi = qsequence.index(pairseq[1]) +k

                        #Append series to remove from accessory
                        denseranges.append((qi-k, qi+k+1))
                        denseranges_ref.append((ri-k, ri+k+1))

                        snpdiffs = list(np.diff(pmindex))
                        rb = sequence[ri-1]
                        qb = qsequence[qi-1]
                        rdict[ri] = rb
                        qdict[ri] = qb
                        for x in snpdiffs:
                            rb=''
                            qb=''
                            ri+=x
                            qi+=x
                            rb = sequence[ri-1]
                            qb = qsequence[qi-1]
                            rdict[ri] = rb
                            qdict[ri] = qb

                        denseSNP_L += (Lseqlen+1)

                        break

    return(SNPs, rdict, qdict, denseSNP_L, denseranges, denseranges_ref)

def get_indels(kmer2ranges, k, kmers2, kmers1):
    '''
    Match kmers flanking series in alternate genome(s)
    if they're contiguous in alternate
    '''

    totalins = 0
    insdict = {}
    insranges = []
    totalinsbp = 0

    for pair2 in kmer2ranges:
        rangediff2 = pair2[1] - pair2[0]
        if (rangediff2 >= (k+2)) and (rangediff2 <= (k+49)):
            kmer1 = kmers2[pair2[0]-1]
            kmer2 = kmers2[pair2[1]]
            refpos = get_indices([kmer1, kmer2], kmers1)
            if (refpos[0] + k +1) == (refpos[1]):
                indelsize = (rangediff2-k)+1
                totalinsbp += indelsize
                totalins+=1
                insranges.append(pair2)
    return totalins, insdict, insranges, totalinsbp


def get_SNPs(middlekinter, klist1pos, klist2pos, sequence, qfilename, refdict):
    '''
    Find and extract SNPs estimated by KmerAperture
    '''

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

        refdict[refpos] = refbase
        querydict[refpos] = SNP

    return(querydict, refdict)

def outputs(refdict, querynamedict, rdict, qdict, queryaccdict, notcoresites, positions_dict):
    '''
    Generate kmeraperture outputs
    '''
    print('\n\nGenerating output...\n')
    #Output matrix of SNP sites and refpos
    df = pd.DataFrame(list(refdict.items()), columns = ['refpos','refbase'])
    df.columns = ['refpos','refbase']
    querynamedict['reference'] = refdict
    dft = pd.DataFrame.from_dict(querynamedict, orient='columns')
    dft = dft.mask(dft.isnull() | (dft == '') | (dft.isna()), dft['reference'], axis=0)
    dft.sort_index(inplace=True)
    dft.to_csv('SNPmatrix.polymorphic.csv')

    #Output SNP sites
    concatenated = dft.apply(lambda x: ''.join(x.astype(str)))
    with open('polymorphicsites.fasta', 'w') as f:
        for col_name, sequence in concatenated.items():
            f.write(f'>{col_name.split("./")[-1]}\n{sequence}\n')

    #Output core genome
    querynamedict['reference'] = positions_dict
    df_core = pd.DataFrame.from_dict(querynamedict, orient='columns')
    df_core = df_core.mask(df_core.isnull() | (df_core == '') | (df_core.isna()), df_core['reference'], axis=0)
    df_core.drop(index=notcoresites, inplace=True)
    df_core.sort_index(inplace=True)
    concatenatedcore = df_core.apply(lambda x: ''.join(x.astype(str)))
    with open('core_alignment.fasta', 'w') as f:
        for col_name, sequence in concatenatedcore.items():
            f.write(f'>{col_name.split("./")[-1]}\n{sequence}\n')

    #Output Accessory coordinates
    with open('accessory_coords.json', 'w') as f:
        json.dump(queryaccdict, f)

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
    output.write('gID,SNP,indels,acc1,acc2\n')
    querynamedict = {}
    refdict = {}

    refacclists = []

    queryaccdict = {}

    for genome2 in gList:
        #Generate canonical kmers for query genome
        print(f'Reading in query genome {genome2}')
        if pyonly:
            kmers2 = read_kmers_from_file(genome2, ksize)
        else:
            proc = subprocess.Popen([ocamlparser, genome2, str(ksize)], stdout=subprocess.PIPE, encoding='utf8')
            kmers2 = proc.stdout.read().split()
        print(f'Processing {genome2} k-mers')
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

        newdense,rdict,qdict,denseSNP_L,denseranges,denseranges_ref = find_dense_SNP2(kmer2ranges_,kmer1ranges_,ksize,kmers2,kmers1,reference,genome2)

        refdict.update(rdict)
        querydict.update(qdict)
        querynamedict[genome2] = querydict

        numinsertions, insertionsdict, insranges, totalinsbp = get_indels(kmer2ranges_, ksize, kmers2, kmers1)

        #Adjust accessory to not include dense SNPs and indels
        acclength1_=0
        acclength2_=0
        acc2only, acclength2_ = filter_accessory_positions(accranges2, insranges, denseranges, ksize)
        acc1only, acclength1_ = filter_accessory_positions(accranges1, insranges, denseranges_ref, ksize)

        queryaccdict[genome2] = acc2only
        refacclists.append(acc1only)

        #Total SNPs
        SNPs = matchedSNPs+newdense
        result =f"{genome2},{SNPs},{numinsertions},{acclength1_},{acclength2_}\n"
        output.write(result)

    refsequencefull=''
    for record in screed.open(reference):
        refsequencefull += record.sequence
    positions_dict = {i: char for i, char in enumerate(refsequencefull)}
    notcoresites = get_core(refacclists)

    outputs(
    refdict,
    querynamedict,
    rdict,
    qdict,
    queryaccdict,
    notcoresites,
    positions_dict
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
