# KmerAperture

![Python](https://badges.aleen42.com/src/python.svg) ![conda](https://img.shields.io/badge/%E2%80%8B-conda-%2344A833.svg?style=flat&logo=anaconda&logoColor=44A833)

Alignment-free estimation of core and accessory distances in closely related genomes.



## Dependencies

KmerAperture.py is written in python3 and requires the following python packages:

> - numpy
> - screed
> - pandas

## Usage

First set up a conda environment with the appropriate dependencies:

```console
conda env create -f environment.yml
conda activate env-KmerAperture
```

Then run the main script `KmerAperture.py` referencing your fasta directory and reference genome:

```shell
python KmerAperture.py --fastas <fasta dir> --reference <ref file fasta>
```

## Input

A reference genome and a directory of assembled query genomes (fasta format)

## Output

A .csv file with genomeID, estimated SNPs and estimated sequence presence in the reference (acc1) and query genome (acc2)

Optionally output a SNP matrix of polymorphic sites estimated by KmerAperture

...more to come


## Options

Flag &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Short flag | Description | Required | Default val
--------------|------------|-------------|----------|--------------
`--fastas` |  `-f` |  Provide path to query fastas directory | ✅
`--reference` |     `-r` |  Path and file name of fasta reference genome | ✅ | 
`--kmersize` |      `-k` |  k size |                             | 21
`--polySNP` |     `-p` |  Optionally output a matrix of estimated polymorphic SNPs sites |   | False

<br />
<br />
<br />

## Licence

Please note that the code for KmerAperture is distributed under the terms of the GNU GPL v3 license, for more details see https://www.gnu.org/copyleft/gpl.html
