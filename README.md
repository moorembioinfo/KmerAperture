# KmerAperture

![Python](https://badges.aleen42.com/src/python.svg) 

Alignment-free estimation of SNP distances in closely related genomes



## Dependencies

KmerAperture.py is written in python3 and requires the following python packages:

> - numpy
> - screed

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



## Options

Flag &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Short flag | Description | Required | Default val
--------------|------------|-------------|----------|--------------
`--fastas` |  `-f` |  Provide path to query fastas directory | ✅
`--reference` |     `-r` |  Path and file name of fasta reference genome | ✅ | 
`--kmersize` |      `-k` |  k size |                             | 21
`--sensitive` |        |  Run in sensitive mode | | False

<br />
<br />
<br />
