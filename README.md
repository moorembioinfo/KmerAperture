# KmerAperture

![Python](https://badges.aleen42.com/src/python.svg) ![conda](https://img.shields.io/badge/%E2%80%8B-conda-%2344A833.svg?style=flat&logo=anaconda&logoColor=44A833)

Alignment-free estimation of core and accessory distances in closely related genomes.



## Dependencies

KmerAperture.py is written in python3 and requires the following python packages:

> - numpy
> - screed
> - pandas


And for preclustering:

> - matplotlib
> - scipy
> - networkX  


And sourmash: https://github.com/sourmash-bio/sourmash


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

You may also precluster your genomes if you suspect them of being relatively diverse (such as species-wide)

## Input

A reference genome and a directory of assembled query genomes (fasta format)

<br />
<br />


## Options

Flag &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Short flag | Description | Required | Default val
--------------|------------|-------------|----------|--------------
`--fastas` |  `-f` |  Provide path to query fastas directory | ✅
`--reference` |     `-r` |  Path and file name of fasta reference genome | ✅ | 
`--kmersize` |      `-k` |  k size |                             | 21
`--polySNP` |     `-p` |  Optionally output a matrix of estimated polymorphic SNPs sites |   | False

<br />
<br />

## Output

A .csv file with genomeID, estimated SNPs and estimated sequence presence in the reference (acc1) and query genome (acc2)

Optionally output a SNP matrix of polymorphic sites estimated by KmerAperture

...more to come

<br />
<br />

## Pre-cluster (suspected) diverse genomes prior to KmerAperture

Run:

```shell
python sourmash_precluster.py --fastas <fasta dir> 
```

The output is a heirarchical dendrogram (UPGMA) of your genomes MinHash (Jaccard) distances. Select a clustering threshold that cuts across long branches (major lineages). Then cluster based on this threshold with:

```shell
python sourmash_precluster.py --fastas <fasta dir> --refine <threshold>
```
<br />
<br />

## Cite

If you find KmerAperture useful please cite:

KmerAperture: Retaining k-mer synteny for alignment-free estimation of within-lineage core and accessory differences
https://www.biorxiv.org/content/10.1101/2022.10.12.511870v1

```console
@article {Moore2022.10.12.511870,
	author = {Moore, Matthew P and Laager, Mirjam and Didelot, Xavier},
	title = {KmerAperture: Retaining k-mer synteny for alignment-free estimation of within-lineage core and accessory differences},
	elocation-id = {2022.10.12.511870},
	year = {2022},
	doi = {10.1101/2022.10.12.511870},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/10/16/2022.10.12.511870},
	eprint = {https://www.biorxiv.org/content/early/2022/10/16/2022.10.12.511870.full.pdf},
	journal = {bioRxiv}
}
```

If you use preclustering, please also cite sourmash:
https://joss.theoj.org/papers/10.21105/joss.00027#

```console
@article {Brown2016,  
         doi = {10.21105/joss.00027},  
	 url = {https://doi.org/10.21105/joss.00027},  
	 year = {2016},  
	 publisher = {The Open Journal},  
	 volume = {1}, number = {5},  
	 pages = {27},  
	 author = {C. Titus Brown and Luiz Irber},  
	 title = {sourmash: a library for MinHash sketching of DNA}, journal = {Journal of Open Source Software} }
```

<br />
<br />
<br />

## Licence

Please note that the code for KmerAperture is distributed under the terms of the GNU GPL v3 license, for more details see https://www.gnu.org/copyleft/gpl.html
