# KmerAperture

![Python](https://badges.aleen42.com/src/python.svg) ![conda](https://img.shields.io/badge/%E2%80%8B-conda-%2344A833.svg?style=flat&logo=anaconda&logoColor=44A833)

Alignment-free estimation of core and accessory distances in closely related genomes.



## Dependencies

KmerAperture.py is written in python3 and requires the following python packages:

> - numpy
> - screed
> - pandas 
> - biopython

## Setup

First set up a conda environment with the appropriate dependencies:

```shell
conda env create -f environment.yml
conda activate env-KmerAperture
```

Compile the kmer parser. This is recommended for speed
```shell
conda install ocaml -c conda-forge
ocamlopt.opt -O3 -o parser/KmerApertureParser parser/KmerApertureParser.ml -ccopt -static
```
Alternatively KmerAperture can be run with ```--pyonly``` and will take longer to read the genomes

## Usage

Run the main script `KmerAperture` referencing your fasta directory and reference genome (required):

```shell
KmerAperture --fastas <fasta dir> --reference <ref file fasta>
```

The following options are also available:

Flag &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Short flag | Description | Required | Default val
--------------|------------|-------------|----------|--------------
`--fastas` |  `-f` |  Provide path to query fastas directory | ✅
`--reference` |     `-r` |  Path and file name of fasta reference genome | ✅ | 
`--kmersize` |      `-k` |  k size |                             | 21
`--pyonly` |     `-py` |  Run without fast ocaml kmer parser |   | False



You may also precluster your genomes (cf below) if you suspect them of being relatively diverse (such as species-wide)


## Input

A reference genome and a directory of assembled query genomes (fasta format)


## Output

Filename | Description | 
--------------|-------|
`SNPmatrix.polymorphic.csv` |  Matrix of polymorphic sites for all genomes against the reference |
`polymorphicsites.fasta` |     Multifasta file of polymorphic sites compared with the reference  |
`core_alignment.fasta` |      Core genome SNP alignment. For sequence present in ≥95% of genomes including the reference |
`accessory_coords.json` |     Coordinates of sequence accessory to each query genome compared with the reference |
`{referencename}_{k}.csv` | Comma separated results for each genomes SNP and indel count compared with reference and accessory size |

To use the core genome SNPs for a phylogeny with branch lengths corrected for recombination use [iqtree](https://github.com/Cibiv/IQ-TREE) and [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML):

```shell
conda install -c conda-forge -c bioconda -c defaults clonalframeml
conda install -c bioconda iqtree
iqtree -s core_alignment.fasta
ClonalFrameML core_alignment.fasta.treefile core_alignment.fasta kmeraperture
```

<br />
<br />

## Pre-cluster (suspected) diverse genomes prior to KmerAperture

Install sourmash (https://github.com/sourmash-bio/sourmash):

```shell
conda install -c conda-forge -c bioconda sourmash
```

Run:

```shell
python sourmash_precluster.py --fastas <fasta dir> 
```

The output is a hierarchical dendrogram (UPGMA) of your genomes MinHash (Jaccard) distances. Select a clustering threshold that cuts across long branches (major lineages). Then cluster based on this threshold with:

```shell
python sourmash_precluster.py --fastas <fasta dir> --threshold <threshold>
```

Obtain the threshold for clustering by viewing the dendrogram and selecting a distance that cuts horizontally through blue lines only. You may also assess the reference placement with preclustering and consider a more appropriate reference for genome clusters if needed.


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
