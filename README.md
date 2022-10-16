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
	abstract = {By decomposing genome sequences into k-mers, it is possible to estimate genome similarity without alignment. Dimension-reduction techniques such as k-mer minimisers (MinHash), have been developed and are often accurate approximations of distances based on full k-mer sets. These and other alignment-free methods avoid the enormous temporal and computational expense of alignment or mapping. However, these k-mer set comparisons are not entirely accurate within-species and can be completely inaccurate within-lineage. This is due, in part, to their inability to distinguish core polymorphism from accessory differences. KmerAperture takes the relative complements of a pair of whole genome k-mer sets and matches back to the original enumerated k-mer lists to gain positional information. SNPs are expected to result in contiguous series of unique k-mers of length L= k. On the other hand, series of length L \&gt; 2k are likely to be caused by accessory differences of length L-k+1; when the start and end of the sequence are contiguous with homologous sequence. KmerAperture was benchmarked against Jaccard similarity and {\textquoteright}split k-mer analysis{\textquoteright} (SKA) using datasets including a diverse lineage (Clostridioides difficile RT017), a lower core diversity sub-lineage with a large accessory genome (Escherichia coli ST1193) and a very low core diversity simulated population with accessory content not associated with number of SNPs. We present a new algorithm that demonstrates that with the few available axioms of how core and accessory sequence diversity is represented in k-mers, we can accurately distinguish them and estimate both. By matching unique k-mers back to their original lists we regain their synteny and may make inferences about their likely cause.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2022/10/16/2022.10.12.511870},
	eprint = {https://www.biorxiv.org/content/early/2022/10/16/2022.10.12.511870.full.pdf},
	journal = {bioRxiv}
}
```


<br />
<br />
<br />

## Licence

Please note that the code for KmerAperture is distributed under the terms of the GNU GPL v3 license, for more details see https://www.gnu.org/copyleft/gpl.html
