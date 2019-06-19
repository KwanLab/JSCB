# JSCB

This is a fork of an implementation of the JSCB algorithm developed by Azad and Lawrence ([Nucleic Acids Research **2007**, *35*, 4629-4639](https://academic.oup.com/nar/article/35/14/4629/1012074)) and further validated by Mehul *et al.* ([Open Biology **2017**, *7*, 170094](https://royalsocietypublishing.org/doi/10.1098/rsob.170094)). Please cite these papers if you use this code. We have added a wrapper script that can handle draft quality genomes (i.e. with multiple contigs).

Installation
------------
First download the repository using `git clone`, then add the repository's directory to your `$PATH` variable. The `jscb` file is a Fortran executable compiled for Linux x86. If you have some other system, you will have to recompile the `jscb.f` file. On Ubuntu or other Debian-derived linux distributions you will have to install a few modules before running:

```bash
sudo apt-get install build-essential libgfortran3
```

Usage
-----
First annotate your genome (we recommend the [Prokka](https://github.com/tseemann/prokka) pipeline). Then feed the resulting genbank file into the `run_jscb.py` script:

```bash
run_jscb.py --genbank input_genome.gbk --output_dir directory/path
```

In the output directory, there will be a tab-delimited table called `genomic_islands_summary.tsv`, which has the following columns:

* Genomic island ID (GI_ID)
* Contig
* Start coordinate
* End coordinate
* Genes (comma-delimited list)

There will also be a file called `hgt_genes_summary.tsv`, which has the following columns:

* Locus tag
* Cluster ID
* Cluster size (number of genes)
* Gene length

Original README.md contents below
---------------------------------

Usage:
```
python jscb.py inputfile.gb
```

Input:
JS-CB takes genbank(full) as input file

Output:
JS-CB outputs three tab separated files:
1) JSCB_output.tsv: Has three columns- Genomic Island (GI)-ID, GI-start co-ordinate and GI-end co-ordinate.
2) JSCB_ouput.gi: Contains four columns- Gene number, Cluster id, Cluster size and Gene length.
3) JSCB_output.clus: Contains four columns- Cluster id, Cluster size, Gene number, and Gene length.

Genomic island prediction can be found in JSCB_output.tsv file. JSCB_output.gi/.clus files contains cluster configuration information.

Requirements:
The program requires Biopython. The binary JSCB file is compiled using F95 (Fortran95) compiler
