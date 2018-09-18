# BIO694 - Next generation sequencing (NGS) II. Transcriptomes, Variant Calling and Biological Interpretation
## September 17-18th 2018
## University of Zürich (UZH) 

![alt text](https://github.com/carlalbc/URPP_tutorials/blob/master/img/Logo_URPP_kl2.png)

## URPP Evolution in action

# Part II.- RNA sequencing: Transcriptomes and differential gene expression analyses

## Step 1: Read mapping with STAR aligner:

It is recommended using the STAR aligner for all genomes where there are no alternative alleles. For genomes such as hg38 that have alt alleles, hisat2 should be used as it handles the alts correctly and STAR does not yet. Use Tophat2 only if you do not have enough RAM available to run STAR (about 30 GB). The documentation for STAR is available [here](https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf)

Today we will work with data from the Zebrafish at different stages of differentiation. The data files are contained in the subdirectory called data and are the following:

- 2cells_1.fastq and 2cells_2.fastq: these files are based on RNA-seq data of a 2-cell zebrafish embryo, and
- 6h_1.fastq and 6h_2.fastq: these files are based on RNA-seq data of zebrafish embryos 6h post fertilization.

## Step 2: Install STAR aligner

```
sudo apt install rna-star
```

- The documentation for STAR is available [here](https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf)

## Step 3: Prepare all the directories and download the data

- Remember you can go from one directory to the next using ***cd***. Now let's create a new directory called STARGenome.

```
mkdir danRer 
mkdir danRer/reference danRer/data
cd danRer
wget ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_1.fastq ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_2.fastq ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_1.fastq ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_2.fastq -P data
wget ftp.ensembl.org/pub/release-93/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz -P reference
wget ftp://ftp.ensembl.org/pub/release-93/gtf/danio_rerio/Danio_rerio.GRCz11.93.gtf.gz -P reference
```

- If you feel adventurous you can run FastQC on all the reads in the **data** folder we just created using **cd** and find out the read length and quality:

```
fastqc *.gz
firefox *.html
```

## Step 4: Generate the genome index with STAR

- Run STAR in "genomeGenerate" mode

```
 STAR  [options]... --genomeDir REFERENCE   --readFilesIn R1.fq R2.fq
```

```
mkdir STARindex
gunzip reference/*.gz
STAR --runMode genomeGenerate --genomeDir STARindex --genomeFastaFiles reference/Danio_rerio.GRCz11.dna.toplevel.fa.gz --sjdbGTFfile reference/Danio_rerio.GRCz11.93.gtf.gz --sjdbOverhang 75 --runThreadN 1 
```



14 -- runThreadN 1 \ # can be used to define more processors
A survey of best practices for RNA-seq data analysis
Ana ConesaEmail author, Pedro MadrigalEmail author, Sonia Tarazona, David Gomez-Cabrero, Alejandra Cervera, Andrew McPherson, Michał Wojciech Szcześniak, Daniel J. Gaffney, Laura L. Elo, Xuegong Zhang and Ali MortazaviEmail author
Genome Biology201617:13
https://doi.org/10.1186/s13059-016-0881-8
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8

