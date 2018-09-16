# BIO694 - Next generation sequencing (NGS) II. Genomes, Variant Calling and Biological Interpretation
## September 17-18th 2018, University of Zürich (UZH)
## URPP Evolution in action


# Part I.- DataQC, pre-processing and mapping genomes 

Before we get started --remember, linux is your friend :penguin::penguin::penguin:

Just to remind you some useful commands for Today's workflow:

```
ls                           (listing files)
more                         (view the files)
cd or cd ..                  (access or exit folders)
mv                           (moving files)
cp                           (copying files)
mkdir                        (make directories)
wget                         (download files from the web)  
```

## I) Quality Assesment: Pre-processing the reads

We will be working with data from the [Lenski lab](http://myxo.css.msu.edu/ecoli/genomicsdat.html) and Desai. They are known for making long-term evolution experiments in *E. coli* since the early 00’s. The strain we will work with today is an *E. coli* from a long-term evolution experiment (LTEE). The twelve LTEE populations have been serially propagated
in the same medium for more than 60,000 generations, with samples preserved every 500 generations.  

The sequencing data is from: [Good, B. H., McDonald, M. J., Barrick, J. E., Lenski, R. E., & Desai, M. M. (2017). The dynamics of molecular evolution over 60,000 generations. Nature, 551(7678), 45–50. doi:10.1038/nature24287](https://www.nature.com/articles/nature24287
).

Let's get started! 

## a) Downloading raw sequencing reads from a database.

There are two main databases, the **Sequence Read Archive** (SRA, US based) and the **European nucleotide archive**(ENA, EU based). 

Today we will download the raw reads from the ENA. The project accession name is [PRJNA380528](https://www.ebi.ac.uk/ena/data/view/PRJNA380528). Please, follow the next steps before downloading the reads:
 
1. Create the following directories and go to the main directory (remember it's good to keep things tidy!):

``` 
mkdir fastq                  (creates a folder called fastq)
mkdir fastq/SRR6170103       (creates a subdirectory for the fastq files)
cd fastq/SRR6170103          (goes to the SRR6170103 inside the fastq directory)
mkdir FastQC                 (we will store the results from the QC here later)
``` 

2. Download the fastq files from ENA - we will work with paired-end reads from E. coli:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/003/SRR6170103/SRR6170103_1.fastq.gz -P SRR6170103_1 (gets the first fastq paired-end read from ENA and stores it in the subdirectory we just created)

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/003/SRR6170103/SRR6170103_2.fastq.gz -P SRR6170103_2 (gets second fastq paired-end read file from ENA and stores it in the subdirectory we just created)
```
Let’s check the fastq files:

```
less SRR6170103_1/SRR6170103_1.fastq.gz           (exit with Ctrl+Z)
head SRR6170103_1/SRR6170103_1.fastq.gz           (shows the first 10 lines)
``` 

*You can also use tail to see the end of the file.*

Now that you have seen the files, continue with the rest of the workflow.


b) **Run FASTQC:** when running new software it is always useful to understand it first. A quick glimpse to different options can be obtained by looking at the in-built help:

´´´ 
fastqc --help
´´´
**Now let’s assess the quality of our fastq files containing paired-end sequencing reads by running FastQC:**

We can do this either through the command-line(recommended, jump to option two) or by directly running FastQC and opening the files from the window by running the following in the terminal:
 
1) Graphical Interface option:
```
fastqc &
```

That will open FastQC and you will be able to open the fastq files directly with the software. If you prefer to use the command-line (recommended) do the following:

2) Command-line option:

´´´  
fastqc SRR6170103_1.fastq.gz SRR6170103_2.fastq.gz (wait till it’s done running)
mv *.zip *.html FastQC (move the files to the FastQC folder)
cd FastQC
´´´

**Reminder:** You can always check where you are in the terminal using:

```
pwd          (shows the location where you are in the terminal)
```


* Open the FastQC results with your favorite html visualizer (i.e firefox, chrome, etc.) or if you prefer it, you can open the file through your graphical interface by directly clicking on it.

```
firefox SRR6170103_1_fastqc.html
```
You will see the following:
  



Pretty good quality! :octocat:


You could get *really bad quality reads* like the following:


  

Those are really bad quality reads! :octocat:



* Go through this manual https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf to understand each of the results you have gotten in your report.

## Questions:


1. What’s the total number of sequences in each of the paired-end fastq files?  (Hint: look at the Basic statistics module)

SRR6170103_1 Total number of sequences  _________

SRR6170103_2 Total number of sequences  _________


1. What’s the type of encoding used?
2. What’s the length of the reads? 
3. What are the warnings we get for each of the fastq files? 
4. Which sequence seems to be overrepresented? 

There is an expected drop in quality at the 3’ end of the sequences and also to get an overrepresentation of adaptor sequences. We will trim the low quality ends and remove adaptors next.

**Do not close the FastQC window as we will compare the original files to the ones we will produce after adapter removal and quality filtering.**


c) Trimming, removing adaptors and low quality reads with Trimmomatic: a java tool for performing a range of trimming tasks on Illumina paired end and single end read data. The manual can be found [here](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

1. We use Trimmomatic to remove adapter sequences:

```
java -jar trimmomatic-0.35.jar PE -phred33  SRR6170103_1.fastq.gz SRR6170103_2.fastq.gz SRR6170103_1_paired.fastq.gz SRR6170103_1_unpaired.fastq.gz SRR6170103_2_paired.fastq.gz SRR6170103_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

The parameters used for Trimmomatic are defined as follows:
1) **-PE**
data is paired end
2) **-phred33**
Quality scores are 33 offset
3) **-threads 1**
number of threads to use
4) **-trimlog logfile**
name of logfile for summary information
5)**paired_end1.fastq**
name of input fastq file for left reads
**paired_end2.fastq**
name of input fastq file for right reads
**Left_paired.fastq**
paired trimmed output fastq file for left reads
**Left_unpaired.fastq**
unpaired trimmed output fastq file for left reads
**Right_paired.fastq**
paired trimmed output fastq file for right reads
**Right_unpaired.fastq**
unpaired trimmed output fastq file for right reads
**ILLUMINACLIP**
parameters for the adapter clipping
**TruSeq3-PE-2.fa** 
text file of adapter sequences to search for
**:2:40:15**
adapter-read alignment settings – see manual
**MINLEN:36**
delete reads trimmed below length MINLEN

## Questions: 

1. While you wait, what does :2:40:15 means? (Check Trimmomatic’s manual)
2. How many reads survive after Trimmomatic? (Hint: Check the messages left in the terminal after Trimmomatic)

b) Using Trimmomatic to filter low quality reads

Because we saw that some of the reads were low quality at the beginning-end we will remove them

## II) Mapping the reads to a reference genome of E.coli using BWA. link to BWA aligner.


