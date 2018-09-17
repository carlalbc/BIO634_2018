# BIO694 - Next generation sequencing (NGS) II. Genomes, Variant Calling and Biological Interpretation
## September 17-18th 2018, University of Zürich (UZH)![alt text]

![alt text](https://github.com/carlalbc/URPP_tutorials/blob/master/img/Logo_URPP_kl2.png)

## URPP Evolution in action

# Part I.- DataQC, pre-processing and mapping genomes 

Before we get started --remember, linux is your friend :penguin::penguin::penguin:

- Just to remind you some useful commands for Today's workflow:

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

- The sequencing data is from: [Good, B. H., McDonald, M. J., Barrick, J. E., Lenski, R. E., & Desai, M. M. (2017). The dynamics of molecular evolution over 60,000 generations. Nature, 551(7678), 45–50. doi:10.1038/nature24287](https://www.nature.com/articles/nature24287
).

Let's get started! 

## Step 1: Downloading raw sequencing reads from a database.

- There are two main databases, the **Sequence Read Archive** (SRA, US based) and the **European nucleotide archive** (ENA, EU based). 

- Today we will download the raw reads from the ENA. The project accession name is [PRJNA380528](https://www.ebi.ac.uk/ena/data/view/PRJNA380528). Please, follow the next steps before downloading the fastq files:
 
1. Create the following directories and go to the main directory (remember it's good to keep things tidy!):

``` 
mkdir fastq                  (creates a folder called fastq)
mkdir fastq/SRR6170103       (creates a subdirectory for the fastq files)
cd fastq/SRR6170103          (goes to the SRR6170103 inside the fastq directory)
mkdir FastQC                 (we will store the results from the QC here later)
``` 

2. Download the fastq files from ENA - we will work with paired-end reads from E. coli:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/003/SRR6170103/SRR6170103_1.fastq.gz  (gets the first fastq paired-end read from ENA and stores it in the subdirectory we just created)

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/003/SRR6170103/SRR6170103_2.fastq.gz  (gets the second fastq paired-end read file from ENA and stores it in the subdirectory we just created)
```
- Let’s check the fastq files:

```
less SRR6170103/SRR6170103_1.fastq.gz           (exit with Ctrl+Z)
head SRR6170103/SRR6170103_1.fastq.gz           (shows the first 10 lines)
``` 

- You can also use ***tail*** to see the end of the file.

Now that you have seen the files, continue with the rest of the workflow.


## Step 2: Run FASTQC

When running new software it is always useful to understand it first. A quick glimpse to different options can be obtained by looking at the in-built help:

```  
fastqc --help	

Usage: fastqc seqfile1 seqfile2 .. seqfileN
```

- Now let’s assess the quality of our fastq files containing paired-end sequencing reads by running FastQC:

- We can do this either through the command-line (recommended, jump to option two) or by directly running FastQC and opening the files from the window by running the following in the terminal:
 
- Graphical Interface option 1:

```
fastqc &
```

That will open FastQC and you will be able to open the fastq files directly with the program. If you prefer to use the command-line (recommended) do the following:

- Command-line option 2:

```  
fastqc SRR6170103_1.fastq.gz SRR6170103_2.fastq.gz                       (wait till it’s done running)
mv *.zip *.html FastQC                                                   (move the files to the FastQC folder)
cd FastQC
``` 

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


2. What’s the type of encoding used?
3. What’s the length of the reads? 
4. What are the warnings we get for each of the fastq files? 
5. Which sequence seems to be overrepresented? 

There is an expected drop in quality at the 3’ end of the sequences and also to get an overrepresentation of adaptor sequences. We will trim the low quality ends and remove adaptors next.

**Do not close the FastQC window as we will compare the original files to the ones we will produce after adapter removal and quality filtering.**


c) Trimming, removing adaptors and low quality reads with Trimmomatic: a java tool for performing a range of trimming tasks on Illumina paired end and single end read data. The manual can be found [here](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

## Step 3: Use Trimmomatic to remove adapter sequences

```
java -jar trimmomatic-0.35.jar PE -phred33  SRR6170103_1.fastq.gz SRR6170103_2.fastq.gz SRR6170103_1_paired.fastq.gz SRR6170103_1_unpaired.fastq.gz SRR6170103_2_paired.fastq.gz SRR6170103_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

The parameters used for Trimmomatic are defined as follows:
1) **PE**
data is paired end
2) **phred33**
Quality scores are 33 offset
3) **threads 1**
number of threads to use
4) **trimlog logfile**
name of logfile for summary information
5) **paired_end1.fastq**
name of input fastq file for left reads
6) **paired_end2.fastq**
name of input fastq file for right reads
7) **Left_paired.fastq**
paired trimmed output fastq file for left reads
8) **Left_unpaired.fastq**
unpaired trimmed output fastq file for left reads
9) **Right_paired.fastq**
paired trimmed output fastq file for right reads
10) **Right_unpaired.fastq**
unpaired trimmed output fastq file for right reads
11) **ILLUMINACLIP**
parameters for the adapter clipping
12) **TruSeq3-PE-2.fa** 
text file of adapter sequences to search for
13) **:2:40:15**
adapter-read alignment settings – see manual
14) **MINLEN:36**
delete reads trimmed below length MINLEN

## Questions: 

1. What does :2:40:15 means?      (Check Trimmomatic’s manual)
2. How many reads survive after Trimmomatic? (Hint: Check the messages left in the terminal after Trimmomatic)


**NOTE: Remember that trimmomatic only deletes reads if the length after trimming of adapter sequences is less than MINLEN (which we set to 36bp)**

## Step 4: Use Trimmomatic to filter low quality reads

Next, we can remove low quality reads of the sequences by trimming the bases at the 3' end of the reads with the following command:

```
java -jar ../../../software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -threads 1 -trimlog logfile2 SRR6170103_1_paired.fastq.gz SRR6170103_2_paired.fastq.gz SRR6170103_1_trim_paired.fastq SRR6170103_1_unpaired.fastq SRR6170103_2_trim_paired.fastq SRR6170103_2_trim_unpaired.fastq SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
```
The parameters used for Trimmomatic are defined as follows:

1) **PE**
data is paired end
2) **phred33**
Quality scores are 33 offset
3) **threads 1** 
number of threads to use
4)**trimlog logfile2**
name of logfile for summary information
5)**Left_paired.fastq**
name of input adapter trimmed left fastq file
6)**Right_paired.fastq**
name of input adapter trimmed right fastq file
7)**Left_trim_paired.fastq**
paired trimmed output fastq file for left reads
8)**Left_unpaired.fastq**
unpaired trimmed output fastq file for left reads
9)**Right_paired.fastq**
paired trimmed output fastq file for right reads
10)**Right_unpaired.fastq**
unpaired trimmed output fastq file for right reads
11)**LEADING:3**
Trim 5’ bases with quality score < 3
12)**TRAILING:3**
Trim 3’ bases with quality score < 3
13)**SLIDINGWINDOW:4:15** see manual for explanation
14)**MINLEN:36**
delete reads trimmed below length MINLEN


#Questions:

- What was the screenlog of Trimmomatic this time? How many reads were removed? Whas there any?
- What does the SLIDINGWINDOW:4:15 means? Check the manual.

## Step 5: Use SOAPec to correct sequencing errors

If not installed already, download SOAPec

```
mkdir software
wget http://sourceforge.net/projects/soapdenovo2/files/ErrorCorrection/SOAPec_v2.01.tar.gz/download -P software
cd software
tar -zxf SOAPec_v2.01.tar.gz
cd ..
```

Now we need to create a file with the fastq files location for SOAPec to be able to work. Use your favorite text editor (ie. gedit) and write the path of our fastq files. If you can't remember, use ***pwd**.

```
gedit files.txt &
```

It will open up gedit where you have to write the path of the files that you created in Step 4. 
In **my** case it will look like this:

```
/home/vega/URPP_2018/BIO634-2018/PartI/fastq/SRR6170103/SRR6170103_1_trim_paired.fastq
/home/vega/URPP_2018/BIO634-2018/PartI/fastq/SRR6170103/SRR6170103_2_trim_paired.fastq
```
Save the file (Ctrl+S) and exit (Alt+F4). Then, run the KmerFreq_AR command below and when it finishes run the Corrector_AR command.

Depending on your location and where you downloaded and installed SOAPEc the path will change slightly. Remember to use **pwd** to know where you are.

```
/home/vega/URPP_2018/BIO634-2018/software/SOAPec_v2.01/bin/KmerFreq_AR -k 16 -t 1 -q 33 -p Error_Corr files.txt > kmerfreq16.log 2> kmerfreq16.err
/home/vega/URPP_2018/BIO634-2018/software/SOAPec_v2.01/bin/Corrector_AR -k 16 -Q 33 -t 1 -o 3 Error_Corr.freq.cz Error_Corr.freq.cz.len files.txt > Corr16.log 2>Corr16.err
```
## Questions

1. Inspect the output of Corrector_AR output file SRR6170103_1_trim_paired.fastq.cor.stat (Hint, use ***more*** / ***less***)

For SRR6170103_1:

- How many bases were corrected by the "Fast method"?  ______________

- How many bases were corrected by the "BBtree method"?  ______________

What about for SRR6170103_2?

- How many bases were corrected by the "Fast method"?  ______________

- How many bases were corrected by the "BBtree method"?  ______________

- **Are there any differences between the two? What conclusions can you make?**

_____________________________________________________________________________________


**If you closed FastQC before, open it and run it for the trimmed paired-end reads, do they look very different?**


Congratulations you now can continue with the next step of our practice! :octocat:



II) Mapping the reads to a reference genome of *E. coli* using Burrows Wheeler Aligner (BWA).

Here's the [manual](http://bio-bwa.sourceforge.net/bwa.shtml)




