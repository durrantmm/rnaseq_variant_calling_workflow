# rnaseq_variant_calling_workflow

## Download
Clone the git repository with the command:

    git clone git@github.com:durrantmm/rnaseq_variant_calling_workflow.git

## Installation

This workflow uses a conda environment to satisfy all the necessary dependencies.
 
Make sure you have anaconda installed. Download it [here](https://www.continuum.io/downloads).

You should be able to install the workflow simply by entering:

    bash install.sh

In your console. Now it's time to move on to configuration.

## Configuration
This is a very important step in running the workflow properly.

Open the provided `config.yaml` file to get started

### Set GATK and Picard execution paths
The first two parameters of the `config.yaml` file are

    gatk_path: "java -jar /path/to/GenomeAnalysisTK.jar"
    picard_path: "java -jar /path/to/Picard.jar"

These are two strings that allow the workflow to execute GATK and Picard.
You can download these from 
https://broadinstitute.github.io/picard/ 
and 
https://software.broadinstitute.org/gatk/download/

GATK requires you to make an account to download.

Once they are both installed, you can enter their execution paths in the `config.yaml` file.

You may want to add an addition flag, `-Xmx50G`, to make sure that Java allocates enough memory. The final
`config.yaml` file would look something like this:

    gatk_path: "java -jar -Xmx50G /path/to/GenomeAnalysisTK.jar"
    picard_path: "java -jar -Xmx50G /path/to/Picard.jar"
    
### Set working directory
The next parameter of the `config.yaml` file is `wd: test`. This is the default working directory, and it contains
some test data to try out. You can give this a try by typing
    
    snakemake --cores 2
    
It should take a few minutes to run, but once it runs correctly you know you are all set.

## Setting up a workflow
It is essential that you set up your workflow properly before you run it. Follow these directions exactly.

#### IMPORTANT RULES AND TIPS TO REMEMBER:
* When naming files, periods must only be used where specified. They cannot be used as a general name delimiter.
    * `myfile.R1.fq.gz` - CORRECT - this is the right way to name a fastq file
    * `myfile.Run1.Replicate2.R1.fq.gz` - INCORRECT - this is the WRONG way to name a fastq file.
    * `myfile_Run1_Replicate2.R1.fq.gz` - CORRECT - Replacing the improper periods with underscores is appropriate.
* Files can be be symbolic links to save disk space.

### Overview
Navigate to your working directory. Make sure you have chosen a file system that has lots of free disk space available.

From within your working directory, create 5 directories as follows:

    mkdir 0.fastq;
    mkdir 0.reference_genome_fasta;
    mkdir 0.gencode;
    mkdir 0.dbsnp;
    mkdir 0.adar_sites;
    
Now we'll go through each of the directories and describe how to populate them.

### `0.fastq`
This is the directory that contains all of the paired-end FASTQ files to be analyzed. Single-end analysis is possible,
but hasn't been built into this workflow automatically. You can figure that one out on your own.

These fastq files must follow this precise format.

    <sample1>.R1.fq.gz
    <sample1>.R2.fq.gz
    <sample2>.R1.fq.gz
    <sample2>.R2.fq.gz
    ...
    
This folder contains as many samples as you want. Make sure that the only periods are found in the `.R#.fq.gz` file extension.

### `0.reference_genome_fasta`
This is the reference genome fasta file to be used in this analysis. You can have as many reference genomes as you want.
Each genome must also have matching files in `0.gencode`, `0.dbsnp`, and `0.adar_sites`. The identifying name must
be the same between these directories. 

They must be named as follows:

    <genome>.fa

They cannot be gzipped. 

You can download the hg19 human genome into this directory by typing:

    wget https://s3-us-west-1.amazonaws.com/mdurrant/biodb/genomes/human/hg19/hg19.fa.gz

And then uncompress it:

    gunzip hg19.fa.gz

### `0.gencode`
This folder contains the gencode GTF file to be used by the STAR aligner.
 
It must follow the naming convention:

    <genome>.gtf
    
It must have the same `<genome>` name as the corresponding reference genome in 
`0.reference_genome_fasta`

A version corresponding to hg19 can be downloaded from http://www.gencodegenes.org/releases/19.html

You can also download one with wget as:

    wget https://s3-us-west-1.amazonaws.com/mdurrant/biodb/gencode/gtf/hg19.gtf.gz

And uncompress it:

    gunzip hg19.fa.gz
    
### `0.dbsnp`
This folder contains the dbsnp VCF file.

It must follow the naming convention

    <genome>.vcf

It must have the same `<genome>` name as the corresponding reference genome in 
`0.reference_genome_fasta` and `0.gencode`.

A version of dbsnp corresponding to hg19 can be downloaded from 
 
    wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz

or

    wget https://s3-us-west-1.amazonaws.com/mdurrant/biodb/dbsnp/gtf/hg19.gtf.gz
    

### `0.adar_sites`
This folder contains known ADAR A-to-I RNA editing sites. They are excluded outright from the
resulting RNA-seq variant calls.

It must follow the naming convention

    <genome>.txt

It must have the same `<genome>` name as the corresponding reference genome in 
`0.reference_genome_fasta`, `0.gencode`, and `0.dbsnp`.

This was originally obtained from the website http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19_v2.txt 

You can download one from 
 
    wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz
    
## Run Snakemake
You should now be able to run the snakemake workflow.

You can run this with the command

    snakemake
    
From the `rnaseq_variant_calling_workflow` directory.

You can set the number of cores used with

    snakemake --cores 6
