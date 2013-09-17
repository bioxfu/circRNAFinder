#! /usr/bin/bash

# This tutorial helps you get started with circRNAFinder by demonstrating how to detect, annotate and visualize circRNAs using public RNA-Seq data. It also could be used as a template to build a pipeline for your own data.

# Step1: Create a working directory
mkdir example
cd example

# Step2: Download and build human genome/transcript sequences and gene annotations. (You only have to do it once.)
# The genome/transcript sequences are downloaded from USCS Genome Browser or Ensembl database in FASTA format and built into Bowtie indices. The latest version of gene annotations are downloaded from Ensembl database in GTF format. The longest transcript for each gene are extracted and converted into a custom format. Currently, SeqAnnoDownloadBuild.py supports the genome of human(hg19), mouse(mm10) and Arabidopsis(TAIR10). 
mkdir build
cd build
SeqAnnoDownloadBuild.py hg19 --download --genome
SeqAnnoDownloadBuild.py hg19 --build --genome
SeqAnnoDownloadBuild.py hg19 --download --transcriptome
SeqAnnoDownloadBuild.py hg19 --build --transcriptome
SeqAnnoDownloadBuild.py hg19 --download --annotation
SeqAnnoDownloadBuild.py hg19 --build --annotation
cd ..

## Step3: Create fasta folder, download the public RNA-Seq data from SRA database
mkdir fasta
cd fasta
## Download .sra files manually at 
# http://www.ncbi.nlm.nih.gov/sra/SRX218203 (RIBOMINUS_HEK293)
# http://www.ncbi.nlm.nih.gov/sra/SRX105931 (CD19)
# http://www.ncbi.nlm.nih.gov/sra/SRX105932 (CD34)
# http://www.ncbi.nlm.nih.gov/sra/SRX105933 (neutrophils)

# Or if you've installed the Aspera, you can download them using the following commands in terminal:

# ~/.aspera/connect/bin/ascp -QTr -k1 -l640M -L log --retry-timeout=10 -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR650/SRR650317/SRR650317.sra . 

# ~/.aspera/connect/bin/ascp -QTr -k1 -l640M -L log --retry-timeout=10 -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR364/SRR364679/SRR364679.sra . 
# ~/.aspera/connect/bin/ascp -QTr -k1 -l640M -L log --retry-timeout=10 -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR364/SRR364680/SRR364680.sra . 
# ~/.aspera/connect/bin/ascp -QTr -k1 -l640M -L log --retry-timeout=10 -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR364/SRR364681/SRR364681.sra . 

# ~/.aspera/connect/bin/ascp -QTr -k1 -l640M -L log --retry-timeout=10 -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR384/SRR384963/SRR384963.sra . 
# ~/.aspera/connect/bin/ascp -QTr -k1 -l640M -L log --retry-timeout=10 -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR384/SRR384964/SRR384964.sra . 
# ~/.aspera/connect/bin/ascp -QTr -k1 -l640M -L log --retry-timeout=10 -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR384/SRR384965/SRR384965.sra . 


## Step4: Convert .sra files into FASTA format, collapse the identical reads
# HEK293 (Paired-end reads) #
sra2fasta.py -p SRR650317.sra
fastx_collapser.py -f SRR650317_1.fasta.gz,SRR650317_2.fasta.gz HEK293.fa

# CD19 (Single-end reads) #
sra2fasta.py -s SRR364679.sra
sra2fasta.py -s SRR384963.sra
fastx_collapser.py -f SRR364679.fasta.gz,SRR384963.fasta.gz CD19.fa

# CD34 (Single-end reads) #
sra2fasta.py -s SRR364680.sra
sra2fasta.py -s SRR384964.sra
fastx_collapser.py -f SRR364680.fasta.gz,SRR384964.fasta.gz CD34.fa

# neutrophils (Single-end reads) #
sra2fasta.py -s SRR364681.sra
sra2fasta.py -s SRR384965.sra
fastx_collapser.py -f SRR364681.fasta.gz,SRR384965.fasta.gz neutrophils.fa
cd ..

## Step5: Edit configure files

# Example circRNAFind.cfg file
#[lib1]
#sample_name = HEK293
#reads_file = ./fasta/HEK293.fa
#genome_index = ./build/hg19_dna
#trans_index = ./build/hg19_trans
#genome_seq = ./build/hg19_dna.fa

#[lib2]
#sample_name = CD34
#reads_file = ./fasta/CD34.fa
#genome_index = ./build/hg19_dna
#trans_index = ./build/hg19_trans
#genome_seq = ./build/hg19_dna.fa

#[lib3]
#sample_name = CD19
#reads_file = ./fasta/CD19.fa
#genome_index = ./build/hg19_dna
#trans_index = ./build/hg19_trans
#genome_seq = ./build/hg19_dna.fa

#[lib4]
#sample_name = neutrophils
#reads_file = ./fasta/neutrophils.fa
#genome_index = ./build/hg19_dna
#trans_index = ./build/hg19_trans
#genome_seq = ./build/hg19_dna.fa

# Example circRNAAnno.cfg file
#[group1]
#sample_names = CD19,CD34,neutrophils,HEK293
#genome_index = ./build/hg19_dna
#genome_seq = ./build/hg19_dna.fa
#gene_bed = ./build/hg19.gtf.build
#id_prefix = hsa

# Example circRNAView.cfg file
#[group1]
#gff_file =  ./group1/table/combLib_annoTab.gff 
#png_folder = ./group1/IGV/
#Genome = hg19

## Step6: Run circRNAFinder
circRNAFind.py circRNAFind.cfg -s 0
circRNAAnno.py circRNAAnno.cfg -s 0
igv.sh &
circRNAView.py circRNAView.cfg

