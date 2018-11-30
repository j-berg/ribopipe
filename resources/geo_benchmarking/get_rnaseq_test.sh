#!/bin/sh

#RiboPipe
#An assembly and analysis pipeline for sequencing data
#alias: ripopipe
#Copyright (C) 2018  Jordan A. Berg
#jordan <dot> berg <at> biochem <dot> utah <dot> edu
#This program is free software: you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your option) any later
#version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with
#this program.  If not, see <https://www.gnu.org/licenses/>.

#Commands to download and run short read SE RNAseq data through RiboPipe (used for figure 4)

USER=$1
LOCATION=$2

#Download SRA Toolkit
conda install -c bioconda sra-tools

#Download Dou, 2017 data, GSE99028
#Single-ended, 75 bp reads were mildly trimmed using Trimmomatic (version 0.32)
#to remove leading or trailing nucleotides whose sequencing quality was below 3.
#Reads whose length fell below 30 bp after trimming were also removed from
#downstream analysis. STAR (version 2.3.0e) was used for mapping reads to
#reference genome (hg19), requiring a minimum alignment score of 10. The
#expression level of RefSeq genes was quantified using featureCounts
#(version 1.5.0) and normalized using DESeq2.
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99028

mkdir $LOCATION/rnaseq_fastq
cd $LOCATION/rnaseq_fastq

#sh-cGAS etoposide-replicate 2; Homo sapiens; RNA-Seq
prefetch -v SRR5574056
fastq-dump --outdir $LOCATION/rnaseq_fastq $USER/ncbi/public/sra/SRR5574056.sra
mv SRR5574056.fastq rnaseq_b_cGAS_2_RNA.fastq
gzip rnaseq_b_cGAS_2_RNA.fastq

#sh-cGAS etoposide-replicate 1; Homo sapiens; RNA-Seq
prefetch -v SRR5574055
fastq-dump --outdir $LOCATION/rnaseq_fastq $USER/ncbi/public/sra/SRR5574055.sra
mv SRR5574055.fastq rnaseq_b_cGAS_1_RNA.fastq
gzip rnaseq_b_cGAS_1_RNA.fastq

#GSM2630605: sh-NTC etoposide-replicate 2; Homo sapiens; RNA-Seq
prefetch -v SRR5574054
fastq-dump --outdir $LOCATION/rnaseq_fastq $USER/ncbi/public/sra/SRR5574054.sra
mv SRR5574054.fastq rnaseq_c_NTC_2_RNA.fastq
gzip rnaseq_c_NTC_2_RNA.fastq

#GSM2630604: sh-NTC etoposide-replicate 1; Homo sapiens; RNA-Seq
prefetch -v SRR5574053
fastq-dump --outdir $LOCATION/rnaseq_fastq $USER/ncbi/public/sra/SRR5574053.sra
mv SRR5574053.fastq rnaseq_c_NTC_1_RNA.fastq
gzip rnaseq_c_NTC_1_RNA.fastq

#GSM2630603: Proliferating control (PD34); Homo sapiens; RNA-Seq
prefetch -v SRR5574052
fastq-dump --outdir $LOCATION/rnaseq_fastq $USER/ncbi/public/sra/SRR5574052.sra
mv SRR5574052.fastq rnaseq_a_control_PD34_RNA.fastq
gzip rnaseq_a_control_2_RNA.fastq

#GSM2630602: Proliferating control (PD29); Homo sapiens; RNA-Seq
prefetch -v SRR5574051
fastq-dump --outdir $LOCATION/rnaseq_fastq $USER/ncbi/public/sra/SRR5574051.sra
mv SRR5574051.fastq rnaseq_a_control_PD29_RNA.fastq
gzip rnaseq_a_control_1_RNA.fastq
