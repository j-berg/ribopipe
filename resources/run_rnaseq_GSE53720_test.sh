#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o slurmjob-%j
#SBATCH --partition=kingspeak

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

#Path to ncbi/public/sra where sra-tools prefetch saves files
USER=$1

#Download SRA Toolkit
#conda install -c bioconda sra-tools

#Need to have set up whatever package we need to run with before
#Ensure only running locally installed conda cohort
module purge

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

#import ripopipe
mkdir $SCRDIR/raw_dang_berger
mkdir $SCRDIR/out_dang_berger

#catenate raw data if needed
cd $SCRDIR/raw_dang_berger

prefetch -v
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR1066660.sra
mv SRR1066660.fastq SRR1066660_WT_CR_B.fastq
rm $USER/ncbi/public/sra/SRR1066660.sra

prefetch -v
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR1066659.sra
mv SRR1066659.fastq SRR1066659_WT_CR_A.fastq
rm $USER/ncbi/public/sra/SRR1066659.sra

prefetch -v
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR1066658.sra
mv SRR1066658.fastq SRR1066658_WT_NR_B.fastq
rm $USER/ncbi/public/sra/SRR1066658.sra

prefetch -v
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR1066657.sra
mv SRR1066657.fastq SRR1066657_WT_NR_A.fastq
rm $USER/ncbi/public/sra/SRR1066657.sra

#run script --- add directory of raw data
#Trimming Universal miRNA cloning linker (New England Biolabs) as specified in linked protocol (Ingolia,2010)
ribopipe rnaseq -i $SCRDIR/raw_dang_berger/ -o $SCRDIR/out_dang_berger/ -r yeast -e dang_berger_2014 \
  -s WT_NR_A WT_NR_B WT_CR_A WT_CR_B \
  -p HISAT2 -a TGGAATTCTCGGGTGCCAAGG --platform ILLUMINA --full_genome --replicates
