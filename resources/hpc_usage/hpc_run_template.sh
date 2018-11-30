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

#Set up whatever package we need to run with
#Requires system or local installs of these programs
module load hisat2 samtools python3 fastx_toolkit fastqc picard plastid star bedtools deeptools plastid
#This does not install htseq, needed for alignment and must manually install if HPC does not provide or ask them to add

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

#import ripopipe
mkdir $SCRDIR/raw_calvin
mkdir $SCRDIR/out_calvin
cp /uufs/chpc.utah.edu/common/home/$USER/raw_calvin/* $SCRDIR/raw_calvin/.


#catenate raw data if needed
cd $SCRDIR/.

#run script --- add directory of raw data
#Using Qiagen Illumina 3' adaptor, but shouldn't need adaptor trimming anyways, just a placeholder for now
ribopipe riboseq -i $SCRDIR/raw_calvin/ -o $SCRDIR/out_calvin/ -r yeast -e calvin_test_1 -s cccp_00m cccp_05m cccp_10m cccp_15m cccp_30m cccp_60m -p HISAT2 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --platform ILLUMINA
