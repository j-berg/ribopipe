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

#Commands to download and run Ribosome Profiling data through RiboPipe (used for figure 3)

#Download SRA Toolkit
#conda install -c bioconda sra-tools

#Need to have set up whatever package we need to run with before
#Ensure only running locally installed conda cohort
module purge

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

#import ripopipe
mkdir $SCRDIR/raw_ingolia
mkdir $SCRDIR/out_ingolia

#catenate raw data if needed
cp /scratch/general/lustre/$USER/sra-files/ingolia_raw/* $SCRDIR/raw_ingolia

cd $SCRDIR/

#run script --- add directory of raw data
#Trimming Universal miRNA cloning linker (New England Biolabs) as specified in linked protocol (Ingolia,2010)
ribopipe riboseq -i $SCRDIR/raw_ingolia/ -o $SCRDIR/out_ingolia/ -r yeast -e ingolia_2015 \
  -s a_WT_DED1_1_15deg b_WT_DED1_2_15deg c_ded1_cs_1_15deg d_ded1_cs_2_15deg \
  e_WT_DED1_1_37deg f_WT_DED1_2_37deg g_ded1_ts_1_37deg h_ded_ts_2_37deg \
  i_WT_TIF1_1_30deg j_WT_TIF1_2_30deg k_tif1_ts_1_30deg l_tif1_ts_2_30deg \
  m_WT_TIF1_1_37deg n_WT_TIF1_2_37deg o_tif1_ts_1_37deg p_tif1_ts_2_37deg \
  -p HISAT2 -a CTGTAGGCACCATCAAT --platform ILLUMINA --count_cutoff 128
