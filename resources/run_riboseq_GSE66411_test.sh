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

#Path to ncbi/public/sra where sra-tools prefetch saves files
NCBI=$1

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
cd $SCRDIR/raw_ingolia


prefetch -v SRR1822476
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822476.sra
mv SRR1822476.fastq GSE66411_p_tif1-ts_replicate_2_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822476.sra

prefetch -v SRR1822475
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822475.sra
mv SRR1822475.fastq GSE66411_o_tif1-ts_replicate_1_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822475.sra

prefetch -v SRR1822474
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822474.sra
mv SRR1822474.fastq GSE66411_n_wild-type_TIF1_replicate_2_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822474.sra

prefetch -v SRR1822473
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822473.sra
mv SRR1822473.fastq GSE66411_m_wild-type_TIF1_replicate_1_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822473.sra

prefetch -v SRR1822472
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822472.sra
mv SRR1822472.fastq GSE66411_p_tif1-ts_replicate_2_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822472.sra

prefetch -v SRR1822471
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822471.sra
mv SRR1822471.fastq GSE66411_o_tif1-ts_replicate_1_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822471.sra

prefetch -v SRR1822470
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822470.sra
mv SRR1822470.fastq GSE66411_n_wild-type_TIF1_replicate_2_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822470.sra

prefetch -v SRR1822469
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822469.sra
mv SRR1822469.fastq GSE66411_m_wild-type_TIF1_replicate_1_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822469.sra


prefetch -v SRR1822468
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822468.sra
mv SRR1822468.fastq GSE66411_l_tif1-ts_replicate_2_30_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822468.sra

prefetch -v SRR1822467
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822467.sra
mv SRR1822467.fastq GSE66411_k_tif1-ts_replicate_1_30_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822467.sra

prefetch -v SRR1822466
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822466.sra
mv SRR1822466.fastq GSE66411_j_wild-type_TIF1_replicate_2_30_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822466.sra

prefetch -v SRR1822465
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822465.sra
mv SRR1822465.fastq GSE66411_i_wild-type_TIF1_replicate_1_30_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822465.sra

prefetch -v SRR1822464
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822464.sra
mv SRR1822464.fastq GSE66411_l_tif1-ts_replicate_2_30_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822464.sra

prefetch -v SRR1822463
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822463.sra
mv SRR1822463.fastq GSE66411_k_tif1-ts_replicate_1_30_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822463.sra

prefetch -v SRR1822462
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822462.sra
mv SRR1822462.fastq GSE66411_j_wild-type_TIF1_replicate_2_30_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822462.sra

prefetch -v SRR1822461
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822461.sra
mv SRR1822461.fastq GSE66411_i_wild-type_TIF1_replicate_1_30_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822461.sra



prefetch -v SRR1822460
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822460.sra
mv SRR1822460.fastq GSE66411_h_ded1-ts_replicate_2_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822460.sra

prefetch -v SRR1822459
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822459.sra
mv SRR1822459.fastq GSE66411_g_ded1-ts_replicate_1_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822459.sra

prefetch -v SRR1822458
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822458.sra
mv SRR1822458.fastq GSE66411_f_wild-type_DED1_replicate_2_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822458.sra

prefetch -v SRR1822457
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822457.sra
mv SRR1822457.fastq GSE66411_e_wild-type_DED1_replicate_1_37_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822457.sra

prefetch -v SRR1822456
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822456.sra
mv SRR1822456.fastq GSE66411_h_ded1-ts_replicate_2_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822456.sra

prefetch -v SRR1822455
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822455.sra
mv SRR1822455.fastq GSE66411_g_ded1-ts_replicate_1_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822455.sra

prefetch -v SRR1822454
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822454.sra
mv SRR1822454.fastq GSE66411_f_wild-type_DED1_replicate_2_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822454.sra

prefetch -v SRR1822453
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822453.sra
mv SRR1822453.fastq GSE66411_e_wild-type_DED1_replicate_1_37_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822453.sra



prefetch -v SRR1822452
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822452.sra
mv SRR1822452.fastq GSE66411_d_ded1-cs_replicate_2_15_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822452.sra

prefetch -v SRR1822451
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822451.sra
mv SRR1822451.fastq GSE66411_c_ded1-cs_replicate_1_15_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822451.sra

prefetch -v SRR1822450
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822450.sra
mv SRR1822450.fastq GSE66411_b_wild-type_DED1_replicate_2_15_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822450.sra

prefetch -v SRR1822449
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822449.sra
mv SRR1822449.fastq GSE66411_a_wild-type_DED1_replicate_1_15_deg_RNA.fastq
rm $NCBI/ncbi/public/sra/SRR1822449.sra

prefetch -v SRR1822448
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822448.sra
mv SRR1822448.fastq GSE66411_d_ded1-cs_replicate_2_15_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822448.sra

prefetch -v SRR1822447
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822447.sra
mv SRR1822447.fastq GSE66411_c_ded1-cs_replicate_1_15_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822447.sra

prefetch -v SRR1822446
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822446.sra
mv SRR1822446.fastq GSE66411_b_wild-type_DED1_replicate_2_15_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822446.sra

prefetch -v SRR1822445
fastq-dump --outdir ./ $NCBI/ncbi/public/sra/SRR1822445.sra
mv SRR1822445.fastq GSE66411_a_wild-type_DED1_replicate_1_15_deg_FP.fastq
rm $NCBI/ncbi/public/sra/SRR1822445.sra



cd $SCRDIR/

#run script --- add directory of raw data
#Trimming Universal miRNA cloning linker (New England Biolabs) as specified in linked protocol (Ingolia,2010)
ribopipe riboseq -i $SCRDIR/raw_ingolia/ -o $SCRDIR/out_ingolia/ -r yeast -e ingolia_2015 \
  -s a_WT_DED1_1_15deg b_WT_DED1_2_15deg c_ded1_cs_1_15deg d_ded1_cs_2_15deg \
  e_WT_DED1_1_37deg f_WT_DED1_2_37deg g_ded1_ts_1_37deg h_ded_ts_2_37deg \
  i_WT_TIF1_1_30deg j_WT_TIF1_2_30deg k_tif1_ts_1_30deg l_tif1_ts_2_30deg \
  m_WT_TIF1_1_37deg n_WT_TIF1_2_37deg o_tif1_ts_1_37deg p_tif1_ts_2_37deg \
  -p HISAT2 -a CTGTAGGCACCATCAAT --platform ILLUMINA --count_cutoff 128
