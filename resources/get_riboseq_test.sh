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

#Commands to download and run Ribosome Profiling data through RiboPipe (used for figure 3)

USER=$1
LOCATION=$2

#Download SRA Toolkit
conda install -c bioconda sra-tools

#Download Ingolia, 2009 data, GSE13750
#Already trimmed?
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13750

mkdir $LOCATION/ingolia_fastq
cd $LOCATION/ingolia_fastq

#mRNA_rich_2
prefetch -v SRR014386
fastq-dump --outdir $LOCATION/ingolia_fastq $USER/ncbi/public/sra/SRR014386.sra
mv SRR014386.fastq SRR014386_mRNA_rich_2_1.fastq

prefetch -v SRR014387
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014387.sra
mv SRR014387.fastq SRR014387_mRNA_rich_2_2.fastq

prefetch -v SRR028774
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR028774.sra
mv SRR028774.fastq SRR028774_mRNA_rich_2_3.fastq

cat SRR014386_mRNA_rich_2_1.fastq SRR014387_mRNA_rich_2_2.fastq SRR028774_mRNA_rich_2_3.fastq > ingolia_a_rich_2_RNA.fastq
gzip ingolia_a_rich_2_RNA.fastq

#mRNA_rich_1
prefetch -v SRR014385
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014385.sra
mv SRR014385.fastq ingolia_a_rich_1_RNA.fastq

gzip ingolia_a_rich_1_RNA.fastq

#mRNA_starved_2
prefetch -v SRR014383
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014383.sra
mv SRR014383.fastq SRR014383_mRNA_starved_2_1.fastq

prefetch -v SRR014384
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014384.sra
mv SRR014384.fastq SRR014384_mRNA_starved_2_1.fastq

cat SRR014383_mRNA_starved_2_1.fastq SRR014384_mRNA_starved_2_1.fastq > ingolia_b_starved_2_RNA.fastq
gzip ingolia_b_starved_2_RNA.fastq

#mRNA_starved_1
prefetch -v SRR014382
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014382.sra
mv SRR014382.fastq ingolia_b_starved_1_RNA.fastq

gzip ingolia_b_starved_1_RNA.fastq

#footprints_rich_2
prefetch -v SRR014377
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014377.sra
mv SRR014377.fastq SRR014377_footprints_rich_2_1.fastq

prefetch -v SRR014378
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014378.sra
mv SRR014378.fastq SRR014378_footprints_rich_2_2.fastq

prefetch -v SRR014379
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014379.sra
mv SRR014379.fastq SRR014379_footprints_rich_2_3.fastq

prefetch -v SRR014380
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014380.sra
mv SRR014380.fastq SRR014380_footprints_rich_2_4.fastq

prefetch -v SRR014381
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014381.sra
mv SRR014381.fastq SRR014381_footprints_rich_2_5.fastq

cat SRR014377_footprints_rich_2_1.fastq SRR014378_footprints_rich_2_2.fastq SRR014379_footprints_rich_2_3.fastq SRR014380_footprints_rich_2_4.fastq SRR014381_footprints_rich_2_5.fastq > ingolia_a_rich_2_FP.fastq
gzip ingolia_a_rich_2_FP.fastq

#footprints_rich_1
prefetch -v SRR014374
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014374.sra
mv SRR014374.fastq SRR014374_footprints_rich_1_1.fastq

prefetch -v SRR014375
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014375.sra
mv SRR014375.fastq SRR014375_footprints_rich_1_2.fastq

prefetch -v SRR014376
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014376.sra
mv SRR014376.fastq SRR014376_footprints_rich_1_3.fastq

cat SRR014374_footprints_rich_1_1.fastq SRR014375_footprints_rich_1_2.fastq SRR014376_footprints_rich_1_3.fastq > ingolia_a_rich_1_FP.fastq
gzip ingolia_a_rich_1_FP.fastq

#footprints_starved_2
prefetch -v SRR014370
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014370.sra
mv SRR014370.fastq SRR014370_footprints_starved_2_1.fastq

prefetch -v SRR014371
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014371.sra
mv SRR014371.fastq SRR014371_footprints_starved_2_2.fastq

prefetch -v SRR014372
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014372.sra
mv SRR014372.fastq SRR014372_footprints_starved_2_3.fastq

prefetch -v SRR014373
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014373.sra
mv SRR014373.fastq SRR014373_footprints_starved_2_4.fastq

cat SRR014370_footprints_starved_2_1.fastq SRR014371_footprints_starved_2_2.fastq SRR014372_footprints_starved_2_3.fastq SRR014373_footprints_starved_2_4.fastq > ingolia_b_starved_2_FP.fastq
gzip ingolia_b_starved_2_FP.fastq

#footprints_starved_1
prefetch -v SRR014368
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014368.sra
mv SRR014368.fastq SRR014368_footprints_starved_1_1.fastq

prefetch -v SRR014369
fastq-dump --outdir ./ $USER/ncbi/public/sra/SRR014369.sra
mv SRR014369.fastq SRR014369_footprints_starved_1_2.fastq

cat SRR014368_footprints_starved_1_1.fastq SRR014369_footprints_starved_1_2.fastq > ingolia_b_starved_1_FP.fastq
gzip ingolia_b_starved_1_FP.fastq
