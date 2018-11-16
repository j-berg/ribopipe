#!/bin/bash

#Create reference folder
mkdir yeast_reference_star
cd yeast_reference_star
mkdir source_files
cd source_files
#Download chromosome fasta files
for CHR in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito; do
    curl -O ftp://ftp.ensembl.org/pub/release-94/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.${CHR}.fa.gz; done

#Download GTF
wget ftp://ftp.ensembl.org/pub/release-94/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.94.gtf.gz

#Download ncrna fasta
wget ftp://ftp.ensembl.org/pub/release-94/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz

#Unzip files
gzip -d *.gz

#make copy of gtf in main folder
cp Saccharomyces_cerevisiae.R64-1-1.94.gtf ../transcripts.gtf
cd ../

#Make refFlat file
conda install ucsc-gtftogenepred
gtftogenepred transcripts.gtf yeast_refFlat.txt

#metagene prep
metagene generate transcripts_cds_start --landmark cds_start --annotation_files transcripts.gtf --downstream 200

#Build
mkdir genome 
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles source_files/*dna.chromosome.*.fa --sjdbGTFfile Saccharomyces_cerevisiae.R64-1-1.94.gtf --genomeSAindexNbases 11

mkdir ncrna
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ncrna/ --genomeFastaFiles source_files/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa --genomeSAindexNbases 7

cd ../
zip -r yeast_reference_star.zip yeast_reference_star/
