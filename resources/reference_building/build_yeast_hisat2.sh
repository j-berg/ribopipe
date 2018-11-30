#!/bin/bash

#Create reference folder
mkdir yeast_reference_hisat2
cd yeast_reference_hisat2
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

#Remove any ncrna from gtf file, truncate
cp ../build_truncate_gtf.py ./
python build_truncate_gtf.py
rm build_truncate_gtf.py

#Make refFlat file
conda install -y ucsc-gtftogenepred
gtftogenepred transcripts.gtf transcripts_refFlat.txt
gtftogenepred transcripts_coding.gtf transcripts_coding_refFlat.txt
gtftogenepred transcripts_truncated.gtf transcripts_truncated_refFlat.txt

#remove this and update scipy again to be able to run plastid
conda uninstall -y ucsc-gtftogenepred
conda install -y scipy

#metagene prep
metagene generate transcripts_cds_start --landmark cds_start --annotation_files transcripts.gtf --downstream 200

#Build
cd source_files
hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.II.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.III.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.IV.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.V.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.VI.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.VII.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.VIII.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.IX.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.X.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.XI.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.XII.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.XIII.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.XIV.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.XV.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.XVI.fa,Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.Mito.fa genome

hisat2-build Saccharomyces_cerevisiae.R64-1-1.ncrna.fa ncrna

mv genome.* ../
mv ncrna.* ../

cd ../../
zip -r yeast_reference_hisat2.zip yeast_reference_hisat2/
