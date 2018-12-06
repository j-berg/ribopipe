"""
RiboPipe
An assembly and analysis pipeline for sequencing data
alias: ripopipe

Copyright (C) 2018  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""

"""
IMPORT DEPENDENCIES
"""
import os
import sys
from .utils import file_list, compile_size_distribution
import concurrent.futures

"""
FUNCTIONS
"""
#Run assembly file-wise
def assemble(args):

    #Unpackage file data and check which transcript ref files to use
    file, dir_dict, args_dict, transcripts = args[0], args[1], args[2], args[3]

    #Map to full genome, do ncrna depletion
    if 'full_genome' in args_dict and args_dict['full_genome'] == True:
        if str(args_dict['program']).upper() == 'STAR':
            #STAR -- align curated reads to reference, uses default parameters
            os.system("STAR --runThreadN 1 --genomeDir " + str(dir_dict['reference']) + "genome --readFilesIn " + str(dir_dict['trimdir']) + file + " --outFileNamePrefix " + str(dir_dict['aligndir']) + file[:-6] + "_star_out")

        elif str(args_dict['program']).upper() == 'HISAT2':
            #hisat2 -- align curated reads to references
            os.system("hisat2 --quiet -x " + str(dir_dict['reference']) + "genome -U " + str(dir_dict['trimdir']) + file + " -S " + str(dir_dict['aligndir']) + file[:-6] + "_hisat2_out.sam 2> " + str(dir_dict['aligndir']) + file[:-6] + "_genome_report.txt")

        else:
            sys.exit(1)

    #Map to coding genome, do ncrna depletion
    else:
        if str(args_dict['program']).upper() == 'STAR':
            #STAR -- in silico rRNA removal, uses default parameters
            os.system("STAR --runThreadN 1 --genomeDir " + str(dir_dict['reference']) + "ncrna --readFilesIn " + str(dir_dict['trimdir']) + file + " --outReadsUnmapped Fastx --outFileNamePrefix " + str(dir_dict['aligndir']) + file[:-6] + "_norrna_")
            os.system("fastqc -q " + str(dir_dict['aligndir']) + file[:-6] + "_norrna_Unmapped.out.mate1 -o " + str(dir_dict['postncdir']))

            #STAR -- align curated reads to reference, uses default parameters
            os.system("STAR --runThreadN 1 --genomeDir " + str(dir_dict['reference']) + "genome --readFilesIn " + str(dir_dict['aligndir']) + file[:-6] + "_norrna_Unmapped.out.mate1 --outFileNamePrefix " + str(dir_dict['aligndir']) + file[:-6] + "_star_out")

        elif str(args_dict['program']).upper() == 'HISAT2':
            #hisat2 -- in silico rRNA removal
            os.system("hisat2 --quiet -x " + str(dir_dict['reference']) + "ncrna --un=" + str(dir_dict['aligndir']) + file[:-6] + "_norrna.fastq -U " + str(dir_dict['trimdir']) + file + " 2> " + str(dir_dict['aligndir']) + file[:-6] + "_ncrna_report.txt")
            os.system("fastqc -q " + str(dir_dict['aligndir']) + file[:-6] + "_norrna.fastq -o " + str(dir_dict['postncdir']))

            #hisat2 -- align curated reads to references
            os.system("hisat2 --quiet -x " + str(dir_dict['reference']) + "genome -U " + str(dir_dict['aligndir']) + file[:-6] + "_norrna.fastq -S " + str(dir_dict['aligndir']) + file[:-6]+"_hisat2_out.sam 2> "+str(dir_dict['aligndir'])+file[:-6]+"_genome_report.txt")

        else:
            sys.exit(1)

    #samtools -- sort, count, tabulate alignment output
    if str(args_dict['program']).upper() == 'STAR':
        os.system("samtools sort " + dir_dict['aligndir'] + file[:-6] + "_star_outAligned.out.sam -o " + dir_dict['aligndir'] + file[:-6] + "_star_sorted.sam")
        os.system("htseq-count " + dir_dict['aligndir'] + file[:-6] + "_star_sorted.sam " + dir_dict['reference'] + transcripts + " >> " + dir_dict['aligndir'] + file[:-6] + "_pre.csv")
        os.system("cat " + dir_dict['aligndir'] + file[:-6] + "_pre.csv | tr -s '[:blank:]' ',' > " + dir_dict['countsdir'] + file[:-6] + ".csv")

        #create a sorted bam, bed, and bigwig file for meta analyses
        os.system("samtools view -S -b " + dir_dict['aligndir'] + file[:-6] + "_star_sorted.sam > " + dir_dict['bamdir'] + file[:-6] + "_star_sorted.bam")
        os.system("samtools index " + dir_dict['bamdir'] + "/" + file[:-6] + "_star_sorted.bam")
        os.system("bedtools bamtobed -i " + dir_dict['bamdir'] + file[:-6] + "_star_sorted.bam > " + dir_dict['beddir'] + file[:-6] + "_star_sorted.bed")
        os.system("bamCoverage -b " + dir_dict['bamdir'] + file[:-6] + "_star_sorted.bam -o " + dir_dict['bwdir'] + file[:-6] + "_star_sorted.bw")

    elif str(args_dict['program']).upper() == 'HISAT2':
        os.system("samtools sort " + dir_dict['aligndir'] + file[:-6] + "_hisat2_out.sam -o " + dir_dict['aligndir'] + file[:-6] + "_hisat2_sorted.sam")
        os.system("htseq-count " + dir_dict['aligndir'] + file[:-6] + "_hisat2_sorted.sam " + dir_dict['reference'] + transcripts + " >> " + dir_dict['aligndir'] + file[:-6] + "_pre.csv")
        os.system("cat " + dir_dict['aligndir'] + file[:-6] + "_pre.csv | tr -s '[:blank:]' ',' > " + dir_dict['countsdir'] + file[:-6] + ".csv")

        #create a sorted bam, bed, and bigwig file for meta analyses
        os.system("samtools view -S -b " + dir_dict['aligndir'] + file[:-6] + "_hisat2_sorted.sam > " + dir_dict['bamdir'] + file[:-6] + "_hisat2_sorted.bam")
        os.system("samtools index " + dir_dict['bamdir'] + "/" + file[:-6] + "_hisat2_sorted.bam")
        os.system("bedtools bamtobed -i " + dir_dict['bamdir'] + file[:-6] + "_hisat2_sorted.bam > " + dir_dict['beddir'] + file[:-6] + "_hisat2_sorted.bed")
        os.system("bamCoverage -b " + dir_dict['bamdir'] + file[:-6] + "_hisat2_sorted.bam -o " + dir_dict['bwdir'] + file[:-6] + "_sorted_coverage.bw")
    else:
        sys.exit(1)


"""
MAIN
"""
def align(args_dict, dir_dict, directory, transcripts):

    #Make list of files to process alignment and package input data for multiprocessing
    align_list = file_list(directory)
    args_iter = ([file, dir_dict, args_dict, transcripts] for file in align_list)

    #Execute multiprocessed alignment of files
    with concurrent.futures.ProcessPoolExecutor(max_workers=args_dict['max_processors']) as executor:
        for file in zip(args_iter, executor.map(assemble, args_iter)):
            print(file, "has been aligned and counted.")

    if 'full_genome' in args_dict and args_dict['full_genome'] == False:
        compile_size_distribution(dir_dict, 'postncdir', 'fastqc_post_ncRNA_depletion_distribution_summary.pdf')
