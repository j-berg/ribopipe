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
import os, sys
from .utils import check_directory

"""
FUNCTIONS
"""
#Create directories for 'trim' submodule
def create_trim_directories(inputDir, outputDir):

    inputDir = check_directory(inputDir)
    outputDir = check_directory(outputDir)

    #create quality control output directories
    os.system('mkdir ' + outputDir + 'quality_control')
    os.system('mkdir ' + outputDir + 'quality_control/pre-trim_fastqc/')
    os.system('mkdir ' + outputDir + 'quality_control/post-trim_fastqc/')

    #create dictionary with output directory paths
    directories = {
        'preqcdir': str(outputDir + 'quality_control/pre-trim_fastqc/'),
        'postqcdir': str(outputDir + 'quality_control/post-trim_fastqc/'),
    }

    return inputDir, outputDir, directories

#Create directories for 'align' submodule
def create_assemble_directories(inputDir, outputDir):

    inputDir = check_directory(inputDir)
    outputDir = check_directory(outputDir)

    #create assembly output directories
    os.system('mkdir ' + outputDir + 'assembly')
    os.system("""mkdir """ + outputDir + """assembly/alignment_output |
                mkdir """ + outputDir + """assembly/counts |
                mkdir """ + outputDir + """assembly/post-ncRNA_depletion |
                mkdir """ + outputDir + """assembly/bam_files |
                mkdir """ + outputDir + """assembly/bed_files |
                mkdir """ + outputDir + """assembly/bigwig_files"""
    )

    #create dictionary with output directory paths
    directories = {
        'inputDir': inputDir,
        'aligndir': str(outputDir + 'assembly/alignment_output/'),
        'countsdir': str(outputDir + 'assembly/counts/'),
        'postncdir': str(outputDir + 'assembly/post-ncRNA_depletion/'),
        'bamdir': str(outputDir + 'assembly/bam_files/'),
        'beddir': str(outputDir + 'assembly/bed_files/'),
        'bwdir': str(outputDir + 'assembly/bigwig_files/')
    }

    return directories

#Create directories for full pipeline ('riboseq' or 'rnaseq'
def create_directories(inputDir, outputDir):

    #create output directories
    os.system('mkdir ' + outputDir + 'assembly')
    os.system('mkdir ' + outputDir + 'analysis')
    os.system('mkdir ' + outputDir + 'highlights')
    os.system("""mkdir """ + outputDir + """assembly/trimmed_output |
                mkdir """ + outputDir + """assembly/alignment_output |
                mkdir """ + outputDir + """assembly/counts |
                mkdir """ + outputDir + """assembly/bam_files |
                mkdir """ + outputDir + """assembly/bed_files |
                mkdir """ + outputDir + """assembly/bigwig_files |
                mkdir """ + outputDir + """analysis/periodicity |
                mkdir """ + outputDir + """analysis/quality_control |
                mkdir """ + outputDir + """analysis/metaORF |
                mkdir """ + outputDir + """highlights/_metrics"""
    )
    os.system("""mkdir """ + outputDir + """analysis/quality_control/pre-trim_fastqc |
                mkdir """ + outputDir + """analysis/quality_control/post-trim_fastqc |
                mkdir """ + outputDir + """analysis/quality_control/post-ncRNA_depletion"""
    )

    #create dictionary with output directory paths
    directories = {
        'inputDir': inputDir,
        'trimdir': str(outputDir + 'assembly/trimmed_output/'),
        'aligndir': str(outputDir + 'assembly/alignment_output/'),
        'countsdir': str(outputDir + 'assembly/counts/'),
        'bamdir': str(outputDir + 'assembly/bam_files/'),
        'beddir': str(outputDir + 'assembly/bed_files/'),
        'bwdir': str(outputDir + 'assembly/bigwig_files/'),
        'plastdir': str(outputDir + 'analysis/periodicity/'),
        'qualdir': str(outputDir + 'analysis/quality_control/'),
        'preqcdir': str(outputDir + 'analysis/quality_control/pre-trim_fastqc/'),
        'postqcdir': str(outputDir + 'analysis/quality_control/post-trim_fastqc/'),
        'postncdir': str(outputDir + 'analysis/quality_control/post-ncRNA_depletion/'),
        'picarddir': str(outputDir + 'analysis/metaORF/'),
        'highlights': str(outputDir + 'highlights/'),
        'metrics': str(outputDir + 'highlights/_metrics/')
    }

    return directories

"""
MAIN
"""
def make_directories(inputDir, outputDir):

    inputDir = check_directory(inputDir)
    outputDir = check_directory(outputDir)
    directories = create_directories(inputDir, outputDir)

    return inputDir, outputDir, directories
