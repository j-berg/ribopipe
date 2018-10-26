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
import pandas as pd
import argparse

"""
FUNCTIONS
"""
#Read in a simple csv as a pandas dataframe
def read_df(file_name):

    #Read in csv to pandas dataframe, use first row as column headers
    df_name = pd.read_csv(file_name, sep=",", index_col=0)

    return df_name

#Check transcript file formatting based on user input
def check_transcript(reference_status):

    #Initialize variables
    file_ref = ''
    file_flat = ''

    #Check if using normal or truncated transcripts and set file variables
    if reference_status == True:
        file_ref = 'transcripts.gtf'
        file_flat = 'transcripts_refFlat.txt'
    else:
        file_ref = 'transcripts_50.gtf'
        file_flat = 'transcripts_refFlat_50.txt'

    return file_ref, file_flat

#Check directory formatting
def check_directory(directory):

    #Check input directory name is formatted correctly and fix if necessary
    if directory.endswith('/'):
        pass
    else:
        directory += '/'

    return directory

#Make a list of the files in a given directory
def file_list(directory):

    #Initialize blank file list to fill
    file_list = []

    #Walk through raw data files within given directory
    for subdir, dirs, files in os.walk(directory):
        for f in files:
            if f.startswith('.') or f.endswith('.bai') or f.endswith('.png') or f.endswith('.pdf'): #ignore hidden files and other uninterested files (particular for some submodules)
                pass
            else:
                file_list.append(f)

    #Sort files in alphabetical order (helps in formatting the count tables correctly)
    file_list = sorted(file_list)

    return tuple(file_list)

#Gather reference information based on user input
def prep_reference(args_dict, dir_dict, __path__):

    #Specify reference directory path based on given model organism name
    if str(args_dict['program']).upper() == 'HISAT2':
        if str(args_dict['reference']).upper() == 'YEAST':
            dir_dict['reference'] = str(__path__) + '/references/yeast_reference_hisat2/'
        elif str(args_dict['reference']).upper() == 'MOUSE':
            dir_dict['reference'] = str(__path__) + '/references/mouse_reference_hisat2/'
        elif str(args_dict['reference']).upper() == 'HUMAN':
            dir_dict['reference'] = str(__path__) + '/references/human_reference_hisat2/'
        else:
            pass
    elif str(args_dict['program']).upper() == 'STAR':
        if str(args_dict['reference']).upper() == 'YEAST':
            dir_dict['reference'] = str(__path__) + '/references/yeast_reference_star/'
        elif str(args_dict['reference']).upper() == 'MOUSE':
            dir_dict['reference'] = str(__path__) + '/references/mouse_reference_star/'
        elif str(args_dict['reference']).upper() == 'HUMAN':
            dir_dict['reference'] = str(__path__) + '/references/human_reference_star/'
        else:
            pass
    else:
        pass

    return dir_dict['reference']

#Check that read lengths are consistent
def check_length(args_dict):

    #Check read length min is not greater than max
    if int(args_dict['read_length_min']) > int(args_dict['read_length_max']):
        print("Invalid read lengths -- min read length cannot be greater than max read length")
        sys.exit(1)

#Check if file is comma separated
def check_csv(file):

    #Make sure file is a csv (for use with converting count tables to gene names and RPKM values)
    if str(file).endswith('.csv'):
        pass
    else:
        raise argparse.ArgumentTypeError("Input file type incorrect")
        sys.exit(1)

    return str(file)

#Check files are in proper format and modify if necessary
def check_files(args_dict):


    for subdir, dirs, files in os.walk(args_dict['input']): #walk through raw data files within main data directory

        for x in files:

            #Check that raw files are unzipped
            if x.startswith('.'): #ignore hidden files
                continue
            if x[-3:] == '.gz':
                os.system("gzip -d " + args_dict['input'] + x)
            elif x[-4:] == '.zip':
                os.system("unzip -q " + args_dict['input'] + x)
            else:
                pass

            #Check file type is correct for trimming/alignment submodules
            if 'txt' in x or 'fastq' in x:
                pass
            else:
                print("Incorrect file type: " + x + '\n.fastq or .txt file required.')
                sys.exit(1)
