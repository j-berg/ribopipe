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
from math import ceil
import matplotlib
matplotlib.use('agg') #remove need for -X server connect
import matplotlib.pyplot as plt

"""
FUNCTIONS
"""
#Read in a simple csv as a pandas dataframe
def read_df(file_name):

    #Read in csv to pandas dataframe, use first row as column headers
    df_name = pd.read_csv(file_name, sep=",", index_col=0)

    return df_name

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
            if f.startswith('.') or f.endswith('.bai') or f.endswith('.html') or f.endswith('.json') or f.endswith('.png') or f.endswith('.pdf'): #ignore hidden files and other uninterested files (particular for some submodules)
                pass
            else:
                file_list.append(f)

    #Sort files in alphabetical order (helps in formatting the count tables correctly)
    file_list = sorted(file_list)

    return tuple(file_list)

def curate_reference(args_dict, __path__):

    print('Creating custom curated reference for RiboPipe...')

    if str(args_dict['program']).upper() == 'HISAT2':
        if str(args_dict['reference']).upper() == 'YEAST':
            print('Option not currently available...')
        elif str(args_dict['reference']).upper() == 'MOUSE':
            print('Option not currently available...')
        elif str(args_dict['reference']).upper() == 'HUMAN':
            print('Option not currently available...')
        else:
            pass
    elif str(args_dict['program']).upper() == 'STAR':
        if str(args_dict['reference']).upper() == 'YEAST':
            print('Option not currently available...')
        elif str(args_dict['reference']).upper() == 'MOUSE':
            print('Option not currently available...')
        elif str(args_dict['reference']).upper() == 'HUMAN':
            os.system('sh ' + str(__path__) + '/references/building/build_human_star.sh ' + str(args_dict['location']) + ' ' + str(__path__) + ' ' + str(args_dict['cores']))
        else:
            pass

#Gather reference information based on user input
def prep_reference(args_dict, dir_dict, __path__):

    #Allow for custom path to reference file
    if args_dict['custom'] == True:
        dir_dict['reference'] = str(args_dict['reference'])

    #Specify reference directory path based on given model organism name
    else:
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

#Plot summary PDF of all size distributions of all files post-quality trimming
def compile_size_distribution(dir_dict, location, title):

    files = file_list(dir_dict[location])
    zip_files = []

    for x in files:
        if x.endswith('.zip'):
            zip_files.append(x)

        #Keep axes happy to avoid 'IndexError: too many indices for array' error
    if len(zip_files)/2 < 2:
        plot_rows = 2
        fig_size = (15,16)
    else:
        plot_rows = ceil(len(zip_files)/2)
        fig_size = (15,(8*(int(len(zip_files)/2))))

    fig, axes = plt.subplots(nrows=plot_rows, ncols=2, figsize=fig_size)
    plt.subplots_adjust(bottom = .3)

    file_number = 0
    ax_y = 0

    for file in zip_files:
        if file.endswith('.zip'):
            os.system('unzip -q ' + dir_dict[location] + file + ' -d ' + dir_dict[location])
            print(dir_dict[location])

            with open(dir_dict[location] + file[:-4] + '/fastqc_data.txt', 'r') as f:
                for line in f:
                    if '#Length' in line:

                        df = pd.DataFrame(pd.DataFrame(columns=['length','counts']))
                        x = 0
                        for line in f: # now you are at the lines you want
                            if '>>END_MODULE' in line:
                                break
                            else:
                                data1, data2 = line.split("\t")
                                if data2.endswith('\n'):
                                    data2 = data2.strip()
                                try:
                                    df.loc[x] = [int(data1),float(data2)]
                                except:
                                    pass
                                x += 1

        #prepare subplots
        if file_number % 2 == 0:
            ax_x = 0
        else:
            ax_x = 1

        if file_number != 0:
            if file_number % 2 == 0:
                ax_y += 1

        df.plot.line(x='length', y='counts', title=file[8:-11], ax=axes[ax_y,ax_x])
        file_number += 1

    #Save catenated figure
    fig.savefig(dir_dict['highlights'] + title, dpi=600)
