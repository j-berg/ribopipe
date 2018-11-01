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
import pandas as pd
from math import ceil
import matplotlib.pyplot as plt
from .utils import file_list
import concurrent.futures

"""
FUNTIONS
"""
#Plot summary PDF of all size distributions of all files post-quality trimming
def compile_size_distribution(dir_dict):

    files = file_list(dir_dict['postqcdir'])
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
            os.system('unzip -q ' + dir_dict['postqcdir'] + file + ' -d ' + dir_dict['postqcdir'])

            with open(dir_dict['postqcdir'] + file[:-4] + '/fastqc_data.txt', 'r') as f:
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
    fig.savefig(dir_dict['highlights'] + 'fastqc_posttrim_distribution_summary.pdf', dpi=600)

#Trim adaptors file-wise if an adaptor sequence is given by user
def trim_adaptor(args):

    file, dir_dict, args_dict = args[0], args[1], args[2]

    os.system("fastqc -q " + str(args_dict['input']) + file + " -o " + str(dir_dict['preqcdir']))
    try:
        os.system("fastx_clipper -a " + str(args_dict['adaptor']) + " -l " + str(args_dict['read_length_min']) + " -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "pre_" + file)
        os.system("fastq_quality_filter -q " + str(args_dict['read_quality']) + " -i " + str(args_dict['output']) + "pre_" + file + " -o " + str(args_dict['output']) + "trimmed_" + file)
    except:
        #add -Q33 flag for illumina encoded quality scoring
        os.system("fastx_clipper -Q33 -a " + str(args_dict['adaptor']) + " -l " + str(args_dict['read_length_min']) + " -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "pre_" + file)
        os.system("fastq_quality_filter -Q33 -q " + str(args_dict['read_quality']) + " -i " + str(args_dict['output']) + "pre_" + file + " -o " + str(args_dict['output']) + "trimmed_" + file)
    os.system("fastqc -q " + str(args_dict['output']) + "trimmed_" + file + " -o " + str(dir_dict['postqcdir']))
    os.system("rm " + str(args_dict['output']) + "pre_" + file)

#Trim adaptors file-wise if an adaptor sequence is not given by user
def trim_adaptorless(args): #not complete

    file, dir_dict, args_dict = args[0], args[1], args[2]

    print('feature coming soon')

"""
MAIN
"""
def trim(args_dict, dir_dict, directory):

    trim_list = file_list(directory)

    args_iter = ([file, dir_dict, args_dict] for file in trim_list)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args_dict['max_processors']) as executor:

        if args_dict['adaptor'] is not None:

            for file in zip(args_iter, executor.map(trim_adaptor, args_iter)):
                print(file, "has been trimmed.")

        else:

            for file in zip(args_iter, executor.map(trim_adaptorless, args_iter)):
                print(file, "has been trimmed.")

    compile_size_distribution(dir_dict)
