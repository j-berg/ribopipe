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
import matplotlib
matplotlib.use('agg') #remove need for -X server connect
import matplotlib.pyplot as plt
from .utils import file_list, compile_size_distribution
import concurrent.futures

"""
FUNTIONS
"""
#Trim adaptors file-wise if an adaptor sequence is given by user
def trim_adaptor(args):

    file, dir_dict, args_dict = args[0], args[1], args[2]

    #Run initial quality control
    os.system("fastqc -q " + str(args_dict['input']) + file + " -o " + str(dir_dict['preqcdir']))

    #Trim files
    if args_dict['platform'].upper() == 'SANGER':
        if len(args_dict['adaptor']) == 1 and args_dict['adaptor'][0].upper() == "NONE":
            os.system("fastq_quality_filter -q " + str(args_dict['read_quality']) + " -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "trimmed_" + file)
        else:
            counter = 1
            for x in args_dict['adaptor']:
                os.system("fastx_clipper -a " + str(args_dict['adaptor']) + " -l " + str(args_dict['read_length_min']) + " -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "pre_" + file)
                if counter != len(args_dict['adaptor']):
                    os.system('mv ' + str(args_dict['output']) + "pre_" + file + ' ' + str(args_dict['input']) + file)
                    counter += 1
                else:
                    os.system("fastq_quality_filter -q " + str(args_dict['read_quality']) + " -i " + str(args_dict['output']) + "pre_" + file + " -o " + str(args_dict['output']) + "trimmed_" + file)

    elif args_dict['platform'].upper() == 'ILLUMINA':
        #add -Q33 flag for illumina encoded quality scoring
        if len(args_dict['adaptor']) == 1 and args_dict['adaptor'][0].upper() == "NONE":
            os.system("fastq_quality_filter -Q33 -q " + str(args_dict['read_quality']) + " -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "trimmed_" + file)
        else:
            counter = 1

            for x in args_dict['adaptor']:
                os.system("fastx_clipper -Q33 -a " + str(args_dict['adaptor']) + " -l " + str(args_dict['read_length_min']) + " -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "pre_" + file)
                if counter != len(args_dict['adaptor']):
                    os.system('mv ' + str(args_dict['output']) + "pre_" + file + ' ' + str(args_dict['input']) + file)
                    counter += 1
                else:
                    os.system("fastq_quality_filter -Q33 -q " + str(args_dict['read_quality']) + " -i " + str(args_dict['output']) + "pre_" + file + " -o " + str(args_dict['output']) + "trimmed_" + file)
    else:
        pass

    #Run post-trim quality control
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

    compile_size_distribution(dir_dict, 'postqcdir', 'fastqc_posttrim_distribution_summary.pdf')
