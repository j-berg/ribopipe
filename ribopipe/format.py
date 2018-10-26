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
import numpy as np
from math import floor
from .utils import file_list
from .normalize import *

"""
FUNCTIONS
"""
#Format sample names for count table
def add_long_names(samples, count_list):

    #Give table column names input sample names
    list_long_names = []

    #Check for fresh output directory
    if (len(count_list) / 2) == len(samples):
        #initialize counters
        is_fp = 0
        is_sample = 0

        for x in count_list:
            if is_fp % 2 == 0:
                list_long_names.append(samples[int(floor(is_sample))] + '_FP')
            else:
                list_long_names.append(samples[int(floor(is_sample))] + '_RNA')

            is_fp += 1
            is_sample += 0.5

    elif len(count_list) == len(samples):
        #initialize counters
        is_sample = 0

        for x in count_list:
            list_long_names.append(samples[int(floor(is_sample))])

            is_sample += 1

    else:
        print('Something went wrong. You might be re-running this program with a previously used output directory.')
        sys.exit(1)

    return list_long_names

#Save pandas dataframe to csv file
def save_table(df, countsdir, args_dict, table_type):

    #Output raw table with experiment name in file title if given
    if args_dict['experiment'] != None:
        df.to_csv(countsdir + args_dict['experiment'] + "_" + table_type + ".csv", sep=",")

    else:
        df.to_csv(countsdir + "_" + table_type + ".csv", sep=",")

"""
MAIN
"""
def format(args_dict, countsdir):

    #Retrieve list of count tables to be concatenated
    count_list = file_list(countsdir)

    #get gene list from the first, as well as length(#rows)
    with open(str(countsdir) + str(count_list[0])) as f:
        gene_names = pd.read_csv(f, header=None, usecols=[0], dtype=str)
        row_num = len(gene_names) - 5
        gene_names = gene_names[:-5]

    #populate dataframe with expression values
    df = pd.DataFrame(index=range(row_num))
    df['gene_names'] = gene_names

    for x in count_list:

        with open(countsdir + x) as f:
            reader = pd.read_csv(f, header=None, usecols=[1])
            df[x] = pd.Series(df.index)
            length = len(reader)
            df[x] = reader

    #Output raw counts table compiled, normalize data
    df = df.set_index('gene_names') #has to go after for this frame, but has to go before data population of normed frame, not sure why

    #For Ribosome Profiling data
    df.columns = add_long_names(args_dict['samples'], count_list)
    save_table(df, countsdir, args_dict, 'raw_counts_compiled')

    #Threshold counts if requested
    if args_dict['count_cutoff'] is not None:
        df = cutoff(df, int(args_dict['count_cutoff']))

    #Normalize counts to RPM values
    df_rpm = rpm_convert(df)
    save_table(df_rpm, countsdir, args_dict, 'rpm_normalized')

    if 'type' in args_dict:
        if args_dict['type'] == 'riboseq' and args_dict['footprints_only'] == False:
            #Normalize translation efficiency
            df_te = te_convert(df_rpm, args_dict['samples'], log2=True)
            save_table(df_te, countsdir, args_dict, 'translation_normalized')

    return df
