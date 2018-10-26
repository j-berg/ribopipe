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
import numpy as np
import pandas as pd
import sys
import csv

"""
FUNCTIONS
"""
# RPM normalization
def rpm_convert(df):

    #Normalize dataframe by sample RPM metric
    df_rpm = df / (df.sum() / 1e6)

    return df_rpm

# Translation efficiency normalization
def te_convert(df, samples, log2=True):

    #Scale all DataFrame values for division where x = x + 1e-7
    df += 1e-7

    #Initialize counters division scheme for fp/rna
    y = 0
    z = 1

    #Perform translation efficiency calculations
    for x in samples:
        df[x] = df[df.columns[y]]/df[df.columns[z]]
        #Move counters for next set of samples
        y = y + 2
        z = z + 2

    #Get TE_normalized columns
    df_te = df[samples]

    #Perform log2 scaling of data
    if log2 == True:
        df_log = np.log2(df_te)

    return df_log

# Count cut-off trimmming
def cutoff(df, min_counts):

    #Trim dataframe of rows with any value less than the min_counts value
    df = df.T[df.T.columns[df.T.min() > min_counts]]
    df = df.T

    return df

def name_convert(df, dict):
    df['gene_names'] = df['gene_names'].map(dict).fillna(df['gene_names']) #Convert names to common if available
    df = df.set_index('gene_names', drop=True) #Set gene_names column to index and drop column

    return df

def rpkm_convert(df, dict):
    df['length'] = df['gene_names'].map(dict) #map lengths to system names
    df = df.dropna() #drop if no length available
    df = df.set_index('gene_names', drop=True) #set system names as index and drop column
    df_numeric = df.apply(pd.to_numeric) #Make sample columns numeric for math
    df_numeric = df_numeric.div(df_numeric.length, axis=0) #Divide all columns by length
    del df_numeric['length'] #drop length column

    df_rpkm = df_numeric / (df_numeric.sum() / 1e6) #get PM metric for RPKs
    df_rpkm['gene_names'] =  df_rpkm.index #set index

    return df_rpkm

def de_convert(args_dict):
    system2common = {}

    #make dictionaries
    with open (args_dict['reference'], mode='r') as infile:
        dict_raw = csv.reader(infile)

        for col in dict_raw:
            system2common[col[0]] = col[1]

def run_dictionary(args_dict):

    #populate dictionaries
    system2common = {}
    system2length = {}

    #make dictionaries
    with open (args_dict['reference'], mode='r') as infile:
        dict_raw = csv.reader(infile)

        for col in dict_raw:
            system2common[col[0]] = col[1]
            system2length[col[0]] = col[2]

    #process file of interest
    with open (args_dict['input'], mode='r') as infile:
        df = pd.read_csv(infile)

    #Convert table system names to common if available
    if args_dict['conversion'] == 'common':

        df = name_convert(df, system2common)
        df.to_csv(args_dict['output'] + '_system2common.csv',sep=",")

    #Convert table system names to common if available
    elif args_dict['conversion'] == 'common_rpkm':

        df_rpkm = rpkm_convert(df, system2length)
        del df_rpkm.index.name #delete index
        df = name_convert(df_rpkm, system2common)
        if 'te_convert' in args_dict and args_dict['te_convert'] == True:
            df = te_convert(df, args_dict['samples'])
            df.to_csv(args_dict['output'] + '_system2common_rpkm_te.csv',sep=",")
        else:
            df.to_csv(args_dict['output'] + '_system2common_rpkm.csv',sep=",")

    elif args_dict['conversion'] == 'rpkm':

        df = rpkm_convert(df, system2length)
        del df['gene_names'] #delete index
        if 'te_convert' in args_dict and args_dict['te_convert'] == True:
            df = te_convert(df, args_dict['samples'])
            df.to_csv(args_dict['output'] + '_system2rpkm_te.csv',sep=",")
        else:
            df.to_csv(args_dict['output'] + '_system2rpkm.csv',sep=",")

    else:
        sys.exit(1)
