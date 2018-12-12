import pandas as pd
import numpy as np
import csv
import re
pd.options.mode.chained_assignment = None
from multiprocessing import cpu_count, Pool

"""
FUNCTIONS
"""
#Truncate 45 nt
def truncate_chunk(df_copy):
    df_copy['plus'] = df_copy[[2,3,4,6,8]].apply(lambda x:
        (x[3] + 45) if x[2] == "exon" and x[3] + 45 <= x[4] and x[6] == "+" and "exon_number \"1\"" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] + 45 > x[4] and x[6] == "+" and "exon_number \"1\"" in x[8] else x[3]),axis=1)

    df_copy['minus'] = df_copy[[2,3,4,6,8]].apply(lambda x:
        (x[4] - 45) if x[2] == "exon" and x[3] <= x[4] - 45 and x[6] == "-" and "exon_number \"1\"" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] > x[4] - 45 and x[6] == "-" and "exon_number \"1\"" in x[8] else x[4]),axis=1)

    #remove exon1s that are too short
    df_copy = df_copy[~df_copy['plus'].isin(['delete_this'])]
    df_copy = df_copy[~df_copy['minus'].isin(['delete_this'])]

    #copy new coordinates back to original columns
    df_copy[3] = df_copy['plus']
    df_copy[4] = df_copy['minus']

    #remove placeholder columns
    df_copy = df_copy.drop(columns=['plus','minus'])
    return df_copy

def parallelize(data, func):

    cores = cpu_count() #Number of CPU cores on your system
    partitions = cpu_count() #Define as many partitions as you want

    data_split = np.array_split(data, partitions)
    pool = Pool(cores)
    data = pd.concat(pool.map(func, data_split))
    pool.close()
    pool.join()
    return data

"""
MAIN
"""
def truncate(args_dict):

    #Import gtf reference file to
    print("reading in reference")
    if str(args_dict['input']).endswith('.gtf'):
        df = pd.read_csv(args_dict['input'], sep="\t", header=None, comment='#', low_memory=False)
    else:
        sys.exit(1)

    #Remove rows where "gene_id \"Y" is not present
    print("gathering coding genes")
    df_coding = df[df.iloc[:, 8].str.contains('protein_coding') == True]

    #Save to .gtf file (tsv)
    print("saving coding genes to reference file")
    df_coding.to_csv(str(args_dict['input'])[:-4] + "_coding.gtf",sep="\t",header=None, index=False, quoting=csv.QUOTE_NONE)

    #Truncate first 45 nt of mRNA fasta (yeast) -- but with just the mRNA as reference, will that lose data since ribosome will need 30 nt of match to bind?
    print("deep copying coding genes")
    df_coding_copy = df_coding.copy()

    print("multiprocessing reference chunks -- this may take a while...")
    data = parallelize(df_coding_copy, truncate_chunk)

    print("saving truncated reference")
    data.to_csv(str(args_dict['input'])[:-4] + "_truncated.gtf",sep="\t",header=None, index=False, quoting=csv.QUOTE_NONE)
