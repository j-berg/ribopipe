import pandas as pd
import numpy as np
import csv
import re
pd.options.mode.chained_assignment = None

def truncate(args_dict):
    #Import gtf reference file to
    if str(args_dict['input']).endswith('.gtf') or str(args_dict['input']).endswith('.gff'):
        df = pd.read_csv(args_dict['input'], sep="\t", header=None, comment='#')
    else:
        sys.exit(1)

    #Check that every df.iloc[:,8] contains 'gene_id'
    df_id = df.iloc[:, 8].str.contains('gene_id')

    #Not all are for coding sequences -- remove
    df_id = df.iloc[:, 8].str.contains('gene_id \"Y')

    #Remove rows where "gene_id \"Y" is not present
    df_coding = df[df.iloc[:, 8].str.contains('gene_id \"Y') == True]

    #See what fell out
    df_noncoding = df[df.iloc[:, 8].str.contains('gene_id \"Y') != True]

    #Save to .gtf file (tsv)
    df_coding.to_csv(str(args_dict['input'])[:-4] + "_coding.gtf",sep="\t",header=None, index=False, quoting=csv.QUOTE_NONE)


    #Truncate first 45 nt of mRNA fasta (yeast) -- but with just the mRNA as reference, will that lose data since ribosome will need 30 nt of match to bind?
    df_copy = df_coding.copy() #Make a deep copy where original won't be mutated in next step

    #Truncate 45 nt
    naughty_genes = []

    for index, row in df_copy.iterrows():
        if "exon_number \"1\"" in row[8] and "exon_number \"11\"" not in row[8]:
            if row[6] == "+":
                if row[2] == "exon" and (row[3] + 45) <= row[4]:
                    df_copy.loc[index, 3] = row[3] + 45
                elif row[2] == "exon" and (row[3] + 45) > row[4]:
                    naughty_genes.append(index)
                else:
                    pass


            elif row[6] == "-":
                if row[2] == "exon" and (row[4] - 45) >= row[3]:
                    df_copy.loc[index, 4] = row[4] - 45
                elif row[2] == "exon" and (row[4] - 45) < row[3]:
                    naughty_genes.append(index)
                else:
                    pass

    #delete all rows with that id in row[8]
    df1 = df_copy.drop(naughty_genes)

    df1.to_csv(str(args_dict['input'])[:-4] + "_coding_truncated.gtf",sep="\t",header=None, index=False, quoting=csv.QUOTE_NONE)
