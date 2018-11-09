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
import seaborn as sns
import math
import matplotlib
matplotlib.use('agg') #remove need for -X server connect
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import linregress
from .utils import file_list
import concurrent.futures

"""
FUNCTIONS
"""
#Compile all images of a given analysis into one PDF
def compile_images(directory, plotdir, summary_name):

    #Get files for plotting
    image_list = file_list(directory)

    #Keep axes happy to avoid 'IndexError: too many indices for array' error
    if len(image_list)/2 < 2:
        plot_rows = 2
        fig_size = (15,16)
    else:
        plot_rows = math.ceil(len(image_list)/2)
        fig_size = (15,(8*(int(len(image_list)/2))))

    #Initiate plotting necessities
    fig, axes = plt.subplots(nrows=plot_rows, ncols=2, figsize=fig_size)
    plt.subplots_adjust(bottom = .3)

    #Initiate counters
    file_number = 0
    ax_y = 0
    plot_type = ''

    for file in image_list:

        with open(directory + file, 'r') as f:

            if file.endswith('metrics'):

                lab_x = 'meta_position'
                lab_y = 'normalized_coverage(all_reads)'

                df = pd.DataFrame(columns=[lab_x,lab_y])

                for line in f:

                    if 'normalized_position' in line:
                        plot_type = 'metaORF'
                        x = 0

                        for line in f:
                            #get relevant data and store in dataframe
                            data1, data2 = line.split("\t")

                            if data2.endswith('\n'):
                                data2 = data2.strip()

                            try:
                                df.loc[x] = [int(data1),float(data2)]
                            except:
                                df.loc[x] = [int(data1),0]

                            x += 1

                            if data1 == '100':
                                break

            elif file.endswith('profile.txt'):
                lab_x = 'transcript_position(nt)'
                lab_y = 'metagene_average(all_reads)'

                df = pd.DataFrame(columns=[lab_x,lab_y])

                for line in f:

                    if 'metagene_average' in line:
                        plot_type = 'periodicity'
                        x = 0

                        for line in f:
                            #get relevant data and store in dataframe
                            data1, data2, data3 = line.split("\t")
                            if data2 == 'nan':
                                data2 = 0

                            df.loc[x] = [int(data1),float(data2)]
                            x += 1

                            if data1 == '199': #this requires plastid reference to be made the same way each time
                                break

            else:
                continue

        #prepare subplots
        if file_number % 2 == 0:
            ax_x = 0
        else:
            ax_x = 1

        if file_number != 0:
            if file_number % 2 == 0:
                ax_y += 1

        #format plot name
        if summary_name == 'meta_ORF_analysis':
            name = file[8:-26]
        elif summary_name == '3_nt_periodicity_analysis':
            name = file[8:-47]
        else:
            pass

        #plot
        try:
            df.plot.line(x=lab_x, y=lab_y, title=name, ax=axes[ax_y,ax_x])
        except:
            print(name + 'did not contain any meta data, skipping...')
            pass
        file_number += 1

    #Save catenated figure
    fig.savefig(plotdir + summary_name + '_summary.pdf', dpi=600)

#Plot quality graphs (RPF vs RNA scatters plots with coefficient of determination, R^2)
def quality(df, plotdir, type):

    #Initialize data
    sample_list = list(df)
    df = df.T

    #Initialize counters
    x = 0
    y = 1
    file_number = 0
    ax_y = 0

    #Initiate plotting necessities
    if len(sample_list)/2 < 2:
        plot_rows = 2
        fig_size = (15,16)
    else:
        plot_rows = math.ceil(len(sample_list)/2)
        fig_size = (15,(7.5*(int(len(sample_list)/2))))

    #Initiate plotting necessities
    fig, axes = plt.subplots(nrows=plot_rows, ncols=2, figsize=fig_size)
    plt.subplots_adjust(bottom = .3)

    while x < len(sample_list):
        #to prevent the divide by 0 error, loss of R2 value
        df += 1e-01

        #Pass samples to lists
        sampl1 = df.loc[sample_list[x]].values.tolist()
        sampl2 = df.loc[sample_list[y]].values.tolist()

        #Scale samples on log10 if user-required <----- Maybe do this after LM?
        sampl1 = np.log10(sampl1)
        sampl2 = np.log10(sampl2)

        #Run linear regression on user input samples and collect stats
        slope, intercept, r_value, p_value, std_err = linregress(sampl1, sampl2)

        #Calculate r_sq with conserved sign for anticorrelation
        if r_value < 0:
            r_sq = -1 * (r_value**2)
        else:
            r_sq = r_value**2

        #Create dataframe of user input sample lists for graphing
        df_graph = pd.DataFrame({sample_list[x]:sampl1, sample_list[y]:sampl2})

        #prepare subplots
        if file_number % 2 == 0:
            ax_x = 0
        else:
            ax_x = 1

        if file_number != 0:
            if file_number % 2 == 0:
                ax_y += 1

        #Plot scatter of two genes of interest
        df_graph.plot.scatter(x=sample_list[x],y=sample_list[y], c='Black', s=0.5, ax=axes[ax_y,ax_x], title="r$^2$ = " + str(round(r_sq, 2)) + " (p = " + str(round(p_value,4)) + ")")

        x += 2
        y += 2
        file_number += 1

    #Save graph
    if type == 'riboseq':
        fig.savefig(plotdir + 'RPFvsRNA_correlation_summary.pdf', dpi=600, bbox_inches="tight")
    elif type == 'rnaseq':
        fig.savefig(plotdir + 'RNA_replicates_correlation_summary.pdf', dpi=600, bbox_inches="tight")
    elif type == 'custom':
        fig.savefig(plotdir + 'custom_correlation_summary.pdf', dpi=600, bbox_inches="tight")
    else:
        sys.exit(1)
#Run Picard and Plastid meta-analyses file-wise
def meta_run(args):

    #Get variables
    file, dir_dict, args_dict, transcripts_flat = args[0], args[1], args[2], args[3]

    #Perform meta ORF analysis for each file
    os.system("picard CollectRnaSeqMetrics QUIET=true REF_FLAT=" + dir_dict['reference'] + transcripts_flat + " STRAND_SPECIFICITY=NONE INPUT=" + dir_dict['bamdir'] + file + " OUTPUT=" + dir_dict['picarddir'] + file[:-4] + "_rna_metrics")

    #Perform periodicity analysis for each file
    if 'type' in args_dict and args_dict['type'] == 'riboseq':
        os.system('metagene count -q ' + dir_dict['reference'] + 'cds_start_200_rois.txt ' + dir_dict['plastdir'] + file[:-4] + '_periodicity --count_files ' + dir_dict['bamdir'] + file + ' --fiveprime --offset 14 --normalize_over 30 200 --min_counts 50 --cmap Blues --title ' + file[:-4])

"""
MAIN
"""
def meta_analysis(args_dict, dir_dict, transcripts_flat):

    #Prepare variables
    bam_files = file_list(dir_dict['bamdir'])
    args_iter = ([file, dir_dict, args_dict, transcripts_flat] for file in bam_files if file.endswith('.bam'))

    #Run each file through the meta-analysis suite/PARALLEL
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for file in zip(args_iter, executor.map(meta_run, args_iter)):
            print(file, "has been analyzed.")

    #Compile meta-analyses outputs into summary sheets
    #MetaORF
    compile_images(dir_dict['picarddir'], dir_dict['highlights'], 'meta_ORF_analysis')

    #Periodicity
    if 'type' in args_dict and args_dict['type'] == 'riboseq':
        compile_images(dir_dict['plastdir'], dir_dict['highlights'],'3_nt_periodicity_analysis')
