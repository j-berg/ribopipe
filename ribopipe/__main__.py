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
from .__init__ import __version__
from .arguments import get_arguments
from .utils import *
from .make_directories import *
from .rrna_prober import rrna_prober
from .trim import trim
from .align import align
from .format import format
from .normalize import *
from .quality import quality, meta_analysis
from .messages import *

"""
INITIALIZATION PARAMETERS
"""
#Retrieve path for scripts used in this pipeline
__path__, ribopipe_main = os.path.split(__file__)

"""
MAIN
"""
def main(args=None):

    """
    PRINT LICENSE INFORMATION
    """
    msg_license()

    """
    INITIALIZE ARGUMENTS
    See argparsing.py script for defaults and other information
    """
    args, args_dict = get_arguments(args, __version__)

    """
    RUN SUBMODULES
    """
    #Install dependencies for supercomputing use
    if args.cluster == True:
        print("Run bash script from resources folder before running ribopipe")
        sys.exit(1)

    #Run sequencing pipelines
    if args.cmd == 'riboseq' or args.cmd == 'rnaseq':

        #Initialize directories for output
        args_dict['input'], args_dict['output'], dir_dict = make_directories(args_dict['input'], args_dict['output'])
        dir_dict['reference'] = prep_reference(args_dict, dir_dict, __path__)

        #CHECK ARGUMENTS FOR EXCEPTIONS AND FORMATTING
        check_files(args_dict)
        check_length(args_dict)

        #Run TRIM
        msg_trim()
        trim(args_dict, dir_dict, args_dict['input'])
        os.system("mv " + args_dict['output'] + "trimmed_* " + dir_dict['trimdir'])

        #Run rRNA prober
        if args.cmd == 'riboseq':

            probe_list = []
            zipfiles = file_list(dir_dict['postqcdir'])

            for x in zipfiles:
                if x.endswith('.zip') == True:
                    if 'footprint' in x.lower() or 'fp' in x.lower():
                        probe_list.append(dir_dict['postqcdir'] + x)

            probe_out = rrna_prober(probe_list, args_dict['min_overlap']) #use inputDir to get FASTQC files and output to outputDir/analysis
            with open(dir_dict['highlights'] + 'rrna_prober_output.txt', "w") as text_file:
                print(probe_out, file=text_file)

        #Run ASSEMBLE
        msg_align()
        #Prep transcripts reference to pull from
        if args.cmd == 'riboseq':
            transcripts = 'transcripts_45.gtf'
            transcripts_flat = 'transcripts_refFlat_45.txt'
        elif args.cmd == 'rnaseq':
            transcripts = 'transcripts.gtf'
            transcripts_flat = 'transcripts_refFlat.txt'
        else:
            sys.exit(1)
        align(args_dict, dir_dict, dir_dict['trimdir'], transcripts)

        #Run CATENATE COUNTS
        msg_count()
        if args.cmd == 'riboseq':
            args_dict['type'] = 'riboseq'
        df = format(args_dict, dir_dict['countsdir'])

        #Run QC ANALYSIS
        msg_checking()
        if (args.cmd == 'riboseq' and 'footprints_only' in args_dict and args_dict['footprints_only'] == False) or (args.cmd == 'rnaseq' and 'replicates' in args_dict and args_dict['replicates'] == True):
            quality(df, dir_dict['highlights'], args.cmd)

        meta_analysis(args_dict, dir_dict, transcripts_flat)

        #Output metrics to csv in outputDir

        #Rezip raw files
        msg_cleaning()
        os.system("gzip " + args_dict['input'] + "*.fastq")
        msg_finish()

    #Run trimming module
    elif args.cmd == 'trim':

        print('Use this module at your own risk, currently untested')

        #CHECK ARGUMENTS FOR EXCEPTIONS AND FORMATTING
        check_files(args_dict)
        check_length(args_dict)
        args_dict['input'] = check_directory(args_dict['input'])

        #Make output directories
        args_dict['input'], args_dict['output'], dir_dict = create_trim_directories(args_dict['input'], args_dict['output'])

        #Run TRIM
        trim(args_dict, dir_dict, args_dict['input'])
        os.system("mv " + args_dict['output'] + "trimmed_* " + dir_dict['trimdir'])

    #Run alignment module
    elif args.cmd == 'align':

        print('Use this module at your own risk, currently untested')

        #CHECK ARGUMENTS FOR EXCEPTIONS AND FORMATTING
        check_files(args_dict)
        args_dict['input'] = check_directory(args_dict['input'])

        #Make output directories
        args_dict['input'], args_dict['output'], dir_dict = create_assemble_directories(args_dict['input'], args_dict['output'])

        #Run ASSEMBLE
        align(args_dict, dir_dict, args_dict['input'])

        #Run CATENATE COUNTS
        msg_count()
        df = format(args_dict, dir_dict['countsdir'])

    #Run quality module
    elif args.cmd == 'quality':

        print('Use this module at your own risk, currently untested')

        #Read in csv to pandas DataFrame
        df = read_df(args['input'])

        #Run quality analyses
        meta_analysis(args_dict, dir_dict)

        if bool(args['replicates']) == True:
            quality(df, dir_dict['highlights'], 'custom')

    #Run rRNA prober module
    elif args.cmd == 'rrna_prober':

        check_files(args_dict)

        #CHECK ARGUMENTS FOR EXCEPTIONS AND FORMATTING
        probe_list = []
        for x in args_dict['input']:
            if x.endswith('.zip') == True:
                if 'footprint' in x.lower() or 'fp' in x.lower():
                    probe_list.append(args_dict['input'] + x)

        #Run rrna_prober, output to outputDir
        probe_out = rrna_prober(probe_list, args_dict['min_overlap']) #use inputDir to get FASTQC files and output to outputDir/analysis
        with open(args_dict['output']+ 'rrna_prober_output.txt', "w") as text_file:
            print(probe_out, file=text_file)

    #Run gene name/RPKM conversion tool
    elif args.cmd == 'gene_dictionary':

        run_dictionary(args_dict)
        sys.exit(1)

    #Run gene name/RPKM conversion tool
    elif args.cmd == 'diffex':

        print('Use this module at your own risk, currently untested')

        os.system('rscript '+str(__path__)+'/run_deseq.r '+str(args_dict['input'])+' '+str(args_dict['output'])+'_deseq2_table.csv '+str(args_dict['descriptions'])+' '+str(args_dict['replicates'])+' '+str(args_dict['type'])+' '+str(args_dict['custom']))

        if args_dict['name'] is not None:
            #add info to csv output by de
            print('Feature coming soon')

        sys.exit(1)

    #Run local install of dependencies
    elif args.cmd == 'local_install':

        #Run bash script
        sys.exit(1)

    #Display submodule error information
    else:
        print('Invalid selection, must specify RiboPipe module to run')
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":

    sys.exit(main() or 0)
