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
import argparse
import multiprocessing
from textwrap import dedent
from .utils import check_csv
#from gooey import Gooey >>> to make a GUI, call 'ribopipe' in command line

"""
INITIALIZATION PARAMETERS
"""
#Set default values for arguments
DEFAULT_LINKER = 'AACTGTAGGCACCATCAAT' #set null linker sequence for use in ligation-free ribosome profiling
DEFAULT_SEQ_TYPE = 'sanger'
DEFAULT_PROGRAM = 'STAR'
DEFAULT_READ_MIN = 11
DEFAULT_READ_MAX = 50
DEFAULT_READ_QUALITY = 28
DEFAULT_NORM_LEN = False #for future implementation of normalization by transcript length
DEFAULT_5PRIME_TRUNCATE = False
DEFAULT_MAX_PROCESSORS = None
DEFAULT_DE_EQUATION = 'False'

descrip = """\
The RiboPipe submodules can be accessed by executing:
    "ribopipe module_name"

<riboseq> [--help]
    Pipeline for handling raw ribosome profiling sequence data
    Run quality and adaptor trimming, alignment, quality control, and output
        formatting on a directory of raw ribosome profiling sequence data

<rnaseq> [--help]
    Pipeline for handling raw single-end short (<= 100 bps) sequence data
    Run quality and adaptor trimming, alignment, quality control, and output
        formatting on a directory of raw single-end short (<= 100 bps) sequence
        data

<trim> [--help]
    Run quality and adaptor trimming on a directory of raw sequence data

<align> [--help]
    Run alignment and output formatting on a directory of sequence data
    User must ensure data is properly trimmed if desired

<quality> [--help]
    Perform quality analyses on a table <.csv> of sequence counts

<local_install> [--help]
    Install RiboPipe dependencies needed for running any of the other submodules

<rrna_prober> [--help]
    Run rrna_prober, a tool that identified over-represented sequences in
    footprint data for rRNA depletion in the Ribosome Profiling protocol
    Input is a directory of _fastqc_data.txt files

<gene_dictionary> [--help]
    Converts systematic gene names in a counts table to the common gene names
    Converts raw count table into an RPKM normalized count table
    Input is a raw count table <.csv> and a reference table <.csv> formatted as
    follows:
        Column1 data -> Systematic names
        Column2 data -> Corresponding common names
        Column3 data -> Transcript length (nt)

<diffex> [--help]
    Runs pairwise DESeq2 differential expression analysis on the raw counts
    table output in the pipeline and a sample description table.
    Please see https://github.com/j-berg/ribopipe/diffex_template.csv for a
    template sample description table.
    If samples input are not replicates, remove replicates column in
    diffex_template.csv
    If samples are not ribosome profiling, remove the type column in
    diffex_template.csv
    If other information is to be added, use the --name flag by providing a
    headless .csv table with information as follows:
    Column1 data -> Systematic names
    Column2 data -> Information to be added (common names, descriptions, etc.)
"""

"""
MAIN
"""
#Gather arguments
#@Gooey(advanced=True) >>> to make a GUI, call 'ribopipe' in command line
def get_arguments(args, __version__):
    if args is None:
        args = sys.argv[1:] #requires user input

    """
    INITIALIZE PARSER
    """
    parser = argparse.ArgumentParser(prog='RiboPipe', description=dedent(descrip), formatter_class=argparse.RawDescriptionHelpFormatter)
    #optional arguments
    parser.add_argument(
        "-v", "--version",
        help="Print installed version to stout",
        action="version",
        version="%(prog)s " + str(__version__)
        )
    parser.add_argument(
        "--cluster",
        help="Select this option if running on a high-performance cluster\ni.e. ribopipe riboseq --cluster ...",
        action='store_true',
        default=False
        )

    """
    PARENT SUBPARSER PROGRAMS
    """
    #Input files
    files_parser = argparse.ArgumentParser(description='Shared file info', add_help=False)
    files_parser.add_argument(
        "-i", "--input",
        help="Specify full PATH to input directory",
        required=True
        )
    files_parser.add_argument(
        "-o", "--output",
        help="Specify full PATH to output directory",
        required=True
        )
    files_parser.add_argument(
        "-m", "--max_processors",
        help="Number of max processors pipeline can use for multiprocessing tasks (default: No limit)",
        metavar="<int>",
        default=DEFAULT_MAX_PROCESSORS,
        required=False
        )

    #Sample information
    samples_parser = argparse.ArgumentParser(description='Shared sample info', add_help=False)
    samples_parser.add_argument(
        "-r", "--reference",
        help="Specifiy model organism used for experiments. Pipeline will align sequence data to a current reference file for the given organism",
        metavar="<yeast>, <human>, <mouse>",
        required=True
        )
    samples_parser.add_argument(
        "-e", "--experiment",
        help="Provide experiment name to prepend to output files",
        metavar="<string>",
        required=True
        )
    samples_parser.add_argument(
        "-s", "--samples",
        help="Space delimited list of samples in order.\nIf replicates are included, indicate number in sample name to delineate.\nFor ribosome profiling, do not differentiate between footprint and RNA samples.",
        metavar="<string>",
        nargs="+",
        required=True
        )
    samples_parser.add_argument(
        "-p", "--program",
        help="Alignment software to be used to align reads to reference (default: %s)" % DEFAULT_PROGRAM,
        metavar="<HISAT2>, <STAR>",
        default=DEFAULT_PROGRAM
        )

    #Read trimming parameters
    read_parser = argparse.ArgumentParser(description='Shared read info', add_help=False)
    read_parser.add_argument(
        "--read_length_min",
        help="Minimum read length threshold to keep for reads (default: %s)" % DEFAULT_READ_MIN,
        metavar="<int>",
        default=DEFAULT_READ_MIN
        )
    read_parser.add_argument(
        "--read_length_max",
        help="Maximum read length threshold to keep for reads (default: %s)" % DEFAULT_READ_MAX,
        metavar="<int>",
        default=DEFAULT_READ_MAX
        )
    read_parser.add_argument(
        "--read_quality",
        help="PHRED read quality threshold (default: %s)" % DEFAULT_READ_QUALITY,
        metavar="<int>",
        default=DEFAULT_READ_QUALITY
        )

    #Count information parameters
    count_parser = argparse.ArgumentParser(description='Shared count info', add_help=False)
    count_parser.add_argument(
        "--full_transcripts",
        help="Select this option to keep first 45 nt of all transcript references (it is recommended to not include this option as ribosome profiling data for this region is often unreliable\nSpecify as False if running [align] with non-Ribosome Profiling data",
        action='store_true'
        )
    count_parser.add_argument(
        "--count_cutoff",
        help="Minimum counts threshold. Will remove any row in the final count tables if any sample does not meet this cutoff threshold",
        metavar="<int>"
        )

    #Adaptor information parameters
    adaptor_parser = argparse.ArgumentParser(description='Shared adaptors info', add_help=False)
    adaptor_parser.add_argument(
        "-a", "--adaptor",
        help="Sequence of 3' linker (only supports one 3' linker currently) (default: %s)" % DEFAULT_LINKER,
        metavar="<string>",
        default=DEFAULT_LINKER,
        required=True
        )

    """
    MODULE SUBPARSER PROGRAMS
    """
    subparser = parser.add_subparsers(title='Sub-modules', description='Choose one of the following:', dest='cmd')

    #RIBOSEQ subparser program
    riboseq_parser = subparser.add_parser('riboseq', description='Ribosome Profiling Pipeline', parents=[files_parser, samples_parser, read_parser, adaptor_parser, count_parser])
    riboseq_parser.add_argument(
        "-f", "--footprints_only",
        help="Select this option if ONLY providing raw footprint sequence data",
        action='store_true',
        default=False
        )

    #RNASEQ subparser program
    rnaseq_parser = subparser.add_parser('rnaseq', description='RNAseq Pipeline', parents=[files_parser, samples_parser, read_parser, adaptor_parser])
    rnaseq_parser.add_argument(
        "--count_cutoff",
        help="Minimum counts threshold. Will remove any row in the final count tables if any sample does not meet this cutoff threshold",
        metavar="<int>"
        )
    rnaseq_parser.add_argument(
        "--replicates",
        help="Select this option if samples are replicates (do not use if 3+ replicates). Make sure when passing the argument to list samples these are sorted so replicates are next to each other in this list",
        action='store_true',
        default=False
        )

    #TRIM subparser program
    trim_parser = subparser.add_parser('trim', description='Quality trimming submodule', parents=[files_parser, read_parser, adaptor_parser])

    #ALIGN subparser program
    align_parser = subparser.add_parser('align', description='Alignment submodule', parents=[files_parser, samples_parser, count_parser])
    align_parser.add_argument(
        "-t", "--type",
        help="Sequencing type -- ribosome profiling <riboseq> or single-end short read sequence data <rnaseq>",
        required=True,
        metavar="<riboseq>, <rnaseq>"
        )

    #ALIGN subparser program
    quality_parser = subparser.add_parser('quality', description='Quality control submodule')
    #required arguments
    quality_parser.add_argument(
        "-i", "--input",
        help="Input table (must be .csv file) of raw counts",
        required=True,
        type=lambda file:check_csv(file)
        )
    quality_parser.add_argument(
        "-o", "--output",
        help="Specify full PATH to output directory",
        required=True
        )
    quality_parser.add_argument(
        "-t", "--type",
        help="Sequencing type -- ribosome profiling <riboseq> or single-end short read sequence data <rnaseq>",
        required=True,
        metavar="<riboseq>, <rnaseq>"
        )
    quality_parser.add_argument(
        "--replicates",
        help="Select this option if samples are replicates (do not use if 3+ replicates). Make sure when passing the argument to list samples these are sorted so replicates are next to each other in this list",
        action='store_true',
        default=False
        )

    #LOCAL INSTALL subparser program
    local_install_parser = subparser.add_parser('local_install', description='Local install submodule. Install or update dependencies required for RiboPipe')

    #RRNA PROBER subparser program
    probe_parser = subparser.add_parser('rrna_prober', description='rRNA prober submodule')
    #required arguments
    probe_parser.add_argument(
        "-i", "--input",
        help="Space delimited list of zipped files",
        metavar="<string>",
        nargs="+",
        required=True
        )

    #GENE/LENGTH DICTIONARY subparser program
    dict_parser = subparser.add_parser('gene_dictionary', description='Gene name/RPKM conversion submodule')
    #required arguments
    dict_parser.add_argument(
        "-i", "--input",
        help="Input file (must be .csv file)",
        required=True,
        type=check_csv
        )
    dict_parser.add_argument(
        "-o", "--output",
        help="Output file prefix",
        required=True
        )
    dict_parser.add_argument(
        "-r", "--reference",
        help="Input dictionary (must be .csv file, no headers, see ribopipe --help for more information)",
        required=True,
        type=check_csv
        )
    dict_parser.add_argument(
        "-c", "--conversion",
        help="Type of conversion to run",
        required=True,
        metavar="<common>, <common_rpkm>, <rpkm>"
        )
    dict_parser.add_argument(
        "--te_convert",
        help="Provide this flag if you want your count table to be normalized for translation efficiency. Only valid with --conversion <system2commonName_rpkm>, <system2rpkm>",
        action='store_true'
        )
    dict_parser.add_argument(
        "-s","--samples",
        help="Provide this flag if you want your count table to be normalized for translation efficiency",
        metavar="<string>",
        nargs="+",
        required=False
        )

    #DESeq2 module
    de_parser = subparser.add_parser('diffex', description='Run DESeq2 (differential expression analysis) with raw count table and sample_description table')
    #required arguments
    de_parser.add_argument(
        "-i", "--input",
        help="Input file -- raw count table (must be .csv file)",
        required=True,
        type=lambda file:check_csv(file)
        )
    de_parser.add_argument(
        "-o", "--output",
        help="Output file prefix",
        required=True
        )
    de_parser.add_argument(
        "-d", "--descriptions",
        help="Add common names, descriptions, or other information to DESeq output table (must be a .csv file, see ribopipe --help for more information)",
        metavar="<str>",
        type=check_csv
        )
    de_parser.add_argument(
        "--replicates",
        help="Select this option if samples are replicates",
        action='store_true'
        )
    de_parser.add_argument(
        "--type",
        help="Select this option if samples are replicates",
        required=True,
        metavar="<riboseq>, <rnaseq>"
        )
    de_parser.add_argument(
        "--custom",
        help="Use this to specify a custom DESeq2 design equation (default: %s)" % DEFAULT_DE_EQUATION,
        required=False,
        metavar="<str>",
        default=DEFAULT_DE_EQUATION
        )

    """
    COLLECT PARSED ARGUMENTS AND PREPARE FOR DOWNSTREAM USE
    """
    #Print help if no arguments/submodules specified
    if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

    #Parse arguments into NameSpace
    args = parser.parse_args(args)

    #Collect subargs and package
    args_dict = vars(args)

    #Check max_processor input
    if 'max_processors' in args_dict and args_dict['max_processors'] != None:
        args_dict['max_processors'] = int(args_dict['max_processors'])
        if multiprocessing.cpu_count() < args_dict['max_processors']:
            print("ERROR: Cannot specify more cores than are available.\nSpecified " + str(args_dict['max_processors']) + " cores, only " + str(multiprocessing.cpu_count()) + " available.")
            sys.exit(1)

    #Check samples given if doing a TE convert with the genename.rpkm_convert
    if 'te_convert' in args_dict and args_dict['te_convert'] == True and args_dict['samples'] is None:
        print('ERROR: Must specify samples with te_convert flag.')
        sys.exit(1)

    if 'te_convert' in args_dict and args_dict['te_convert'] == True and args_dict['conversion'] == 'common':
        print('ERROR: te_convert not valid with --conversion system2commonName.')
        sys.exit(1)

    return args, args_dict
