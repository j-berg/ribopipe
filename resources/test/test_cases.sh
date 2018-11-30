#!/bin/bash

#RiboPipe
#An assembly and analysis pipeline for sequencing data
#alias: ripopipe

#Copyright (C) 2018  Jordan A. Berg
#jordan <dot> berg <at> biochem <dot> utah <dot> edu

#This program is free software: you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your option) any later
#version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with
#this program.  If not, see <https://www.gnu.org/licenses/>.

USER=/Users/name

#riboseq_basic
ribopipe riboseq -i $USER/Desktop/ribo_test/test_small/input/riboseq_basic -o $USER/Desktop/ribo_test/test_small/output/riboseq_basic -a AACTGTAGGCACCATCAAT --samples 00m 05m -r yeast --experiment ribopipe_basic -p HISAT2
"""
rnaseq_parser referenced before called -> should have been renamed to riboseq for riboseq subparser FIXED

gzip: can't stat: /Users/jordanberg/Desktop/ribo_test/test_small/input/riboseq_basicYeast_00m_RNA_1000000.fastq.gz (/Users/jordanberg/Desktop/ribo_test/test_small/input/riboseq_basicYeast_00m_RNA_1000000.fastq.gz.gz): No such file or directory
  -> check_directory script not fixing input and output PATH names, problem was tried unzipping files before checking directory PATHs were formatted correctly FIXED

COMPLETE
"""

#riboseq_fp_only
ribopipe riboseq -i $USER/Desktop/ribo_test/test_small/input/riboseq_fp_only -o $USER/Desktop/ribo_test/test_small/output/riboseq_fp_only -a AACTGTAGGCACCATCAAT --samples 00m 05m -r yeast --experiment riboseq_fp_only -p HISAT2 -f True
"""
Plotting issue in  compile_size_distribution, IndexError: too many indices for array, likely due to too few samples (likes >3)
  -> Made standard row number = 2 if less than 4 samples and rounded up row number if over 4 samples FIXED

COMPLETE
"""

#riboseq_noreps_oddfiles
ribopipe riboseq -i $USER/Desktop/ribo_test/test_small/input/riboseq_noreps_oddfiles -o $USER/Desktop/ribo_test/test_small/output/riboseq_noreps_oddfiles -a AACTGTAGGCACCATCAAT --samples 00m 05m 15m -r yeast --experiment riboseq_noreps_oddfiles -p HISAT2
"""
Plotting issue in  compile_size_distribution, IndexError: too many indices for array, likely due to too few samples (likes >3), probably carry over from not applying the fix from last test elsewhere
  -> Was an issue with not rounding up and getting a nrow call for 1.25, 1.5 rows, etc, fixed up rounding up FIXED

COMPLETE
"""

#riboseq_oneset
ribopipe riboseq -i $USER/Desktop/ribo_test/test_small/input/riboseq_oneset -o $USER/Desktop/ribo_test/test_small/output/riboseq_oneset -a AACTGTAGGCACCATCAAT --samples 00m -r yeast --experiment riboseq_oneset -p HISAT2
"""
COMPLETE
"""

#riboseq_replicates_evenfiles
ribopipe riboseq -i $USER/Desktop/ribo_test/test_small/input/riboseq_replicates_evenfiles -o $USER/Desktop/ribo_test/test_small/output/riboseq_replicates_evenfiles -a AACTGTAGGCACCATCAAT --samples 00m_1 00m_2 05m_1 05m_2 -r yeast --experiment riboseq_replicates_evenfiles -p HISAT2
"""
COMPLETE
"""

#riboseq_replicates_oddfiles
ribopipe riboseq -i $USER/Desktop/ribo_test/test_small/input/riboseq_replicates_oddfiles -o $USER/Desktop/ribo_test/test_small/output/riboseq_replicates_oddfiles -a AACTGTAGGCACCATCAAT --samples 00m_1 00m_2 05m_1 -r yeast --experiment riboseq_replicates_oddfiles -p HISAT2
"""
COMPLETE
"""

#rnaseq_noreps_oddfiles
ribopipe rnaseq -i $USER/Desktop/ribo_test/test_small/input/rnaseq_noreps_oddfiles -o $USER/Desktop/ribo_test/test_small/output/rnaseq_noreps_oddfiles -a AACTGTAGGCACCATCAAT --samples 00m 05m 15m -r yeast --experiment rnaseq_noreps_oddfiles --replicates -p HISAT2
"""
Running quality when it shouldn't
  -> Changed args to store as bool values, deleted  bool forcing since was saying true if any string present FIXED

COMPLETE
"""

#rnaseq_noreps_twofiles
ribopipe rnaseq -i $USER/Desktop/ribo_test/test_small/input/rnaseq_noreps_twofiles -o $USER/Desktop/ribo_test/test_small/output/rnaseq_noreps_twofiles -a AACTGTAGGCACCATCAAT --samples 00m 05m -r yeast --experiment rnaseq_noreps_twofiles -p HISAT2
"""
COMPLETE
"""

#rnaseq_onefile
ribopipe rnaseq -i $USER/Desktop/ribo_test/test_small/input/rnaseq_onefile -o $USER/Desktop/ribo_test/test_small/output/rnaseq_onefile -a AACTGTAGGCACCATCAAT --samples 00m -r yeast --experiment rnaseq_onefile -p HISAT2
"""
COMPLETE
"""

#rnaseq_replicates_fourfiles
ribopipe rnaseq -i $USER/Desktop/ribo_test/test_small/input/rnaseq_replicates_fourfiles -o $USER/Desktop/ribo_test/test_small/output/rnaseq_replicates_fourfiles -a AACTGTAGGCACCATCAAT --samples 00m_1 00m_2 05m_1 05m_2 -r yeast --experiment rnaseq_replicates_fourfiles --replicates -p HISAT2
"""
COMPLETE
"""

#rnaseq_replicates_oddsets
ribopipe rnaseq -i $USER/Desktop/ribo_test/test_small/input/rnaseq_replicates_oddsets -o $USER/Desktop/ribo_test/test_small/output/rnaseq_replicates_oddsets -a AACTGTAGGCACCATCAAT --samples 00m_1 00m_2 05m_1 05m_2 15m_1 15m_2 -r yeast --experiment rnaseq_replicates_oddsets --replicates -p HISAT2
"""
COMPLETE
"""

#rnaseq_replicates_twofiles
ribopipe rnaseq -i $USER/Desktop/ribo_test/test_small/input/rnaseq_replicates_twofiles -o $USER/Desktop/ribo_test/test_small/output/rnaseq_replicates_twofiles -a AACTGTAGGCACCATCAAT --samples 00m_1 00m_2 -r yeast --experiment rnaseq_replicates_twofiles --replicates -p HISAT2
"""

"""
