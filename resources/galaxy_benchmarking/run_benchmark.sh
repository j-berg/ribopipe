#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/general/lustre/u0690617/slurmjob-%j
#SBATCH --partition=kingspeak

#Set up whatever package we need to run with
module purge

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

#import ripopipe
mkdir $SCRDIR/benchmark
mkdir $SCRDIR/output
cp /uufs/chpc.utah.edu/common/home/$USER/benchmark/* $SCRDIR/benchmark/.


#catenate raw data if needed
cd $SCRDIR/.

#run script --- add directory of raw data
#Using Qiagen Illumina 3' adaptor, but shouldn't need adaptor trimming anyways, just a placeholder for now
ribopipe riboseq -i $SCRDIR/benchmark/ -o $SCRDIR/output/ -r yeast -e benchmark_test -s cccp_00m cccp_05m -p STAR -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --platform ILLUMINA
