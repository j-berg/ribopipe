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

"""
FUNCTIONS
"""
#Install dependencies on local computer
def local():
    #INSTALL DEPENDENCIES
    print("THE FOLLOWING MAY ASK YOU FOR YOUR PASSWORD -- THIS IS IN ORDER TO INSTALL THE NECESSARY DEPENDENCIES WITH ADMIN PERMISSIONS.")
    #os.system("sudo conda install -y -c bioconda setuptools=40.4.3 fastqc=0.11.7 fastx_toolkit=0.0.14 htseq=0.9.1 picard=2.18.3 samtools=1.7 star=2.6.1b bedtools=2.27.1 deeptools=3.1.3 scipy=1.0.0 plastid=0.4.8 pandas=0.22.0 numpy=1.8.0 matplotlib=2.2.0 seaborn=0.9.0")
    os.system("sudo conda install -y -c bioconda setuptools fastqc fastx_toolkit htseq picard samtools star bedtools deeptools scipy plastid pandas numpy matplotlib seaborn")
    sys.exit()

#Install dependencies on supercomputer using module loading
def node():
    #INSTALL DEPENDENCIES, change if need to install in an environment for installation and operation
    os.system("module load hisat2 samtools python3 fastx_toolkit fastqc picard plastid star bedtools deeptools scipy plastid pandas numpy matplotlib seaborn")
