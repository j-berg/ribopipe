<p><img src="https://github.com/j-berg/ribopipe/blob/master/ribopipe_logo_v4.png" class="center" width="17%" height="17%" align="right">

<b><u>RiboPipe v0.1.4-beta</u></b>   
<i>A Flexible Sequence Assembly and Analysis Pipeline</i>  

Author: Jordan A Berg   
Affiliation: Department of Biochemistry, University of Utah, Salt Lake City, Utah, USA  

Contact: jordan \<dot\> berg \<at\> biochem \<dot\> utah \<dot\> edu </p>
<br />
Please cite the following any publications where this software was used to process or analyze data:   
```
Berg, JA, ..., Rutter, JP. (XXXX) RiboPipe: A Flexible Sequence Assembly and Analysis Pipeline. Coming soon.
```
<br />
<b><u>WHAT IS RIBOPIPE?</u></b>   
<i><a href="https://en.wikipedia.org/wiki/Ribosome_profiling">Ribosome profiling</a></i> utilizes Next Generation Sequencing to provide a detailed picture of the protein translation landscape within cells. Cells are lysed, translating ribosomes are isolated, and the ribosome protected mRNA fragments are integrated into a sequencing library. The library is then sequenced and raw data (often in the form of <i><a href="http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm">.fastq</a></i> or <i>.txt</i> files) is generated. This pipeline is flexibly designed to be able to process and perform preliminary analyses on SE (single-end) short (<= 100 bp) read raw sequence data.   

See this <a href="https://www.ncbi.nlm.nih.gov/pubmed/28579404">paper</a> for a recent discussion and detailed protocol of the technique.   


<img src="https://github.com/j-berg/ribopipe/blob/master/riboseq_overview.png" class="center">

<b>RiboPipe</b> is a ribosome profiling raw data assembly and preliminary analysis pipeline intended to ease the process of analyzing ribosome profiling data. It alleviates the pain of having to manually pass each raw read file through the appropriate quality trimming and assembly software. Additionally, it mitigates any potential stress by outputting the necessary quality checking analysis so the user can verify the quality of their run. It also offers the benefit of multiprocessing to make full use of computational resources, as well as faster assemblers to speed up this assembly process.   

Watch this <a href=""><b>video</b></a> for a walkthrough of how to use Ribopipe.

<b><u>LOCAL INSTALLATION:</u></b>   
1)  Make sure Python3, git, and wget are installed (we recommend version <a href='https://www.python.org/downloads/release/python-364/'>3.5.0</a> or higher).   
2)  Download <a href='https://www.anaconda.com/download/#macos'>Conda</a>, a package manager, for your operating system. Double click the ```.pkg``` file if on MacOS, the ```.exe``` file on Windows, or follow these <a href='https://conda.io/docs/user-guide/install/linux.html#install-linux-silent'>instructions</a> on Linux.    
3)  Execute the following lines of code in <a href="https://www.imore.com/how-use-terminal-mac-when-you-have-no-idea-where-start">Terminal</a> (on Mac, open Spotlight and type 'Terminal'):    
      ```linux
      #3.1a: to download current repository:
      git clone https://github.com/j-berg/ribopipe.git
      cd ribopipe/ribopipe/references

      #3.1b: to download specific version
      tag='v0.1.4-beta'
      wget https://github.com/j-berg/ribopipe/archive/$tag.zip
      unzip ribopipe-${tag:1}.zip
      cd ribopipe-${tag:1}/ribopipe/references

      #3.2: get reference
      model='yeast'
      program='hisat2'
      wget https://sourceforge.net/projects/ribopipe/files/${program}_references/${model}_reference_${program}.zip
      unzip ${model}_reference_${program}.zip
      rm ${model}_reference_${program}.zip
      cd ../../
      python3 setup.py install --prefix ~/.local

      #3.3: add script installation location given near the end of the setup scripting output to ~/.bashrc or ~/.bash_profile
      #add to .bashrc
      echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bashrc
      #add to .bash_profile
      echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bash_profile

      #3.4: Test by typing the following:
      ribopipe --help

      #3.5: Install conda dependencies:
      ribopipe install
      ```
    See local_install.sh in the <a href="https://github.com/j-berg/ribopipe/">resources</a> folder for interactive script

<b><u>HPC INSTALLATION:</u></b>   
1)  Make sure Python3, git, and wget are installed (we recommend version <a href='https://www.python.org/downloads/release/python-364/'>3.5.0</a> or higher).   
2)  Execute the following lines of code:   
      ```linux
      #3.1a: to download current repository:
      git clone https://github.com/j-berg/ribopipe.git
      cd ribopipe/ribopipe/references

      #3.1b: to download specific version
      tag='v0.1.4-beta'
      wget https://github.com/j-berg/ribopipe/archive/$tag.zip
      unzip ribopipe-${tag:1}.zip
      cd ribopipe-${tag:1}/ribopipe/references

      #3.2: get reference
      model='yeast'
      program='hisat2'
      wget https://sourceforge.net/projects/ribopipe/files/${program}_references/${model}_reference_${program}.zip
      unzip ${model}_reference_${program}.zip
      rm ${model}_reference_${program}.zip
      cd ../../
      module load python3
      python setup.py install --prefix ~/.local

      #3.3: add script installation location given near the end of the setup scripting output to ~/.bashrc or ~/.bash_profile
      #add to .bashrc
      echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bashrc
      #add to .bash_profile
      echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bash_profile

      #3.4: Test by typing the following:
      ribopipe --help
      ```
    See hpc_install.sh in the <a href="https://github.com/j-berg/ribopipe/">resources</a> folder for interactive script
3)  Modify hpc_run_template.sh in the <a href="https://github.com/j-berg/ribopipe/">resources</a> folder for an example script for submitting the pipeline job to the HPC and make sure dependencies listed in this script are on the HPC system, else they need to be locally installed    
4)  Run the script by executing the following:
      ```linux
      sbatch hpc_run_template.sh
      ```
    If you want the slurm output file to be sent to the SLURM directory to avoid storage space issues on your interactive node, then in the
    ```linux
    #SBATCH -o slurmjob-%j
    ```
    line, replace it with the path to your SLURM directory:
    ```linux
    #SBATCH -o /scratch/general/lustre/INPUT_USER_ID_HERE/slurmjob-%j
    ```   

<b><u>RUNNING THE PROGRAM:</u></b>   
1)  Download your raw sequence data and place in a folder -- this folder should contain all the sequence data and nothing else  
2)  Make sure files follow a pattern naming scheme. For example, if you had 3 genetic backgrounds of ribosome profiling data, the naming scheme would go as follows:  
     ```linux
     ExperimentName_BackgroundA_FP.fastq(.qz)  
     ExperimentName_BackgroundA_RNA.fastq(.qz)  
     ExperimentName_BackgroundB_FP.fastq(.qz)  
     ExperimentName_BackgroundB_RNA.fastq(.qz)  
     ExperimentName_BackgroundC_FP.fastq(.qz)  
     ExperimentName_BackgroundC_RNA.fastq(.qz)
     ```
    If the sample names are replicates, their sample number needs to be indicated  
    If you want the final count table to be in a particular order and the samples ordered that way are not alphabetically, append a letter in front of the sample name to force this ordering.  
      ```linux
      ExperimentName_a_WT_FP.fastq(.qz)  
      ExperimentName_a_WT_RNA.fastq(.qz)  
      ExperimentName_b_exType_FP.fastq(.qz)  
      ExperimentName_b_exType_RNA.fastq(.qz)  
      ```
    If you have replicates
      ```linux
      ExperimentName_a_WT_1_FP.fastq(.qz)  
      ExperimentName_a_WT_1_RNA.fastq(.qz)  
      ExperimentName_a_WT_2_FP.fastq(.qz)  
      ExperimentName_a_WT_2_RNA.fastq(.qz)
      ExperimentName_b_exType_1_FP.fastq(.qz)  
      ExperimentName_b_exType_1_RNA.fastq(.qz)  
      ExperimentName_b_exType_2_FP.fastq(.qz)  
      ExperimentName_b_exType_2_RNA.fastq(.qz)
      ```
    If you are just running RNAseq files through the pipeline, you only need the RNA samples in your input directory and specify the <i>rnaseq</i> module:
    ```linux
    ExperimentName_a_WT_1_RNA.fastq(.qz)   
    ExperimentName_a_WT_2_RNA.fastq(.qz)
    ExperimentName_b_exType_1_RNA.fastq(.qz)  
    ExperimentName_b_exType_2_RNA.fastq(.qz)

    ribopipe rnaseq -i input_directory ...
      ```
3)  Create a folder for pipeline output. This folder should be blank.  
4)  In Terminal, run the pipeline:  
      ```linux
      ribopipe riboseq -i /path/to/input/data/folder -o /path/to/output/data/folder -a AACTGTAGGCACCATCAAT --samples 00m 05m -r yeast --experiment ribopipe_basic -p HISAT2 ...   
      ```
      For other customization inputs, in Terminal type:
      ```linux
      ribopipe --help
      ```
5)  After the pipeline is finished processing, the data (in the case of the RIBOSEQ option) can be accessed along the following path tree:   
<img src="https://github.com/j-berg/ribopipe/blob/master/ribopipe_overview.png" align="center">


<b><u>INTERPRETING THE OUTPUT:</u></b>   
<i>Highlighted meta analyses will be output to the "highlights" folder in your indicated output directory.</i>   
<b>RPF vs RNA</b>: This summary plots the RPF counts vs mRNA counts for each gene in the samples. One would expect these metrics to be well-correlated between samples as translation is dependent on mRNA abundance. Super=translated genes are unusual. An r<sup>2</sup> value > 0.70 generally indicates a good library preparation.   
<b>Periodicity</b>: As ribosomes take 3 nt/1 codon steps down the transcript, a periodicity in read location is expected for a good library. The X axis in these figures indicates the start codon region of all transcripts in the organism and the Y axis is relative abundance of reads at that position.   
<b>metaORF</b>: This plot takes a meta view of all transcripts normalized to create a representative transcript (relative distance down the transcript, X axis). The Y axis indicates relative abundance of reads at that position down the representative transcript.   


***  
******  
The current version includes an empty references folder where reference builds can be stored. Model organisms can be specified for use as references within the pipeline why using '-r human', '-r mouse', etc.  
Reference folder after unzipping should be named "yeast_reference_HISAT2" or "human_reference_STAR", etc.
******  
***  
