<b><u>RiboPipe v0.1.4</u></b>   
<i>A Ribosome Profiling Data Handling Pipeline</i>  

Author: Jordan A Berg, Department of Biochemistry, University of Utah, Salt Lake City, Utah, USA  
Contact: jordan \<dot\> berg \<at\> biochem \<dot\> utah \<dot\> edu 

Please cite the following any publications where this software was used to process or analyze data:   
```
Berg, JA, ..., Rutter, JP. (XXXX) Manuscript name. Journal. DOI.
```

<b><u>WHAT IS RIBOPIPE?</u></b>   
<i><a href="https://en.wikipedia.org/wiki/Ribosome_profiling">Ribosome profiling</a></i> utilizes Next Generation Sequencing to provide a detailed picture of the protein translation landscape within cells. Cells are lysed, translating ribosomes are isolated, and the ribosome protected mRNA fragments are integrated into a sequencing library. The library is then sequenced and raw data (often in the form of <i><a href="http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm">.fastq</a></i> or <i>.txt</i> files) is generated. This pipeline is flexibly designed to be able to process and perform preliminary analyses on SE (single-end) short (<= 100 bp) read raw sequence data.   

See this <a href="https://www.ncbi.nlm.nih.gov/pubmed/28579404">paper</a> for a recent discussion and detailed protocol of the technique.   


<img src="https://github.com/j-berg/ribopipe/blob/master/riboseq_overview.png" class="center">

<b>RiboPipe</b> is a ribosome profiling raw data assembly and preliminary analysis pipeline intended to ease the process of analyzing ribosome profiling data. It alleviates the pain of having to manually pass each raw read file through the appropriate quality trimming and assembly software. Additionally, it mitigates any potential stress by outputting the necessary quality checking analysis so the user can verify the quality of their run. It also offers the benefit of multiprocessing to make full use of computational resources, as well as faster assemblers to speed up this assembly process.   


<b><u>MANUAL INSTALLATION:</u></b>   
Watch this <a href=""><b>video</b></a> for a walkthrough of how to use Ribopipe. 
1)  Make sure Python3 is installed (we recommend version <a href='https://www.python.org/downloads/release/python-364/'>3.6.4</a>).   
2)  Download <a href='https://www.anaconda.com/download/#macos'>Conda</a>, a package manager, for your operating system. Double click the ```.pkg``` file if on MacOS, the ```.exe``` file on Windows, or follow these <a href='https://conda.io/docs/user-guide/install/linux.html#install-linux-silent'>instructions</a> on Linux.    
3) Download <a href="https://github.com/j-berg/ribopipe/releases/tag/0.1.2">RiboPipe repository</a> to location of choice and unzip if necessary.  
4)  Download your <a href="https://sourceforge.net/projects/ribopipe/files/references/">reference</a> file to the <i>"references"</i> folder within the <b>RiboPipe</b> package. Unzip and delete the zipped file so only the unzipped folder remains. Make sure that the right alignment program folder is downloaded -- i.e. if -p HISAT2 is used, download the folder with suffix "_HISAT2"
5)  Add path to software to bash profile.  
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Open <a href="https://www.imore.com/how-use-terminal-mac-when-you-have-no-idea-where-start">Terminal</a> (on Mac, open Spotlight and type 'Terminal')  
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In Terminal, type the following:  
     ```linux
     cd ~/  
     vim ./bash_profile
     ```
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Scroll to the bottom of the file and type: 
     ```linux
     export PATH=path/to/ribopipe/directory:$PATH
     <ESCAPE>:wq
     ```     
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The '/path/to/ribopipe/directory' should be the parent ribopipe directory downloaded that contains the setup.py script. ```<ESCAPE>``` means press the escape key.      
6.1) If on a personal computer, you next need to install RiboPipe and any software dependencies (FASTX-Toolkit, STAR, etc.)   
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In Terminal, navigate to the RiboPipe parent directory like below and type: 
     ```linux
     cd ~/scripts/ribopipe-#.#.#
     python setup.py install
     ribopipe local_install
     ```
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Enter your system password if necessary (this requires to to be admin on your computer)  
6.2) If on a supercomputing node, these dependencies need to be loaded by the admin.  
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;When running the ribopipe program, include the 
     ```linux
     --cluster
     ```
     flag with each run of the pipeline.  

<b><u>CONDA INSTALLATION:</u></b>   
1) Download <a href="https://www.anaconda.com/download/#macos">Conda</a>   
2) Download RiboPipe from Conda:
     ```linux
     conda install -c bioconda ribopipe
     ```
3) Check that RiboPipe is in $PATH, add if needed    
     ```linux
     $PATH
     ```

<b><u>RUN SINGULARITY:</u></b>  
1) <a href="https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps"><b>Install Singularity</b></a>.   
If using an OS other than Linux, you can download a Vagrant Linux Virtual Box in Command Line by doing the following:   
Install brew if not already installed, then download and install Virtualbox and Vagrant Linux machine:   
     ```linux
     #install brew
     /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
     
     #install virtualbox and vagrant
     brew cask install virtualbox
     brew cask install vagrant
     brew cask install vagrant-manager
     
     #load Vagrant Linux environment (only the last two steps are required to login if logged out)
     mkdir singularity-vm
     cd singularity-vm
     vagrant init ubuntu/trusty64
     vagrant up
     vagrant ssh (if password required, type "vagrant")
     ```
    
2) Run RiboPipe singularity:   
     ```linux
     singularity exec ribopipe_latest.sif riboseq -i input_dir -o output_dir ...
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
<img src="https://github.com/j-berg/ribopipe/blob/master/ribopipe_overview.png" class="center">


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
