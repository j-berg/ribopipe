#############
General Usage
#############

The purpose of RiboPipe is to automate the alignment, quality control, and initial analysis of ribosome profiling and other short, single-end read data. It is intended that
input data is a directory of .fastq formatted files. However, when using the intermediate submodules, such as :data:`align` or :data:`quality`, input will vary and is explicated
in the :data:`--help` menu for each submodule.

RiboPipe was created with the novice bioinformatician in mind, and the documentation written assuming no background in using Command Line or programming. As such, instructions
have been created with as much detail as possible. Additionally, `walkthrough videos <>`_ are available. In cases where RiboPipe is being used on a cloud computing cluster, or
Linux machine with Singularity installed, a singularity container image can be used to avoid the installation process. Additionally, alignment references have been created for a
variety of model organisms using current builds to help users avoid this step. Details on how to incorporate a newly created reference not already provided within the RiboPipe
infrastructure can be found `here <advanced-usage.html>`_.

RiboPipe can be run essentially from beginning to end as a pipeline, or as individual submodules. We will describe each option in detail below.

:data:`riboseq` module
^^^^^^^^^^^^^^^^^^^^^^
Pipeline for handling raw ribosome profiling sequence data. Runs quality and adaptor trimming, alignment, quality control, and formatting on a directory of raw ribosome profiling
sequence data.

The :data:`riboseq` module can be run as follows:

.. code-block:: shell

  $ ribopipe riboseq -i $DATA_PATH/raw_ingolia/ -o $DATA_PATH/out_ingolia/ -r yeast -e ingolia_2015 \
    -s a_WT_DED1_1_15deg b_WT_DED1_2_15deg c_ded1_cs_1_15deg d_ded1_cs_2_15deg \
    e_WT_DED1_1_37deg f_WT_DED1_2_37deg g_ded1_ts_1_37deg h_ded_ts_2_37deg \
    i_WT_TIF1_1_30deg j_WT_TIF1_2_30deg k_tif1_ts_1_30deg l_tif1_ts_2_30deg \
    m_WT_TIF1_1_37deg n_WT_TIF1_2_37deg o_tif1_ts_1_37deg p_tif1_ts_2_37deg \
    -p HISAT2 -a CTGTAGGCACCATCAAT --platform ILLUMINA --count_cutoff 32

And an explanation of the submodule commands can be found in the following table:

.. list-table:: riboseq commands
   :widths: 20 10 50
   :header-rows: 1

   * - Argument
     - Usage
     - Description
   * - -i INPUT, --input INPUT
     - Required
     - Specify full PATH to input directory
   * - -o OUTPUT, --output OUTPUT
     - Required
     - Specify full PATH to output directory
   * - -r <yeast>, <human>, <mouse>, --reference <yeast>, <human>, <mouse>
     - Required
     - Specify model organism used for experiments. Pipeline will align sequence data to a current reference file for the given organism
   * - -e <string>, --experiment <string>
     - Required
     - Provide experiment name to prepend to output files
   * - -s <string> [<string> ...], --samples <string> [<string> ...]
     - Required
     - Space delimited list of samples in order. If replicates are included, indicate number in sample name to delineate. For ribosome profiling, do not differentiate between footprint and RNA samples.
   * - -a <string> [<string> ...], --adaptor <string> [<string> ...]
     - Required
     - Sequence of 3' linker (only supports one 3' linker currently) (default: :data:`AACTGTAGGCACCATCAAT`). If no adaptor was used, specify :data:`None`'
   * - -m <int>, --max_processors <int>
     - Optional
     - Number of max processors pipeline can use for multiprocessing tasks (default: No limit)
   * -  -f, --footprints_only
     - Optional
     - Select this option if ONLY providing raw footprint sequence data
   * - --min_overlap <integer>
     - Optional
     - Minimum number of bases that must match on a side to combine sequences for :data:`rrna_prober`
   * - -p <HISAT2>, <STAR>, --program <HISAT2>, <STAR>
     - Optional
     - Alignment software to be used to align reads to reference (default: :data:`STAR`)
   * - --read_length_min <int>
     - Optional
     - Minimum read length threshold to keep for reads (default: :data:`11`)
   * - --read_length_max <int>
     - Optional
     - Maximum read length threshold to keep for reads (default: :data:`50`)
   * - --read_quality <int>
     - Optional
     - PHRED read quality threshold (default: :data:`28`)
   * - --platform <SANGER>, <ILLUMINA>
     - Optional
     - Sequencing platform used (default: :data:`ILLUMINA`)
   * - --full_genome
     - Optional
     - Add this option to map reads to full genome. If not given, will not count any read mapping to the first 45 nt of transcripts for ribosome profiling <riboseq>
      (it is recommended to NOT include this option as ribosome profiling data for this region is often unreliable Specify as False if running [align] with non-Ribosome
      Profiling data
   * - --count_cutoff <int>
     - Optional
     - Minimum counts threshold. Will remove any row in the final count tables if any sample does not meet this cutoff threshold

:data:`rnaseq` module
^^^^^^^^^^^^^^^^^^^^^^
Similar to the :data:`riboseq` module, but tuned for small RNAseq (i.e. anything single-end under 100 bps). It is expected that in the future, RiboPipe will be able to perform
automated paired-end and genome sequencing assembly and alignment.

:data:`trim` module
^^^^^^^^^^^^^^^^^^^^^^


:data:`align` module
^^^^^^^^^^^^^^^^^^^^^^


:data:`quality` module
^^^^^^^^^^^^^^^^^^^^^^


:data:`rrna_prober` module
^^^^^^^^^^^^^^^^^^^^^^^^^^


:data:`gene_dictionary` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


:data:`diffex` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


:data:`truncate` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


:data:`gene_dictionary` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
