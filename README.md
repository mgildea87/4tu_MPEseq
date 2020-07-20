# 4-tu_MPE_seq
See paper for more thorough details (DOI:). Genome index was generated with STAR. Parameters can be found in genomeParameters.txt. MPEseq.py is the pipeline for sequencing data processing and feature counting. Annotation files references in the pipeline can be found in Annotation_files/. The pipeline takes as input a list of fastq.gz files from a paired end sequencing run. Fastq files can be found at (GEO). The pipeline outputs several coalated count files. These files are used directly in the rate_modeling/ r scripts. An exception is the "adjusted" files were generated as described in the methods section of the paper using the 'Combined_branched_counts', 'Combined_branch_to3ss_counts', and 'Combined_concordant_unspliced_counts' files. 

Intron names in the annotation files used throughout are of this format "gene;intron". Introns are numbered from the 5' end of the transcript.

The .r scripts in rate_modeling/ were used to generate the data in the paper and were written to run on data from their respective cell type and modeling method. The "coupled_model" scripts output synthesis, total, 1st, and 2nd step rate estimates as well as 90% confidence intervals for each parameter estimate from the coupled model. The "total_splicing_rate_model" scripts output synthesis, and total splicing rate estimates as well as 90% confidence intervals for each parameter estimate from the total splicing rate model.

The script and supporting files for the position weight matrix splice scores are included in PWM_splicesite_scores.
