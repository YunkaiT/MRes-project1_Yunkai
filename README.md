# MRes-project1_Yunkai

## Pipeline
The comprehensive pipeline developed for this study is detailed in the file 'Total Pipeline.sh'. This pipeline, based on bash scripting, requires several additional files, including: WGS .bam files, 'QDNAseq_from_bam_chrX.R','cytoband_adapted_hg19.csv','312shallowHRD_hg19_1.13_QDNAseq_chrX.R', and a custom WGS file name list generated from your one WGS dataset.

## Regarding the Files
It should be noted that 'QDNAseq_from_bam_chrX.R','cytoband_adapted_hg19.csv','312shallowHRD_hg19_1.13_QDNAseq_chrX.R' are adapted from the original shallowHRD GitHub page, with minor modifications to the output format.

Other files, which are used for data organization, transformation, and processing, are not uploaded due to their easy replicability using R packages like tydiverse.

## Regarding the Data
The input data format for this project is WGS data (.bam file). As all whole-genome sequencing (WGS) data utilized in this project are confidential-including both the raw WGS data and the processed data for subsequent analysis-none of the data used in this project is uploaded tp this public GitHub repository. However, users are encouraged to utilize the provided pipeline and R scripts to process and analyze their own WGS dataset.
