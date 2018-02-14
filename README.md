# crispr_analysis
*in progress*

Workflow used to analyze CRISPR knockout screen data

Authors: Maya Shumyatcher, Blanca Himes.

This set of scripts was initially developed to analyze CRISPR knockout screen data associated with the manuscript "Genome-wide CRISPR screen identifies novel suppressors of endoplasmic reticulum stress induced apoptosis". The goal is to take fastq files (sequenced with Illumina HiSeq) associated with a "project" and:

* Perform preliminary QC
* Align reads to GeCKO or CRISPR reference files
* Perform QC on aligned files
* Perform differential expression of reads aligned to transcripts according to a given reference genome

Several freely available software packages are used to perform most of these steps (see below). 

### Dependencies

The Tuxedo suite of tools and various other programs should be installed: bowtie or bowtie2, tophat, cufflinks, cummerbund, fastqc, trimmomatic, samtools, bamtools, picardtools.
Annotation files should be available: reference genome fasta, gtf, refFlat, and index files. We use ERCC spike-ins, so our reference files include ERCC mix transcripts.
For adapter trimming, we include Ilumina sequences 
The Python scripts make use of modules that include subprocess, os, argparse, sys.
R and various libraries should be available, including 

### input files

Before running the pipeline, characteristics of a set of fastq files for samples that are part of a project are described in a tab-delimited txt file containing the following fields:

	customer_ID		| ID given to sample by customer
	gigpad_ID		| A sample ID that may differ from customer_ID and that is associated with library prior to pooling
	lane			| Lane of sequencer where this sample's pool was run, needed to locate raw fastq file
	index			| Six digit sequence of the library index used for this sample, needed to locate raw fastq file and perform adapter trimming
	ercc_mix		| Mix of ERCC spike used for library construction (options: "1", "2", "-")
	file_directory	| Top-level directory where sample's fastq files were written to following Casava filters
	batch			| gigpad batch number associated with sample, needed to locate raw fastq file and used as directory name for pipeline output files
	label			| Biological condition associated with the sample, provided by customer
	ref_genome		| Rerence genome associated with sample. (options: "hg19", "Zv9", "mm10")
	library_type	| Type of library for sample (options: "PE", "SE", "DGE", "SPE",
							corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end")


### Workflow
