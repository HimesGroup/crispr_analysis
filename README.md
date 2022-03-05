# crispr_analysis

Workflow used to analyze CRISPR knockout screen data

Authors: Maya Shumyatcher, Blanca Himes.

This set of scripts was initially developed to analyze CRISPR knockout screen data associated with the manuscript "Genome-wide CRISPR screen identifies novel suppressors of endoplasmic reticulum stress induced apoptosis". The goal is to take fastq files (sequenced with Illumina HiSeq) associated with a "project" and:

* Perform preliminary QC
* Align reads to GeCKO or CRISPR reference files
* Perform QC on aligned files
* Perform differential expression of reads aligned to transcripts according to a given reference file

Several freely available software packages are used to perform most of these steps (see below). 

### Dependencies

Programs that should be installed: ea-utils, bowtie2, cutadapt, samtools, bamtools.
Annotation files should be available: reference crispr libraries. 
For trimming, sequences based on crispr library used
The Python scripts make use of modules that include subprocess, os.
R and various libraries should be available, including stringr, plyr and dplyr.

### Input Files

Before running the pipeline, characteristics of a set of fastq files for samples that are part of a project are described in a tab-delimited txt file containing the following fields:

```
sample_ID		| ID given to sample 
index			| Six digit sequence of the index for this library
file_directory 	        | Directory where sample's fastq files reside
project			| Name for project associated with sample
label			| Biological condition associated with the sample, provided by customer
ref_genome		| Rerence genome associated with sample. (options: "gecko", "gpcr_CRISPR")
 ```
Note that there are two gRNA library options. The "gecko" option refers to the human GeCKO (Genome-Scale CRISPR Knock-Out) lentiviral pooled libraries, available from [Addgene](https://www.addgene.org/pooled-library/zhang-human-gecko-v2/). The GeCKO v2 libraries consist of over 100,000 unique gRNAs for gene knock-out and target early consecutive exons for genome editing. The "gpcr_CRISPR" option refers to a library designed by [Lu Lab](https://www.hsph.harvard.edu/quan-lu/) at the Harvard School of Public Health to target miRNA and GPCRs.


### Workflow

1) Write and execute an lsf job to perform read alignment for RNA-seq samples associated with a project using gecko_analysis.py:

    > python gecko_analysis.py <i>sample_info_file.txt</i>

     Following execution of this script, various output files will be written for each sample in directories structured as:
    > 
    <i>batch_num</i>/<i>sample_name</i>/sample_out <br>

    This script creates a reference file based on the genome selected. Note that reference file names are hard coded in the python script and should be modified if applicable. The script then demultiplexes sequentially after trimming a single 5' base, since some Gecko libraries have barcodes within reads (i.e. in our dataset, we found that barcodes were not always 5' anchored). More details about the demultiplexing methodology and rationale can be found within the python script. Finally, it writes an lsf file to perform read alignment.
 
2) Normalize raw gRNA counts to the total number of mapped reads for each sample using code in gecko_analysis_normalized_reads.R 

3) Obtain per-sgRNA and per-gene results using the scripts in sgRNAs_to_genes.R 

    Our study design used two replicates and one control sample. If desired, these scripts can be modified to accomodate a different study design. 
    
    Because only two replicates were used, differential expression analysis was done simply by computing fold changes. The scripts in sgRNAs_to_genes.R compute fold changes for case vs. control for replicate1, case2 and the mean of case1 and case2. Ultimately, we used the results based on the mean of case1 and case2 for our study. 
    
    The scripts in sgRNAs_to_genes.R tally the number of sgRNAs with fold change > 2 for each gene and provide a few different methods of ranking per-gene results. For instance, per-gene results may be ranked by the number of sgRNAs with fold change > 2 corresponding to each gene, or by first ranking individual sgRNAs by fold change and subsequently matching these to genes.