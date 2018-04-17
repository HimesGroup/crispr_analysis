#!/usr/bin/python
import subprocess
import os
#import argparse


def get_sample_info(fin):
	"""
	Open tab-delimited txt file containing the following columns:
	v0: sample_ID		| ID given to sample 
	v1: index			| Six digit sequence of the index for this library
	v2: file_directory	| Directory where sample's fastq files reside
	v3: project			| Name for project associated with sample
	v4: label			| Biological condition associated with the sample, provided by customer
	v5: ref_genome		| Rerence genome associated with sample. (options: "hg19", "gecko", "gpcr_CRISPR")
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	d = []
	for x in c:
		sample_id = x.split('\t')[0]
		index = x.split('\t')[1]
		label = x.split('\t')[2]
		ref_genome = x.split('\t')[3]
		d.append([sample_id, index, label, ref_genome])
	return d


### Script to create new fa files for alignment when a new gecko library is provided
def  make_gecko_fa(fin):
	#Read in csv file
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	#Create output file, depending on ref_genome listed in info sheet
	if ref_genome == 'gecko':
		outp = open("human_geckov2_sgRNA.fa", 'w')
	if ref_genome == 'gpcr_CRISPR':
		outp = open("/project/bhimeslab/Reference/gpcr_CRISPR/gpcr_CRISPR_sgRNA.fa", 'w')
	for x in c:
		gene = x.split(',')[0]
		id = x.split(',')[1]
		seq = x.split(',')[2]
		outp.write(">"+gene+"_"+id+"\n"+seq+"\n")
	outp.close()

if ref_genome == 'gecko':
	make_gecko_fa("/project/bhimeslab/Reference/GeCKO/human_geckov2_library.csv")
if ref_genome == 'gpcr_CRISPR':
	make_gecko_fa("/project/bhimeslab/gpcr_screen/gpcr_crispr_library.csv")


def get_genome_ref_files(genome):
	"""
	Location of all sgRNA reference files needed for a given gecko library.
	"""
	if genome == "gecko":
		sgRNA_index = "/project/bhimeslab/Reference/GeCKO/geckov2"
	if genome == "gpcr_CRISPR":
		sgRNA_index = "/project/bhimeslab/Reference/gpcr_CRISPR/gpcr_CRISPR_sgRNA.fa"
	else:
		print 'Unknown genome selected: ', genome
	return(sgRNA_index)



### Script that demultiplexes sequentially after trimming a single 5' base since some Gecko libraries have barcodes within reads
#### Rationale for shifting 88 positions:
	# TCTTGTGGAAAGGACGAAACACCG = 5' primer 24bp
	# 123456789012345678901234
	# AACGATCGATAGGTAAGG = longest barcode = 18bp
	# sgRNAs = 20bp
	# whole read = 150bp
	# 150-44-18=88 could have up to 88 shifted positions and still have sgRNA present
def gecko_demultiplex(fastq_in, project, sample_info_file, path_start, bcfile):
	runs = get_sample_info(sample_info_file)

	base_name = fastq_in.split("/")[-1].split(".fastq")[0]
	for i in range(1, 88):
		prev_base_name = base_name+".Round"+str(i)
		new_base_name = base_name+".Round"+str(i+1)
		if i == 1:
			subprocess.call("fastq-multx -b "+bcfile+" "+fastq_in+" -o "+base_name+".%.fastq", shell=True)
			subprocess.call("wc -l "+base_name+".*.fastq >> demultiplexed_counts.txt", shell=True)
			subprocess.call("cutadapt -u 1 -o "+new_base_name+".fastq "+base_name+".unmatched.fastq", shell=True)
		else:
			subprocess.call("cutadapt -u 1 -o "+new_base_name+".fastq "+prev_base_name+".unmatched.fastq", shell=True)
		subprocess.call("fastq-multx -b "+bcfile+" "+new_base_name+".fastq -o "+new_base_name+".%.fastq", shell=True)
		for k in runs:
			#Get sample information
			curr_sample, curr_barcode, label, ref_genome = k
			subprocess.call("cat "+new_base_name+"."+curr_sample+".fastq >> "+base_name+"."+curr_sample+".fastq", shell=True)
			subprocess.call("wc -l "+new_base_name+"."+curr_sample+".fastq >> demultiplexed_counts.txt", shell=True)
			subprocess.call("rm "+new_base_name+"."+curr_sample+".fastq", shell=True)

		

if ref_genome == 'gecko':
	gecko_demultiplex("/project/bhimeslab/gecko/gpcr_NoIndex.R1.fastq", "gpcr", "/project/bhimeslab/gecko/gpcr_Info_Sheet.txt", "/project/bhimeslab/", "/project/bhimeslab/gecko/gpcr_bcfile.txt")
if ref_genome == 'gpcr_CRISPR':
	gecko_demultiplex("/project/bhimeslab/gecko/mir_NoIndex.R1.fastq", "mir", "/project/bhimeslab/gecko/mir_Info_Sheet.txt", "/project/bhimeslab/", "/project/bhimeslab/gecko/mir_bcfile.txt")


def read_demultiplex_log(fin, project, sample_info_file):
	"""
	Read gecko log file from loop of demultiplex runs and extract the count of demultiplexed reads at each cycle
	"""
	runs = get_sample_info(sample_info_file)
	keys = []
	for k in runs:
		#Get sample information
		curr_sample, curr_barcode, label, ref_genome = k
		keys.append(curr_sample)
	f = open(fin, 'r')
	c = f.read().split('Id\tCount\tFile(s)\n')
	d = dict.fromkeys(keys)
	for x in c[1:]:
		curr_demultiplex = x.split('This is cutadapt 1.9 with Python 2.7.5\n')[0]
		curr_hist = curr_demultiplex.split('\n')
		for y in curr_hist:
			k = y.split('\t')[0]
			if k in d.keys():
				entry = y.split('\t')[1]
				if d[k] ==  None:
					d[k] = [entry]
				else:
					d[k] += [entry]
	outp = open(project+"_demultiplex_counts.txt", "w")
	outp.write("\t".join(d.keys())+"\n")
	count_lines = zip(*d.values())
	for l in count_lines:
		outp.write("\t".join(l)+"\n")
	outp.close()

read_demultiplex_log("gecko_demultiplex_gpcr_030916.log", "gpcr", "gpcr_Info_Sheet.txt")

def gecko_main(fastq_in, project, sample_info_file, path_start):

	#Set up project and sample output directories
	if path_start == "./":
		path_start = os.getcwd()
	if path_start[-1] != "/":
		path_start = path_start+"/"		
	project_dir = path_start+project+"/"
	if not os.path.exists(project_dir):
		os.makedirs(project_dir)

	sgRNA_index = get_genome_ref_files(ref_genome)
	
	#Make lsf file		
	#Make bcfile
	#barcode_f = open("/project/bhimeslab/Arsenic_GeCKO/arsenic_tag_bcfile.txt", 'r')
	#barcodes = barcode_f.read().split('\n')
	#outp.write("cat lane4_Undetermined.R1.fastq | fastx_barcode_splitter.pl --bcfile /project/bhimeslab/Arsenic_GeCKO/arsenic_tag_bcfile.txt --prefix /project/bhimeslab/Arsenic_GeCKO/Arsenic_GeCKO_ --bol --suffix \".fastq\"\n")

	runs = get_sample_info(sample_info_file)
	i=0
	for k in runs:

		#Get sample information
		curr_sample, curr_barcode, label, ref_genome = k
		curr_file_name = project+"_NoIndex.R1."+curr_sample
		out_dir = project_dir+curr_sample+"/"
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)		
		
		#Write lsf file	
		job_name = curr_sample
		outp = open(job_name+".lsf", "w")
		outp.write("#!/bin/bash \n")
		outp.write("#BSUB -L /bin/bash\n")
		outp.write("#BSUB -J "+job_name+"\n")
		outp.write("#BSUB -q normal \n")
		outp.write("#BSUB -o "+job_name+"_%J.out\n")
		outp.write("#BSUB -e "+job_name+"_%J.screen\n")	

		R1 = curr_file_name+".fastq"
		#outp.write("fastqc -o "+out_dir+" "+R1+"\n")
		#outp.write("cutadapt -g ^"+curr_barcode+primer_end+" -a GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC -q 10 -n 5 --minimum-length 19 --maximum-length 21 --overlap 10 -o "+out_dir+curr_sample+"_Trimmed.fastq "+R1+"\n")
		outp.write("cutadapt -g ^TCTTGTGGAAAGGACGAAACACCG -a GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTCGCTAGC -q 10 -n 5 --minimum-length 19 --maximum-length 21 --overlap 10 -o "+out_dir+curr_sample+"_Trimmed.fastq "+R1+"\n")
		#outp.write("fastqc -o "+out_dir+" "+out_dir+curr_sample+"_Trimmed.fastq\n")
		outp.write("bowtie2 -D 20 -R 3 -N 0 -L 16 -i L,0,0.05 --local -x "+sgRNA_index+" -U "+out_dir+curr_sample+"_Trimmed.fastq -S "+out_dir+curr_sample+".sam\n")
		outp.write("samtools view -hS "+out_dir+curr_sample+".sam | grep -e \"^@\" -e \"XM:i:[01][^0-9]\" > "+out_dir+curr_sample+".clean.sam\n")
		outp.write("samtools view -bS "+out_dir+curr_sample+".clean.sam > "+out_dir+curr_sample+".bam\n")
		outp.write("samtools sort "+out_dir+curr_sample+".bam "+out_dir+curr_sample+".sorted\n")
		outp.write("samtools index "+out_dir+curr_sample+".sorted.bam\n")
		outp.write("samtools idxstats "+out_dir+curr_sample+".sorted.bam > "+out_dir+curr_sample+".sorted.bam.stats\n")
		outp.write("bamtools stats -in "+out_dir+curr_sample+".sorted.bam > "+out_dir+curr_sample+".sorted.bamstats\n")
		i += 1
	outp.close()	

gecko_main("gpcr_NoIndex.R1.fastq", "gpcr", "gpcr_Info_Sheet.txt", "/project/bhimeslab/gecko/")
