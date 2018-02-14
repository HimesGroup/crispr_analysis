#!/usr/bin/python
import subprocess
import os
import argparse


def get_sample_info(fin):
	"""
	Open tab-delimited txt file containing the following columns:
	v0: sample_ID		| ID given to sample 
	v1: index			| Six digit sequence of the index for this library
	v2: file_directory	| Directory where sample's fastq files reside
	v3: project			| Name for project associated with sample
	v4: label			| Biological condition associated with the sample, provided by customer
	v5: ref_genome		| Rerence genome associated with sample. (options: "hg19", "hg38")
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


def get_genome_ref_files(genome):
	"""
	Location of all reference files needed for a given genome.
	The ERCC gtf files were appended separately to each species own gtf file
	Current choice: "hg38"
	"""
	if genome == "hg19":
		sgRNA_index = "/project/bhimeslab/Reference/GeCKO/geckov2"
	else:
		print 'Unknown genome selected: ', genome
	return(sgRNA_index)


def make_bcfile(sample_info_file, path_start):
	runs = get_sample_info(sample_info_file)
	#outp = open(path_start+"

	#for curr_sample, k in runs.iteritems():
	for k in runs:
		#Get sample information
		curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type = k

#Index sequencing primers
primer_start = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
primer_end = "TCTTGTGGAAAGGACGAAACACCG"


def gecko_main(fastq_in, project, sample_info_file, path_start):

	#Set up project and sample output directories
	if path_start == "./":
		path_start = os.getcwd()
	if path_start[-1] != "/":
		path_start = path_start+"/"		
	project_dir = path_start+project+"/"
	if not os.path.exists(project_dir):
		os.makedirs(project_dir)

	sgRNA_index = "/project/bhimeslab/Reference/GeCKO/geckov2"
	
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
		curr_file_name = project+"_"+curr_sample
		out_dir = project_dir+curr_sample+"/"
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)		
		
		#Write lsf file	
		job_name = curr_sample
		outp = open(job_name+".lsf", "w")
		outp.write("#!/bin/bash \n")
		outp.write("#BSUB -L /bin/bash\n")
		outp.write("#BSUB -J "+job_name+"\n")
		outp.write("#BSUB -q max_mem64 \n")
		outp.write("#BSUB -o "+job_name+"_%J.out\n")
		outp.write("#BSUB -e "+job_name+"_%J.screen\n")	

		R1 = project_dir+curr_file_name+".fastq"
		#outp.write("fastqc -o "+out_dir+" "+R1+"\n")
		#outp.write("cutadapt -g ^"+curr_barcode+primer_end+" -a GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC -q 10 -n 5 --minimum-length 19 --maximum-length 21 --overlap 10 -o "+out_dir+curr_sample+"_Trimmed.fastq "+R1+"\n")
		#outp.write("fastqc -o "+out_dir+" "+out_dir+curr_sample+"_Trimmed.fastq\n")
		#outp.write("bowtie2 -D 20 -R 3 -N 0 -L 16 -i L,0,0.05 --local -x "+sgRNA_index+" -U "+out_dir+curr_sample+"_Trimmed.fastq -S "+out_dir+curr_sample+".sam\n")
		outp.write("samtools view -hS "+out_dir+curr_sample+".sam | grep -e \"^@\" -e \"XM:i:[01][^0-9]\" > "+out_dir+curr_sample+".clean.sam\n")
		outp.write("samtools view -bS "+out_dir+curr_sample+".clean.sam > "+out_dir+curr_sample+".bam\n")
		outp.write("samtools sort "+out_dir+curr_sample+".bam "+out_dir+curr_sample+".sorted\n")
		outp.write("samtools index "+out_dir+curr_sample+".sorted.bam\n")
		outp.write("samtools idxstats "+out_dir+curr_sample+".sorted.bam > "+out_dir+curr_sample+".sorted.bam.stats\n")
		outp.write("bamtools stats -in "+out_dir+curr_sample+".sorted.bam > "+out_dir+curr_sample+".sorted.bamstats\n")
		i += 1
	outp.close()	


gecko_main("/project/bhimeslab/Arsenic_GeCKO/lane4_Undetermined.R1.fastq", "Arsenic_GeCKO", "/project/bhimeslab/Arsenic_GeCKO/Arsenic_GeCKO_Info_Sheet.txt", "/project/bhimeslab/")

def  make_gecko_fa(fin):
	#Read in csv file
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	#Create output file
	outp = open("human_geckov2_sgRNA.fa", 'w')
	for x in c:
		gene = x.split(',')[0]
		id = x.split(',')[1]
		seq = x.split(',')[2]
		outp.write(">"+gene+"_"+id+"\n"+seq+"\n")
	outp.close()

#make_gecko_fa("/project/bhimeslab/Reference/GeCKO/human_geckov2_library.csv")


#Indexes used for PE, SE, or SPE RNA-seq reads
illumina_indexes = \
	{'ATCACG':['Primer1', primer_start+'ATCACG'+primer_end], \
	 'CGATGT':['TruSeqAdapterIndex2', primer_start+'CGATGT'+primer_end], \
	 'TTAGGC':['TruSeqAdapterIndex3', primer_start+'TTAGGC'+primer_end], \
	 'TGACCA':['TruSeqAdapterIndex4', primer_start+'TGACCA'+primer_end], \
	 'ACAGTG':['TruSeqAdapterIndex5', primer_start+'ACAGTG'+primer_end], \
	 'GCCAAT':['TruSeqAdapterIndex6', primer_start+'GCCAAT'+primer_end]}


def make_adapter_fa(index, index_dictionary, out_name, library_type):
	"""
	For each sample, a specific adapter fasta file is created based on its library type and index 
	This file is used by Trimmomatic to perform adapter trimming
	Index sequences were obtained from manufacturer files and are assigned as used at PCPGM
	"""	 
	outp = open(out_name, "w")
	outp.write(">ReversePrimer/1\nCAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTCTACTATTCTTTCCCCTGCACTGT\n")
	if index in index_dictionary:
		outp.write(">"+index_dictionary[index][0]+"\n"+index_dictionary[index][1]+"\n")
	outp.close()


def main(sample_info_file, discovery, standard_trim, path_start):
	"""
	Dispatches an lsf job to locate fastq files that were output by Casava (GIGPAD A1 routine) and then:
	1) Perform adapter trimming
	2) Run fastqc
	3) Get unique reads 
	4) Run tophat to align reads to reference genome
	5) Obtain various QC metrics on aligned files
	6) Run cufflinks to quantify ERCC spike ins (if applicable)
	7) Run cufflinks to quantify all other transcripts in sample
	
	The directory structure below is specific to the UPenn HPC cluster
	"""
	runs = get_sample_info(sample_info_file)
	#for curr_sample, k in runs.iteritems():
	for k in runs:
		#Get sample information
		curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type = k

		#Get genome reference files
		ref_index, fa, gtf, ref, ERCC_gtf = get_genome_ref_files(ref_genome)

		#Set up project and sample output directories
		if path_start == "./":
			path_start = os.getcwd()
		if path_start[-1] != "/":
			path_start = path_start+"/"		
		project_dir = path_start+project+"/"
		if not os.path.exists(project_dir):
			os.makedirs(project_dir)
		out_dir = project_dir+curr_sample+"/"
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		
		job_name = curr_sample
				
		#This directory structure and naming convention is (obviously!) unique to UPenn 
		if library_type in ("SE", "DGE"):
			R1 = project_dir+curr_sample+".fastq.gz"
		else:
			R1 = project_dir+curr_sample+"_1.fastq.gz"
			R2 = project_dir+curr_sample+"_2.fastq.gz"
		
		#Make lsf file		
		outp = open(job_name+".lsf", "w")
		outp.write("#!/bin/bash \n")
		outp.write("#BSUB -L /bin/bash\n")
		outp.write("#BSUB -J "+job_name+"\n")
		outp.write("#BSUB -q max_mem64 \n")
		outp.write("#BSUB -o "+job_name+"_%J.out\n")
		outp.write("#BSUB -e "+job_name+"_%J.screen\n")
		outp.write("#BSUB -R 'rusage[mem=24000]'\n")
		outp.write("#BSUB -n 12\n")
		
		#Check whether unaligned fastq files that were processed by Casava are present
		#Create directory for each sample
		if not os.path.isfile(R1):
			print "R1 file not found ", R1
		if library_type in ("PE", "SPE"):
			if not os.path.isfile(R2):		
				print "R2 file not found", R2
		
		#In case need to unzip at some point:
		#outp.write("zcat "+R1+".gz > "+R1[:-3]+"\n")
		#outp.write("zcat "+R2+".gz > "+R2[:-3]+"\n")
				
		outp.write("cd "+out_dir+"\n")
		
		#Perform adapter trimming with trimmomatic
		#May perform a standard trimming of bases from reads by amount given above if standard_trim variable is greater than 0. Most will use standard_trim=0 
		#Create fa file of adapters specific to file
		if library_type in ["PE", "SE", "SPE"]:
			make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
		elif library_type == "DGE":
			make_adapter_fa(index, nextflex_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
		if standard_trim == 0:
			R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
			R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"		
			if library_type in ["PE", "SPE"]:
				make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
				outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
			elif library_type in ["SE", "DGE"]:				
				outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
		else:
			R1_trim = out_dir+curr_sample+"_R1_Trim"+str(standard_trim)+".fastq"
			R2_trim = out_dir+curr_sample+"_R2_Trim"+str(standard_trim)+".fastq"
			if library_type in ["PE", "SPE"]:
				outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq HEADCROP:"+standard_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
			elif library_type in ["SE", "DGE"]:
				outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" HEADCROP:"+standard_trim+"ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")					
								
		#Run fastqc on trimmed files.  
		#In some cases fastqc should be run on original files, but we have dropped this as a routine practice because the reports haven't changed much after trimming - adapter contamination has been minimal.
		if library_type in ["PE", "SPE"]:
			outp.write("fastqc -o "+out_dir+" "+R1_trim+" "+R2_trim+"\n")
		elif library_type in ("SE", "DGE"):
			outp.write("fastqc -o "+out_dir+" "+R1_trim+"\n")
		
		#Get total number of reads, unique reads, % unique reads from trimmed file(s). 
		outp.write("cat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+curr_sample+"_ReadCount\n")
		if library_type in ["PE", "SPE"]:
			outp.write("cat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+curr_sample+"_ReadCount\n")
		
		#Run TopHat with options specific to library type
		if discovery == "no":
			if library_type == "PE":
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type == "SPE":
				outp.write("tophat --library-type fr-firststrand -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type in ["DGE", "SE"]:
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+"\n")
		
		elif discovery == "yes":
			if library_type == "PE":
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type == "SPE":
				outp.write("tophat --library-type fr-firststrand -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type in ["DGE", "SE"]:
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+"\n")
				
		#Get samtools mapping stats
		outp.write("cd "+out_dir+"/tophat_out/\n")
		#Create sorted bam file:
		outp.write("samtools sort accepted_hits.bam "+curr_sample+"_accepted_hits.sorted\n")
		#Create indexed bam file:
		outp.write("samtools index "+curr_sample+"_accepted_hits.sorted.bam\n")
		#Write out index stats of where reads align to by chr:
		outp.write("samtools idxstats "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.stats\n")
		#Write out bamtools summary stats:
		outp.write("bamtools stats -in "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.bamstats\n")
		#Run CollectRnaSeqMetrics
		if library_type == "SPE":
			outp.write("java -Xmx2g -jar /opt/software/picard/picard-tools-1.96/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		else:
			outp.write("java -Xmx2g -jar /opt/software/picard/picard-tools-1.96/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=NONE INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		
		#Get number of reads spanning junctions by getting "N"s in CIGAR field of bam file
		#Be sure that Junction Spanning Reads are Added First then Unmapped Reads for proper ordering of fields in report
		outp.write("echo \"Junction Spanning Reads: \" $(bamtools filter -in "+curr_sample+"_accepted_hits.sorted.bam -script /home/bhimes/taffeta/cigarN.script | bamtools count ) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")
		#Get number of unmapped reads
		outp.write("echo Unmapped Reads: $(samtools view -c unmapped.bam) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")		

		#Gather metrics unique to paired-end samples using CollectInsertSizeMetrics
		if library_type in ["PE", "SPE"]:
			outp.write("java -Xmx2g -jar /opt/software/picard/picard-tools-1.96/CollectInsertSizeMetrics.jar HISTOGRAM_FILE="+curr_sample+"_InsertSizeHist.pdf INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_InsertSizeMetrics\n")	
		
		#Run cufflinks to count ERCC spike ins
		outp.write("mkdir "+out_dir+"/cufflinks_out_ERCC/\n")
		outp.write("cd "+out_dir+"/cufflinks_out_ERCC/\n")
		if library_type == "DGE":
			outp.write("cufflinks --library-type fr-unstranded --no-length-correction -G "+ERCC_only+" -p 12 ../tophat_out/accepted_hits.bam \n")
		elif library_type == "SPE":
			outp.write("cufflinks --library-type fr-firststrand -G "+ERCC_only+" -p 12 ../tophat_out/accepted_hits.bam \n")
		else:
			outp.write("cufflinks --library-type fr-unstranded -G "+ERCC_only+" -p 12 ../tophat_out/accepted_hits.bam \n")
		
		#Cufflinks to assemble and quantify transcripts
		outp.write("mkdir "+out_dir+"/cufflinks_out/\n")
		outp.write("cd "+out_dir+"/cufflinks_out/\n")
		if library_type == "DGE":
			outp.write("cufflinks --library-type fr-unstranded --no-length-correction -G "+gtf+" -p 12 ../tophat_out/accepted_hits.bam \n")
		elif library_type == "SPE":
			outp.write("cufflinks --library-type fr-firststrand -G "+gtf+" -p 12 ../tophat_out/accepted_hits.bam \n")			
		else:
			outp.write("cufflinks --library-type fr-unstranded -G "+gtf+" -p 12 ../tophat_out/accepted_hits.bam \n")
		outp.close()
	
		#subprocess.call("bsub < "+job_name+".lsf", shell=True)
		#subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


# if __name__ == "__main__":
# 	parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a project.")
# 	parser.add_argument("--standard_trim", default=0, type=int, help="Number of bases to be trimmed from leftmost end of all reads (default=0)")
# 	parser.add_argument("--discovery", default="yes", type=str, help="Should TopHat be run with options to discover novel transcripts (i.e. disable --no-novel-juncs --transcriptome-only)? "
# 		"(options: yes, no; default=yes) "
# 		"Note: the 'no' option only works with hg19 at the moment")
# 	parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
# 	parser.add_argument("samples_in", help="Path to a tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
# 	args = parser.parse_args()
# 	main(args.samples_in, args.discovery, args.standard_trim, args.path_start)

