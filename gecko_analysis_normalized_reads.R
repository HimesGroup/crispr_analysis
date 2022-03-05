
#projects include mir_screen, gpcr_screen, gecko
project = "gpcr_screen"
sample_info = read.table(paste0("/project/bhimeslab/",project,"/",project,"_Info_Sheet.txt"), header=T, as.is=T, sep="\t")


for (i in c(1:length(sample_info$Sample))){
	curr_sample = sample_info$Sample[i]
	curr_reads = read.table(paste0(project, "/", curr_sample, "/", curr_sample, ".sorted.bam.stats"), header=F, as.is=T)
	if (i == 1) {
		raw.reads = data.frame(curr_reads$V1)
		names(raw.reads) = c("sgRNA")
	}
	raw.reads$new = curr_reads$V3
	names(raw.reads)[i+1] = curr_sample
}	

#Normalize reads
normalized.counts = as.matrix(raw.reads[, -1])
for (i in c(1:length(sample_info$Sample))){
	sample_count = sum(normalized.counts[, i])
	normalized.counts[, i] = 1 + 10^6 * normalized.counts[, i] / sample_count
}

write.table(raw.reads, file=paste0(project, "_raw_reads.txt"), row.names=FALSE, quote=TRUE)
write.table(normalized.counts, file=paste0(project,"_normalized_counts.txt"), row.names=FALSE, quote=FALSE)


