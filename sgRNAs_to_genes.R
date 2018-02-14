rm(list=ls())

library("ggplot2")
library("stringr")
library("plyr")
library("dplyr")

############################################################################
############################# By sgRNA #####################################
############################################################################

raw_reads <- read.table(file="C:/Users/mayashum/Dropbox/CRISPR/arsenic/arsenic_raw_reads.txt", header=TRUE)

#read in data -- counts_by_sgRNA is initial data
counts_by_sgRNA <- read.table(file="C:/Users/mayashum/Dropbox/CRISPR/arsenic/arsenic_normalized_counts.txt", header=TRUE, stringsAsFactors=FALSE)
counts_by_sgRNA <- data.frame(raw_reads$sgRNA, counts_by_sgRNA)
colnames(counts_by_sgRNA) <- c("id","control", "arsenic1", "arsenic2")
counts_by_sgRNA$mean_ars1_ars2 <- rowMeans(counts_by_sgRNA[,c("arsenic1", "arsenic2")])
counts_by_sgRNA$gene <- str_split_fixed(counts_by_sgRNA$id,"_", 2)[,1]
counts_by_sgRNA$sgRNA <- str_split_fixed(counts_by_sgRNA$id,"_", 2)[,2]

#split the data frame by gene name
#rank sgRNAs with respect to other sgRNAs for that gene
counts_by_sgRNA_split <- counts_by_sgRNA %>%
  group_by(gene) %>%
  mutate(rank_gene = row_number(-mean_ars1_ars2))

#rank sgRNAs with respect to all sgRNAs 
counts_by_sgRNA_split$rank_univ <- rank(-counts_by_sgRNA_split$mean_ars1_ars2, ties.method = "first")

#count total number of sgRNAs for each gene
max_sgRNA_rank <- aggregate(counts_by_sgRNA_split$rank_gene, by = list(counts_by_sgRNA_split$gene), max)
colnames(max_sgRNA_rank) <- c("gene", "sgRNA_count")
counts_by_sgRNA_split <- merge(counts_by_sgRNA_split, max_sgRNA_rank, by.x="gene", by.y="gene")


#find fold change of arsenic vs control
#do this for arsenic1, arsenic2 and mean_ars1_ars2
#count total number of sgRNAs with fold change greater than 2 for each gene

#fold change arsenic1
counts_by_sgRNA_split$fc_ars_1 <- counts_by_sgRNA_split$arsenic1/counts_by_sgRNA_split$control

#is FC > 2?  yes/no
for (i in 1:nrow(counts_by_sgRNA_split)) {
  if (counts_by_sgRNA_split$fc_ars_1[i] > 2) {
    counts_by_sgRNA_split$diff_exp_1[i] = 1
  }
  else {
    counts_by_sgRNA_split$diff_exp_1[i] = 0
  }
}

#for arsenic2
counts_by_sgRNA_split$fc_ars_2 <- counts_by_sgRNA_split$arsenic2/counts_by_sgRNA_split$control

for (i in 1:nrow(counts_by_sgRNA_split)) {
  if (counts_by_sgRNA_split$fc_ars_2[i] > 2) {
    counts_by_sgRNA_split$diff_exp_2[i] = 1
  }
  else {
    counts_by_sgRNA_split$diff_exp_2[i] = 0
  }
}

#for mean_ars1_ars2
counts_by_sgRNA_split$fc_mean <- counts_by_sgRNA_split$mean_ars1_ars2/counts_by_sgRNA_split$control

for (i in 1:nrow(counts_by_sgRNA_split)) {
  if (counts_by_sgRNA_split$fc_mean[i] > 2) {
    counts_by_sgRNA_split$diff_exp_mean[i] = 1
  }
  else {
    counts_by_sgRNA_split$diff_exp_mean[i] = 0
  }
}

#count number of sgRNAs per gene w/ FC > 2
diff_exp_count_1 <- aggregate(counts_by_sgRNA_split$diff_exp_1, by = list(counts_by_sgRNA_split$gene), sum)
diff_exp_count_2 <- aggregate(counts_by_sgRNA_split$diff_exp_2, by = list(counts_by_sgRNA_split$gene), sum)
diff_exp_count_mean <- aggregate(counts_by_sgRNA_split$diff_exp_mean, by = list(counts_by_sgRNA_split$gene), sum)

colnames(diff_exp_count_1) <- c("gene", "diff_exp_sgRNAs_count_1")
colnames(diff_exp_count_2) <- c("gene", "diff_exp_sgRNAs_count_2")
colnames(diff_exp_count_mean) <- c("gene", "diff_exp_sgRNAs_count_mean")

counts_by_sgRNA_split <- merge(counts_by_sgRNA_split, diff_exp_count_1, by.x="gene", by.y="gene") 
counts_by_sgRNA_split <- merge(counts_by_sgRNA_split, diff_exp_count_2, by.x="gene", by.y="gene") 
counts_by_sgRNA_split <- merge(counts_by_sgRNA_split, diff_exp_count_mean, by.x="gene", by.y="gene")

#rank genes by the value of their highest read sgRNA
max_sgRNA_val <- aggregate(counts_by_sgRNA_split$mean_ars1_ars2, by = list(counts_by_sgRNA_split$gene), max)
colnames(max_sgRNA_val) <- c("gene", "max_val")
counts_by_sgRNA_split <- merge(counts_by_sgRNA_split, max_sgRNA_val, by.x="gene", by.y="gene")
counts_by_sgRNA_split <- counts_by_sgRNA_split[order(-counts_by_sgRNA_split$max_val, counts_by_sgRNA_split$rank_gene),]

#scatterplots of counts_by_sgRNA
pairs(~log10(control)+log10(arsenic1)+log10(arsenic2), data=counts_by_sgRNA, main="Arsenic: log normalized counts by sgRNA")

write.csv(counts_by_sgRNA_split, file="C:/Users/mayashum/Dropbox/CRISPR/Arsenic/counts_by_sgRNA_split.txt", row.names=F, quote=F)


############################################################################
############################## By Gene #####################################
############################################################################

#counts by gene
#average together all sgRNA reads that correspond to a gene
#includes mean by gene
counts_by_gene <- ddply(counts_by_sgRNA_split, "gene", summarize, sd_ctrl=sd(control), ctrl=mean(control), sd_arsenic1=sd(arsenic1), arsenic1=mean(arsenic1), sd_arsenic2=sd(arsenic2), arsenic2=mean(arsenic2), sd_arsenic_mean=sd(mean_ars1_ars2), mean_ars1_ars2=mean(mean_ars1_ars2), sgRNA_count=mean(sgRNA_count), fc_ars_1=mean(fc_ars_1), fc_ars_2=mean(fc_ars_2), fc_mean=mean(fc_mean), diff_exp_sgRNAs_count_1=mean(diff_exp_sgRNAs_count_1), diff_exp_sgRNAs_count_2=mean(diff_exp_sgRNAs_count_2), diff_exp_sgRNAs_count_mean=mean(diff_exp_sgRNAs_count_mean))
pairs(~log10(ctrl)+log10(arsenic1)+log10(arsenic2), data=counts_by_gene, main="Arsenic - log normalized counts by gene")
counts_by_gene <- counts_by_gene[order(-counts_by_gene$mean_ars1_ars2),,drop=FALSE]
write.csv(counts_by_gene, file="C:/Users/mayashum/Dropbox/CRISPR/Arsenic/counts_by_gene_new.txt", row.names=F, quote=F)
