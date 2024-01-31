#########################################
### QA/QC QIIME2 output files
### from MR DNA for pipe type study
### Lou LaMartina, finalized Jan 25, 2024
#########################################


setwd("~/Desktop/Postdoc/16S/PipeType")
library(readxl)


#################
### load data ###
#################

# see 01_PipeType_dataPrep.R

# OTU counts
counts <- read.csv("DATA/PipeType_OTU_counts.csv")
rownames(counts) <- counts$Sample_ID
counts <- counts[-1]


# taxonomy
taxa <- read.csv("DATA/PipeType_OTU_taxonomy.csv")
rownames(taxa) <- taxa$OTU
identical(colnames(counts), rownames(taxa))


# sample metadata
info <- data.frame(read_xlsx("DATA/PipeType_sample_metadata.xlsx"))
rownames(info) <- info$Sample_ID
info <- info[match(rownames(counts), info$Sample_ID),]
identical(rownames(counts), rownames(info))


# convert to relative abundance
relabun <- counts / rowSums(counts)




#####################
### contamination ###
#####################

# subset blanks
blanks <- counts[grepl("blank", rownames(counts), ignore.case = T),]


# reads per sample
cat(as.integer(mean(rowSums(counts))), "±", as.integer(sd(rowSums(counts))))
# 29899 ± 1028


# reads per blank
cat(as.integer(mean(rowSums(blanks))), "±", as.integer(sd(rowSums(blanks))))
# 29344 ± 210


# number of reads in blanks vs samples
contams <- data.frame(OTU = colnames(relabun), blanks = colSums(blanks / rowSums(blanks)),
                      samples = colSums(relabun[! rownames(relabun) %in% rownames(blanks),]),
                      nreads = colSums(counts))


# OTUs that are 3x more abundant in blanks than samples
contams$blanks_x3 <- contams$blanks * 3
contams <- subset(contams, blanks_x3 > samples)
contams <- merge(contams, taxa[c("OTU", "Genus_species", "Accession")], by = "OTU")
write.csv(contams, "DATA/PipeType_contaminants.csv", row.names = F, na = "", quote = F)


# proportion in dataset
cat(signif(sum(counts[contams$OTU]) / sum(counts) * 100, 2), "%")
# 0.056 %


# remove contams
counts <- counts[! rownames(counts) %in% rownames(blanks), ! colnames(counts) %in% contams]
counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
info <- info[match(rownames(counts), info$Sample_ID),]
taxa <- subset(taxa, OTU %in% colnames(counts))
identical(rownames(counts), rownames(info))
identical(colnames(counts), rownames(taxa))




############
### save ###
############

# OTU abundances
counts <- data.frame(Sample_ID = rownames(counts), counts)
write.csv(counts, "DATA/PipeType_OTU_counts.csv", row.names = F, na = "", quote = F)


# OTU taxonomy
write.csv(taxa, "DATA/PipeType_OTU_taxonomy.csv", row.names = F, na = "", quote = F)

