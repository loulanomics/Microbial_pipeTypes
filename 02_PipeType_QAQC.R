#########################################
### QA/QC QIIME2 output files
### from MR DNA for pipe type study
### Lou LaMartina, finalized Jan 11, 2024
#########################################


setwd("~/Desktop/Postdoc/16S/PipeType")
library(vegan)
library(readxl)
library(scales)


#################
### load data ###
#################

# see 01_PipeType_dataPrep.R

# OTU counts
counts <- read.csv("DATA/PipeType_OTU_counts.csv")
rownames(counts) <- counts$ID
counts <- counts[-1]


# taxonomy
taxa <- read.csv("DATA/PipeType_OTU_taxonomy.csv")
rownames(taxa) <- taxa$OTU
identical(colnames(counts), rownames(taxa))


# sample metadata
info <- data.frame(read_xlsx("DATA/PipeType_sample_metadata.xlsx", sheet = "QC"))
rownames(info) <- info$ID
info <- info[match(rownames(counts), info$ID),]
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


# OTUs in blanks
blank_OTUs <- colSums(relabun[1:3,])


# OTUs in samples
sample_OTUs <- colSums(relabun[-1:-3,])


# OTUs that are 3x more abundant in blanks than samples
contams <- names(which((blank_OTUs * 3) > sample_OTUs))


# proportion in dataset
cat(signif(sum(counts[contams]) / sum(counts) * 100, 2), "%")
# 0.056 %


# remove those
counts_nocontam <- counts[! colnames(counts) %in% contams]
contams <- data.frame(OTU = contams, 
                      subset(taxa, OTU %in% contams)[c("Genus_species", "Accession")], 
                      numReads = colSums(counts[contams]),
                      Prop = colSums(counts[contams]) / sum(counts[! colnames(counts) %in% contams]))
#write.csv(contams, "DATA/PipeType_contaminants.csv", row.names = F, na = "", quote = F)




#################
### PERMANOVA ###
#################

### with blanks ###

# variables to test
vars <- info[c(1,4:6,8:10)]
vars[is.na(vars)] <- "blank"
vars <- subset(vars, Age != "blank")
vars <- data.frame(apply(vars, 2, as.character))


# samples to test
perman <- relabun[rownames(relabun) %in% vars$ID,]
identical(rownames(perman), vars$ID)


# calculate distances
bray.stat <- vegdist(perman, method = "bray")


# PERMANOVA
perman <- adonis2(bray.stat ~ Material * Surface * Location * Line * Age * Diameter_in, data = vars)
perman_wBlanks <- data.frame(Variable = rownames(perman), perman)
colnames(perman_wBlanks)[6] <- "p"



### without blanks ###

# variables to test
vars <- vars[-1:-3,]


# samples to test
perman <- relabun[rownames(relabun) %in% vars$ID,]
identical(rownames(perman), vars$ID)


# calculate distances
bray.stat <- vegdist(perman, method = "bray")


# PERMANOVA
perman <- adonis2(bray.stat ~ Material * Surface * Location * Line * Age * Diameter_in, data = vars)
perman_noBlanks <- data.frame(Variable = rownames(perman), perman)
colnames(perman_noBlanks)[6] <- "p"




############
### save ###
############

# counts - no blanks
counts_noblanks <- counts_nocontam[-1:-3,]
counts_noblanks <- data.frame(ID = rownames(counts_noblanks), counts_noblanks)
write.csv(counts_noblanks, "DATA/PipeType_OTU_counts_noBlanks.csv", row.names = F, na = "", quote = F)


# taxa - no contams
taxa <- subset(taxa, ! OTU %in% contams$OTU)
identical(taxa$OTU, colnames(counts_noblanks[-1]))
write.csv(counts_noblanks, "DATA/PipeType_OTU_taxonomy_noBlanks.csv", row.names = F, na = "", quote = F)


