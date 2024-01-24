#########################################
### organize QIIME2 output files
### from MR DNA for pipe type study
### Lou LaMartina, finalized Jan 23, 2024
#########################################


setwd("~/Desktop/Postdoc/16S/PipeType/DATA")
library(readxl)
library(stringr)


#################
### load data ###
#################

# QIIME2 results
otus <- data.frame(read_xlsx("032522MDillcus515F-zotus.fa.bacteria.OTU.xlsx"))


# sample info
info <- data.frame(read_xlsx("PipeType_sample_metadata.xlsx"))


# failed samples
fails <- data.frame(read_xlsx("032522MD negative listb.xlsx"))




##################
### OTU counts ###
##################

# extract counts
counts <- otus[c(1,8:99)]
rownames(counts) <- counts$otu.name
counts <- counts[-1]
dim(counts) # 4262   92


# remove failed & empty samples
gsub("-| ", "\\.", fails$Sample.name) %in% colnames(counts)
counts <- counts[! colnames(counts) %in% gsub("-| ", "\\.", fails$Sample.name)]
counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
dim(counts) # 4088   73


# change OTU ID
rownames(counts) <- paste0("OTU", sprintf("%04d", as.numeric(gsub("Zotu", "", rownames(counts)))))
identical(rownames(counts), sort(rownames(counts)))


# transpose, keep file name
counts <- data.frame(File_name = colnames(counts), t(counts))
identical(colnames(counts), sort(colnames(counts)))


# R changed the hyphens to periods
counts$File_name <- gsub("\\.", "-", counts$File_name)


# merge, keep blanks
counts <- merge(info[c("File_name", "Sample_ID")], counts, by = "File_name")


# change row names
rownames(counts) <- counts$Sample_ID
counts <- counts[-1:-2]
identical(colnames(counts), sort(colnames(counts)))




################
### taxonomy ###
################

# change Zotu to OTU, add leading zeros
taxa <- otus[1:2]
taxa$otu.name <- paste0("OTU", sprintf("%04d", as.numeric(gsub("Zotu", "", taxa$otu.name))))
colnames(taxa)[1] <- "OTU"
rownames(taxa) <- taxa$OTU
taxa <- subset(taxa, OTU %in% colnames(counts))
identical(colnames(counts), rownames(taxa))


# OTU IDs and taxonomy assignment strings, separated by semi colons
taxa.ls <- strsplit(taxa$Taxonomy, ";")


# extract classifications based on starting with k__, p__, etc
splittaxa.ls <- list()
for (i in paste0("^", c("k", "p", "c", "o", "f", "g", "s"), "__")) {
  for (j in 1:length(taxa.ls)) {
    splittaxa.ls[[i]][[j]] <- unlist(taxa.ls[j])[grep(i, unlist(taxa.ls[j]))]
    splittaxa.ls[[i]][[j]] <- gsub(i, "", splittaxa.ls[[i]][[j]])
  }
}


# combine
splittaxa <- data.frame(do.call(cbind, splittaxa.ls))
colnames(splittaxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Genus_species")
splittaxa$Genus_species <- gsub(" ", "_", splittaxa$Genus_species)
splittaxa <- apply(splittaxa, 2, str_to_title)
taxa <- cbind(taxa, splittaxa)
identical(colnames(counts), rownames(taxa))


# extract accession numbers
taxa$Accession <- NA
taxa$Accession[grep("\\.1 ", taxa$Taxonomy)] <- taxa$Taxonomy[grep("\\.1 ", taxa$Taxonomy)]
taxa$Accession <- sapply(strsplit(taxa$Accession, " "), '[', 1)
taxa <- taxa[-2]
identical(colnames(counts), rownames(taxa))


# cleanup
taxa <- data.frame(apply(taxa, 2, function(x) gsub("Candidatus |Candidatus_", "", x)))
taxa$Genus_species <- str_to_sentence(taxa$Genus_species)
taxa$Genus[taxa$Genus_species == "[Ruminococcus]_Torques"] <- "Ruminococcus"
taxa$Genus_species[taxa$Genus_species == "[Ruminococcus]_Torques"] <- "Ruminococcus_torques"
taxa$Genus_species[grep("Psychrosinus", taxa$Genus_species)] <- "Psychrosinus_fermentans"
taxa$Genus_species <- gsub("\\._3", "", taxa$Genus_species)
taxa$Genus_species[grep("genomosp", taxa$Genus_species)] <- NA
taxa$Genus_species <- gsub("\\._3", "", taxa$Genus_species)
taxa$Genus_species <- gsub("_sp$", "", taxa$Genus_species)
taxa$Genus_species[grep("phytoplasma", taxa$Genus_species)] <- NA
taxa$Genus_species[taxa$Genus_species == "Cyanobacterium_ucyn_a"] <- "Atelocyanobacterium_thalassa"




############
### save ###
############

# OTU abundances
counts <- counts[! rownames(counts) %in% names(which(rowSums(counts) == 0)),]
counts <- data.frame(Sample_ID = rownames(counts), counts)
write.csv(counts, "PipeType_OTU_counts.csv", row.names = F, na = "", quote = F)


# OTU taxonomy
write.csv(taxa, "PipeType_OTU_taxonomy.csv", row.names = F, na = "", quote = F)

