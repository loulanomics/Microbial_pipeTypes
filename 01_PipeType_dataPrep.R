#########################################
### organize QIIME2 output files
### from MR DNA for pipe type study
### Lou LaMartina, finalized Jan 11, 2024
#########################################


setwd("~/Desktop/Postdoc/16S/PipeType/DATA")
library(readxl)
library(stringr)


#################
### load data ###
#################

# QIIME2 results
otus <- data.frame(read_xlsx("MRDNA/032522MDillcus515F-zotus.fa.bacteria.OTU.xlsx"))


# sample info
smps <- data.frame(read_xlsx("MRDNA/00000MRdiversitySampleSubmissionForm.xlsx", sheet = "Lou edits"))


# failed assays
fails <- data.frame(read_xlsx("MRDNA/032522MD negative listb.xlsx"))




#####################
### match samples ###
#####################

# custom
smps$File <- smps$Label


# replace spaces, hyphens, and periods with underscores
(smps$File <- gsub(" |-", "_", smps$File))
(colnames(otus) <- gsub("\\.", "_", colnames(otus)))


# changes
smps$File[smps$File == "MKE4_UT3"] <- "MKE_4_UT3"
smps$File[smps$File == "MKE_5_Tub_1"] <- "MKE_5_Tub1_R1"
smps$File[smps$File == "MKE_5_Tub_2"] <- "MKE_5_Tub2_R1"
smps$File[smps$File == "MKE_5_Tub_3"] <- "MKE_5_Tub3_R1"
smps$File[smps$File == "WK2_DLP_19B"] <- "WK2_LP_19_B"
smps$File[smps$File == "WK2_DLP_20B"] <- "WK2_LP_20_B"
smps$File[smps$File == "WK2_DLP_21_B"] <- "WK2_LP_21_B"


# MKE_5_B1 or MKE 5-T3 are mislabeled. i will change the column name
# of the counts data but remember i did this. considering its
# in the right order when it's changed, i have a feeling this is right
colnames(otus)[colnames(otus) == "MKE_5_B1"] <- "MKE_5_T3"


# same here
colnames(otus)[colnames(otus) == "MKE_7_1"] <- "MKE_7_3"
colnames(otus)[colnames(otus) == "WK2_Fe_12D"] <- "WK2_Fe_2"


# removed failed assays
smps <- subset(smps, ! Label %in% fails$Sample.name)




##################
### OTU counts ###
##################

# extract counts
counts <- otus[c(1,8:99)]
rownames(counts) <- counts$otu_name
counts <- counts[-1]


# change OTU ID
rownames(counts) <- paste0("OTU", sprintf("%04d", as.numeric(gsub("Zotu", "", rownames(counts)))))


# remove failed assays
counts <- counts[smps$File]


# remove empty samples and OTUs
dim(counts) # 4262   74
names(which(colSums(counts) == 0)) # "MKE_7_3"
counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
dim(counts) # 4088   73
smps <- subset(smps, File != "MKE_7_3")




################
### taxonomy ###
################

# OTU IDs and taxonomy assignment strings, separated by semi colons
tax <- otus[1:2]
tax.ls <- strsplit(tax$Taxonomy, ";")


# extract classifications based on starting with k__, p__, etc
splittax.ls <- list()
for (i in paste0("^", c("k", "p", "c", "o", "f", "g", "s"), "__")) {
  for (j in 1:length(tax.ls)) {
    splittax.ls[[i]][[j]] <- unlist(tax.ls[j])[grep(i, unlist(tax.ls[j]))]
    splittax.ls[[i]][[j]] <- gsub(i, "", splittax.ls[[i]][[j]])
  }
}


# combine
splittax <- data.frame(do.call(cbind, splittax.ls))
colnames(splittax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Genus_species")
splittax$Genus_species <- gsub(" ", "_", splittax$Genus_species)
splittax <- apply(splittax, 2, str_to_title)
tax <- cbind(tax, splittax)


# change Zotu to OTU, add leading zeros
tax$otu_name <- paste0("OTU", sprintf("%04d", as.numeric(gsub("Zotu", "", tax$otu_name))))
colnames(tax)[1] <- "OTU"
rownames(tax) <- tax$OTU


# extract accession numbers
tax$Accession <- NA
tax$Accession[grep("\\.1 ", tax$Taxonomy)] <- tax$Taxonomy[grep("\\.1 ", tax$Taxonomy)]
tax$Accession <- sapply(strsplit(tax$Accession, " "), '[', 1)
tax <- tax[-2]


# remove empties
tax <- subset(tax, OTU %in% rownames(counts))




############
### save ###
############

# note: added A/B to samples from different pipes but collected same day
# for example, MKE_dFeB_smth_21Aug18, MKE_dFeA_smth_21Aug18.
# A goes to the older pipe


# OTU abundances
counts <- data.frame(File = colnames(counts), t(counts))
counts <- merge(smps[c("File", "ID")], counts, by = "File")[-1]
write.csv(counts, "PipeType_OTU_counts.csv", row.names = F, na = "", quote = F)


# OTU taxonomy
write.csv(tax, "PipeType_OTU_taxonomy.csv", row.names = F, na = "", quote = F)


# sample metadata was organized in excel

