########################################
### lead pipe microbial communities 
### Lou LaMartina, finalized Jan 24 2024
########################################

setwd("~/Desktop/Postdoc/16S/PipeType")

library(readxl)
library(vegan)
library(reshape2)
library(ggplot2)
library(indicspecies)
library(RColorBrewer)
library(OTUtable)



#################
### load data ###
#################

# see 02_PipeType_QAQC.R

# OTU counts
counts <- read.csv("DATA/PipeType_OTU_counts.csv")
rownames(counts) <- counts$Sample_ID
counts <- counts[-1]


# taxonomy
taxa <- read.csv("DATA/PipeType_OTU_taxonomy.csv")
rownames(taxa) <- taxa$OTU


# sample metadata
info <- data.frame(read_xlsx("DATA/PipeType_sample_metadata.xlsx"))
rownames(info) <- info$Sample_ID
info <- info[match(rownames(counts), info$Sample_ID),]


# subset lead
info <- subset(info, Pipe_material == "lead")
counts <- subset(counts, rownames(counts) %in% info$Sample_ID)
counts <- counts[colSums(counts) > 0]
taxa <- subset(taxa, OTU %in% colnames(counts))
table(info$Site_city)
# Milwaukee  Waukesha 
# 7         3 




############
### glom ###
############

# counts to genus
counts_genus <- data.frame(OTU = colnames(counts), t(counts))
counts_genus <- merge(taxa[c("OTU", "Genus")], counts_genus, by = "OTU")
counts_genus <- aggregate(. ~ Genus, sum, data = counts_genus[-1])
rownames(counts_genus) <- counts_genus$Genus
counts_genus <- data.frame(t(counts_genus[-1]))


# relative abundance
relabun_genus <- counts_genus / rowSums(counts_genus)


# 10 most abundant in each sample
tops <- list()
for (i in rownames(relabun_genus)) {
  tops[[i]] <- names(sort(colSums(relabun_genus[i,]), decreasing = T))[1:10]
}
tops <- unique(unlist(tops))
sum(relabun_genus[tops]) / sum(relabun_genus) # 0.7649632




###############
### barplot ###
###############

# all else is "other"
relabun_genus_top <- relabun_genus
relabun_genus_top$Other <- rowSums(relabun_genus[! colnames(relabun_genus) %in% tops])
relabun_genus_top <- relabun_genus_top[colnames(relabun_genus_top) %in% c(tops, "Other")]
rowSums(relabun_genus_top)


# melt
relabun_genus_top <- melt(data.frame(Sample_ID = rownames(relabun_genus_top), relabun_genus_top), 
                          variable.name = "Genus", value.name = "Relabun")


# sum by site
relabun_genus_top <- merge(info[c("Sample_ID", "Site_name", "Site_city")], relabun_genus_top, by = "Sample_ID")
relabun_genus_top <- aggregate(Relabun ~ Genus + Site_name + Site_city, sum, data = relabun_genus_top[-1])


# sort by abundance
genus_order <- melt(data.frame(Genus = tops, t(relabun_genus[tops])))
genus_order <- aggregate(value ~ Genus, sum, data = genus_order[-2])
genus_order <- genus_order[order(genus_order$value, decreasing = T),]
relabun_genus_top$Genus <- factor(relabun_genus_top$Genus, levels = c(genus_order$Genus, "Other"))


# colors
genus_order$Colors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(tops))


# plot
bar.plot <-
  ggplot(relabun_genus_top, aes(x = Site_name, y = Relabun, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(genus_order$Colors, "grey90")) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  facet_grid(~ Site_city, scales = "free", space = "free") +
  theme_light() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(size = 9, color = "black", face = "italic"),
        strip.background = element_rect(fill = "grey90", color = "grey70"),
        strip.text = element_text(color = "black", size = 10, face = "bold")) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Site", y = "Relative abundance")
ggsave("PLOTS/lead_barplot.pdf", plot = bar.plot, device = "pdf", width = 6, height = 6, unit = "in", bg = "white")


