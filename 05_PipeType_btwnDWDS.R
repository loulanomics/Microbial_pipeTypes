###########################################
### impact of DWDS on microbial communities 
### btwn DWDS + within material
### Lou LaMartina, finalized Jan 24 2024
###########################################

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


# subset ductile iron
info <- subset(info, Pipe_material == "ductileiron")
counts <- subset(counts, rownames(counts) %in% info$Sample_ID)
counts <- counts[colSums(counts) > 0]
taxa <- subset(taxa, OTU %in% colnames(counts))
table(info$Site_city)
# Milwaukee Oak Creek  Waukesha 
# 10         3         1 


# simple DWDS variable
info$DWDS <- "MKE"
info$DWDS[info$Site_city == "Waukesha"] <- "WK"
info$DWDS[info$Site_city == "Oak Creek"] <- "OC"




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




########################
### indicator genera ###
########################

########
### DWDS

# stat
identical(rownames(relabun_genus), info$Sample_ID)
indic_dwds.stat <- multipatt(relabun_genus, as.vector(info$DWDS), 
                             control = how(nperm = 999))


# organize
indic_dwds <- data.frame(indic_dwds.stat$sign)
indic_dwds <- data.frame(Genus = rownames(indic_dwds), indic_dwds)
colnames(indic_dwds) <- gsub("^s\\.", "", colnames(indic_dwds))


# can have more than one group
indic_dwds$Variable <- NA
for (i in 1:nrow(indic_dwds)) {
  indic_dwds$Variable[i] <- paste(colnames(indic_dwds)[which(indic_dwds[i,1:4] == 1)], collapse = "_")
}
indic_dwds$Indicator <- paste(indic_dwds$Variable, indic_dwds$Genus)
indic_dwds <- data.frame(Test = "DWDS", indic_dwds[c("Genus", "Indicator", "Variable", "stat", "p.value")])


# A stat
indic_dwdsA <- data.frame(indic_dwds.stat$A)
indic_dwdsA <- data.frame(Genus = rownames(indic_dwdsA), indic_dwdsA)
indic_dwdsA <- melt(indic_dwdsA, value.name = "A", variable.name = "Variable")
indic_dwdsA$Variable <- gsub("\\.", "_", indic_dwdsA$Variable)
indic_dwdsA$Indicator <- paste(indic_dwdsA$Variable, indic_dwdsA$Genus)
indic_dwdsA <- subset(indic_dwdsA, Indicator %in% indic_dwds$Indicator)


# B stat
indic_dwdsB <- data.frame(indic_dwds.stat$B)
indic_dwdsB <- data.frame(Genus = rownames(indic_dwdsB), indic_dwdsB)
indic_dwdsB <- melt(indic_dwdsB, value.name = "B", variable.name = "Variable")
indic_dwdsB$Variable <- gsub("\\.", "_", indic_dwdsB$Variable)
indic_dwdsB$Indicator <- paste(indic_dwdsB$Variable, indic_dwdsB$Genus)
indic_dwdsB <- subset(indic_dwdsB, Indicator %in% indic_dwds$Indicator)


# combine
indic_dwds_all <- merge(merge(indic_dwds, indic_dwdsA[c("Indicator", "A")], by = "Indicator"),
                   indic_dwdsB[c("Indicator", "B")], by = "Indicator")
rm(indic_dwdsA, indic_dwdsB)
indic_dwds <- subset(indic_dwds_all, p.value < 0.05)
table(indic_dwds$Variable)
# MKE    OC OC_WK    WK 
# 1     4    14     8 



###########
### surface

# stat
identical(rownames(relabun_genus), info$Sample_ID)
indic_surface.stat <- multipatt(relabun_genus, as.vector(info$Sample_surface), 
                                max.order = 1, control = how(nperm = 999))


# organize
indic_surface <- data.frame(indic_surface.stat$sign)
indic_surface <- data.frame(Genus = rownames(indic_surface), indic_surface)
colnames(indic_surface) <- gsub("^s\\.", "", colnames(indic_surface))
indic_surface$Variable <- colnames(indic_surface[2:3])[apply(indic_surface[2:3], 1, which.max)]
indic_surface$Indicator <- paste(indic_surface$Variable, indic_surface$Genus)
indic_surface <- data.frame(Test = "Sample_surface", indic_surface[c("Genus", "Indicator", "Variable", "stat", "p.value")])


# A stat
indic_surfaceA <- data.frame(indic_surface.stat$A)
indic_surfaceA <- data.frame(Genus = rownames(indic_surfaceA), indic_surfaceA)
indic_surfaceA <- melt(indic_surfaceA, value.name = "A", variable.name = "Variable")
indic_surfaceA$Indicator <- paste(indic_surfaceA$Variable, indic_surfaceA$Genus)
indic_surfaceA <- subset(indic_surfaceA, Indicator %in% indic_surface$Indicator)


# B stat
indic_surfaceB <- data.frame(indic_surface.stat$B)
indic_surfaceB <- data.frame(Genus = rownames(indic_surfaceB), indic_surfaceB)
indic_surfaceB <- melt(indic_surfaceB, value.name = "B", variable.name = "Variable")
indic_surfaceB$Indicator <- paste(indic_surfaceB$Variable, indic_surfaceB$Genus)
indic_surfaceB <- subset(indic_surfaceB, Indicator %in% indic_surface$Indicator)


# combine
indic_surface_all <- merge(merge(indic_surface, indic_surfaceA[c("Indicator", "A")], by = "Indicator"),
                           indic_surfaceB[c("Indicator", "B")], by = "Indicator")
rm(indic_surfaceA, indic_surfaceB)
indic_surface <- subset(indic_surface_all, p.value < 0.05 & Genus %in% indic_dwds$Genus)


# save
write.csv(rbind(indic_dwds, indic_surface), "DATA/DWDS_indicators.csv", row.names = F, na = "", quote = F)




###############
### heatmap ###
###############

# summed by site
heat_data <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus)
heat_data <- merge(info[c("Sample_ID", "Site_name")], heat_data, by = "Sample_ID")
heat_data <- aggregate(. ~ Site_name + Site_name, sum, data = heat_data[-1])
rownames(heat_data) <- heat_data$Site_name


# proportions of genera
heat_data <- heat_data[indic_dwds$Genus]
heat_data <- data.frame(t(heat_data))
heat_data <- data.frame(Genus = rownames(heat_data), heat_data / rowSums(heat_data))
rowSums(heat_data[-1])
heat_data <- melt(heat_data, variable.name = "Site_name", value.name = "Prop")


# add data
heat_data <- unique(merge(info[c("DWDS", "Site_name")], heat_data, by = "Site_name"))
heat_data <- merge(indic_dwds[c("Genus", "Variable", "Indicator")], heat_data, by = "Genus")


# order
indic_dwds <- indic_dwds[order(indic_dwds$Indicator),]
heat_data$Indicator <- factor(heat_data$Indicator, levels = indic_dwds$Indicator)
heat_data$Genus <- factor(heat_data$Genus, levels = indic_dwds$Genus)
heat_data$Site_name[heat_data$Site_name == "OakCreek3"] <- "Oak Creek"


# plot
heat.plot <-
  ggplot(heat_data, aes(x = Site_name, y = Genus)) +
  geom_tile(aes(fill = Prop)) +
  facet_grid(Variable ~ DWDS, space = "free", scales = "free",
             labeller = labeller(Variable = c(
               "MKE" = "M",
               "OC" = "Oak Creek",
               "OC_WK" = "Oak Creek + Waukesha",
               "WK" = "Waukesha"), DWDS = c(
                 "MKE" = "Milwaukee",
                 "OC" = "Oak Creek",
                 "WK" = "W"))) +
  scale_fill_gradientn(values = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1),
                       colors = brewer.pal(5, "YlGnBu")) +
  theme_light() +
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 7, face = "italic"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "grey70"),
        strip.text = element_text(color = "black", size = 10, face = "bold")) +
  labs(x = "Site", y = "Indicator genus", fill = "Proportion")
heat.plot
ggsave("PLOTS/dwds_heatmap.pdf", plot = heat.plot, device = "pdf", width = 5.5, height = 6.5, unit = "in")




################
### barplots ###
################

# sum by sample surface
bar_data <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus[unique(indic_dwds$Genus)])
bar_data <- merge(info[c("Sample_ID", "Sample_surface")], bar_data, by = "Sample_ID")
bar_data <- aggregate(. ~ Sample_surface, sum, data = bar_data[-1])
rownames(bar_data) <- bar_data$Sample_surface


# proportions of genera
bar_data <- data.frame(t(bar_data[-1]))
bar_data <- data.frame(Genus = rownames(bar_data), bar_data / rowSums(bar_data))
bar_data <- melt(bar_data, variable.name = "Sample_surface", value.name = "Prop")


# add indic results
bar_data <- merge(indic_surface[c("Genus", "Variable", "Indicator")], bar_data, by = "Genus")


# plot
bar.plot <-
  ggplot(bar_data, aes(y = Genus, x = Prop, fill = Sample_surface)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Variable ~ ., scales = "free", space = "free") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = brewer.pal(9, "YlGnBu")[c(3,5,7)]) +
  theme_light() +
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 9, face = "italic"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "grey70"),
        strip.text = element_text(color = "black", size = 10, face = "bold")) +
  labs(x = "Proportion", y = "Indicator genus", fill = "Surface")
bar.plot
ggsave("PLOTS/dwds_barplot.pdf", plot = bar.plot, device = "pdf", width = 6, height = 2.5, unit = "in")


