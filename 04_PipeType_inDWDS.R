#########################################
### impact of pipe material on 
### microbial communities 
### btwn materials + within DWDS
### Lou LaMartina, finalized Jan 23 2024
#########################################

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


# subset Milwaukee
info <- subset(info, Site_city == "Milwaukee")
counts <- subset(counts, rownames(counts) %in% info$Sample_ID)
counts <- counts[colSums(counts) > 0]
taxa <- subset(taxa, OTU %in% colnames(counts))




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


# only this with max rel abundance > 0.1%
relabun_genus_top <- relabun_genus[which(apply(relabun_genus, 2, max) > 0.001)]




########################
### indicator genera ###
########################

##############
### pipe types

# stat
identical(rownames(relabun_genus_top), info$Sample_ID)
indic_type.stat <- multipatt(relabun_genus_top, as.vector(info$Pipe_material), 
                             max.order = 1, control = how(nperm = 999))


# organize
indic_type <- data.frame(indic_type.stat$sign)
indic_type <- data.frame(Genus = rownames(indic_type), indic_type)
colnames(indic_type) <- gsub("^s\\.", "", colnames(indic_type))
indic_type$Variable <- colnames(indic_type[2:4])[apply(indic_type[2:4], 1, which.max)]
indic_type$Indicator <- paste(indic_type$Variable, indic_type$Genus)
indic_type <- data.frame(Test = "Pipe_material", indic_type[c("Genus", "Indicator", "Variable", "stat", "p.value")])


# A stat
indic_typeA <- data.frame(indic_type.stat$A)
indic_typeA <- data.frame(Genus = rownames(indic_typeA), indic_typeA)
indic_typeA <- melt(indic_typeA, value.name = "A", variable.name = "Variable")
indic_typeA$Indicator <- paste(indic_typeA$Variable, indic_typeA$Genus)
indic_typeA <- subset(indic_typeA, Indicator %in% indic_type$Indicator)


# B stat
indic_typeB <- data.frame(indic_type.stat$B)
indic_typeB <- data.frame(Genus = rownames(indic_typeB), indic_typeB)
indic_typeB <- melt(indic_typeB, value.name = "B", variable.name = "Variable")
indic_typeB$Indicator <- paste(indic_typeB$Variable, indic_typeB$Genus)
indic_typeB <- subset(indic_typeB, Indicator %in% indic_type$Indicator)


# combine
indic_type_all <- merge(merge(indic_type, indic_typeA[c("Indicator", "A")], by = "Indicator"),
                   indic_typeB[c("Indicator", "B")], by = "Indicator")
rm(indic_typeA, indic_typeB)
indic_type <- subset(indic_type_all, p.value < 0.05)
table(indic_type$Variable)
# castiron      copper ductileiron 
# 46           7          19



###########
### surface

# stat
identical(rownames(relabun_genus_top), info$Sample_ID)
indic_surface.stat <- multipatt(relabun_genus_top, as.vector(info$Sample_surface), 
                             max.order = 1, control = how(nperm = 999))


# organize
indic_surface <- data.frame(indic_surface.stat$sign)
indic_surface <- data.frame(Genus = rownames(indic_surface), indic_surface)
colnames(indic_surface) <- gsub("^s\\.", "", colnames(indic_surface))
indic_surface$Variable <- colnames(indic_surface[2:4])[apply(indic_surface[2:4], 1, which.max)]
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
indic_surface <- subset(indic_surface_all, p.value < 0.05 & Genus %in% indic_type$Genus)


# both
indic_both <- subset(indic_type, Genus %in% indic_surface$Genus)
indic_both <- rbind(indic_both, indic_surface)
indic_both <- subset(indic_both, Genus != "Alistipes") # cast iron only


# save
write.csv(rbind(indic_type, indic_surface), "DATA/PipeType_indicators.csv", row.names = F, na = "", quote = F)




###############
### heatmap ###
###############

# remove more
removes <- subset(indic_type, Variable == "castiron" & p.value > 0.01)$Genus


# summed by site
heat_data <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus)
heat_data <- merge(info[c("Sample_ID", "Pipe_ID")], heat_data, by = "Sample_ID")
heat_data <- aggregate(. ~ Pipe_ID, sum, data = heat_data[-1])
rownames(heat_data) <- heat_data$Pipe_ID


# proportions of genera
heat_data <- heat_data[indic_type$Genus[! indic_type$Genus %in% removes]]
heat_data <- data.frame(t(heat_data))
heat_data <- data.frame(Genus = rownames(heat_data), heat_data / rowSums(heat_data))
rowSums(heat_data[-1])
heat_data <- melt(heat_data, variable.name = "Pipe_ID", value.name = "Prop")


# add data
heat_data <- merge(info[c("Pipe_ID", "Site_name", "Pipe_material")], heat_data, by = "Pipe_ID")
heat_data <- merge(indic_type[c("Genus", "Variable", "Indicator")], heat_data, by = "Genus")


# order
indic_type <- indic_type[order(indic_type$Indicator),]
heat_data$Indicator <- factor(heat_data$Indicator, levels = indic_type$Indicator)
heat_data$Genus <- factor(heat_data$Genus, levels = indic_type$Genus)


# plot
heat.plot <-
  ggplot(heat_data, aes(x = Site_name, y = Genus)) +
  geom_tile(aes(fill = Prop)) +
  facet_grid(Variable ~ Pipe_material, space = "free", scales = "free",
             labeller = labeller(
               Pipe_material = c(castiron = "Cast iron pipes",
                                 ductileiron = "Ductile iron pipes",
                                 copper = "Copper pipes",
                                 lead = "Lead pipes"),
               Variable = c(castiron = "Cast iron",
                            ductileiron = "Ductile iron",
                            copper = "Copper"))) +
  scale_fill_gradientn(values = c(0, 0.25, 0.5, 0.75, 1), 
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
ggsave("PLOTS/pipetype_heatmap.pdf", plot = heat.plot, device = "pdf", width = 8, height = 8, unit = "in")




################
### barplots ###
################

# sum by sample surface
bar_data <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus[unique(indic_both$Genus)])
bar_data <- merge(info[c("Sample_ID", "Sample_surface")], bar_data, by = "Sample_ID")
bar_data <- aggregate(. ~ Sample_surface, sum, data = bar_data[-1])
rownames(bar_data) <- bar_data$Sample_surface


# proportions of genera
bar_data <- data.frame(t(bar_data[-1]))
bar_data <- data.frame(Genus = rownames(bar_data), bar_data / rowSums(bar_data))
bar_data <- melt(bar_data, variable.name = "Sample_surface", value.name = "Prop")


# add indic results
bar_data <- merge(indic_surface[c("Genus", "Variable", "Indicator")], bar_data, by = "Genus")


# order by abundance in surface
temp1 <- subset(bar_data, Variable == "smooth" & Sample_surface == "smooth")
temp1 <- temp1[order(temp1$Prop, decreasing = T),]
temp2 <- subset(bar_data, Variable == "tubercle" & Sample_surface == "tubercle")
temp2 <- temp2[order(temp2$Prop, decreasing = T),]
temp3 <- subset(bar_data, Variable == "undertubercle" & Sample_surface == "undertubercle")
temp3 <- temp3[order(temp3$Prop, decreasing = T),]
bar_data$Genus <- factor(bar_data$Genus, levels = c(temp1$Genus, temp2$Genus, temp3$Genus))


# fix
bar_data$Variable <- str_to_title(bar_data$Variable)
bar_data$Variable <- gsub("Under", "Under ", bar_data$Variable)


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
ggsave("PLOTS/pipetype_barplot.pdf", plot = bar.plot, device = "pdf", width = 6, height = 6, unit = "in")


