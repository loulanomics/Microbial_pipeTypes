#########################################
### indicator species of pipe material in
### disparate DWDSs
### Lou LaMartina, finalized Jan 31 2024
#########################################

setwd("~/Desktop/Postdoc/16S/PipeType")

library(readxl)
library(vegan)
library(reshape2)
library(ggplot2)
library(indicspecies)
library(RColorBrewer)
library(stringr)



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


# site labels
info$Site_label <- info$Site_city
info$Site_label[info$Site_label == "Milwaukee"] <- "(Mil.)"
info$Site_label[info$Site_label == "Oak Creek"] <- "(Oak.)"
info$Site_label[info$Site_label == "Waukesha"] <- "(Wau.)"
info$Site_label <- paste(info$Site_label, info$Site_name)
info$Site_label <- gsub("OakCreek", "Site ", info$Site_label)
(info$Site_label <- gsub("N51", "N. 51", info$Site_label))




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
# 54           9          39


# save
write.csv(indic_type, "DATA/PipeType_indicators.csv", row.names = F, na = "", quote = F)




###############
### heatmap ###
###############

# summed by site
heat_data <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus)
heat_data <- merge(info[c("Sample_ID", "Pipe_ID")], heat_data, by = "Sample_ID")
heat_data <- aggregate(. ~ Pipe_ID, sum, data = heat_data[-1])
rownames(heat_data) <- heat_data$Pipe_ID


# proportions of most significant genera
heat_data <- heat_data[indic_type$Genus[indic_type$A > 0.5]]
heat_data <- data.frame(t(heat_data))
heat_data <- data.frame(Genus = rownames(heat_data), heat_data / rowSums(heat_data))
rowSums(heat_data[-1])
heat_data <- melt(heat_data, variable.name = "Pipe_ID", value.name = "Prop")


# add data
heat_data <- merge(info[c("Pipe_ID", "Site_name", "Site_label", "Pipe_material")], heat_data, by = "Pipe_ID")
heat_data <- merge(indic_type[c("Genus", "Variable", "Indicator")], heat_data, by = "Genus")




################
### barplot ###
################

# sum by sample surface
bar_data <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus[sort(unique(as.character(heat_data$Genus)))])
bar_data <- merge(info[c("Sample_ID", "Sample_surface")], bar_data, by = "Sample_ID")
bar_data <- aggregate(. ~ Sample_surface, sum, data = bar_data[-1])
rownames(bar_data) <- bar_data$Sample_surface


# proportions of genera
bar_data <- data.frame(t(bar_data[-1]))
bar_data <- data.frame(Genus = rownames(bar_data), bar_data / rowSums(bar_data))
bar_data <- melt(bar_data, variable.name = "Sample_surface", value.name = "Prop")


# add indic results
bar_data <- merge(bar_data, indic_type[c("Genus", "Variable")], by = "Genus", all.x = T)


# order by abundance in surface
temp1 <- subset(bar_data, Variable == "castiron" & Sample_surface == "smooth")
temp1 <- temp1[order(temp1$Prop, decreasing = T),]
temp2 <- subset(bar_data, Variable == "copper" & Sample_surface == "smooth")
temp2 <- temp2[order(temp2$Prop, decreasing = T),]
temp3 <- subset(bar_data, Variable == "ductileiron" & Sample_surface == "smooth")
temp3 <- temp3[order(temp3$Prop, decreasing = T),]
bar_data$Genus <- factor(bar_data$Genus, levels = c(temp1$Genus, temp2$Genus, temp3$Genus))
heat_data$Genus <- factor(heat_data$Genus, levels = c(temp1$Genus, temp2$Genus, temp3$Genus))


# fix
bar_data$Sample_surface <- str_to_title(bar_data$Sample_surface)
bar_data$Sample_surface <- gsub("Under", "Under ", bar_data$Sample_surface)




############
### plot ###
############

# barplot
bar.plot <-
  ggplot(bar_data, aes(y = Genus, x = Prop, fill = Sample_surface)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Variable ~ ., scales = "free", space = "free",
             labeller = labeller(Variable = c(castiron = "Cast iron",
                                                ductileiron = "Ductile iron",
                                                copper = "Copper"))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = brewer.pal(9, "YlGnBu")[c(3,5,7)]) +
  theme_light() +
  theme(axis.text.x = element_text(color = "black", size = 10, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.position = c(0.6, 1.05),
        plot.margin = unit(c(0.6,0.3,0.1,0), "in"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(fill = guide_legend(title.position = "top", direction = "horizontal",
                             label.hjust = 1, title.hjust = 1)) +
  labs(x = "Proportion", y = "Indicator genus", fill = "Surface")
bar.plot
ggsave("PLOTS/pipetype_barplot.pdf", plot = bar.plot, device = "pdf", width = 4, height = 7.7, unit = "in")


# heatmap
heat.plot <-
  ggplot(heat_data, aes(x = Site_label, y = Genus)) +
  geom_tile(aes(fill = Prop)) +
  facet_grid(Variable ~ Pipe_material, space = "free", scales = "free", switch = "y",
             labeller = labeller(
               Pipe_material = c(castiron = "Cast iron",
                                 ductileiron = "Ductile iron",
                                 copper = "Copper",
                                 lead = "Lead"),
               Variable = c(castiron = "Cast iron",
                            ductileiron = "Ductile iron",
                            copper = "Copper"))) +
  scale_fill_gradientn(values = c(0, 0.25, 0.5, 0.75, 1), 
                       colors = brewer.pal(5, "YlGnBu")) +
  scale_y_discrete(position = "right") +
  theme_light() +
  theme(axis.text.x = element_text(color = "black", size = 9, angle = 90, hjust = 1),
        axis.text.y.right = element_text(color = "black", size = 8, face = "italic", hjust = 0.5),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.position = c(0.28, 1.1),
        plot.margin = unit(c(0.8,0.1,0.1,0.2), "in"),
        strip.background = element_rect(fill = "grey90", color = "grey70"),
        strip.text = element_text(color = "black", size = 10, face = "bold")) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.75, title.position = "top", direction = "horizontal",
                               label.hjust = 0)) +
  labs(x = "Site", y = "Indicator genus", fill = "Proportion")
heat.plot
ggsave("PLOTS/pipetype_heatmap.pdf", plot = heat.plot, device = "pdf", width = 7, height = 9, unit = "in")




#########################
### common in samples ###
#########################

# glom
common_smp <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus)
common_smp <- merge(info[c("Sample_ID", "Site_city", "Pipe_material")], common_smp, by = "Sample_ID")
common_smp$Sample <- paste0(common_smp$Site_city, "_", common_smp$Pipe_material)
common_smp <- aggregate(. ~ Sample, sum, data = common_smp[-c(1:3)])
rownames(common_smp) <- common_smp$Sample
common_smp <- common_smp[-1]
common_smp <- common_smp / rowSums(common_smp)


# 10 most abundant in each sample
tops_smp.ls <- list()
for (i in rownames(common_smp)) {
  tops_smp.ls[[i]] <- names(sort(colSums(common_smp[i,]), decreasing = T))[1:8]
}
tops_smp <- unique(unlist(tops_smp.ls))


# compare
smp_data <- data.frame(t(common_smp[tops_smp])) * 100
smp_data$mean <- apply(smp_data[1:9], 1, mean)
smp_data$sd <- apply(smp_data[1:9], 1, sd)
smp_data$CV <- smp_data$sd / smp_data$mean


# save
write.csv(smp_data, "DATA/common_samples.csv", row.names = F, na = "", quote = F)




########################
### common in cities ###
########################

# glom
common_city <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus)
common_city <- merge(info[c("Sample_ID", "Site_city")], common_city, by = "Sample_ID")
common_city <- aggregate(. ~ Site_city, sum, data = common_city[-1])
rownames(common_city) <- common_city$Site_city
common_city <- common_city[-1]
common_city <- common_city / rowSums(common_city)
rowSums(common_city)


# 10 most abundant
tops_city.ls <- list()
for (i in rownames(common_city)) {
  tops_city.ls[[i]] <- names(sort(colSums(common_city[i,]), decreasing = T))[1:10]
}
tops_city <- unique(unlist(tops_city.ls))


# compare
city_data <- data.frame(t(common_city[tops_city])) * 100
city_data$mean <- apply(city_data[1:3], 1, mean)
city_data$sd <- apply(city_data[1:3], 1, sd)
city_data$CV <- city_data$sd / city_data$mean
city_data$MvW <- city_data$Milwaukee / city_data$Waukesha
city_data$MvO <- city_data$Milwaukee / city_data$Oak.Creek


# save
write.csv(city_data, "DATA/common_cities.csv", row.names = F, na = "", quote = F)




###########################
### common in materials ###
###########################

# glom
common_type <- data.frame(Sample_ID = rownames(relabun_genus), relabun_genus)
common_type <- merge(info[c("Sample_ID", "Pipe_material")], common_type, by = "Sample_ID")
common_type <- aggregate(. ~ Pipe_material, sum, data = common_type[-1])
rownames(common_type) <- common_type$Pipe_material
common_type <- common_type[-1]
common_type <- common_type / rowSums(common_type)
rowSums(common_type)


# 10 most abundant
tops_type.ls <- list()
for (i in rownames(common_type)) {
  tops_type.ls[[i]] <- names(sort(colSums(common_type[i,]), decreasing = T))[1:10]
}
tops_type <- unique(unlist(tops_type.ls))


# compare
type_data <- data.frame(t(common_type[tops_type])) * 100
type_data$mean <- apply(type_data[1:4], 1, mean)
type_data$sd <- apply(type_data[1:4], 1, sd)
type_data$CV <- type_data$sd / type_data$mean
type_data$LvCi <- type_data$lead / type_data$castiron
type_data$LvDi <- type_data$lead / type_data$ductileiron
type_data$LvCo <- type_data$lead / type_data$copper
type_data$CivDi <- type_data$castiron / type_data$ductileiron
type_data$CivCo <- type_data$castiron / type_data$copper
type_data$DivCo <- type_data$ductileiron / type_data$copper


# save
write.csv(type_data, "DATA/common_types.csv", row.names = F, na = "", quote = F)




############
### plot ###
############

# for plotting
smp.m <- common_smp
smp.m$Other <- rowSums(smp.m[! colnames(smp.m) %in% tops_smp])
smp.m <- smp.m[colnames(smp.m) %in% c("Other", tops_smp)]
min(smp.m$Other) # 0.306013
smp.m$Other <- smp.m$Other - 0.3


# melt
smp.m <- melt(data.frame(Sample = rownames(smp.m), smp.m), variable.name = "Genus", value.name = "Relabun")
smp.m$City <- sapply(strsplit(smp.m$Sample, "_"), '[', 1)
smp.m$Sample <- sapply(strsplit(smp.m$Sample, "_"), '[', 2)


# colors
colors <- data.frame(Genus = c(names(sort(colSums(relabun_genus[tops_smp]), decreasing = T)), "Other"),
                     Color = c(colorRampPalette(brewer.pal(11, "Spectral"))(length(tops_smp)), "grey90"))
smp.m$Genus <- factor(smp.m$Genus, levels = colors$Genus)
cbind(seq(from = 30, to = 100, length.out = 100),
      scales::rescale(seq(from = 30, to = 100, length.out = 100)))
smp.m$Sample <- str_to_title(smp.m$Sample)
smp.m$Sample <- gsub("iron", "\niron", smp.m$Sample)


# plot
stacked.plot <- 
  ggplot(smp.m, aes(x = Sample, y = Relabun, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ City, scales = "free", space = "free") +
  scale_fill_manual(values = colors$Color) +
  scale_y_continuous(expand = c(0.01, 0.01), breaks = c(0, 0.3, 0.65, 1), labels = c("0-30%\n(Other)", "50%", "75%", "100%")) +
  scale_x_discrete(labels = c("Milwaukee" = "Mil.", "Waukesha" = "Wau.", "Oak Creek" = "Oak")) +
  theme_light() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 9, face = "italic"),
        strip.background = element_rect(fill = "grey90", color = "grey70"),
        strip.text = element_text(color = "black", size = 10, face = "bold")) +
  labs(y = "Percent genus")
stacked.plot
ggsave("PLOTS/stacked.pdf", plot = stacked.plot, device = "pdf", width = 8, height = 5, unit = "in")

