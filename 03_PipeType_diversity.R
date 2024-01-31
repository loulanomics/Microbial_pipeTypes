########################################
### microbial communities in DWDS
### Lou LaMartina, finalized Jan 24 2024
########################################

setwd("~/Desktop/Postdoc/16S/PipeType")

library(readxl)
library(vegan)
library(ape)
library(reshape2)
library(ggplot2)
library(RColorBrewer)



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


# relative abundance
relabun <- counts / rowSums(counts)




##################
### ordination ###
##################

# calculate distances
bray.stat <- vegdist(relabun, method = "bray")


# do ordination
pcoa.stat <- pcoa(bray.stat)


# extract values (first 2 eigenvectors = axis 1 & 2)
pcoa_data <- data.frame(pcoa.stat$vectors[,1:2])
pcoa_data$Sample_ID <- rownames(pcoa_data)
pcoa_data <- unique(merge(pcoa_data, info, by = "Sample_ID"))
round(pcoa.stat$values$Relative_eig[1:2] * 100, digits = 1)
# 24.8 21.2


# plot
pcoa.plot <-
  ggplot(pcoa_data, aes(x = Axis.1, y = Axis.2, color = Pipe_material, shape = Site_city)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = brewer.pal(11, "Spectral")[c(9,4,11,2)],
                     labels = c("cast iron", "copper", "ductile iron", "lead")) +
  theme_classic() +
  coord_fixed() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        legend.title =  element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 1)) +
  labs(x = "Axis 1 (24.8%)", y = "Axis 2 (21.2%)", shape = "Location", color = "Pipe\nmaterial")
pcoa.plot
ggsave("PLOTS/pcoa.pdf", plot = pcoa.plot, device = "pdf", width = 9, height = 6, units = "in")


# cluster 1
mini1 <- subset(pcoa_data, Axis.1 > 0.05 & Axis.2 < -0.21)
mini1.plot <-
  ggplot(mini1, aes(x = Axis.1, y = Axis.2, color = Pipe_material, shape = Site_city)) +
  geom_point(size = 3) +
  scale_color_manual(values = brewer.pal(11, "Spectral")[c(9,4,11,2)],
                     labels = c("cast iron", "copper", "ductile iron", "lead")) +
  scale_y_continuous(limits = c(-0.234, -0.217), breaks = c(-0.22, -0.23), position = "right") +
  theme_classic() +
  coord_fixed() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(color = "grey80", linewidth = 0.25),
        axis.line.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(color = "grey80"),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 1))
ggsave("PLOTS/mini1.pdf", plot = mini1.plot, device = "pdf", width = 4, height = 2, units = "in")


# cluster 2
mini2 <- subset(pcoa_data, Axis.1 > 0.2 & Axis.2 > 0.2)
mini2.plot <-
  ggplot(mini2, aes(x = Axis.1, y = Axis.2, color = Pipe_material, shape = Site_city)) +
  geom_point(size = 3) +
  scale_color_manual(values = brewer.pal(11, "Spectral")[c(9,4,2)],
                     labels = c("cast iron", "copper", "ductile iron", "lead")) +
  scale_x_continuous(limits = c(0.264, 0.282), breaks = c(0.27, 0.28), position = "top") +
  scale_y_continuous(limits = c(0.27, 0.288), breaks = c(0.27, 0.28), position = "right") +
  theme_classic() +
  coord_fixed() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "grey80", linewidth = 0.25),
        axis.ticks = element_line(color = "grey80", linewidth = 1),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 1))
ggsave("PLOTS/mini2.pdf", plot = mini2.plot, device = "pdf", width = 2, height = 2, units = "in")


# ranges
signif(range(subset(pcoa_data, Site_city == "Milwaukee")$Axis.1), 2) # -0.56  0.28
signif(range(subset(pcoa_data, Site_city == "Milwaukee")$Axis.2), 2) # -0.23  0.29

signif(range(subset(pcoa_data, Site_city == "Waukesha")$Axis.1), 2) # 0.13 0.21
signif(range(subset(pcoa_data, Site_city == "Waukesha")$Axis.2), 2) # 0.13 0.22

signif(range(subset(pcoa_data, Site_city == "Oak Creek" & Sample_ID != "OC_1_HS_H")$Axis.1), 2) # 0.066 0.096
signif(range(subset(pcoa_data, Site_city == "Oak Creek" & Sample_ID != "OC_1_HS_H")$Axis.2), 2) # -0.23 -0.18
signif(subset(pcoa_data, Sample_ID == "OC_1_HS_H")[c("Axis.1", "Axis.2")], 2) # -0.18   0.27




#################
### distances ###
#################

# extract distances
dist_data <- data.frame(Sample_ID_x = rownames(relabun), scores(bray.stat))


# get comparisons
dist_data <- melt(dist_data, variable.name = "Sample_ID_y", value.name = "Bray", id.vars = "Sample_ID_x")
dist_data$Sample_ID_y <- as.character(dist_data$Sample_ID_y)
dist_data$Compare <- NA
for (i in 1:nrow(dist_data)) {
  dist_data$Compare[i] <- paste(sort(unlist(dist_data[i, c("Sample_ID_x", "Sample_ID_y")])), collapse = "__")
}
dist_data <- unique(dist_data[dist_data$Bray > 0, c("Compare", "Bray")])


# add info - x
dist_data$Sample_ID <- sapply(strsplit(dist_data$Compare, "__"), '[', 1)
dist_data <- merge(dist_data, info[c("Sample_ID", "Sample_name", "Site_city", "Pipe_material", "Pipe_ID", "Sample_surface")], by = "Sample_ID")
colnames(dist_data)[colnames(dist_data) %in% colnames(info)] <- 
  paste0(colnames(dist_data)[colnames(dist_data) %in% colnames(info)], "_x")


# add info - y
dist_data$Sample_ID <- sapply(strsplit(dist_data$Compare, "__"), '[', 2)
dist_data <- merge(dist_data, info[c("Sample_ID", "Sample_name", "Site_city", "Pipe_material", "Pipe_ID", "Sample_surface")], by = "Sample_ID")
colnames(dist_data)[colnames(dist_data) %in% colnames(info)] <- 
  paste0(colnames(dist_data)[colnames(dist_data) %in% colnames(info)], "_y")
write.csv(dist_data, "DATA/BrayCurtis_distances.csv", row.names = F, na = "", quote = F)




########################
### median abs diffs ###
########################

# within city
dist_city <- subset(dist_data, Site_city_x == Site_city_y)
dist_city$Within <- "1 city"
nrow(dist_city) # 1296


# within surface
dist_surface <- subset(dist_data, Sample_surface_x == Sample_surface_y)
dist_surface$Within <- "2 surface"
nrow(dist_surface) # 1103


# within type
dist_type <- subset(dist_data, Pipe_material_x == Pipe_material_y)
dist_type$Within <- "3 type"
nrow(dist_type) # 811


# type + city
dist_city_type <- subset(dist_data, Site_city_x == Site_city_y & Pipe_material_x == Pipe_material_y)
dist_city_type$Within <- "4 city+type"
nrow(dist_city_type) # 451
dist_city_type$Var <- paste(dist_city_type$Site_city_x, dist_city_type$Pipe_material_x)


# type + surface
dist_surface_type <- subset(dist_data, Sample_surface_x == Sample_surface_y & Pipe_material_x == Pipe_material_y)
dist_surface_type$Within <- "5 surface+type"
nrow(dist_surface_type) # 355
dist_surface_type$Var <- paste(dist_surface_type$Sample_surface_x, dist_surface_type$Pipe_material_x)


# within pipe
dist_pipe <- subset(dist_data, Pipe_ID_x == Pipe_ID_y)
dist_pipe$Within <- "6 pipe"
nrow(dist_pipe) # 153


# within sample
dist_sample <- subset(dist_data, Sample_name_x == Sample_name_y)
dist_sample$Within <- "7 sample"
nrow(dist_sample) # 53


# median absolute differences
mads <- rbind(
  data.frame(Test = "1 City",
             Variable = aggregate(Bray ~ Site_city_x, mad, data = dist_city[c("Bray", "Site_city_x")])[,1],
             MAD = aggregate(Bray ~ Site_city_x, mad, data = dist_city[c("Bray", "Site_city_x")])[,2],
             n = data.frame(table(info$Site_city))$Freq),
  data.frame(Test = "2 Surface",
             Variable = aggregate(Bray ~ Sample_surface_x, mad, data = dist_surface[c("Bray", "Sample_surface_x")])[,1],
             MAD = aggregate(Bray ~ Sample_surface_x, mad, data = dist_surface[c("Bray", "Sample_surface_x")])[,2],
             n = data.frame(table(info$Sample_surface))$Freq),
  data.frame(Test = "3 Type",
           Variable = aggregate(Bray ~ Pipe_material_x, mad, data = dist_type[c("Bray", "Pipe_material_x")])[,1],
             MAD = aggregate(Bray ~ Pipe_material_x, mad, data = dist_type[c("Bray", "Pipe_material_x")])[,2],
             n = data.frame(table(info$Pipe_material))$Freq),
  data.frame(Test = "4 City+Type",
             Variable = aggregate(Bray ~ Var, mad, data = dist_city_type[c("Bray", "Var")])[,1],
             MAD = aggregate(Bray ~ Var, mad, data = dist_city_type[c("Bray", "Var")])[,2],
             n = subset(data.frame(table(paste(info$Site_city, info$Pipe_material))), Freq > 1)$Freq),
  data.frame(Test = "5 Surface+Type",
             Variable = aggregate(Bray ~ Var, mad, data = dist_surface_type[c("Bray", "Var")])[,1],
             MAD = aggregate(Bray ~ Var, mad, data = dist_surface_type[c("Bray", "Var")])[,2],
             n = subset(data.frame(table(paste(info$Sample_surface, info$Pipe_material))), Freq > 1)$Freq))
write.csv(mads, "DATA/MAD_scores.csv", row.names = F, na = "", quote = F)


# significantly different variances?
signif(leveneTest(Bray ~ as.factor(Site_city_x), data = dist_city)$`Pr(>F)`[1], 3)         # 3.44e-32
signif(leveneTest(Bray ~ as.factor(Sample_surface_x), data = dist_surface)$`Pr(>F)`[1], 3) # 0.0212
signif(leveneTest(Bray ~ as.factor(Pipe_material_x), data = dist_type)$`Pr(>F)`[1], 3)     # 8.94e-17
signif(leveneTest(Bray ~ as.factor(Var), data = dist_city_type)$`Pr(>F)`[1], 3)            # 1.48e-17
signif(leveneTest(Bray ~ as.factor(Var), data = dist_surface_type)$`Pr(>F)`[1], 3).        # 8.37e-11


# plot
mads.plot <-
  ggplot(mads, aes(x = Variable, y = MAD, fill = Test)) +
  geom_col() +
  geom_text(data = mads, aes(x = Variable, y = 0, label = n), nudge_y = -0.01, size = 2) +
  facet_grid(~ Test, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("PLOTS/mads.pdf", plot = mads.plot, device = "pdf", width = 9, height = 6, units = "in")




#################
### PERMANOVA ###
#################

# permanova
identical(rownames(relabun), info$Sample_ID)
perman <- adonis2(bray.stat ~ Site_city + Pipe_material + Sample_surface,
                  data = info, na.action = na.exclude)


# extract results
perman <- data.frame(Variable = rownames(perman), perman)
colnames(perman)[6] <- "p"
perman$R2 <- signif(perman$R2, 3)


# post-hoc pairwise multilevel comparison
posthoc <- pairwise.adonis2(relabun ~ Site_city + Pipe_material + Sample_surface, 
                            data = info)
posthoc <- do.call(rbind, posthoc[-1])
posthoc <- data.frame(Comparison = sapply(strsplit(rownames(posthoc), "\\."), '[', 1),
                      Variable = sapply(strsplit(rownames(posthoc), "\\."), '[', 2),
                      posthoc)


# keep p values
posthoc$Pr..F. <- signif(posthoc$Pr..F., 3)
posthoc.p <- dcast(Variable ~ Comparison, value.var = "Pr..F.", data = posthoc)
colnames(posthoc.p)[-1] <- paste0(colnames(posthoc.p)[-1], ".p")


# keep R squared values
posthoc$R2 <- signif(posthoc$R2, 3)
posthoc.R2 <- dcast(Variable ~ Comparison, value.var = "R2", data = posthoc)
colnames(posthoc.R2)[-1] <- paste0(colnames(posthoc.p)[-1], ".R2")


# combine p and R2
posthoc.p <- merge(posthoc.p, posthoc.R2, by = "Variable")
posthoc.p <- posthoc.p[c("Variable", sort(colnames(posthoc.p)[-1]))]


# combine permanova
perman <- merge(perman, posthoc.p, by = "Variable")
write.csv(perman, "DATA/PipeType_permanova.csv", row.names = F, na = "", quote = F)


# confidence intervals on the differences between the mean distance-to-centroid
plot(TukeyHSD(betadisper(bray.stat, info$Site_city)))
plot(TukeyHSD(betadisper(bray.stat, info$Pipe_material)))
plot(TukeyHSD(betadisper(bray.stat, info$Sample_surface)))
