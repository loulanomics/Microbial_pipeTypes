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




#################
### PERMANOVA ###
#################

# permanova
identical(rownames(relabun), info$Sample_ID)
perman <- adonis2(bray.stat ~ Pipe_line * Pipe_material * Pipe_age_years * 
                    Pipe_diameter_in * Sample_surface * Sample_local, 
                  data = info, na.action = na.exclude)


# extract results
perman <- data.frame(Variable = rownames(perman), perman)
colnames(perman)[6] <- "p"
perman$R2 <- signif(perman$R2, 3)


# save
write.csv(perman, "DATA/permanova.csv", row.names = F, na = "", quote = F)
