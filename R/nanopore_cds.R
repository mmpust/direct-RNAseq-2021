
# title: "CDS analysis"
# author: "Marie-Madlen Pust"
# Last updated: "23 10 2021"

############################################################################################################
# Clean environment
rm(list = ls())

# set working directory
setwd("R/")

# function for installing/importing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

packages <- c("Rsubread", "ggplot2","edgeR","pheatmap","plyr","dplyr","tidyr","ggpubr","vegan","lemon","factoextra","stringr","ggplotify","tidyverse","RColorBrewer","purrr")

############################################################################################################
# load packages
ipak(packages)

# load input
input_bam <- list.files(path = "b_minimap2/", pattern = ".bam", full.names = TRUE)
UserDefinedAnnotationRef <- "NN2_SG17M_CS_curated.gtf"
Ref_genome <- "NN2_SG17M_CS.fna"

# import annotation file
import_gtf <- read.table(UserDefinedAnnotationRef, sep = '\t', header = FALSE)
colnames(import_gtf) <- c("Isolate", "database", "featureType", "start", "end", "V6", "strandType", "V8", "V9")
import_gtf$V6 <- NULL
import_gtf$V8 <- NULL
import_gtf$V9 <- str_replace_all(import_gtf$V9, "transcript_id ", "")
import_gtf$V9 <- str_replace_all(import_gtf$V9, "gene_name ", "")

import_gtf <- separate(data = import_gtf, col = V9, into = c("transcript_id", "gene_name"), sep = "\\;")
import_gtf$gene_name2 <- ifelse(import_gtf$database == "Aragorn:001002", "tRNA", 
                               ifelse(import_gtf$database == "barrnap:0.9", "rRNA", 
                                      ifelse(import_gtf$database == "	Infernal:001001", "ncRNA", import_gtf$gene_name)))
import_gtf$gene_name2 <- ifelse(import_gtf$gene_name2 == "", import_gtf$start, import_gtf$gene_name2)
import_gtf$gene_name <- NULL
import_gtf$length <- import_gtf$end - import_gtf$start

transcript_id_list <- import_gtf$transcript_id
gene_name_list <- import_gtf$gene_name2
transcriptID_geneName <- data.frame(cbind(transcript_id_list, gene_name_list))
gene_name_list_2 <- import_gtf$gene_name2
names(gene_name_list_2) <- import_gtf$transcript_id

# import reads and count with featureCounts, antisense
featureCounts_PA_antisense <- featureCounts(input_bam, annot.ext = UserDefinedAnnotationRef, isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                                            GTF.featureType = c("CDS"), strandSpecific = 2, isLongRead = TRUE, useMetaFeatures = FALSE, allowMultiOverlap = TRUE, countMultiMappingReads = TRUE)

# import reads and count with featureCounts, sense
featureCounts_PA_sense <- featureCounts(input_bam, annot.ext = UserDefinedAnnotationRef, isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                                        GTF.featureType = c("CDS"), strandSpecific = 1, isLongRead = TRUE, useMetaFeatures = FALSE, allowMultiOverlap = TRUE, countMultiMappingReads = TRUE)

# extract count table, antisense
featureCounts_PA_antisense_df <- featureCounts_PA_antisense$counts
# extract count table, sense
featureCounts_PA_sense_df <- featureCounts_PA_sense$counts
# set column names in antisense table
colnames(featureCounts_PA_antisense_df) <- c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3","NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3", "SG17M_4h_BR1", "SG17M_4h_BR2", "SG17M_4h_BR3", "SG17M_8h_BR1","SG17M_8h_BR2","SG17M_8h_BR3")
# set column names in sense table
colnames(featureCounts_PA_sense_df) <- colnames(featureCounts_PA_antisense_df)
# convert data matrix to dataframe
featureCounts_PA_antisense_df <- data.frame(featureCounts_PA_antisense_df)
featureCounts_PA_sense_df <- data.frame(featureCounts_PA_sense_df)

mylist_antisense <- c()
for(i in rownames(featureCounts_PA_antisense_df)){
    newelem <- ifelse(i %in% names(gene_name_list_2), gene_name_list_2[[i]], i)
    mylist_antisense <- append(mylist_antisense, newelem)}
featureCounts_PA_antisense_df$name2 <- mylist_antisense
featureCounts_PA_antisense_df$name <- map(strsplit(featureCounts_PA_antisense_df$name2, "_"),1)
featureCounts_PA_antisense_df$name2 <- NULL
rownames(featureCounts_PA_antisense_df) <- NULL
featureCounts_PA_antisense_df <- featureCounts_PA_antisense_df %>% group_by(name) %>% summarise_all(funs(sum))
featureCounts_PA_antisense_df <- data.frame(featureCounts_PA_antisense_df)
rownames(featureCounts_PA_antisense_df) <- featureCounts_PA_antisense_df$name
featureCounts_PA_antisense_df$name <- NULL

mylist_sense <- c()
for(i in rownames(featureCounts_PA_sense_df)){
    newelem <- ifelse(i %in% names(gene_name_list_2), gene_name_list_2[[i]], i)
    mylist_sense <- append(mylist_sense, newelem)}
featureCounts_PA_sense_df$name2 <- mylist_sense
featureCounts_PA_sense_df$name <- map(strsplit(featureCounts_PA_sense_df$name2, "_"),1)
featureCounts_PA_sense_df$name2 <- NULL
rownames(featureCounts_PA_sense_df) <- NULL
featureCounts_PA_sense_df <- featureCounts_PA_sense_df %>% group_by(name) %>% summarise_all(funs(sum))
featureCounts_PA_sense_df <- data.frame(featureCounts_PA_sense_df)
rownames(featureCounts_PA_sense_df) <- featureCounts_PA_sense_df$name
featureCounts_PA_sense_df$name <- NULL

# Generate DeSEQ2 objects
y_antisense <- DGEList(counts=featureCounts_PA_antisense_df)
y_sense <- DGEList(counts=featureCounts_PA_sense_df)
# Grouping factor is added
group <- c("NN2_4h", "NN2_4h", "NN2_4h", "NN2_8h", "NN2_8h", "NN2_8h", "SG17M_4h", "SG17M_4h", "SG17M_4h", "SG17M_8h", "SG17M_8h", "SG17M_8h")
y_antisense <- DGEList(counts=y_antisense, group=group)
y_sense <- DGEList(counts=y_sense, group=group)

# Filter out lowly expressed genes using the following commands:
keep_antisense <- filterByExpr(y_antisense)
keep_sense <- filterByExpr(y_sense)
y_antisense <- y_antisense[keep_antisense, keep.lib.sizes=FALSE]
y_sense <- y_sense[keep_sense, keep.lib.sizes=FALSE]

# TMM read count normalisation, antisense
y_antisense <- calcNormFactors(y_antisense, method="TMM")
y_sense <- calcNormFactors(y_sense, method="TMM")

# Obtain copies per million
logcpm_antisense <- cpm(y_antisense, log=TRUE)
logcpm_sense <- cpm(y_sense, log=TRUE)
# Convert final output to data frames
logcpm_antisense_df <- data.frame(logcpm_antisense)
logcpm_sense_df <- data.frame(logcpm_sense)

############################################################################################################
# Start data analysis
# Ordination 
set.seed(2)
# samples as row
logcpm_antisense_saw <- t(logcpm_antisense)
logcpm_sense_saw <- t(logcpm_sense)

# non-metric multidimensional scaling with Bray-Curtis distances, two dimensions
mds_rna_antisense <- metaMDS(logcpm_antisense_saw, distance="bray", k=2, trymax=100, autotransform=FALSE)
mds_rna_sense <- metaMDS(logcpm_sense_saw, distance="bray", k=2, trymax=100, autotransform=FALSE)

# generate stressplots
stressplot(mds_rna_antisense)
stressplot(mds_rna_sense)

# prepare nmds dataframe for ggplot plotting
data.scores_antisense = as.data.frame(scores(mds_rna_antisense))
data.scores_sense = as.data.frame(scores(mds_rna_sense))
rownames(data.scores_antisense)
data.scores_antisense$time <- c("4h", "4h", "4h", "8h", "8h", "8h", "4h", "4h", "4h", "8h", "8h", "8h")
data.scores_antisense$isolate <- c("NN2", "NN2", "NN2", "NN2", "NN2", "NN2", "SG17M", "SG17M", "SG17M", "SG17M", "SG17M", "SG17M")
data.scores_sense$time <- data.scores_antisense$time
data.scores_sense$isolate <- data.scores_antisense$isolate

# plot nmds results for sense RNA
sense <- ggplot() + geom_point(data=data.scores_sense, aes(x=NMDS1, y=NMDS2, color=isolate, shape=time), size = 4, alpha=1) +
  theme_pubr(border=TRUE, legend="bottom") + scale_color_manual(values=c("darkred", "darkblue", "orange")) + xlim(-0.1,0.1) + ylim(-0.1,0.1) +
  geom_hline(yintercept = 0.00, alpha=0.5) + geom_vline(xintercept = 0.00, alpha=0.5) +
  theme(legend.title = element_blank()) + geom_label(aes(x=0.08, y = 0.1, label = "sense"), fill="white")

# plot nmds results for antisense RNA
antisense <- ggplot() + geom_point(data=data.scores_antisense, aes(x=NMDS1, y=NMDS2, color=isolate, shape=time), size = 4, alpha=1) +
  theme_pubr(border=TRUE, legend="bottom") + scale_color_manual(values=c("darkred", "darkblue", "orange")) + xlim(-0.1,0.1) + ylim(-0.1,0.1) +
  geom_hline(yintercept = 0.00, alpha=0.5) + geom_vline(xintercept = 0.00, alpha=0.5) +
  theme(legend.title = element_blank()) + geom_label(aes(x=0.08, y = 0.1, label = "antisense"), fill="white")

# merge sense and antisense NMDS plots
nmds_all <- ggarrange(sense, antisense, nrow=1, ncol=2, labels=c("a", "b"), common.legend = TRUE, legend = "right")
# Save nmds plot
ggsave(nmds_all, file="nmds_ont_plot.tif", device="tiff", dpi=300, width=8.69, height=3.96)

############################################################################################################
# Principal component analysis (PCA)
set.seed(2)
# Run PCA with and without center and scaling of data, respectively
pca_patternRec_antisense <- prcomp(logcpm_antisense_saw, center = TRUE, scale. = FALSE)
pca_patternRec_sense <- prcomp(logcpm_sense_saw, center = TRUE, scale. = FALSE)

# Visualise eigenvalues
eig_var_antisense <- fviz_eig(pca_patternRec_antisense)
eig_var_antisense$data$type = "antisense"
eig_var_sense <- fviz_eig(pca_patternRec_sense)
eig_var_sense$data$type = "sense"
# combine data frames
eig_var_both <- data.frame(rbind(eig_var_antisense$data, eig_var_sense$data))
# extract information of first and second dimension
eig_var_both_dim2 <- subset(eig_var_both, dim == "2" | dim == "1")
eig_var_both_dim2$dimension <- ifelse(eig_var_both_dim2$dim == "1", "Dim1", "Dim2") 
# generate time variable
time <- c("4h", "4h", "4h", "8h", "8h", "8h", "4h", "4h", "4h", "8h", "8h", "8h")
# generate strain variable
isolate <- c("NN2", "NN2", "NN2", "NN2", "NN2", "NN2", "SG17M", "SG17M", "SG17M", "SG17M", "SG17M", "SG17M")
# construct PCA plot of individuals, antisense
pcaInd_time_antisense <- fviz_pca_ind(pca_patternRec_antisense) 
pcaInd_time_antisense$data$isolate <- isolate
pcaInd_time_antisense$data$time <- time
# construct PCA plot of individuals, sense
pcaInd_time_sense <- fviz_pca_ind(pca_patternRec_sense) 
pcaInd_time_sense$data$isolate <- isolate
pcaInd_time_sense$data$time <- time
# plot dimension information
eig_var_antisense_sense_plot <- ggplot(eig_var_both_dim2) +
  geom_col(aes(x=type, y=eig, fill=type), alpha = 0.8) + xlab(" ") + ylab("Explained variance (in %)") +
  scale_fill_manual(values=c("black", "grey60")) +
  theme_pubr(border=TRUE, legend = "none", base_size = 9) + facet_wrap(~dimension) +
  coord_flip() + theme(axis.text = element_text(size=7), axis.title = element_text(size=7))
# create PCA plot, sense
sense_pca_plot <- ggplot(pcaInd_time_sense$data) +
  geom_point(aes(x=x, y=y, color=isolate, shape=time), size=3, alpha=1) +
  xlab("Dim1 (34 %)") + ylab("Dim2 (21 %)") +
  theme_pubr(border=TRUE, legend="bottom", base_size = 9) +
  xlim(-40, 60) + ylim(-40, 55) +
  scale_color_manual(values=c("darkred", "darkblue", "orange")) +
  geom_hline(yintercept = 0.00, alpha=0.5) + geom_vline(xintercept = 0.00, alpha=0.5) +
  theme(legend.title = element_blank(), axis.text = element_text(size=7), axis.title = element_text(size=7), legend.text = element_text(size=6)) +
  geom_label(aes(x=-20, y = 50, label = "sense"), fill="white", size=2.3)
# create PCA plot, antisense
antisense_pca_plot <- ggplot(pcaInd_time_antisense$data) +
  geom_point(aes(x=x, y=y, color=isolate, shape=time), size=3, alpha=1) +
  xlab("Dim1 (25 %)") + ylab("Dim2 (21 %)") +
  theme_pubr(border=TRUE, legend="bottom", base_size = 9) +
  xlim(-40, 60) + ylim(-40, 55) +
  scale_color_manual(values=c("darkred", "darkblue", "orange")) +
  geom_hline(yintercept = 0.00, alpha=0.5) + geom_vline(xintercept = 0.00, alpha=0.5) +
  theme(legend.title = element_blank(), axis.text = element_text(size=7), axis.title = element_text(size=7), legend.text = element_text(size=6)) + 
  geom_label(aes(x=-20, y = 50, label = "antisense"), fill="white", size=2.3)


pcaInd_time_sense$data$rna_type <- "sense"
pcaInd_time_antisense$data$rna_type <- "antisense"
pcaInd_time_sense_antisense <- data.frame(rbind(pcaInd_time_sense$data, pcaInd_time_antisense$data))
pcaInd_time_sense_antisense$colour_isolate <- with(pcaInd_time_sense_antisense,
                                                ifelse(rna_type == "sense" & isolate == "NN2", "NN2_sense",
                                                       ifelse(rna_type == "antisense" & isolate == "NN2", "NN2_antisense",
                                                              ifelse(rna_type == "sense" & isolate == "SG17M", 
                                                                     "SG17M_sense", "SG17M_antisense"))))
# generate merged PCA plot
merged_pca_plot <-
  ggplot(pcaInd_time_sense_antisense, aes(x=x, y=y)) +
  geom_point(aes(colour=colour_isolate, shape=time),size=2, alpha=1) +
  xlab("Dim1 (34 %)") + ylab("Dim2 (21 %)") +
  theme_pubr(border=TRUE, base_size = 9, legend = "bottom") +
  scale_color_manual(values=c("brown2", "brown4", "aquamarine2", "aquamarine4"),
                     labels=c("NN2_antisense", "NN2_sense", "SG17M_antisense", "SG17M_sense")) +
  geom_hline(yintercept = 0.00, alpha=0.5, linetype="dotted", colour="grey") + 
  geom_vline(xintercept = 0.00, alpha=0.5, linetype="dotted", colour="grey") +
  theme(legend.title = element_blank(), axis.text = element_text(size=7), axis.title = element_text(size=7), 
        legend.text = element_text(size=6)) + guides(colour=guide_legend(ncol=3))

# make multiplot panel 1
a <- ggarrange(sense_pca_plot, antisense_pca_plot, nrow=1, ncol=2, labels=c("A", "B"), common.legend = TRUE, font.label = list(color = "black", size=8), legend = "bottom")
# make multiplot panel 2
a2 <- ggarrange(eig_var_antisense_sense_plot, nrow = 1, ncol = 1, labels = "C", font.label = list(color = "black", size=8))
# make multiplot panel 3
a3 <- ggarrange(a, a2, nrow=2, heights = c(1,0.4))
# save multiplot figure
ggsave(a3, file="save_figures/Figure_08.tif", device="tiff", width=12, height=12, dpi=300, units = "cm", bg="white")

############################################################################################################
# Extract 30 sense variables explaining most of the variance
res.var.sense <- get_pca_var(pca_patternRec_sense)
res.var.sense <- data.frame(res.var.sense$contrib)        
res.var.sense <- res.var.sense[order(-res.var.sense$Dim.2),]
res.var.sense_50 <- data.frame(res.var.sense[1:30,])
list_sense_50 <- rownames(res.var.sense_50)
list_sense_50 <- str_replace(list_sense_50, " ", "")
rownames(logcpm_sense_df) <- str_replace(rownames(logcpm_sense_df), " ", "")
logcpm_sense_df_h <- logcpm_sense_df[rownames(logcpm_sense_df) %in% list_sense_50,]

rownames(logcpm_sense_df_h) <- c("rihB", "prpL", "920454", "merP", "lgrE", "dat", "acrE", "4224471", "inh", "4787153", "5067956", "5071406", "5126084", "5171225", "5204149", "5606105",  "ftsA", "5611644", "5612036", "hin", "virD4", "5628279", "5643301", "1145152", "2702075",  "2710069", "3055598", "3056251", "5244612","5856913") 
breaksList = seq(0, 15, by = 1)
# make sense heatmap
heatmap_genes_4h_sense <- pheatmap(logcpm_sense_df_h, angle_col = "45", fontsize = 8, cellwidth = 20, cellheight = 19, scale = "none",
                                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList,
                                   show_colnames = T, cutree_cols = 2, legend = F, treeheight_col = 5, treeheight_row = 0)
# store sense heatmap as ggplot object
heatmap_genes_4h_sense_gg <- as.ggplot(heatmap_genes_4h_sense, scale=1)

# Extract 30 antisense variables explaining most of the variance
res.var.antisense <- get_pca_var(pca_patternRec_antisense)
res.var.antisense <- data.frame(res.var.antisense$contrib)        
res.var.antisense <- res.var.antisense[order(-res.var.antisense$Dim.2),]
res.var.antisense_50 <- data.frame(res.var.antisense[1:30,])
list_antisense_50 <- rownames(res.var.antisense_50)
list_antisense_50 <- str_replace(list_antisense_50, " ", "")
rownames(logcpm_antisense_df) <- str_replace(rownames(logcpm_antisense_df), " ", "")
logcpm_antisense_df_h <- logcpm_antisense_df[rownames(logcpm_antisense_df) %in% list_antisense_50,]

rownames(logcpm_antisense_df_h) <- c("rihB", "pntAA", "aam", "madA", "yycB", "mexR", "fecI", "628030","721083", "phzS", "hpcH", "lnt", "xerD", "folE", "eptA", "bluF", "virS","ugpQ", "lasR", "gloC", "slyB", "4786592", "ftsZ", "adh1", "ppx","mcbR", "852731", "945784", "merR", "5244051")
# make antisense heatmap
heatmap_genes_4h_antisense <- pheatmap(logcpm_antisense_df_h, angle_col = "45", fontsize = 8, 
                                       cellwidth = 20, cellheight = 19, scale = "none",
                                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, 
                                       show_colnames = T, cutree_cols = 2, legend = T, treeheight_col = 5, treeheight_row = 0)
# store antisense heatmap as ggplot object
heatmap_genes_4h_antisense_gg <- as.ggplot(heatmap_genes_4h_antisense, scale=1)

# merge sense and antisense heatmaps
heatmaps <- ggarrange(heatmap_genes_4h_sense_gg, heatmap_genes_4h_antisense_gg, nrow=1, labels=c("A", "B"), common.legend = TRUE, legend = "top")
# store sense-antisense heatmap
ggsave(heatmaps, file="save_figures/Figure_09.tif", device="tiff", dpi=300, width=25, height=25, units = "cm", bg = "white")
