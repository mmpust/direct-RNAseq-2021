# title: "most_expressed"
# author: "Marie Pust"
# date: "17 6 2021"

# Clean environment
rm(list = ls())

library(rlist)
library(ggplot2)
library(readr)
library(readxl)
library(Rsubread)
library(edgeR)
library(stringr)
library(tidyr)
library(ggpubr)
library(circlize)
library(dplyr)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(vegan)
library(rcompanion)

setwd("/R")

input_bam_NN2<- list.files(
  path = "NN2_run/d_raw_bam",
  pattern = ".bam", full.names = TRUE)

input_bam_sg17m <- list.files(
  path = "SG17M_run/d_raw_bam",
  pattern = ".bam", full.names = TRUE)

UserDefinedAnnotationRef <- "NN2_ENO_SG17M_curated.gtf"
UserDefinedAnnotationRef_NN2 <- "NN2_ENO_curated.gtf"
UserDefinedAnnotationRef_sg17m <- "SG17M_ENO_curated.gtf"

import_gtf_NN2 <- read.table(UserDefinedAnnotationRef_NN2, sep = '\t', header = FALSE)
colnames(import_gtf_NN2) <- c("Isolate", "database", "featureType", "start", "end", "V6", "strandType", "V8", "V9")
import_gtf_NN2$V6 <- NULL
import_gtf_NN2$V8 <- NULL
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "transcript_id ", "")
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "gene_name ", "")
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "tRNA_type ", "")
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "rRNA_type ", "")
import_gtf_NN2 <- separate(data = import_gtf_NN2, col = V9, into = c("transcript_id", "gene_name2"), sep = "\\;")

import_gtf_NN2$length <- import_gtf_NN2$end - import_gtf_NN2$start
import_gtf_NN2$gene_name3 <-
  with(import_gtf_NN2, 
       ifelse(database == "barrnap:0.9", "rRNA",
              ifelse(database == "Infernal:001001", "ncRNA",
                     ifelse(database == "Aragorn:001002", "tRNA",
                            ifelse(database == "IslandViewer4", "GI",
                                   ifelse(database == "Prodigal:002006", "CDS",database))))))

nn2_accessory <- subset(import_gtf_NN2, gene_name3 == "GI")
nn2_accessory_bed <- select(nn2_accessory, c(start, end))
nn2_tRNA <- subset(import_gtf_NN2, gene_name3 == "tRNA")
nn2_tRNA_bed <- select(nn2_accessory, c(start, end))


import_gtf_sg17m <- read.table(UserDefinedAnnotationRef_sg17m, sep = '\t', header = FALSE)
colnames(import_gtf_sg17m) <- c("Isolate", "database", "featureType", "start", "end", "V6", "strandType", "V8", "V9")
import_gtf_sg17m$V6 <- NULL
import_gtf_sg17m$V8 <- NULL
import_gtf_sg17m$V9 <- str_replace_all(import_gtf_sg17m$V9, "transcript_id ", "")
import_gtf_sg17m$V9 <- str_replace_all(import_gtf_sg17m$V9, "gene_name ", "")
import_gtf_sg17m$V9 <- str_replace_all(import_gtf_sg17m$V9, "tRNA_type ", "")
import_gtf_sg17m$V9 <- str_replace_all(import_gtf_sg17m$V9, "rRNA_type ", "")
import_gtf_sg17m <- separate(data = import_gtf_sg17m, col = V9, into = c("transcript_id", "gene_name2"), sep = "\\;")
import_gtf_sg17m$length <- import_gtf_sg17m$end - import_gtf_sg17m$start
import_gtf_sg17m$gene_name3 <-
  with(import_gtf_sg17m, 
       ifelse(database == "barrnap:0.9", "rRNA",
              ifelse(database == "Infernal:001001", "ncRNA",
                     ifelse(database == "Aragorn:001002", "tRNA",
                            ifelse(database == "IslandViewer4", "GI",
                                   ifelse(database == "Prodigal:002006", "CDS",database))))))

sg17m_accessory <- subset(import_gtf_sg17m, gene_name3 == "GI")
sg17m_accessory_bed <- select(sg17m_accessory, c(start, end))
sg17m_tRNA <- subset(import_gtf_sg17m, gene_name3 == "tRNA")
sg17m_tRNA_bed <- select(sg17m_accessory, c(start, end))

featureCounts_PA_antisenseNN2 <- 
  featureCounts(input_bam_NN2, annot.ext = UserDefinedAnnotationRef_NN2, 
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 2, isLongRead = TRUE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

featureCounts_PA_senseNN2 <- 
  featureCounts(input_bam_NN2, annot.ext = UserDefinedAnnotationRef_NN2,
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 1, isLongRead = TRUE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

featureCounts_PA_antisensesg17m <- 
  featureCounts(input_bam_sg17m, annot.ext = UserDefinedAnnotationRef_sg17m, 
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 2, isLongRead = TRUE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

featureCounts_PA_sensesg17m <- 
  featureCounts(input_bam_sg17m, annot.ext = UserDefinedAnnotationRef_sg17m,
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 1, isLongRead = TRUE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

featureCounts_PA_antisense_df_NN2 <- featureCounts_PA_antisenseNN2$counts
featureCounts_PA_antisense_df_NN2 <- data.frame(featureCounts_PA_antisense_df_NN2)
featureCounts_PA_sense_df_NN2 <- featureCounts_PA_senseNN2$counts
featureCounts_PA_sense_df_NN2 <- data.frame(featureCounts_PA_sense_df_NN2)
colnames(featureCounts_PA_antisense_df_NN2) <- c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3", 
                                                 "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3")
colnames(featureCounts_PA_sense_df_NN2) <- c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3", 
                                             "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3")

featureCounts_PA_antisense_df_sg17m <- featureCounts_PA_antisensesg17m$counts
featureCounts_PA_antisense_df_sg17m <- data.frame(featureCounts_PA_antisense_df_sg17m)
featureCounts_PA_sense_df_sg17m <- featureCounts_PA_sensesg17m$counts
featureCounts_PA_sense_df_sg17m <- data.frame(featureCounts_PA_sense_df_sg17m)
colnames(featureCounts_PA_antisense_df_sg17m) <- c("sg17m_4h_BR1", "sg17m_4h_BR2", "sg17m_4h_BR3", 
                                                 "sg17m_8h_BR1", "sg17m_8h_BR2", "sg17m_8h_BR3")
colnames(featureCounts_PA_sense_df_sg17m) <- c("sg17m_4h_BR1", "sg17m_4h_BR2", "sg17m_4h_BR3", 
                                             "sg17m_8h_BR1", "sg17m_8h_BR2", "sg17m_8h_BR3")

featureCounts_PA_sense_df_NN2_norm <- featureCounts_PA_sense_df_NN2
featureCounts_PA_sense_df_NN2_norm_4h <- featureCounts_PA_sense_df_NN2_norm[,1:3]
featureCounts_PA_sense_df_NN2_norm_4h_1 <- DGEList(counts=featureCounts_PA_sense_df_NN2_norm_4h)
featureCounts_PA_sense_df_NN2_norm_4h_2 <- calcNormFactors(featureCounts_PA_sense_df_NN2_norm_4h_1, 
                                                           method="TMM")
featureCounts_PA_sense_df_NN2_norm_4h_3 <- cpm(featureCounts_PA_sense_df_NN2_norm_4h_2, log=FALSE)
featureCounts_PA_sense_df_NN2_norm_4h_4 <- data.frame(featureCounts_PA_sense_df_NN2_norm_4h_3)

featureCounts_PA_sense_df_NN2_norm_8h <- featureCounts_PA_sense_df_NN2_norm[,4:6]
featureCounts_PA_sense_df_NN2_norm_8h_1 <- DGEList(counts=featureCounts_PA_sense_df_NN2_norm_8h)
featureCounts_PA_sense_df_NN2_norm_8h_2 <- calcNormFactors(featureCounts_PA_sense_df_NN2_norm_8h_1, 
                                                           method="TMM")
featureCounts_PA_sense_df_NN2_norm_8h_3 <- cpm(featureCounts_PA_sense_df_NN2_norm_8h_2, log=FALSE)
featureCounts_PA_sense_df_NN2_norm_8h_4 <- data.frame(featureCounts_PA_sense_df_NN2_norm_8h_3)

featureCounts_PA_antisense_df_NN2_norm <- featureCounts_PA_antisense_df_NN2
featureCounts_PA_antisense_df_NN2_norm_4h <- featureCounts_PA_antisense_df_NN2_norm[,1:3]
featureCounts_PA_antisense_df_NN2_norm_4h_1 <- DGEList(counts=featureCounts_PA_antisense_df_NN2_norm_4h)
featureCounts_PA_antisense_df_NN2_norm_4h_2 <- calcNormFactors(featureCounts_PA_antisense_df_NN2_norm_4h_1,
                                                               method="TMM")
featureCounts_PA_antisense_df_NN2_norm_4h_3 <- cpm(featureCounts_PA_antisense_df_NN2_norm_4h_2, log=FALSE)
featureCounts_PA_antisense_df_NN2_norm_4h_4 <- data.frame(featureCounts_PA_antisense_df_NN2_norm_4h_3)
featureCounts_PA_antisense_df_NN2_norm_8h <- featureCounts_PA_antisense_df_NN2_norm[,4:6]
featureCounts_PA_antisense_df_NN2_norm_8h_1 <- DGEList(counts=featureCounts_PA_antisense_df_NN2_norm_8h)
featureCounts_PA_antisense_df_NN2_norm_8h_2 <- calcNormFactors(featureCounts_PA_antisense_df_NN2_norm_8h_1,
                                                               method="TMM")
featureCounts_PA_antisense_df_NN2_norm_8h_3 <- cpm(featureCounts_PA_antisense_df_NN2_norm_8h_2, log=FALSE)
featureCounts_PA_antisense_df_NN2_norm_8h_4 <- data.frame(featureCounts_PA_antisense_df_NN2_norm_8h_3)

featureCounts_PA_sense_df_sg17m_norm <- featureCounts_PA_sense_df_sg17m
featureCounts_PA_sense_df_sg17m_norm_4h <- featureCounts_PA_sense_df_sg17m_norm[,1:3]
featureCounts_PA_sense_df_sg17m_norm_4h_1 <- DGEList(counts=featureCounts_PA_sense_df_sg17m_norm_4h)
featureCounts_PA_sense_df_sg17m_norm_4h_2 <- calcNormFactors(featureCounts_PA_sense_df_sg17m_norm_4h_1, 
                                                             method="TMM")
featureCounts_PA_sense_df_sg17m_norm_4h_3 <- cpm(featureCounts_PA_sense_df_sg17m_norm_4h_2, log=FALSE)
featureCounts_PA_sense_df_sg17m_norm_4h_4 <- data.frame(featureCounts_PA_sense_df_sg17m_norm_4h_3)
featureCounts_PA_sense_df_sg17m_norm_8h <- featureCounts_PA_sense_df_sg17m_norm[,4:6]
featureCounts_PA_sense_df_sg17m_norm_8h_1 <- DGEList(counts=featureCounts_PA_sense_df_sg17m_norm_8h)
featureCounts_PA_sense_df_sg17m_norm_8h_2 <- calcNormFactors(featureCounts_PA_sense_df_sg17m_norm_8h_1, 
                                                             method="TMM")
featureCounts_PA_sense_df_sg17m_norm_8h_3 <- cpm(featureCounts_PA_sense_df_sg17m_norm_8h_2, log=FALSE)
featureCounts_PA_sense_df_sg17m_norm_8h_4 <- data.frame(featureCounts_PA_sense_df_sg17m_norm_8h_3)

featureCounts_PA_antisense_df_sg17m_norm <- featureCounts_PA_antisense_df_sg17m
featureCounts_PA_antisense_df_sg17m_norm_4h <- featureCounts_PA_antisense_df_sg17m_norm[,1:3]
featureCounts_PA_antisense_df_sg17m_norm_4h_1 <- DGEList(counts=featureCounts_PA_antisense_df_sg17m_norm_4h)
featureCounts_PA_antisense_df_sg17m_norm_4h_2 <- calcNormFactors(featureCounts_PA_antisense_df_sg17m_norm_4h_1,
                                                                 method="TMM")
featureCounts_PA_antisense_df_sg17m_norm_4h_3 <- cpm(featureCounts_PA_antisense_df_sg17m_norm_4h_2, log=FALSE)
featureCounts_PA_antisense_df_sg17m_norm_4h_4 <- data.frame(featureCounts_PA_antisense_df_sg17m_norm_4h_3)
featureCounts_PA_antisense_df_sg17m_norm_8h <- featureCounts_PA_antisense_df_sg17m_norm[,4:6]
featureCounts_PA_antisense_df_sg17m_norm_8h_1 <- DGEList(counts=featureCounts_PA_antisense_df_sg17m_norm_8h)
featureCounts_PA_antisense_df_sg17m_norm_8h_2 <- calcNormFactors(featureCounts_PA_antisense_df_sg17m_norm_8h_1,
                                                                 method="TMM")
featureCounts_PA_antisense_df_sg17m_norm_8h_3 <- cpm(featureCounts_PA_antisense_df_sg17m_norm_8h_2, log=FALSE)
featureCounts_PA_antisense_df_sg17m_norm_8h_4 <- data.frame(featureCounts_PA_antisense_df_sg17m_norm_8h_3)


# dev.off()
# circular plots, NN2, 4h
antisense_NN2_norm_4h_circ <- cpm(featureCounts_PA_antisense_df_NN2, log=TRUE)
antisense_NN2_norm_4h_circ <- data.frame(antisense_NN2_norm_4h_circ)[,c(1:3)]
antisense_NN2_norm_4h_circ_gtf <- subset(import_gtf_NN2, transcript_id %in% rownames(antisense_NN2_norm_4h_circ))
antisense_NN2_norm_4h_circ_gtf <- antisense_NN2_norm_4h_circ_gtf[order(antisense_NN2_norm_4h_circ_gtf$transcript_id),]
antisense_NN2_norm_4h_circ <- antisense_NN2_norm_4h_circ[order(rownames(antisense_NN2_norm_4h_circ)),]

antisense_NN2_norm_4h_circ$start <- antisense_NN2_norm_4h_circ_gtf$start
antisense_NN2_norm_4h_circ$end <- antisense_NN2_norm_4h_circ_gtf$end
antisense_NN2_norm_4h_circ$isolate <- "NN2"

antisense_NN2_norm_4h_circ_BR1 <- select(antisense_NN2_norm_4h_circ, c(isolate,start, end, NN2_4h_BR1))
color_NN2_BR1_4h <- with(antisense_NN2_norm_4h_circ_BR1,
                                             ifelse(NN2_4h_BR1 >= 12, "palegreen4", 
                                               ifelse(NN2_4h_BR1 < 12 & NN2_4h_BR1 >= 10, "palegreen3", 
                                                      ifelse(NN2_4h_BR1 < 10 & NN2_4h_BR1 >= 9, "palegreen2",
                                                             ifelse(NN2_4h_BR1 < 9 & NN2_4h_BR1 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(NN2_4h_BR1 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_NN2_norm_4h_circ_BR2 <- select(antisense_NN2_norm_4h_circ, c(isolate,start, end, NN2_4h_BR2))
color_NN2_BR2_4h <- with(antisense_NN2_norm_4h_circ_BR2,
                                             ifelse(NN2_4h_BR2 >= 12, "palegreen4", 
                                               ifelse(NN2_4h_BR2 < 12 & NN2_4h_BR2 >= 10, "palegreen3", 
                                                      ifelse(NN2_4h_BR2 < 10 & NN2_4h_BR2 >= 9, "palegreen2",
                                                             ifelse(NN2_4h_BR2 < 9 & NN2_4h_BR2 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(NN2_4h_BR2 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_NN2_norm_4h_circ_BR3 <- select(antisense_NN2_norm_4h_circ, c(isolate, start, end, NN2_4h_BR3))
color_NN2_BR3_4h <- with(antisense_NN2_norm_4h_circ_BR3,
                                             ifelse(NN2_4h_BR3 >= 12, "palegreen4", 
                                               ifelse(NN2_4h_BR3 < 12 & NN2_4h_BR3 >= 10, "palegreen3", 
                                                      ifelse(NN2_4h_BR3 < 10 & NN2_4h_BR3 >= 9, "palegreen2",
                                                             ifelse(NN2_4h_BR3 < 9 & NN2_4h_BR3 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(NN2_4h_BR3 < 6, "palegoldenrod",
                                                                           "NA"))))))

# circular plots, NN2, 8h
antisense_NN2_norm_8h_circ <- cpm(featureCounts_PA_antisense_df_NN2, log=TRUE)
antisense_NN2_norm_8h_circ <- data.frame(antisense_NN2_norm_8h_circ)[,c(4:6)]
antisense_NN2_norm_8h_circ_gtf <- subset(import_gtf_NN2, transcript_id %in% rownames(antisense_NN2_norm_8h_circ))
antisense_NN2_norm_8h_circ_gtf <- antisense_NN2_norm_8h_circ_gtf[order(antisense_NN2_norm_8h_circ_gtf$transcript_id),]
antisense_NN2_norm_8h_circ <- antisense_NN2_norm_8h_circ[order(rownames(antisense_NN2_norm_8h_circ)),]

antisense_NN2_norm_8h_circ$start <- antisense_NN2_norm_8h_circ_gtf$start
antisense_NN2_norm_8h_circ$end <- antisense_NN2_norm_8h_circ_gtf$end
antisense_NN2_norm_8h_circ$isolate <- "NN2"

antisense_NN2_norm_8h_circ_BR1 <- select(antisense_NN2_norm_8h_circ, c(isolate,start, end, NN2_8h_BR1))
color_NN2_BR1_8h <- with(antisense_NN2_norm_8h_circ_BR1,
                                             ifelse(NN2_8h_BR1 >= 12, "palegreen4", 
                                               ifelse(NN2_8h_BR1 < 12 & NN2_8h_BR1 >= 10, "palegreen3", 
                                                      ifelse(NN2_8h_BR1 < 10 & NN2_8h_BR1 >= 9, "palegreen2",
                                                             ifelse(NN2_8h_BR1 < 9 & NN2_8h_BR1 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(NN2_8h_BR1 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_NN2_norm_8h_circ_BR2 <- select(antisense_NN2_norm_8h_circ, c(isolate,start, end, NN2_8h_BR2))
color_NN2_BR2_8h <- with(antisense_NN2_norm_8h_circ_BR2,
                                             ifelse(NN2_8h_BR2 >= 12, "palegreen4", 
                                               ifelse(NN2_8h_BR2 < 12 & NN2_8h_BR2 >= 10, "palegreen3", 
                                                      ifelse(NN2_8h_BR2 < 10 & NN2_8h_BR2 >= 9, "palegreen2",
                                                             ifelse(NN2_8h_BR2 < 9 & NN2_8h_BR2 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(NN2_8h_BR2 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_NN2_norm_8h_circ_BR3 <- select(antisense_NN2_norm_8h_circ, c(isolate,start, end, NN2_8h_BR3))
color_NN2_BR3_8h <- with(antisense_NN2_norm_8h_circ_BR3,
                                             ifelse(NN2_8h_BR3 >= 12, "palegreen4", 
                                               ifelse(NN2_8h_BR3 < 12 & NN2_8h_BR3 >= 10, "palegreen3", 
                                                      ifelse(NN2_8h_BR3 < 10 & NN2_8h_BR3 >= 9, "palegreen2",
                                                             ifelse(NN2_8h_BR3 < 9 & NN2_8h_BR3 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(NN2_8h_BR3 < 6, "palegoldenrod",
                                                                           "NA"))))))


tiff(filename = "save_figures/circular_plot_nn2_4h_800dpi.tif", units="cm", res=800, width=18, height=16)
set.seed(123)
circos.clear()
par(mar=c(0, 0, 0, 0))
circos.par(start.degree = 90,
           cell.padding = c(0, 0, 0, 0))

# Initialize genome (bed file with genome sizes)
genome <- data.frame(chr=c("NN2"), start = c(1), end = c(7000000))
circos.genomicInitialize(genome, plotType = c("axis"), major.by = 1000000,
                         axis.labels.cex = 0.7, labels.cex = 0.7)

# Add track with annotation
feature <- data.frame(chr = c(nn2_accessory$Isolate), 
                      start = c(nn2_accessory$start), 
                      end = c(nn2_accessory$end))

circos.genomicTrack(feature, ylim=c(0,1), track.height = 0.1, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, col="gold", border="gold")}, bg.col = "dodgerblue4")

feature_tRNA <- data.frame(chr = c(nn2_tRNA$Isolate), 
                      start = c(nn2_tRNA$start), 
                      end = c(nn2_tRNA$end),
                      value = 2)

circos.genomicDensity(feature_tRNA, count_by = "number", col="black", area=TRUE)


antisense_cov1 <- data.frame(chr = c(antisense_NN2_norm_4h_circ_BR1$isolate),
                            start = c(antisense_NN2_norm_4h_circ_BR1$start),
                            end = c(antisense_NN2_norm_4h_circ_BR1$end),
                            depth = c(antisense_NN2_norm_4h_circ_BR1$NN2_4h_BR1))

circos.genomicTrack(antisense_cov1, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_NN2_BR1_4h, ...)})

antisense_cov2 <- data.frame(chr = c(antisense_NN2_norm_4h_circ_BR2$isolate),
                            start = c(antisense_NN2_norm_4h_circ_BR2$start),
                            end = c(antisense_NN2_norm_4h_circ_BR2$end),
                            depth = c(antisense_NN2_norm_4h_circ_BR2$NN2_4h_BR2))

circos.genomicTrack(antisense_cov2, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_NN2_BR2_4h, ...)})

antisense_cov3 <- data.frame(chr = c(antisense_NN2_norm_4h_circ_BR3$isolate),
                            start = c(antisense_NN2_norm_4h_circ_BR3$start),
                            end = c(antisense_NN2_norm_4h_circ_BR3$end),
                            depth = c(antisense_NN2_norm_4h_circ_BR3$NN2_4h_BR3))

circos.genomicTrack(antisense_cov3, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_NN2_BR3_4h, ...)})
circos.text(6200000,10, "A", cex=1.5, facing = "downward", font=2)
text(0,0, "NN2-4h", cex=0.8)
dev.off()

tiff(filename = "save_figures/circular_plot_nn2_8h_800dpi.tif", units="cm", res=800, width=18, height=16)
circos.clear()
par(mar=c(0, 0, 0, 0))
circos.par(start.degree = 90,
           cell.padding = c(0, 0, 0, 0))

# Initialize genome (bed file with genome sizes)
genome <- data.frame(chr=c("NN2"), start = c(1), end = c(7000000))
circos.genomicInitialize(genome, plotType = c("axis"), major.by = 1000000,
                         axis.labels.cex = 0.7, labels.cex = 0.7)

# Add track with annotation
circos.genomicTrack(feature, ylim=c(0,1), track.height = 0.1, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, col="gold", border="gold")}, bg.col = "dodgerblue4")

circos.genomicDensity(feature_tRNA, count_by = "number", col="black", area=TRUE)


antisense_cov4 <- data.frame(chr = c(antisense_NN2_norm_8h_circ_BR1$isolate),
                            start = c(antisense_NN2_norm_8h_circ_BR1$start),
                            end = c(antisense_NN2_norm_8h_circ_BR1$end),
                            depth = c(antisense_NN2_norm_8h_circ_BR1$NN2_8h_BR1))

circos.genomicTrack(antisense_cov4, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_NN2_BR1_8h, ...)})

antisense_cov5 <- data.frame(chr = c(antisense_NN2_norm_8h_circ_BR2$isolate),
                            start = c(antisense_NN2_norm_8h_circ_BR2$start),
                            end = c(antisense_NN2_norm_8h_circ_BR2$end),
                            depth = c(antisense_NN2_norm_8h_circ_BR2$NN2_8h_BR2))

circos.genomicTrack(antisense_cov5, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_NN2_BR2_8h, ...)})

antisense_cov6 <- data.frame(chr = c(antisense_NN2_norm_8h_circ_BR3$isolate),
                            start = c(antisense_NN2_norm_8h_circ_BR3$start),
                            end = c(antisense_NN2_norm_8h_circ_BR3$end),
                            depth = c(antisense_NN2_norm_8h_circ_BR3$NN2_8h_BR3))

circos.genomicTrack(antisense_cov6, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_NN2_BR3_8h, ...)})
circos.text(6200000,10, "B", cex=1.5, facing = "downward", font=2)
text(0,0, "NN2-8h", cex=0.8)
dev.off()

# circular plots, SG17M, 4h
antisense_sg17m_norm_4h_circ <- cpm(featureCounts_PA_antisense_df_sg17m, log=TRUE)
antisense_sg17m_norm_4h_circ <- data.frame(antisense_sg17m_norm_4h_circ)[,c(1:3)]
antisense_sg17m_norm_4h_circ_gtf <- subset(import_gtf_sg17m, transcript_id %in% rownames(antisense_sg17m_norm_4h_circ))
antisense_sg17m_norm_4h_circ_gtf <- antisense_sg17m_norm_4h_circ_gtf[order(antisense_sg17m_norm_4h_circ_gtf$transcript_id),]
antisense_sg17m_norm_4h_circ <- antisense_sg17m_norm_4h_circ[order(rownames(antisense_sg17m_norm_4h_circ)),]

antisense_sg17m_norm_4h_circ$start <- antisense_sg17m_norm_4h_circ_gtf$start
antisense_sg17m_norm_4h_circ$end <- antisense_sg17m_norm_4h_circ_gtf$end
antisense_sg17m_norm_4h_circ$isolate <- "SG17M"
antisense_sg17m_norm_4h_circ_BR1 <- select(antisense_sg17m_norm_4h_circ, c(isolate,start, end, sg17m_4h_BR1))
color_sg17m_BR1_4h <- with(antisense_sg17m_norm_4h_circ_BR1,
                                             ifelse(sg17m_4h_BR1 >= 12, "palegreen4", 
                                               ifelse(sg17m_4h_BR1 < 12 & sg17m_4h_BR1 >= 10, "palegreen3", 
                                                      ifelse(sg17m_4h_BR1 < 10 & sg17m_4h_BR1 >= 9,
                                                             "palegreen2",
                                                             ifelse(sg17m_4h_BR1 < 9 & sg17m_4h_BR1 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(sg17m_4h_BR1 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_sg17m_norm_4h_circ_BR2 <- select(antisense_sg17m_norm_4h_circ, c(isolate,start, end, sg17m_4h_BR2))
color_sg17m_BR2_4h <- with(antisense_sg17m_norm_4h_circ_BR2,
                                             ifelse(sg17m_4h_BR2 >= 12, "palegreen4", 
                                               ifelse(sg17m_4h_BR2 < 12 & sg17m_4h_BR2 >= 10, "palegreen3", 
                                                      ifelse(sg17m_4h_BR2 < 10 & sg17m_4h_BR2 >= 9,
                                                             "palegreen2",
                                                             ifelse(sg17m_4h_BR2 < 9 & sg17m_4h_BR2 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(sg17m_4h_BR2 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_sg17m_norm_4h_circ_BR3 <- select(antisense_sg17m_norm_4h_circ, c(isolate,start, end, sg17m_4h_BR3))
color_sg17m_BR3_4h <- with(antisense_sg17m_norm_4h_circ_BR3,
                                             ifelse(sg17m_4h_BR3 >= 12, "palegreen4", 
                                               ifelse(sg17m_4h_BR3 < 12 & sg17m_4h_BR3 >= 10, "palegreen3", 
                                                      ifelse(sg17m_4h_BR3 < 10 & sg17m_4h_BR3 >= 9,
                                                             "palegreen2",
                                                             ifelse(sg17m_4h_BR3 < 9 & sg17m_4h_BR3 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(sg17m_4h_BR3 < 6, "palegoldenrod",
                                                                           "NA"))))))

tiff(filename = "save_figures/circular_plot_sg17m_4h_800dpi.tif", units="cm", res=800, width=18, height=16)
set.seed(123)
circos.clear()
par(mar=c(0, 0, 0, 0))
circos.par(start.degree = 90,
           cell.padding = c(0, 0, 0, 0))

# Initialize genome (bed file with genome sizes)
genome <- data.frame(chr=c("SG17M"), start = c(1), end = c(7000000))
circos.genomicInitialize(genome, plotType = c("axis"), major.by = 1000000,
                         axis.labels.cex = 0.7, labels.cex = 0.7)

# Add track with annotation
feature_sg17m <- data.frame(chr = c(sg17m_accessory$Isolate), 
                      start = c(sg17m_accessory$start), 
                      end = c(sg17m_accessory$end))

circos.genomicTrack(feature_sg17m, ylim=c(0,1), track.height = 0.1, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, col="gold", border="gold")}, bg.col = "dodgerblue4")

feature_tRNA_sg17m <- data.frame(chr = c(sg17m_tRNA$Isolate), 
                      start = c(sg17m_tRNA$start), 
                      end = c(sg17m_tRNA$end),
                      value = 2)

circos.genomicDensity(feature_tRNA_sg17m, count_by = "number", col="black", area=TRUE)

antisense_cov1_sg17m <- data.frame(chr = c(antisense_sg17m_norm_4h_circ_BR1$isolate),
                            start = c(antisense_sg17m_norm_4h_circ_BR1$start),
                            end = c(antisense_sg17m_norm_4h_circ_BR1$end),
                            depth = c(antisense_sg17m_norm_4h_circ_BR1$sg17m_4h_BR1))

circos.genomicTrack(antisense_cov1_sg17m, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_sg17m_BR1_4h, ...)})

antisense_cov2_sg17m <- data.frame(chr = c(antisense_sg17m_norm_4h_circ_BR2$isolate),
                            start = c(antisense_sg17m_norm_4h_circ_BR2$start),
                            end = c(antisense_sg17m_norm_4h_circ_BR2$end),
                            depth = c(antisense_sg17m_norm_4h_circ_BR2$sg17m_4h_BR2))

circos.genomicTrack(antisense_cov2_sg17m, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_sg17m_BR2_4h, ...)})

antisense_cov3_sg17m <- data.frame(chr = c(antisense_sg17m_norm_4h_circ_BR3$isolate),
                            start = c(antisense_sg17m_norm_4h_circ_BR3$start),
                            end = c(antisense_sg17m_norm_4h_circ_BR3$end),
                            depth = c(antisense_sg17m_norm_4h_circ_BR3$sg17m_4h_BR3))

circos.genomicTrack(antisense_cov3_sg17m, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_sg17m_BR3_4h, ...)})
circos.text(6200000,10, "C", cex=1.5, facing = "downward", font=2)
text(0,0, "SG17M-4h", cex=0.8)
dev.off()

# circular plots, SG17M, 8h
antisense_sg17m_norm_8h_circ <- cpm(featureCounts_PA_antisense_df_sg17m, log=TRUE)
antisense_sg17m_norm_8h_circ <- data.frame(antisense_sg17m_norm_8h_circ)[,c(4:6)]
antisense_sg17m_norm_8h_circ_gtf <- subset(import_gtf_sg17m, transcript_id %in% rownames(antisense_sg17m_norm_8h_circ))
antisense_sg17m_norm_8h_circ_gtf <- antisense_sg17m_norm_8h_circ_gtf[order(antisense_sg17m_norm_8h_circ_gtf$transcript_id),]
antisense_sg17m_norm_8h_circ <- antisense_sg17m_norm_8h_circ[order(rownames(antisense_sg17m_norm_8h_circ)),]

antisense_sg17m_norm_8h_circ$start <- antisense_sg17m_norm_8h_circ_gtf$start
antisense_sg17m_norm_8h_circ$end <- antisense_sg17m_norm_8h_circ_gtf$end
antisense_sg17m_norm_8h_circ$isolate <- "SG17M"
antisense_sg17m_norm_8h_circ_BR1 <- select(antisense_sg17m_norm_8h_circ, c(isolate,start, end, sg17m_8h_BR1))
color_sg17m_BR1_8h <- with(antisense_sg17m_norm_8h_circ_BR1,
                                             ifelse(sg17m_8h_BR1 >= 12, "palegreen4", 
                                               ifelse(sg17m_8h_BR1 < 12 & sg17m_8h_BR1 >= 10, "palegreen3", 
                                                      ifelse(sg17m_8h_BR1 < 10 & sg17m_8h_BR1 >= 9,
                                                             "palegreen2",
                                                             ifelse(sg17m_8h_BR1 < 9 & sg17m_8h_BR1 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(sg17m_8h_BR1 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_sg17m_norm_8h_circ_BR2 <- select(antisense_sg17m_norm_8h_circ, c(isolate,start, end, sg17m_8h_BR2))
color_sg17m_BR2_8h <- with(antisense_sg17m_norm_8h_circ_BR2,
                                             ifelse(sg17m_8h_BR2 >= 12, "palegreen4", 
                                               ifelse(sg17m_8h_BR2 < 12 & sg17m_8h_BR2 >= 10, "palegreen3", 
                                                      ifelse(sg17m_8h_BR2 < 10 & sg17m_8h_BR2 >= 9,
                                                             "palegreen2",
                                                             ifelse(sg17m_8h_BR2 < 9 & sg17m_8h_BR2 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(sg17m_8h_BR2 < 6, "palegoldenrod",
                                                                           "NA"))))))

antisense_sg17m_norm_8h_circ_BR3 <- select(antisense_sg17m_norm_8h_circ, c(isolate,start, end, sg17m_8h_BR3))
color_sg17m_BR3_8h <- with(antisense_sg17m_norm_8h_circ_BR3,
                                             ifelse(sg17m_8h_BR3 >= 12, "palegreen4", 
                                               ifelse(sg17m_8h_BR3 < 12 & sg17m_8h_BR3 >= 10, "palegreen3", 
                                                      ifelse(sg17m_8h_BR3 < 10 & sg17m_8h_BR3 >= 9,
                                                             "palegreen2",
                                                             ifelse(sg17m_8h_BR3 < 9 & sg17m_8h_BR3 >= 6,
                                                                    "lightgreen",
                                                                    ifelse(sg17m_8h_BR3 < 6, "palegoldenrod",
                                                                           "NA"))))))

tiff(filename = "save_figures/circular_plot_sg17m_8h_800dpi.tif", units="cm", res=800, width=18, height=16)
set.seed(123)
circos.clear()
par(mar=c(0, 0, 0, 0))
circos.par(start.degree = 90,
           cell.padding = c(0, 0, 0, 0))

# Initialize genome (bed file with genome sizes)
genome <- data.frame(chr=c("SG17M"), start = c(1), end = c(7000000))
circos.genomicInitialize(genome, plotType = c("axis"), major.by = 1000000,
                         axis.labels.cex = 0.7, labels.cex = 0.7)

# Add track with annotation
circos.genomicTrack(feature_sg17m, ylim=c(0,1), track.height = 0.1, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, col="gold", border="gold")}, bg.col = "dodgerblue4")

circos.genomicDensity(feature_tRNA_sg17m, count_by = "number", col="black", area=TRUE)

antisense_cov4_sg17m <- data.frame(chr = c(antisense_sg17m_norm_8h_circ_BR1$isolate),
                            start = c(antisense_sg17m_norm_8h_circ_BR1$start),
                            end = c(antisense_sg17m_norm_8h_circ_BR1$end),
                            depth = c(antisense_sg17m_norm_8h_circ_BR1$sg17m_8h_BR1))

circos.genomicTrack(antisense_cov4_sg17m, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_sg17m_BR1_8h, ...)})

antisense_cov5_sg17m <- data.frame(chr = c(antisense_sg17m_norm_8h_circ_BR2$isolate),
                            start = c(antisense_sg17m_norm_8h_circ_BR2$start),
                            end = c(antisense_sg17m_norm_8h_circ_BR2$end),
                            depth = c(antisense_sg17m_norm_8h_circ_BR2$sg17m_8h_BR2))

circos.genomicTrack(antisense_cov5_sg17m, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_sg17m_BR2_8h, ...)})

antisense_cov6_sg17m <- data.frame(chr = c(antisense_sg17m_norm_8h_circ_BR3$isolate),
                            start = c(antisense_sg17m_norm_8h_circ_BR3$start),
                            end = c(antisense_sg17m_norm_8h_circ_BR3$end),
                            depth = c(antisense_sg17m_norm_8h_circ_BR3$sg17m_8h_BR3))

circos.genomicTrack(antisense_cov6_sg17m, track.height = 0.1,
                    stack = TRUE, panel.fun = function(region, value, ...) {
    cex = scale(value[[1]]) / 5
    i = getI(...)
    circos.genomicPoints(region, value, cex = cex, pch = 16, col = color_sg17m_BR3_8h, ...)})
circos.text(6200000,10, "D", cex=1.5, facing = "downward", font=2)
text(0,0, "SG17M-8h", cex=0.8)
dev.off()


# extract regions with most alignments in each replicate
n_features = 3
feature_length = 10000

# NN2, antisense, 4h, three most expressed
nn2_antisense_n10high_list_4h <- list()
nn2_antisense_n10high_4h_1 <- featureCounts_PA_antisense_df_NN2_norm_4h_4[order(
  featureCounts_PA_antisense_df_NN2_norm_4h_4$NN2_4h_BR1, 
  featureCounts_PA_antisense_df_NN2_norm_4h_4$NN2_4h_BR2,
  featureCounts_PA_antisense_df_NN2_norm_4h_4$NN2_4h_BR3, decreasing=TRUE),]
nn2_antisense_n10high_list_4h = list.append(nn2_antisense_n10high_list_4h, rownames(nn2_antisense_n10high_4h_1)[1:n_features])
nn2_antisense_n10high_list_4h <- unlist(nn2_antisense_n10high_list_4h)
featureCounts_PA_antisense_df_NN2_4h_high <- subset(featureCounts_PA_antisense_df_NN2,
                                                    rownames(featureCounts_PA_antisense_df_NN2) %in%
                                                      nn2_antisense_n10high_list_4h)
gtf_nn2_antisense_4h_high <- subset(import_gtf_NN2, 
                                    transcript_id %in% rownames(featureCounts_PA_antisense_df_NN2_4h_high))
gtf_nn2_antisense_4h_high$extract_start <- gtf_nn2_antisense_4h_high$start - feature_length
gtf_nn2_antisense_4h_high$extract_end <- gtf_nn2_antisense_4h_high$end + feature_length


# NN2, antisense, 4h, three percent most expressed
nn2_antisense_3per_list_4h <- list()
nn2_antisense_3per_4h_1 <- featureCounts_PA_antisense_df_NN2_norm_4h_4[order(
  featureCounts_PA_antisense_df_NN2_norm_4h_4$NN2_4h_BR1, 
  featureCounts_PA_antisense_df_NN2_norm_4h_4$NN2_4h_BR2,
  featureCounts_PA_antisense_df_NN2_norm_4h_4$NN2_4h_BR3, decreasing=TRUE),]

features_3per_NN2_4h <- (nrow(nn2_antisense_3per_4h_1) * 3) / 100
nn2_antisense_3per_list_4h = list.append(nn2_antisense_3per_list_4h, rownames(nn2_antisense_3per_4h_1)[1:features_3per_NN2_4h])

nn2_antisense_3per_list_4h <- unlist(nn2_antisense_3per_list_4h)
featureCounts_NN2_4h_3per <- subset(featureCounts_PA_antisense_df_NN2,
                                                    rownames(featureCounts_PA_antisense_df_NN2) %in%
                                                      nn2_antisense_3per_list_4h)

gtf_nn2_antisense_4h_3per <- subset(import_gtf_NN2, transcript_id %in% rownames(featureCounts_NN2_4h_3per))
gtf_nn2_antisense_4h_3per$time <- "4h"

# NN2, antisense, 8h
nn2_antisense_n10high_list_8h <- list()
nn2_antisense_n10high_8h_1 <- featureCounts_PA_antisense_df_NN2_norm_8h_4[order(
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR1, 
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR2,
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR3, decreasing=TRUE),]

features_3per_NN2_8h <- (nrow(nn2_antisense_n10high_8h_1) * 3) / 100

nn2_antisense_n10high_list_8h = list.append(nn2_antisense_n10high_list_8h,
                                            rownames(nn2_antisense_n10high_8h_1)[1:n_features])
nn2_antisense_n10high_list_8h <- unlist(nn2_antisense_n10high_list_8h)
featureCounts_PA_antisense_df_NN2_8h_high <- subset(featureCounts_PA_antisense_df_NN2,
                                                    rownames(featureCounts_PA_antisense_df_NN2) %in%
                                                      nn2_antisense_n10high_list_8h)
gtf_nn2_antisense_8h_high <- subset(import_gtf_NN2, 
                                    transcript_id %in% rownames(featureCounts_PA_antisense_df_NN2_8h_high))
gtf_nn2_antisense_8h_high$extract_start <- gtf_nn2_antisense_8h_high$start - feature_length
gtf_nn2_antisense_8h_high$extract_end <- gtf_nn2_antisense_8h_high$end + feature_length


# NN2, antisense, 8h, three percent most expressed
nn2_antisense_3per_list_8h <- list()
nn2_antisense_3per_8h_1 <- featureCounts_PA_antisense_df_NN2_norm_8h_4[order(
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR1, 
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR2,
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR3, decreasing=TRUE),]

features_3per_NN2_8h <- (nrow(nn2_antisense_3per_8h_1) * 3) / 100
nn2_antisense_3per_list_8h = list.append(nn2_antisense_3per_list_8h, rownames(nn2_antisense_3per_8h_1)[1:features_3per_NN2_8h])

nn2_antisense_3per_list_8h <- unlist(nn2_antisense_3per_list_8h)
featureCounts_NN2_8h_3per <- subset(featureCounts_PA_antisense_df_NN2,
                                                    rownames(featureCounts_PA_antisense_df_NN2) %in%
                                                      nn2_antisense_3per_list_8h)

gtf_nn2_antisense_8h_3per <- subset(import_gtf_NN2, transcript_id %in% rownames(featureCounts_NN2_8h_3per))
gtf_nn2_antisense_8h_3per$time <- "8h"

gtf_nn2_antisense_3per <- data.frame(rbind(gtf_nn2_antisense_4h_3per, gtf_nn2_antisense_8h_3per))
gtf_nn2_antisense_3per$gene_name3 <- NULL
gtf_nn2_antisense_3per_unknown <- subset(gtf_nn2_antisense_3per, gene_name2 == "")
gtf_nn2_antisense_3per_unknown <- select(gtf_nn2_antisense_3per_unknown, c(Isolate, start, end))

#write.csv(gtf_nn2_antisense_3per_unknown, file = "nn2_antisense_threePercentMost_hypothetical.csv",
 #         row.names = FALSE)
#write.csv(gtf_nn2_antisense_3per, file = "nn2_antisense_threePercentMost.csv",
        #  row.names = FALSE)

# SG17M, antisense, 4h
sg17m_antisense_n10high_list_4h <- list()
sg17m_antisense_n10high_4h_1 <- featureCounts_PA_antisense_df_sg17m_norm_4h_4[order(
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR1, 
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR2,
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR3, decreasing=TRUE),]

features_3per_sg17m_4h <- (nrow(sg17m_antisense_n10high_4h_1) * 3) / 100

sg17m_antisense_n10high_list_4h = list.append(sg17m_antisense_n10high_list_4h, rownames(sg17m_antisense_n10high_4h_1)[1:n_features])
sg17m_antisense_n10high_list_4h <- unlist(sg17m_antisense_n10high_list_4h)
featureCounts_PA_antisense_df_sg17m_4h_high <- subset(featureCounts_PA_antisense_df_sg17m,
                                                    rownames(featureCounts_PA_antisense_df_sg17m) %in%
                                                      sg17m_antisense_n10high_list_4h)
gtf_sg17m_antisense_4h_high <- subset(import_gtf_sg17m, 
                                    transcript_id %in% rownames(featureCounts_PA_antisense_df_sg17m_4h_high))
gtf_sg17m_antisense_4h_high$extract_start <- gtf_sg17m_antisense_4h_high$start - feature_length
gtf_sg17m_antisense_4h_high$extract_end <- gtf_sg17m_antisense_4h_high$end + feature_length


# SG17M, antisense, 4h, three percent most expressed
sg17m_antisense_3per_list_4h <- list()
sg17m_antisense_3per_4h_1 <- featureCounts_PA_antisense_df_sg17m_norm_4h_4[order(
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR1, 
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR2,
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR3, decreasing=TRUE),]

features_3per_sg17m_4h <- (nrow(sg17m_antisense_3per_4h_1) * 3) / 100
sg17m_antisense_3per_list_4h = list.append(sg17m_antisense_3per_list_4h, rownames(sg17m_antisense_3per_4h_1)[1:features_3per_sg17m_4h])

sg17m_antisense_3per_list_4h <- unlist(sg17m_antisense_3per_list_4h)
featureCounts_sg17m_4h_3per <- subset(featureCounts_PA_antisense_df_sg17m,
                                                    rownames(featureCounts_PA_antisense_df_sg17m) %in%
                                                      sg17m_antisense_3per_list_4h)

gtf_sg17m_antisense_4h_3per <- subset(import_gtf_sg17m, transcript_id %in% rownames(featureCounts_sg17m_4h_3per))
gtf_sg17m_antisense_4h_3per$time <- "4h"

# SG17M, antisense, 8h
sg17m_antisense_n10high_list_8h <- list()
sg17m_antisense_n10high_8h_1 <- featureCounts_PA_antisense_df_sg17m_norm_8h_4[order(
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR1, 
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR2,
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR3, decreasing=TRUE),]

features_3per_SG17M_8h <- (nrow(sg17m_antisense_n10high_8h_1) * 3) / 100

sg17m_antisense_n10high_list_8h = list.append(sg17m_antisense_n10high_list_8h,
                                            rownames(sg17m_antisense_n10high_8h_1)[1:n_features])
sg17m_antisense_n10high_list_8h <- unlist(sg17m_antisense_n10high_list_8h)
featureCounts_PA_antisense_df_sg17m_8h_high <- subset(featureCounts_PA_antisense_df_sg17m,
                                                    rownames(featureCounts_PA_antisense_df_sg17m) %in%
                                                      sg17m_antisense_n10high_list_8h)
gtf_sg17m_antisense_8h_high <- subset(import_gtf_sg17m, 
                                    transcript_id %in% rownames(featureCounts_PA_antisense_df_sg17m_8h_high))
gtf_sg17m_antisense_8h_high$extract_start <- gtf_sg17m_antisense_8h_high$start - feature_length
gtf_sg17m_antisense_8h_high$extract_end <- gtf_sg17m_antisense_8h_high$end + feature_length


# SG17M, antisense, 8h, three percent most expressed
sg17m_antisense_3per_list_8h <- list()
sg17m_antisense_3per_8h_1 <- featureCounts_PA_antisense_df_sg17m_norm_8h_4[order(
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR1, 
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR2,
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR3, decreasing=TRUE),]

features_3per_sg17m_8h <- (nrow(sg17m_antisense_3per_8h_1) * 3) / 100
sg17m_antisense_3per_list_8h = list.append(sg17m_antisense_3per_list_8h, 
                                           rownames(sg17m_antisense_3per_8h_1)[1:features_3per_sg17m_8h])

sg17m_antisense_3per_list_8h <- unlist(sg17m_antisense_3per_list_8h)
featureCounts_sg17m_8h_3per <- subset(featureCounts_PA_antisense_df_sg17m, 
                                      rownames(featureCounts_PA_antisense_df_sg17m) %in% sg17m_antisense_3per_list_8h)
gtf_sg17m_antisense_8h_3per <- subset(import_gtf_sg17m, transcript_id %in% rownames(featureCounts_sg17m_8h_3per))
gtf_sg17m_antisense_8h_3per$time <- "8h"


gtf_sg17m_antisense_3per <- data.frame(rbind(gtf_sg17m_antisense_4h_3per, gtf_sg17m_antisense_8h_3per))
gtf_sg17m_antisense_3per$gene_name3 <- NULL
gtf_sg17m_antisense_3per_unknown <- subset(gtf_sg17m_antisense_3per, gene_name2 == "")
gtf_sg17m_antisense_3per_unknown <- select(gtf_sg17m_antisense_3per_unknown, c(Isolate, start, end))

#write.csv(gtf_sg17m_antisense_3per_unknown, file = "sg17m_antisense_threePercentMost_hypothetical.csv",
  #        row.names = FALSE)
#write.csv(gtf_sg17m_antisense_3per, file = "sg17m_antisense_threePercentMost.csv",
  #        row.names = FALSE)

# NN2 ####
NN2_4h_BR1_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR1.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR1_depth_fa <- data.frame(NN2_4h_BR1_depth_fa)
NN2_4h_BR1_depth_fa$side <- "forward_antisense"
NN2_4h_BR1_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR1.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR1_depth_ra <- data.frame(NN2_4h_BR1_depth_ra)
NN2_4h_BR1_depth_ra$side <- "reverse_antisense"
NN2_4h_BR1_depth <- data.frame(rbind(NN2_4h_BR1_depth_fa, NN2_4h_BR1_depth_ra))
NN2_4h_BR1_depth$X6 <- NULL
NN2_4h_BR1_depth$replicate <- "BR1"


NN2_4h_BR1_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR1.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR1_depth_fs <- data.frame(NN2_4h_BR1_depth_fs)
NN2_4h_BR1_depth_fs$side <- "forward_sense"
NN2_4h_BR1_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR1.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR1_depth_rs <- data.frame(NN2_4h_BR1_depth_rs)
NN2_4h_BR1_depth_rs$side <- "reverse_sense"
NN2_4h_BR1_depth_sense <- data.frame(rbind(NN2_4h_BR1_depth_fs, NN2_4h_BR1_depth_rs))
NN2_4h_BR1_depth_sense$X6 <- NULL
NN2_4h_BR1_depth_sense$replicate <- "BR1"




NN2_4h_BR2_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR2.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR2_depth_fa <- data.frame(NN2_4h_BR2_depth_fa)
NN2_4h_BR2_depth_fa$side <- "forward_antisense"
NN2_4h_BR2_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR2.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR2_depth_ra <- data.frame(NN2_4h_BR2_depth_ra)
NN2_4h_BR2_depth_ra$side <- "reverse_antisense"
NN2_4h_BR2_depth <- data.frame(rbind(NN2_4h_BR2_depth_fa, NN2_4h_BR2_depth_ra))
NN2_4h_BR2_depth$X6 <- NULL
NN2_4h_BR2_depth$replicate <- "BR2"



NN2_4h_BR2_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR2.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR2_depth_fs <- data.frame(NN2_4h_BR2_depth_fs)
NN2_4h_BR2_depth_fs$side <- "forward_sense"
NN2_4h_BR2_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR2.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR2_depth_rs <- data.frame(NN2_4h_BR2_depth_rs)
NN2_4h_BR2_depth_rs$side <- "reverse_sense"
NN2_4h_BR2_depth_sense <- data.frame(rbind(NN2_4h_BR2_depth_fs, NN2_4h_BR2_depth_rs))
NN2_4h_BR2_depth_sense$X6 <- NULL
NN2_4h_BR2_depth_sense$replicate <- "BR2"



NN2_4h_BR3_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR3.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR3_depth_fa <- data.frame(NN2_4h_BR3_depth_fa)
NN2_4h_BR3_depth_fa$side <- "forward_antisense"
NN2_4h_BR3_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR3.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR3_depth_ra <- data.frame(NN2_4h_BR3_depth_ra)
NN2_4h_BR3_depth_ra$side <- "reverse_antisense"
NN2_4h_BR3_depth <- data.frame(rbind(NN2_4h_BR3_depth_fa, NN2_4h_BR3_depth_ra))
NN2_4h_BR3_depth$X6 <- NULL
NN2_4h_BR3_depth$replicate <- "BR3"


NN2_4h_BR3_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR3.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR3_depth_fs <- data.frame(NN2_4h_BR3_depth_fs)
NN2_4h_BR3_depth_fs$side <- "forward_sense"
NN2_4h_BR3_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_4h_BR3.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR3_depth_rs <- data.frame(NN2_4h_BR3_depth_rs)
NN2_4h_BR3_depth_rs$side <- "reverse_sense"
NN2_4h_BR3_depth_sense <- data.frame(rbind(NN2_4h_BR3_depth_fs, NN2_4h_BR3_depth_rs))
NN2_4h_BR3_depth_sense$X6 <- NULL
NN2_4h_BR3_depth_sense$replicate <- "BR3"

NN2_4h_depth <- data.frame(rbind(NN2_4h_BR1_depth, NN2_4h_BR2_depth, NN2_4h_BR3_depth))
colnames(NN2_4h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

NN2_4h_depth_sense <- data.frame(rbind(NN2_4h_BR1_depth_sense, NN2_4h_BR2_depth_sense, NN2_4h_BR3_depth_sense))
colnames(NN2_4h_depth_sense) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

NN2_8h_BR1_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR1.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR1_depth_fa <- data.frame(NN2_8h_BR1_depth_fa)
NN2_8h_BR1_depth_fa$side <- "forward_antisense"
NN2_8h_BR1_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR1.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR1_depth_ra <- data.frame(NN2_8h_BR1_depth_ra)
NN2_8h_BR1_depth_ra$side <- "reverse_antisense"
NN2_8h_BR1_depth <- data.frame(rbind(NN2_8h_BR1_depth_fa, NN2_8h_BR1_depth_ra))
NN2_8h_BR1_depth$X6 <- NULL
NN2_8h_BR1_depth$replicate <- "BR1"



NN2_8h_BR1_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR1.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR1_depth_fs <- data.frame(NN2_8h_BR1_depth_fs)
NN2_8h_BR1_depth_fs$side <- "forward_sense"
NN2_8h_BR1_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR1.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR1_depth_rs <- data.frame(NN2_8h_BR1_depth_rs)
NN2_8h_BR1_depth_rs$side <- "reverse_sense"
NN2_8h_BR1_depth_sense <- data.frame(rbind(NN2_8h_BR1_depth_fs, NN2_8h_BR1_depth_rs))
NN2_8h_BR1_depth_sense$X6 <- NULL
NN2_8h_BR1_depth_sense$replicate <- "BR1"


NN2_8h_BR2_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR2.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR2_depth_fa <- data.frame(NN2_8h_BR2_depth_fa)
NN2_8h_BR2_depth_fa$side <- "forward_antisense"
NN2_8h_BR2_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR2.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR2_depth_ra <- data.frame(NN2_8h_BR2_depth_ra)
NN2_8h_BR2_depth_ra$side <- "reverse_antisense"
NN2_8h_BR2_depth <- data.frame(rbind(NN2_8h_BR2_depth_fa, NN2_8h_BR2_depth_ra))
NN2_8h_BR2_depth$X6 <- NULL
NN2_8h_BR2_depth$replicate <- "BR2"


NN2_8h_BR2_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR2.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR2_depth_fs <- data.frame(NN2_8h_BR2_depth_fs)
NN2_8h_BR2_depth_fs$side <- "forward_sense"
NN2_8h_BR2_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR2.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR2_depth_rs <- data.frame(NN2_8h_BR2_depth_rs)
NN2_8h_BR2_depth_rs$side <- "reverse_sense"
NN2_8h_BR2_depth_sense <- data.frame(rbind(NN2_8h_BR2_depth_fs, NN2_8h_BR2_depth_rs))
NN2_8h_BR2_depth_sense$X6 <- NULL
NN2_8h_BR2_depth_sense$replicate <- "BR2"


NN2_8h_BR3_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR3.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR3_depth_fa <- data.frame(NN2_8h_BR3_depth_fa)
NN2_8h_BR3_depth_fa$side <- "forward_antisense"
NN2_8h_BR3_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR3.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR3_depth_ra <- data.frame(NN2_8h_BR3_depth_ra)
NN2_8h_BR3_depth_ra$side <- "reverse_antisense"
NN2_8h_BR3_depth <- data.frame(rbind(NN2_8h_BR3_depth_fa, NN2_8h_BR3_depth_ra))
NN2_8h_BR3_depth$X6 <- NULL
NN2_8h_BR3_depth$replicate <- "BR3"

NN2_8h_BR3_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR3.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR3_depth_fs <- data.frame(NN2_8h_BR3_depth_fs)
NN2_8h_BR3_depth_fs$side <- "forward_sense"
NN2_8h_BR3_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/test/NN2_CS_8h_BR3.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR3_depth_rs <- data.frame(NN2_8h_BR3_depth_rs)
NN2_8h_BR3_depth_rs$side <- "reverse_sense"
NN2_8h_BR3_depth_sense <- data.frame(rbind(NN2_8h_BR3_depth_fs, NN2_8h_BR3_depth_rs))
NN2_8h_BR3_depth_sense$X6 <- NULL
NN2_8h_BR3_depth_sense$replicate <- "BR3"

NN2_8h_depth <- data.frame(rbind(NN2_8h_BR1_depth, NN2_8h_BR2_depth, NN2_8h_BR3_depth))
colnames(NN2_8h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

NN2_8h_depth_sense <- data.frame(rbind(NN2_8h_BR1_depth_sense, NN2_8h_BR2_depth_sense, NN2_8h_BR3_depth_sense))
colnames(NN2_8h_depth_sense) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")


# Antisense hotspot 1
NN2_4h_depth_hsp1 <- subset(NN2_4h_depth, start > 764600 & end < 768903)
NN2_4h_gtf_hsp1 <- subset(import_gtf_NN2, start > 764600 & end < 768903)
NN2_4h_depth_hsp1 <- NN2_4h_depth_hsp1[order(NN2_4h_depth_hsp1$replicate, NN2_4h_depth_hsp1$start),]
rownames(NN2_4h_depth_hsp1) <- NULL
NN2_4h_depth_hsp1$num <- rownames(NN2_4h_depth_hsp1)
NN2_4h_depth_hsp1$num <- as.numeric(as.character(NN2_4h_depth_hsp1$num))
NN2_4h_gtf_hsp1$gene_name3 <- c("tRNA", "tRNA", "tRNA", "CDS", "tRNA","unknown", "CDS", "CDS", "CDS")

nn2_hsp1_plot <- ggplot() + geom_segment(data=NN2_4h_depth_hsp1, 
                                         aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "darkorange2", "blue", "forestgreen")) +
  xlab("\n") + ylab("Read IDs (NN2-4h antisense transcripts)") + 
  geom_segment(data=NN2_4h_gtf_hsp1, aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-10, yend=-10, color=gene_name3), arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-18, x = 764865), label="Tyr-Anticodon-GTA", size=2.5, color="blue") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-29, x = 764896), label="Gly-Anticodon-TCC", size=2.5, color="blue") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-40, x = 765379), label="Thr-Anticodon-GGT", size=2.5, color="blue") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-18, x = 765929), label="tufA", size=2.5, color="darkorange2") +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-29, x = 766518), label="Trp-Anticodon-CCA", size=2.5, color="blue") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-40, x = 766990), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-18, x = 767190), label="nusG", size=2.5, color="darkorange2") +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-18, x = 767827), label="rplK", size=2.5, color="darkorange2") +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-18, x = 768508), label="rplA", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(764300, 768903))


# Antisense hotspot 2
NN2_4h_depth_hsp2 <- subset(NN2_4h_depth, start > 5577511 & end < 5578609)
NN2_4h_gtf_hsp2 <- subset(import_gtf_NN2, start > 5577511 & end < 5578609)

NN2_4h_depth_hsp2 <- NN2_4h_depth_hsp2[order(NN2_4h_depth_hsp2$replicate, NN2_4h_depth_hsp2$start),]
rownames(NN2_4h_depth_hsp2) <- NULL
NN2_4h_depth_hsp2$num <- rownames(NN2_4h_depth_hsp2)
NN2_4h_depth_hsp2$num <- as.numeric(as.character(NN2_4h_depth_hsp2$num))

nn2_hsp2_plot <- ggplot() + geom_segment(data=NN2_4h_depth_hsp2, aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3","darkorange2", "blue")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_4h_gtf_hsp2, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-6, yend=-6, color=gene_name3), arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=NN2_4h_gtf_hsp2, aes(y=-18, x = 5577760), label="Lys-Anticodon-TTT (pKLC102)", size=2.5, color="blue") + 
  geom_label(data=NN2_4h_gtf_hsp2, aes(y=-18, x = 5578100), label="soj", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(5577511, 5578609))


# Antisense hotspot 3
NN2_4h_depth_hsp3 <- subset(NN2_4h_depth, start > 5600220 & end < 5606993)
NN2_4h_gtf_hsp3 <- subset(import_gtf_NN2, start > 5600220 & end < 5606993)
NN2_4h_depth_hsp3$length <- NN2_4h_depth_hsp3$end - NN2_4h_depth_hsp3$start

NN2_4h_gtf_hsp3$gene_name3 <- c("CDS (hypothetical)", "rapA", "CDS (hypothetical)",
                                "CDS (hypothetical)", "CDS (hypothetical)")
NN2_4h_depth_hsp3 <- NN2_4h_depth_hsp3[order(NN2_4h_depth_hsp3$replicate, NN2_4h_depth_hsp3$start),]
rownames(NN2_4h_depth_hsp3) <- NULL
NN2_4h_depth_hsp3$num <- rownames(NN2_4h_depth_hsp3)
NN2_4h_depth_hsp3$num <- as.numeric(as.character(NN2_4h_depth_hsp3$num))

nn2_hsp3_plot <- ggplot() + geom_segment(data=NN2_4h_depth_hsp3, 
                                         aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "forestgreen", 
                              "darkorange2", "forestgreen", "forestgreen", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_4h_gtf_hsp3, aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)), y=-10, yend=-10, color=gene_name3), 
               arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=NN2_4h_gtf_hsp3, aes(y=-65, x = 5601720), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=NN2_4h_gtf_hsp3, aes(y=-65, x = 5603720), label="rapA", size=2.5, color="darkorange2") + 
  geom_label(data=NN2_4h_gtf_hsp3, aes(y=-65, x = 5605820), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=NN2_4h_gtf_hsp3, aes(y=-46, x = 5605920), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=NN2_4h_gtf_hsp3, aes(y=-25, x = 5606020), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(5600220, 5606993))

# Antisense hotspot 4
NN2_4h_depth_hsp4 <- subset(NN2_4h_depth, start > 325534 & end < 332094)
NN2_4h_gtf_hsp4 <- subset(import_gtf_NN2, start > 325534 & end < 332094)
NN2_4h_depth_hsp4$length <- NN2_4h_depth_hsp4$end - NN2_4h_depth_hsp4$start

NN2_4h_depth_hsp4 <- NN2_4h_depth_hsp4[order(NN2_4h_depth_hsp4$replicate, NN2_4h_depth_hsp4$start),]
rownames(NN2_4h_depth_hsp4) <- NULL
NN2_4h_depth_hsp4$num <- rownames(NN2_4h_depth_hsp4)
NN2_4h_depth_hsp4$num <- as.numeric(as.character(NN2_4h_depth_hsp4$num))
NN2_4h_gtf_hsp4$gene_name4 <- c("CDSu", "CDSu", "CDS", "CDSu", "CDS")

nn2_hsp4_plot <- ggplot() + geom_segment(data=NN2_4h_depth_hsp4, 
                                         aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "darkorange2", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_4h_gtf_hsp4, aes(x=ifelse(strandType == "+", start, 
                                                  ifelse(strandType == "\\.", start, end)), 
                                         xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                                         y=-8, yend=-8, color=gene_name4), arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=NN2_4h_gtf_hsp4, aes(y=-14, x = 327590), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=NN2_4h_gtf_hsp4, aes(y=-31, x = 327908), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=NN2_4h_gtf_hsp4, aes(y=-14, x = 328704), label="intQ", size=2.5, color="darkorange2") +
  geom_label(data=NN2_4h_gtf_hsp4, aes(y=-31, x = 329606), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=NN2_4h_gtf_hsp4, aes(y=-14, x = 330415), label="qseF", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(326934, 331794))

nn2_antisense_hsp_4h <- ggarrange(nn2_hsp1_plot, nn2_hsp2_plot, nn2_hsp4_plot, nrow = 1, labels=c("A", "B", "C"))


# Antisense hotspot 1
NN2_8h_depth_hsp1 <- subset(NN2_8h_depth, start > 764600 & end < 768903)
NN2_8h_gtf_hsp1 <- subset(import_gtf_NN2, start > 764600 & end < 768903)
NN2_8h_depth_hsp1 <- NN2_8h_depth_hsp1[order(NN2_8h_depth_hsp1$replicate, NN2_8h_depth_hsp1$start),]
rownames(NN2_8h_depth_hsp1) <- NULL
NN2_8h_depth_hsp1$num <- rownames(NN2_8h_depth_hsp1)
NN2_8h_depth_hsp1$num <- as.numeric(as.character(NN2_8h_depth_hsp1$num))
NN2_8h_gtf_hsp1$gene_name3 <- c("tRNA", "tRNA", "tRNA", "CDS", "tRNA","unknown", "CDS", "CDS", "CDS")

nn2_hsp1_plot_8h <- ggplot() + geom_segment(data=NN2_8h_depth_hsp1, aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3",  "darkorange2", "blue", "forestgreen")) +
  xlab("\n") + ylab("Read IDs (NN2-8h antisense transcripts)") + 
  geom_segment(data=NN2_8h_gtf_hsp1, aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-3, yend=-3, color=gene_name3), arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-8, x = 764865), label="Tyr-Anticodon-GTA", size=2.5, color="blue") + 
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-14, x = 764896), label="Gly-Anticodon-TCC", size=2.5, color="blue") + 
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-20, x = 765379), label="Thr-Anticodon-GGT", size=2.5, color="blue") + 
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-8, x = 765929), label="tufA", size=2.5, color="darkorange2") +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-14, x = 766518), label="Trp-Anticodon-CCA", size=2.5, color="blue") + 
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-20, x = 766990), label="CDS (hypothetical)", size=2, color="forestgreen") +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-8, x = 767190), label="nusG", size=2.5, color="darkorange2") +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-8, x = 767827), label="rplK", size=2.5, color="darkorange2") +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-8, x = 768508), label="rplA", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(label=scales::comma, limits=c(764300, 768903))

# Antisense hotspot 2
NN2_8h_depth_hsp2 <- subset(NN2_8h_depth, start > 5577511 & end < 5578609)
NN2_8h_gtf_hsp2 <- subset(import_gtf_NN2, start > 5577511 & end < 5578609)
NN2_8h_depth_hsp2 <- NN2_8h_depth_hsp2[order(NN2_8h_depth_hsp2$replicate, NN2_8h_depth_hsp2$start),]
rownames(NN2_8h_depth_hsp2) <- NULL
NN2_8h_depth_hsp2$num <- rownames(NN2_8h_depth_hsp2)
NN2_8h_depth_hsp2$num <- as.numeric(as.character(NN2_8h_depth_hsp2$num))

nn2_hsp2_plot_8h <- ggplot() + geom_segment(data=NN2_8h_depth_hsp2, 
                                            aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "darkorange2", "blue")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_8h_gtf_hsp2, aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-5, yend=-5, color=gene_name3), arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=NN2_8h_gtf_hsp2, aes(y=-23, x = 5577760), label="Lys-Anticodon-TTT (pKLC102)", size=2.5, color="blue") + 
  geom_label(data=NN2_8h_gtf_hsp2, aes(y=-23, x = 5578100), label="soj", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(label=scales::comma, limits=c(5577511, 5578609))


# Antisense hotspot 3
NN2_8h_depth_hsp3 <- subset(NN2_8h_depth, start > 5600220 & end < 5606993)
NN2_8h_gtf_hsp3 <- subset(import_gtf_NN2, start > 5600220 & end < 5606993)
NN2_8h_depth_hsp3$length <- NN2_8h_depth_hsp3$end - NN2_8h_depth_hsp3$start

NN2_8h_gtf_hsp3$gene_name3 <- c("CDS (hypothetical)", "rapA", "CDS (hypothetical)",
                                "CDS (hypothetical)", "CDS (hypothetical)")
NN2_8h_depth_hsp3 <- NN2_8h_depth_hsp3[order(NN2_8h_depth_hsp3$replicate, NN2_8h_depth_hsp3$start),]
rownames(NN2_8h_depth_hsp3) <- NULL
NN2_8h_depth_hsp3$num <- rownames(NN2_8h_depth_hsp3)
NN2_8h_depth_hsp3$num <- as.numeric(as.character(NN2_8h_depth_hsp3$num))

nn2_hsp3_plot_8h <- ggplot() + geom_segment(data=NN2_8h_depth_hsp3, 
                                            aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "forestgreen", 
                              "darkorange2", "forestgreen", "forestgreen", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_8h_gtf_hsp3, aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-10, yend=-10, color=gene_name3), arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-69, x = 5601720), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-69, x = 5603720), label="rapA", size=2.5, color="darkorange2") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-69, x = 5605820), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-50, x = 5605920), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-29, x = 5606020), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(label=scales::comma, limits=c(5600220, 5606993))


# Antisense hotspot 4
NN2_8h_depth_hsp4 <- subset(NN2_8h_depth, start > 325534 & end < 332094)
NN2_8h_gtf_hsp4 <- subset(import_gtf_NN2, start > 325534 & end < 332094)
NN2_8h_depth_hsp4$length <- NN2_8h_depth_hsp4$end - NN2_8h_depth_hsp4$start

NN2_8h_depth_hsp4 <- NN2_8h_depth_hsp4[order(NN2_8h_depth_hsp4$replicate, NN2_8h_depth_hsp4$start),]
rownames(NN2_8h_depth_hsp4) <- NULL
NN2_8h_depth_hsp4$num <- rownames(NN2_8h_depth_hsp4)
NN2_8h_depth_hsp4$num <- as.numeric(as.character(NN2_8h_depth_hsp4$num))
NN2_8h_gtf_hsp4$gene_name4 <- c("CDSu", "CDSu", "CDS", "CDSu", "CDS")

nn2_hsp4_plot_8h <- ggplot() + geom_segment(data=NN2_8h_depth_hsp4, 
                                            aes(x=start, xend=end, y=num, yend=num, colour=replicate), 
                                            size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "darkorange2", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_8h_gtf_hsp4, aes(x=ifelse(strandType == "+", start, 
                            ifelse(strandType == "\\.", start, end)), xend=ifelse(strandType == "+", end, 
                               ifelse(strandType == "\\.", end, start)), y=-12, yend=-12, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=NN2_8h_gtf_hsp4, aes(y=-23, x = 327590), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=NN2_8h_gtf_hsp4, aes(y=-54, x = 327908), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=NN2_8h_gtf_hsp4, aes(y=-23, x = 328704), label="intQ", size=2.5, color="darkorange2") +
  geom_label(data=NN2_8h_gtf_hsp4, aes(y=-54, x = 329606), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=NN2_8h_gtf_hsp4, aes(y=-23, x = 330415), label="qseF", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(label=scales::comma, limits=c(326934, 331794))


nn2_antisense_hsp_8h <-
  ggarrange(nn2_hsp1_plot_8h, nn2_hsp2_plot_8h, nn2_hsp4_plot_8h, nrow = 1,
          labels=c("D", "E", "F"))

antisense_hotspots_nn2 <- ggarrange(nn2_antisense_hsp_4h, nn2_antisense_hsp_8h, nrow=2)

ggsave(antisense_hotspots_nn2, filename="save_figures/Figure_06.tif",
      dpi=300, device="tiff", units="cm", width=28,height=23)









# SG17M ####
SG17M_4h_BR1_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR1.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR1_depth_fa <- data.frame(SG17M_4h_BR1_depth_fa)
SG17M_4h_BR1_depth_fa$side <- "forward_antisense"
SG17M_4h_BR1_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR1.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR1_depth_ra <- data.frame(SG17M_4h_BR1_depth_ra)
SG17M_4h_BR1_depth_ra$side <- "reverse_antisense"
SG17M_4h_BR1_depth <- data.frame(rbind(SG17M_4h_BR1_depth_fa, SG17M_4h_BR1_depth_ra))
SG17M_4h_BR1_depth$X6 <- NULL
SG17M_4h_BR1_depth$replicate <- "BR1"


SG17M_4h_BR1_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR1.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR1_depth_fs <- data.frame(SG17M_4h_BR1_depth_fs)
SG17M_4h_BR1_depth_fs$side <- "forward_sense"
SG17M_4h_BR1_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR1.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR1_depth_rs <- data.frame(SG17M_4h_BR1_depth_rs)
SG17M_4h_BR1_depth_rs$side <- "reverse_sense"
SG17M_4h_BR1_depth_sense <- data.frame(rbind(SG17M_4h_BR1_depth_fs, SG17M_4h_BR1_depth_rs))
SG17M_4h_BR1_depth_sense$X6 <- NULL
SG17M_4h_BR1_depth_sense$replicate <- "BR1"


SG17M_4h_BR2_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR2.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR2_depth_fa <- data.frame(SG17M_4h_BR2_depth_fa)
SG17M_4h_BR2_depth_fa$side <- "forward_antisense"
SG17M_4h_BR2_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR2.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR2_depth_ra <- data.frame(SG17M_4h_BR2_depth_ra)
SG17M_4h_BR2_depth_ra$side <- "reverse_antisense"
SG17M_4h_BR2_depth <- data.frame(rbind(SG17M_4h_BR2_depth_fa, SG17M_4h_BR2_depth_ra))
SG17M_4h_BR2_depth$X6 <- NULL
SG17M_4h_BR2_depth$replicate <- "BR2"


SG17M_4h_BR2_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR2.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR2_depth_fs <- data.frame(SG17M_4h_BR2_depth_fs)
SG17M_4h_BR2_depth_fs$side <- "forward_sense"
SG17M_4h_BR2_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR2.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR2_depth_rs <- data.frame(SG17M_4h_BR2_depth_rs)
SG17M_4h_BR2_depth_rs$side <- "reverse_sense"
SG17M_4h_BR2_depth_sense <- data.frame(rbind(SG17M_4h_BR2_depth_fs, SG17M_4h_BR2_depth_rs))
SG17M_4h_BR2_depth_sense$X6 <- NULL
SG17M_4h_BR2_depth_sense$replicate <- "BR2"


SG17M_4h_BR3_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR3.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR3_depth_fa <- data.frame(SG17M_4h_BR3_depth_fa)
SG17M_4h_BR3_depth_fa$side <- "forward_antisense"
SG17M_4h_BR3_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR3.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR3_depth_ra <- data.frame(SG17M_4h_BR3_depth_ra)
SG17M_4h_BR3_depth_ra$side <- "reverse_antisense"
SG17M_4h_BR3_depth <- data.frame(rbind(SG17M_4h_BR3_depth_fa, SG17M_4h_BR3_depth_ra))
SG17M_4h_BR3_depth$X6 <- NULL
SG17M_4h_BR3_depth$replicate <- "BR3"


SG17M_4h_BR3_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR3.sorted.forward_sense.bed",  "\t", 
                                    escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR3_depth_fs <- data.frame(SG17M_4h_BR3_depth_fs)
SG17M_4h_BR3_depth_fs$side <- "forward_sense"
SG17M_4h_BR3_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_4h_BR3.sorted.reverse_antisense.bed",  "\t", 
                                    escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR3_depth_rs <- data.frame(SG17M_4h_BR3_depth_rs)
SG17M_4h_BR3_depth_rs$side <- "reverse_sense"
SG17M_4h_BR3_depth_sense <- data.frame(rbind(SG17M_4h_BR3_depth_fs, SG17M_4h_BR3_depth_rs))
SG17M_4h_BR3_depth_sense$X6 <- NULL
SG17M_4h_BR3_depth_sense$replicate <- "BR3"

SG17M_4h_depth <- data.frame(rbind(SG17M_4h_BR1_depth, SG17M_4h_BR2_depth, SG17M_4h_BR3_depth))
colnames(SG17M_4h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

SG17M_4h_depth_sense <- data.frame(rbind(SG17M_4h_BR1_depth_sense, SG17M_4h_BR2_depth_sense, SG17M_4h_BR3_depth_sense))
colnames(SG17M_4h_depth_sense) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

SG17M_8h_BR1_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR1.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR1_depth_fa <- data.frame(SG17M_8h_BR1_depth_fa)
SG17M_8h_BR1_depth_fa$side <- "forward_antisense"
SG17M_8h_BR1_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR1.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR1_depth_ra <- data.frame(SG17M_8h_BR1_depth_ra)
SG17M_8h_BR1_depth_ra$side <- "reverse_antisense"
SG17M_8h_BR1_depth <- data.frame(rbind(SG17M_8h_BR1_depth_fa, SG17M_8h_BR1_depth_ra))
SG17M_8h_BR1_depth$X6 <- NULL
SG17M_8h_BR1_depth$replicate <- "BR1"


SG17M_8h_BR1_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR1.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR1_depth_fs <- data.frame(SG17M_8h_BR1_depth_fs)
SG17M_8h_BR1_depth_fs$side <- "forward_sense"
SG17M_8h_BR1_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR1.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR1_depth_rs <- data.frame(SG17M_8h_BR1_depth_rs)
SG17M_8h_BR1_depth_rs$side <- "reverse_sense"
SG17M_8h_BR1_depth_sense <- data.frame(rbind(SG17M_8h_BR1_depth_fs, SG17M_8h_BR1_depth_rs))
SG17M_8h_BR1_depth_sense$X6 <- NULL
SG17M_8h_BR1_depth_sense$replicate <- "BR1"


SG17M_8h_BR2_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR2.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR2_depth_fa <- data.frame(SG17M_8h_BR2_depth_fa)
SG17M_8h_BR2_depth_fa$side <- "forward_antisense"
SG17M_8h_BR2_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR2.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR2_depth_ra <- data.frame(SG17M_8h_BR2_depth_ra)
SG17M_8h_BR2_depth_ra$side <- "reverse_antisense"
SG17M_8h_BR2_depth <- data.frame(rbind(SG17M_8h_BR2_depth_fa, SG17M_8h_BR2_depth_ra))
SG17M_8h_BR2_depth$X6 <- NULL
SG17M_8h_BR2_depth$replicate <- "BR2"


SG17M_8h_BR2_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR2.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR2_depth_fs <- data.frame(SG17M_8h_BR2_depth_fs)
SG17M_8h_BR2_depth_fs$side <- "forward_sense"
SG17M_8h_BR2_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR2.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR2_depth_rs <- data.frame(SG17M_8h_BR2_depth_rs)
SG17M_8h_BR2_depth_rs$side <- "reverse_sense"
SG17M_8h_BR2_depth_sense <- data.frame(rbind(SG17M_8h_BR2_depth_fs, SG17M_8h_BR2_depth_rs))
SG17M_8h_BR2_depth_sense$X6 <- NULL
SG17M_8h_BR2_depth_sense$replicate <- "BR2"


SG17M_8h_BR3_depth_fa <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR3.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR3_depth_fa <- data.frame(SG17M_8h_BR3_depth_fa)
SG17M_8h_BR3_depth_fa$side <- "forward_antisense"
SG17M_8h_BR3_depth_ra <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR3.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR3_depth_ra <- data.frame(SG17M_8h_BR3_depth_ra)
SG17M_8h_BR3_depth_ra$side <- "reverse_antisense"
SG17M_8h_BR3_depth <- data.frame(rbind(SG17M_8h_BR3_depth_fa, SG17M_8h_BR3_depth_ra))
SG17M_8h_BR3_depth$X6 <- NULL
SG17M_8h_BR3_depth$replicate <- "BR3"


SG17M_8h_BR3_depth_fs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR3.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR3_depth_fs <- data.frame(SG17M_8h_BR3_depth_fs)
SG17M_8h_BR3_depth_fs$side <- "forward_sense"
SG17M_8h_BR3_depth_rs <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/test/SG17M_CS_8h_BR3.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR3_depth_rs <- data.frame(SG17M_8h_BR3_depth_rs)
SG17M_8h_BR3_depth_rs$side <- "reverse_sense"
SG17M_8h_BR3_depth_sense <- data.frame(rbind(SG17M_8h_BR3_depth_fs, SG17M_8h_BR3_depth_rs))
SG17M_8h_BR3_depth_sense$X6 <- NULL
SG17M_8h_BR3_depth_sense$replicate <- "BR3"

SG17M_8h_depth <- data.frame(rbind(SG17M_8h_BR1_depth, SG17M_8h_BR2_depth, SG17M_8h_BR3_depth))
colnames(SG17M_8h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

SG17M_8h_depth_sense <- data.frame(rbind(SG17M_8h_BR1_depth_sense, SG17M_8h_BR2_depth_sense, SG17M_8h_BR3_depth_sense))
colnames(SG17M_8h_depth_sense) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

#SG17M_4h_depth$length <- SG17M_4h_depth$end - SG17M_4h_depth$start
#SG17M_8h_depth$length <- SG17M_8h_depth$end - SG17M_8h_depth$start

#write.csv(SG17M_4h_depth, file="SG17M_4h_length.csv", row.names = FALSE)
#write.csv(SG17M_8h_depth, file="SG17M_8h_length.csv", row.names = FALSE)


# Antisense hotspot 1, 4h
SG17M_4h_depth_hsp1 <- subset(SG17M_4h_depth, start > 989053 & end < 993602)
SG17M_4h_gtf_hsp1 <- subset(import_gtf_sg17m, start > 989053 & end < 993602)
SG17M_4h_depth_hsp1 <- SG17M_4h_depth_hsp1[order(SG17M_4h_depth_hsp1$replicate,SG17M_4h_depth_hsp1$start),]
rownames(SG17M_4h_depth_hsp1) <- NULL
SG17M_4h_depth_hsp1$num <- rownames(SG17M_4h_depth_hsp1)
SG17M_4h_depth_hsp1$num <- as.numeric(as.character(SG17M_4h_depth_hsp1$num))
SG17M_4h_gtf_hsp1$gene_name3 <- c("tRNA", "tRNA", "tRNA", "CDS", "tRNA","unknown", "CDS", "CDS", "CDS")

SG17M_hsp1_plot <- ggplot() + geom_segment(data=SG17M_4h_depth_hsp1, 
               aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "darkorange2", "blue", "forestgreen")) +
  xlab("\n") + ylab("Read IDs (SG17M-4h antisense transcripts)") + 
  geom_segment(data=SG17M_4h_gtf_hsp1, aes(x=ifelse(strandType == "+", start, 
                            ifelse(strandType == "\\.", start, end)), xend=ifelse(strandType == "+", end, 
                               ifelse(strandType == "\\.", end, start)),
                   y=-30, yend=-30, color=gene_name3), arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-72, x = 989644), label="Tyr-Anticodon-GTA", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-130, x = 989695), label="Gly-Anticodon-TCC", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-185, x = 989828), label="Thr-Anticodon-GGT", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-72, x = 990728), label="tufA", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-130, x = 991377), label="Trp-Anticodon-CCA", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-185, x = 991625), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-72, x = 992000), label="nusG", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-72, x = 992526), label="rplK", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-72, x = 993157), label="rplA", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(989053, 993602))


# Antisense hotspot 2
SG17M_4h_depth_hsp2 <- subset(SG17M_4h_depth, start > 5814083 & end < 5815379)
SG17M_4h_gtf_hsp2 <- subset(import_gtf_sg17m, start > 5814083 & end < 5815379)
SG17M_4h_depth_hsp2 <- SG17M_4h_depth_hsp2[order(SG17M_4h_depth_hsp2$replicate,
                                             SG17M_4h_depth_hsp2$start),]
SG17M_4h_gtf_hsp2$gene_name4 <- c("Lys-Anticodon-TTT (pKLC102)", "soj")
rownames(SG17M_4h_depth_hsp2) <- NULL
SG17M_4h_depth_hsp2$num <- rownames(SG17M_4h_depth_hsp2)
SG17M_4h_depth_hsp2$num <- as.numeric(as.character(SG17M_4h_depth_hsp2$num))

SG17M_hsp2_plot <- ggplot() + geom_segment(data=SG17M_4h_depth_hsp2, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), 
               size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "blue", "darkorange2")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=SG17M_4h_gtf_hsp2, aes(x=ifelse(strandType == "+", start, 
                            ifelse(strandType == "\\.", start, end)),  xend=ifelse(strandType == "+", end, 
                               ifelse(strandType == "\\.", end, start)), y=-25, yend=-25, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-158, x = 5814393), label="Lys-Anticodon-TTT (pKLC102)", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-158, x = 5814923), label="soj", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(5814083, 5815379))



# Antisense hotspot 3
SG17M_4h_depth_hsp3 <- subset(SG17M_4h_depth, start > 5836292 & end < 5844865)
SG17M_4h_gtf_hsp3 <- subset(import_gtf_sg17m, start > 5836292 & end < 5844865)
SG17M_4h_depth_hsp3 <- SG17M_4h_depth_hsp3[order(SG17M_4h_depth_hsp3$replicate,SG17M_4h_depth_hsp3$start),]
SG17M_4h_gtf_hsp3$gene_name4 <- c("CDS (hypothetical)","rapA", "CDS (hypothetical)", 
                                  "CDS (hypothetical)", "CDS (hypothetical)")

rownames(SG17M_4h_depth_hsp3) <- NULL
SG17M_4h_depth_hsp3$num <- rownames(SG17M_4h_depth_hsp3)
SG17M_4h_depth_hsp3$num <- as.numeric(as.character(SG17M_4h_depth_hsp3$num))

SG17M_hsp3_plot <- ggplot() + geom_segment(data=SG17M_4h_depth_hsp3, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "forestgreen", "darkorange2", 
                              "forestgreen", "forestgreen", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=SG17M_4h_gtf_hsp3, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-50, yend=-50, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-385, x = 5838292), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-385, x = 5840292), label="rapA", size=2.5, color="darkorange2") + 
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-385, x = 5842392), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-265, x = 5843392), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-145, x = 5843872), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma)


# Antisense hotspot 4
SG17M_4h_depth_hsp4 <- subset(SG17M_4h_depth, start > 551623 & end < 556783)
SG17M_4h_gtf_hsp4 <- subset(import_gtf_sg17m, start > 551623 & end < 556783)
SG17M_4h_depth_hsp4 <- SG17M_4h_depth_hsp4[order(SG17M_4h_depth_hsp4$replicate, SG17M_4h_depth_hsp4$start),]

SG17M_4h_gtf_hsp4$gene_name4 <- c("CDS (hypothetical)","CDS (hypothetical)", 
                                  "CDS", "CDS (hypothetical)", "CDS")

rownames(SG17M_4h_depth_hsp4) <- NULL
SG17M_4h_depth_hsp4$num <- rownames(SG17M_4h_depth_hsp4)
SG17M_4h_depth_hsp4$num <- as.numeric(as.character(SG17M_4h_depth_hsp4$num))

SG17M_hsp4_plot <- ggplot() + geom_segment(data=SG17M_4h_depth_hsp4, 
                                           aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "darkorange2", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=SG17M_4h_gtf_hsp4, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-24, yend=-24, color=gene_name4), arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=SG17M_4h_gtf_hsp4, aes(y=-80, x = 552417), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=SG17M_4h_gtf_hsp4, aes(y=-156, x = 552597), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=SG17M_4h_gtf_hsp4, aes(y=-80, x = 553623), label="intQ", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_hsp4, aes(y=-156, x = 554595), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp4, aes(y=-80, x = 555204), label="qseF", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma)


# Antisense hotspot 1, 8h
SG17M_8h_depth_hsp1 <- subset(SG17M_8h_depth, start > 989053 & end < 993602)
SG17M_8h_gtf_hsp1 <- subset(import_gtf_sg17m, start > 989053 & end < 993602)
SG17M_8h_depth_hsp1 <- SG17M_8h_depth_hsp1[order(SG17M_8h_depth_hsp1$replicate, SG17M_8h_depth_hsp1$start),]
rownames(SG17M_8h_depth_hsp1) <- NULL
SG17M_8h_depth_hsp1$num <- rownames(SG17M_8h_depth_hsp1)
SG17M_8h_depth_hsp1$num <- as.numeric(as.character(SG17M_8h_depth_hsp1$num))
SG17M_8h_gtf_hsp1$gene_name4 <- c("tRNA", "tRNA", "tRNA", 
                                "CDS", "tRNA","unknown", "CDS", "CDS", "CDS")

SG17M_hsp1_plot_8h <- ggplot() + geom_segment(data=SG17M_8h_depth_hsp1, 
               aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3",
                              "darkorange2", "blue", "forestgreen")) +
  xlab("\n  ") + ylab("Read IDs (SG17M-8h antisense transcripts)") + 
  geom_segment(data=SG17M_8h_gtf_hsp1, 
               aes(x=ifelse(strandType == "+", start, 
                            ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, 
                               ifelse(strandType == "\\.", end, start)),
                   y=-8, yend=-8, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-18, x = 989644), label="Tyr-Anticodon-GTA", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-29, x = 989695), label="Gly-Anticodon-TCC", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-41, x = 989828), label="Thr-Anticodon-GGT", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-18, x = 990728), label="tufA", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-29, x = 991377), label="Trp-Anticodon-CCA", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-41, x = 991625), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-18, x = 992000), label="nusG", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-18, x = 992526), label="rplK", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-18, x = 993157), label="rplA", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(989053, 993602))


# Antisense hotspot 2
SG17M_8h_depth_hsp2 <- subset(SG17M_8h_depth, start > 5814083 & end < 5815379)
SG17M_8h_gtf_hsp2 <- subset(import_gtf_sg17m, start > 5814083 & end < 5815379)
SG17M_8h_depth_hsp2 <- SG17M_8h_depth_hsp2[order(SG17M_8h_depth_hsp2$replicate,
                                             SG17M_8h_depth_hsp2$start),]
rownames(SG17M_8h_depth_hsp2) <- NULL
SG17M_8h_gtf_hsp2$gene_name4 <- c("Lys-Anticodon-TTT (pKLC102)", "soj")
SG17M_8h_depth_hsp2$num <- rownames(SG17M_8h_depth_hsp2)
SG17M_8h_depth_hsp2$num <- as.numeric(as.character(SG17M_8h_depth_hsp2$num))

SG17M_hsp2_plot_8h <- ggplot() + geom_segment(data=SG17M_8h_depth_hsp2, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3","blue", "darkorange2")) + 
  xlab("Genome position") + ylab(" ") + 
  geom_segment(data=SG17M_4h_gtf_hsp2, 
               aes(x=ifelse(strandType == "+", start, 
                            ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, 
                               ifelse(strandType == "\\.", end, start)),
                   y=-5, yend=-5, color=gene_name4), arrow = arrow(length = unit(0.06, "inches")), size = 1) +
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-22, x = 5814393), label="Lys-Anticodon-TTT (pKLC102)", size=2.5, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-22, x = 5814923), label="soj", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(5814083, 5815379))

# Antisense hotspot 3
SG17M_8h_depth_hsp3 <- subset(SG17M_8h_depth, start > 5836292 & end < 5844865)
SG17M_8h_gtf_hsp3 <- subset(import_gtf_sg17m, start > 5836292 & end < 5844865)
SG17M_8h_depth_hsp3$length <- SG17M_8h_depth_hsp3$end - SG17M_8h_depth_hsp3$start

SG17M_8h_gtf_hsp3$gene_name4 <- c("CDS (hypothetical)","rapA", 
                                  "CDS (hypothetical)", "CDS (hypothetical)", "CDS (hypothetical)")
SG17M_8h_depth_hsp3 <- SG17M_8h_depth_hsp3[order(SG17M_8h_depth_hsp3$replicate, SG17M_8h_depth_hsp3$start),]
rownames(SG17M_8h_depth_hsp3) <- NULL
SG17M_8h_depth_hsp3$num <- rownames(SG17M_8h_depth_hsp3)
SG17M_8h_depth_hsp3$num <- as.numeric(as.character(SG17M_8h_depth_hsp3$num))

SG17M_hsp3_plot_8h <- ggplot() + geom_segment(data=SG17M_8h_depth_hsp3, 
                                              aes(x=start, xend=end, y=num, yend=num, 
                                             colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3", "forestgreen", "darkorange2", 
                              "forestgreen", "forestgreen", "forestgreen")) +
  xlab(" ") + ylab(" ") + 
  geom_segment(data=SG17M_8h_gtf_hsp3, aes(x=ifelse(strandType == "+", start, 
                            ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, 
                               ifelse(strandType == "\\.", end, start)), y=-20, yend=-20, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=SG17M_8h_gtf_hsp3, aes(y=-95, x = 5838292), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=SG17M_8h_gtf_hsp3, aes(y=-95, x = 5840292), label="rapA", size=2.5, color="darkorange2") + 
  geom_label(data=SG17M_8h_gtf_hsp3, aes(y=-95, x = 5842392), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-70, x = 5843392), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-45, x = 5843872), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(angle=45,hjust=1,size=8)) +
  scale_x_continuous(label=scales::comma)

# Antisense hotspot 4
SG17M_8h_depth_hsp4 <- subset(SG17M_8h_depth, start > 551623 & end < 556783)
SG17M_8h_gtf_hsp4 <- subset(import_gtf_sg17m, start > 551623 & end < 556783)
SG17M_8h_depth_hsp4 <- SG17M_8h_depth_hsp4[order(SG17M_8h_depth_hsp4$replicate, SG17M_8h_depth_hsp4$start),]
SG17M_8h_gtf_hsp4$gene_name4 <- c("CDS (hypothetical)","CDS (hypothetical)", 
                                  "CDS", "CDS (hypothetical)", "CDS")

rownames(SG17M_8h_depth_hsp4) <- NULL
SG17M_8h_depth_hsp4$num <- rownames(SG17M_8h_depth_hsp4)
SG17M_8h_depth_hsp4$num <- as.numeric(as.character(SG17M_8h_depth_hsp4$num))

SG17M_hsp4_plot_8h <- ggplot() +
  geom_segment(data=SG17M_8h_depth_hsp4, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "gold2", "lightcyan3",
                              "darkorange2", "forestgreen")) + xlab("\n") + ylab(" ") + 
  geom_segment(data=SG17M_8h_gtf_hsp4, aes(x=ifelse(strandType == "+", start, 
                            ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, 
                               ifelse(strandType == "\\.", end, start)), y=-24, yend=-24, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 1, alpha=1) +
  geom_label(data=SG17M_8h_gtf_hsp4, aes(y=-80, x = 552417), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=SG17M_8h_gtf_hsp4, aes(y=-147, x = 552597), label="CDS (hypothetical)", size=2.5, color="forestgreen") + 
  geom_label(data=SG17M_8h_gtf_hsp4, aes(y=-80, x = 553623), label="intQ", size=2.5, color="darkorange2") +
  geom_label(data=SG17M_8h_gtf_hsp4, aes(y=-147, x = 554595), label="CDS (hypothetical)", size=2.5, color="forestgreen") +
  geom_label(data=SG17M_8h_gtf_hsp4, aes(y=-80, x = 555204), label="qseF", size=2.5, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,size=8)) +
  scale_x_continuous(label=scales::comma)


sg17m_4h_antisense_hotspots <- ggarrange(SG17M_hsp1_plot, SG17M_hsp2_plot, SG17M_hsp4_plot, nrow= 1,
          labels=c("A", "B", "C"))

sg17m_8h_antisense_hotspots <- ggarrange(SG17M_hsp1_plot_8h, SG17M_hsp2_plot_8h, SG17M_hsp4_plot_8h, nrow= 1,
          labels=c("D", "E", "F"))

sg17m_antisense_hsp <- ggarrange(sg17m_4h_antisense_hotspots, sg17m_8h_antisense_hotspots, nrow=2)

ggsave(sg17m_antisense_hsp, filename="save_figures/Figure_07.tif", dpi=300, device="tiff", units="cm", width=28,height=23)
