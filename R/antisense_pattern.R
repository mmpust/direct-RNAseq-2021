
# title: "Antisense pattern"
# author: "Marie Pust"
# date: "05 7 2021"


# clean environment
rm(list = ls())

# load libraries
library(rlist)
library(ggplot2)
library(readr)
library(readxl)
library(Rsubread)
library(stringr)
library(tidyr)
library(ggpubr)
library(circlize)
library(dplyr)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)

# load input data
setwd("/d_Ranalysis")
input_bam_NN2<- list.files(
  path = "/NN2_run/d_raw_bam",
  pattern = ".bam", full.names = TRUE)


input_bam_sg17m <- list.files(
  path = "/SG17M_run/d_raw_bam",
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
featureCounts_PA_sense_df_NN2_norm_4h_2 <- calcNormFactors(featureCounts_PA_sense_df_NN2_norm_4h_1, method="TMM")
featureCounts_PA_sense_df_NN2_norm_4h_3 <- cpm(featureCounts_PA_sense_df_NN2_norm_4h_2, log=FALSE)
featureCounts_PA_sense_df_NN2_norm_4h_4 <- data.frame(featureCounts_PA_sense_df_NN2_norm_4h_3)

featureCounts_PA_sense_df_NN2_norm_8h <- featureCounts_PA_sense_df_NN2_norm[,4:6]
featureCounts_PA_sense_df_NN2_norm_8h_1 <- DGEList(counts=featureCounts_PA_sense_df_NN2_norm_8h)
featureCounts_PA_sense_df_NN2_norm_8h_2 <- calcNormFactors(featureCounts_PA_sense_df_NN2_norm_8h_1, method="TMM")
featureCounts_PA_sense_df_NN2_norm_8h_3 <- cpm(featureCounts_PA_sense_df_NN2_norm_8h_2, log=FALSE)
featureCounts_PA_sense_df_NN2_norm_8h_4 <- data.frame(featureCounts_PA_sense_df_NN2_norm_8h_3)

featureCounts_PA_antisense_df_NN2_norm <- featureCounts_PA_antisense_df_NN2
featureCounts_PA_antisense_df_NN2_norm_4h <- featureCounts_PA_antisense_df_NN2_norm[,1:3]
featureCounts_PA_antisense_df_NN2_norm_4h_1 <- DGEList(counts=featureCounts_PA_antisense_df_NN2_norm_4h)
featureCounts_PA_antisense_df_NN2_norm_4h_2 <- calcNormFactors(featureCounts_PA_antisense_df_NN2_norm_4h_1, method="TMM")
featureCounts_PA_antisense_df_NN2_norm_4h_3 <- cpm(featureCounts_PA_antisense_df_NN2_norm_4h_2, log=FALSE)
featureCounts_PA_antisense_df_NN2_norm_4h_4 <- data.frame(featureCounts_PA_antisense_df_NN2_norm_4h_3)
featureCounts_PA_antisense_df_NN2_norm_8h <- featureCounts_PA_antisense_df_NN2_norm[,4:6]
featureCounts_PA_antisense_df_NN2_norm_8h_1 <- DGEList(counts=featureCounts_PA_antisense_df_NN2_norm_8h)
featureCounts_PA_antisense_df_NN2_norm_8h_2 <- calcNormFactors(featureCounts_PA_antisense_df_NN2_norm_8h_1, method="TMM")
featureCounts_PA_antisense_df_NN2_norm_8h_3 <- cpm(featureCounts_PA_antisense_df_NN2_norm_8h_2, log=FALSE)
featureCounts_PA_antisense_df_NN2_norm_8h_4 <- data.frame(featureCounts_PA_antisense_df_NN2_norm_8h_3)


featureCounts_PA_sense_df_sg17m_norm <- featureCounts_PA_sense_df_sg17m
featureCounts_PA_sense_df_sg17m_norm_4h <- featureCounts_PA_sense_df_sg17m_norm[,1:3]
featureCounts_PA_sense_df_sg17m_norm_4h_1 <- DGEList(counts=featureCounts_PA_sense_df_sg17m_norm_4h)
featureCounts_PA_sense_df_sg17m_norm_4h_2 <- calcNormFactors(featureCounts_PA_sense_df_sg17m_norm_4h_1, method="TMM")
featureCounts_PA_sense_df_sg17m_norm_4h_3 <- cpm(featureCounts_PA_sense_df_sg17m_norm_4h_2, log=FALSE)
featureCounts_PA_sense_df_sg17m_norm_4h_4 <- data.frame(featureCounts_PA_sense_df_sg17m_norm_4h_3)
featureCounts_PA_sense_df_sg17m_norm_8h <- featureCounts_PA_sense_df_sg17m_norm[,4:6]
featureCounts_PA_sense_df_sg17m_norm_8h_1 <- DGEList(counts=featureCounts_PA_sense_df_sg17m_norm_8h)
featureCounts_PA_sense_df_sg17m_norm_8h_2 <- calcNormFactors(featureCounts_PA_sense_df_sg17m_norm_8h_1, method="TMM")
featureCounts_PA_sense_df_sg17m_norm_8h_3 <- cpm(featureCounts_PA_sense_df_sg17m_norm_8h_2, log=FALSE)
featureCounts_PA_sense_df_sg17m_norm_8h_4 <- data.frame(featureCounts_PA_sense_df_sg17m_norm_8h_3)

featureCounts_PA_antisense_df_sg17m_norm <- featureCounts_PA_antisense_df_sg17m
featureCounts_PA_antisense_df_sg17m_norm_4h <- featureCounts_PA_antisense_df_sg17m_norm[,1:3]
featureCounts_PA_antisense_df_sg17m_norm_4h_1 <- DGEList(counts=featureCounts_PA_antisense_df_sg17m_norm_4h)
featureCounts_PA_antisense_df_sg17m_norm_4h_2 <- calcNormFactors(featureCounts_PA_antisense_df_sg17m_norm_4h_1, method="TMM")
featureCounts_PA_antisense_df_sg17m_norm_4h_3 <- cpm(featureCounts_PA_antisense_df_sg17m_norm_4h_2, log=FALSE)
featureCounts_PA_antisense_df_sg17m_norm_4h_4 <- data.frame(featureCounts_PA_antisense_df_sg17m_norm_4h_3)
featureCounts_PA_antisense_df_sg17m_norm_8h <- featureCounts_PA_antisense_df_sg17m_norm[,4:6]
featureCounts_PA_antisense_df_sg17m_norm_8h_1 <- DGEList(counts=featureCounts_PA_antisense_df_sg17m_norm_8h)
featureCounts_PA_antisense_df_sg17m_norm_8h_2 <- calcNormFactors(featureCounts_PA_antisense_df_sg17m_norm_8h_1, method="TMM")
featureCounts_PA_antisense_df_sg17m_norm_8h_3 <- cpm(featureCounts_PA_antisense_df_sg17m_norm_8h_2, log=FALSE)
featureCounts_PA_antisense_df_sg17m_norm_8h_4 <- data.frame(featureCounts_PA_antisense_df_sg17m_norm_8h_3)


# circular plots, NN2, 4h
antisense_NN2_norm_4h_circ <- cpm(featureCounts_PA_antisense_df_NN2_norm_4h_2, log=TRUE)
antisense_NN2_norm_4h_circ <- data.frame(antisense_NN2_norm_4h_circ)
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
                                                             ifelse(NN2_4h_BR1 < 9 & NN2_4h_BR1 >= 6, "lightgreen",
                                                                    ifelse(NN2_4h_BR1 < 6, "palegoldenrod", "NA"))))))



antisense_NN2_norm_4h_circ_BR2 <- select(antisense_NN2_norm_4h_circ, c(isolate,start, end, NN2_4h_BR2))
color_NN2_BR2_4h <- with(antisense_NN2_norm_4h_circ_BR2,
                                             ifelse(NN2_4h_BR2 >= 12, "palegreen4", 
                                               ifelse(NN2_4h_BR2 < 12 & NN2_4h_BR2 >= 10, "palegreen3", 
                                                      ifelse(NN2_4h_BR2 < 10 & NN2_4h_BR2 >= 9, "palegreen2",
                                                             ifelse(NN2_4h_BR2 < 9 & NN2_4h_BR2 >= 6, "lightgreen",
                                                                    ifelse(NN2_4h_BR2 < 6, "palegoldenrod", "NA"))))))

antisense_NN2_norm_4h_circ_BR3 <- select(antisense_NN2_norm_4h_circ, c(isolate,start, end, NN2_4h_BR3))
color_NN2_BR3_4h <- with(antisense_NN2_norm_4h_circ_BR3,
                                             ifelse(NN2_4h_BR3 >= 12, "palegreen4", 
                                               ifelse(NN2_4h_BR3 < 12 & NN2_4h_BR3 >= 10, "palegreen3", 
                                                      ifelse(NN2_4h_BR3 < 10 & NN2_4h_BR3 >= 9, "palegreen2",
                                                             ifelse(NN2_4h_BR3 < 9 & NN2_4h_BR3 >= 6, "lightgreen",
                                                                    ifelse(NN2_4h_BR3 < 6, "palegoldenrod", "NA"))))))

# tiff(filename = "save_figures/circular_plot_nn2_4h.tif", units="cm", res=800, width=18, height=16)
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
                      circos.genomicRect(region, value, col="gold")}, bg.col = "dodgerblue4")

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

# circular plots, NN2, 8h
antisense_NN2_norm_8h_circ <- cpm(featureCounts_PA_antisense_df_NN2_norm_8h_2, log=TRUE)
antisense_NN2_norm_8h_circ <- data.frame(antisense_NN2_norm_8h_circ)
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
                                                             ifelse(NN2_8h_BR1 < 9 & NN2_8h_BR1 >= 6, "lightgreen",
                                                                    ifelse(NN2_8h_BR1 < 6, "palegoldenrod", "NA"))))))

antisense_NN2_norm_8h_circ_BR2 <- select(antisense_NN2_norm_8h_circ, c(isolate,start, end, NN2_8h_BR2))
color_NN2_BR2_8h <- with(antisense_NN2_norm_8h_circ_BR2,
                                             ifelse(NN2_8h_BR2 >= 12, "palegreen4", 
                                               ifelse(NN2_8h_BR2 < 12 & NN2_8h_BR2 >= 10, "palegreen3", 
                                                      ifelse(NN2_8h_BR2 < 10 & NN2_8h_BR2 >= 9, "palegreen2",
                                                             ifelse(NN2_8h_BR2 < 9 & NN2_8h_BR2 >= 6, "lightgreen",
                                                                    ifelse(NN2_8h_BR2 < 6, "palegoldenrod", "NA"))))))

antisense_NN2_norm_8h_circ_BR3 <- select(antisense_NN2_norm_8h_circ, c(isolate,start, end, NN2_8h_BR3))
color_NN2_BR3_8h <- with(antisense_NN2_norm_8h_circ_BR3,
                                             ifelse(NN2_8h_BR3 >= 12, "palegreen4", 
                                               ifelse(NN2_8h_BR3 < 12 & NN2_8h_BR3 >= 10, "palegreen3", 
                                                      ifelse(NN2_8h_BR3 < 10 & NN2_8h_BR3 >= 9, "palegreen2",
                                                             ifelse(NN2_8h_BR3 < 9 & NN2_8h_BR3 >= 6, "lightgreen",
                                                                    ifelse(NN2_8h_BR3 < 6, "palegoldenrod", "NA"))))))

# tiff(filename = "save_figures/circular_plot_nn2_8h.tif", units="cm", res=800, width=18, height=16)
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
                      circos.genomicRect(region, value, col="gold")}, bg.col = "dodgerblue4")

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
antisense_sg17m_norm_4h_circ <- cpm(featureCounts_PA_antisense_df_sg17m_norm_4h_2, log=TRUE)
antisense_sg17m_norm_4h_circ <- data.frame(antisense_sg17m_norm_4h_circ)
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
                                                      ifelse(sg17m_4h_BR1 < 10 & sg17m_4h_BR1 >= 9, "palegreen2",
                                                             ifelse(sg17m_4h_BR1 < 9 & sg17m_4h_BR1 >= 6, "lightgreen",
                                                                    ifelse(sg17m_4h_BR1 < 6, "palegoldenrod", "NA"))))))



antisense_sg17m_norm_4h_circ_BR2 <- select(antisense_sg17m_norm_4h_circ, c(isolate,start, end, sg17m_4h_BR2))
color_sg17m_BR2_4h <- with(antisense_sg17m_norm_4h_circ_BR2,
                                             ifelse(sg17m_4h_BR2 >= 12, "palegreen4", 
                                               ifelse(sg17m_4h_BR2 < 12 & sg17m_4h_BR2 >= 10, "palegreen3", 
                                                      ifelse(sg17m_4h_BR2 < 10 & sg17m_4h_BR2 >= 9, "palegreen2",
                                                             ifelse(sg17m_4h_BR2 < 9 & sg17m_4h_BR2 >= 6, "lightgreen",
                                                                    ifelse(sg17m_4h_BR2 < 6, "palegoldenrod", "NA"))))))

antisense_sg17m_norm_4h_circ_BR3 <- select(antisense_sg17m_norm_4h_circ, c(isolate,start, end, sg17m_4h_BR3))
color_sg17m_BR3_4h <- with(antisense_sg17m_norm_4h_circ_BR3,
                                             ifelse(sg17m_4h_BR3 >= 12, "palegreen4", 
                                               ifelse(sg17m_4h_BR3 < 12 & sg17m_4h_BR3 >= 10, "palegreen3", 
                                                      ifelse(sg17m_4h_BR3 < 10 & sg17m_4h_BR3 >= 9, "palegreen2",
                                                             ifelse(sg17m_4h_BR3 < 9 & sg17m_4h_BR3 >= 6, "lightgreen",
                                                                    ifelse(sg17m_4h_BR3 < 6, "palegoldenrod", "NA"))))))
# tiff(filename = "save_figures/circular_plot_sg17m_4h.tif", units="cm", res=800, width=18, height=16)
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
                      circos.genomicRect(region, value, col="gold")}, bg.col = "dodgerblue4")

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
antisense_sg17m_norm_8h_circ <- cpm(featureCounts_PA_antisense_df_sg17m_norm_8h_2, log=TRUE)
antisense_sg17m_norm_8h_circ <- data.frame(antisense_sg17m_norm_8h_circ)
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
                                                      ifelse(sg17m_8h_BR1 < 10 & sg17m_8h_BR1 >= 9, "palegreen2",
                                                             ifelse(sg17m_8h_BR1 < 9 & sg17m_8h_BR1 >= 6, "lightgreen",
                                                                    ifelse(sg17m_8h_BR1 < 6, "palegoldenrod", "NA"))))))

antisense_sg17m_norm_8h_circ_BR2 <- select(antisense_sg17m_norm_8h_circ, c(isolate,start, end, sg17m_8h_BR2))
color_sg17m_BR2_8h <- with(antisense_sg17m_norm_8h_circ_BR2,
                                             ifelse(sg17m_8h_BR2 >= 12, "palegreen4", 
                                               ifelse(sg17m_8h_BR2 < 12 & sg17m_8h_BR2 >= 10, "palegreen3", 
                                                      ifelse(sg17m_8h_BR2 < 10 & sg17m_8h_BR2 >= 9, "palegreen2",
                                                             ifelse(sg17m_8h_BR2 < 9 & sg17m_8h_BR2 >= 6, "lightgreen",
                                                                    ifelse(sg17m_8h_BR2 < 6, "palegoldenrod", "NA"))))))

antisense_sg17m_norm_8h_circ_BR3 <- select(antisense_sg17m_norm_8h_circ, c(isolate,start, end, sg17m_8h_BR3))
color_sg17m_BR3_8h <- with(antisense_sg17m_norm_8h_circ_BR3,
                                             ifelse(sg17m_8h_BR3 >= 12, "palegreen4", 
                                               ifelse(sg17m_8h_BR3 < 12 & sg17m_8h_BR3 >= 10, "palegreen3", 
                                                      ifelse(sg17m_8h_BR3 < 10 & sg17m_8h_BR3 >= 9, "palegreen2",
                                                             ifelse(sg17m_8h_BR3 < 9 & sg17m_8h_BR3 >= 6, "lightgreen",
                                                                    ifelse(sg17m_8h_BR3 < 6, "palegoldenrod", "NA"))))))
# tiff(filename = "save_figures/circular_plot_sg17m_8h.tif", units="cm", res=800, width=18, height=16)
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
                      circos.genomicRect(region, value, col="gold")}, bg.col = "dodgerblue4")

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

# NN2, antisense, 4h
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
gtf_nn2_antisense_8h_high

# NN2, antisense, 8h
nn2_antisense_n10high_list_8h <- list()
nn2_antisense_n10high_8h_1 <- featureCounts_PA_antisense_df_NN2_norm_8h_4[order(
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR1, 
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR2,
  featureCounts_PA_antisense_df_NN2_norm_8h_4$NN2_8h_BR3, decreasing=TRUE),]
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


# SG17M, antisense, 4h
sg17m_antisense_n10high_list_4h <- list()
sg17m_antisense_n10high_4h_1 <- featureCounts_PA_antisense_df_sg17m_norm_4h_4[order(
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR1, 
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR2,
  featureCounts_PA_antisense_df_sg17m_norm_4h_4$sg17m_4h_BR3, decreasing=TRUE),]
sg17m_antisense_n10high_list_4h = list.append(sg17m_antisense_n10high_list_4h, rownames(sg17m_antisense_n10high_4h_1)[1:n_features])
sg17m_antisense_n10high_list_4h <- unlist(sg17m_antisense_n10high_list_4h)
featureCounts_PA_antisense_df_sg17m_4h_high <- subset(featureCounts_PA_antisense_df_sg17m,
                                                    rownames(featureCounts_PA_antisense_df_sg17m) %in%
                                                      sg17m_antisense_n10high_list_4h)
gtf_sg17m_antisense_4h_high <- subset(import_gtf_sg17m, 
                                    transcript_id %in% rownames(featureCounts_PA_antisense_df_sg17m_4h_high))
gtf_sg17m_antisense_4h_high$extract_start <- gtf_sg17m_antisense_4h_high$start - feature_length
gtf_sg17m_antisense_4h_high$extract_end <- gtf_sg17m_antisense_4h_high$end + feature_length


# SG17M, antisense, 8h
sg17m_antisense_n10high_list_8h <- list()
sg17m_antisense_n10high_8h_1 <- featureCounts_PA_antisense_df_sg17m_norm_8h_4[order(
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR1, 
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR2,
  featureCounts_PA_antisense_df_sg17m_norm_8h_4$sg17m_8h_BR3, decreasing=TRUE),]
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

# NN2 ####
NN2_4h_BR1_depth_fa <- read_delim("NN2_CS_4h_BR1.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR1_depth_fa <- data.frame(NN2_4h_BR1_depth_fa)
NN2_4h_BR1_depth_fa$side <- "forward_antisense"
NN2_4h_BR1_depth_ra <- read_delim("NN2_CS_4h_BR1.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR1_depth_ra <- data.frame(NN2_4h_BR1_depth_ra)
NN2_4h_BR1_depth_ra$side <- "reverse_antisense"
NN2_4h_BR1_depth <- data.frame(rbind(NN2_4h_BR1_depth_fa, NN2_4h_BR1_depth_ra))
NN2_4h_BR1_depth$X6 <- NULL
NN2_4h_BR1_depth$replicate <- "BR1"

NN2_4h_BR2_depth_fa <- read_delim("NN2_CS_4h_BR2.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR2_depth_fa <- data.frame(NN2_4h_BR2_depth_fa)
NN2_4h_BR2_depth_fa$side <- "forward_antisense"
NN2_4h_BR2_depth_ra <- read_delim("NN2_CS_4h_BR2.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR2_depth_ra <- data.frame(NN2_4h_BR2_depth_ra)
NN2_4h_BR2_depth_ra$side <- "reverse_antisense"
NN2_4h_BR2_depth <- data.frame(rbind(NN2_4h_BR2_depth_fa, NN2_4h_BR2_depth_ra))
NN2_4h_BR2_depth$X6 <- NULL
NN2_4h_BR2_depth$replicate <- "BR2"


NN2_4h_BR3_depth_fa <- read_delim("NN2_CS_4h_BR3.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR3_depth_fa <- data.frame(NN2_4h_BR3_depth_fa)
NN2_4h_BR3_depth_fa$side <- "forward_antisense"
NN2_4h_BR3_depth_ra <- read_delim("NN2_CS_4h_BR3.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_4h_BR3_depth_ra <- data.frame(NN2_4h_BR3_depth_ra)
NN2_4h_BR3_depth_ra$side <- "reverse_antisense"
NN2_4h_BR3_depth <- data.frame(rbind(NN2_4h_BR3_depth_fa, NN2_4h_BR3_depth_ra))
NN2_4h_BR3_depth$X6 <- NULL
NN2_4h_BR3_depth$replicate <- "BR3"
NN2_4h_BR3_depth

NN2_4h_depth <- data.frame(rbind(NN2_4h_BR1_depth, NN2_4h_BR2_depth, NN2_4h_BR3_depth))
colnames(NN2_4h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

# Antisense hotspot 1
NN2_4h_depth_hsp1 <- subset(NN2_4h_depth, start > 766300 & end < 767903)
NN2_4h_gtf_hsp1 <- subset(import_gtf_NN2, start > 766300 & end < 767903)
NN2_4h_gtf_hsp1$gene_name4 <- c("tRNA-Trp-Anticodon-CCA", "CDS (undefined)", "nusG")
NN2_4h_depth_hsp1 <- NN2_4h_depth_hsp1[order(NN2_4h_depth_hsp1$start),]
rownames(NN2_4h_depth_hsp1) <- NULL
NN2_4h_depth_hsp1$num <- rownames(NN2_4h_depth_hsp1)
NN2_4h_depth_hsp1$num <- as.numeric(as.character(NN2_4h_depth_hsp1$num))

nn2_hsp1_plot <-
  ggplot() +
  geom_segment(data=NN2_4h_depth_hsp1, 
               aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "forestgreen", "darkorange2", "blue")) +
  xlab("\n") + ylab("Read IDs (NN2-4h antisense transcripts)") + 
  geom_segment(data=NN2_4h_gtf_hsp1, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-20, yend=-20, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-40, x = 766488), 
             label="Trp-Anticodon-CCA", size=2, color="blue") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-40, x = 766910), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-40, x = 767260), 
             label="nusG", size=2, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(766300, 767903))


# Antisense hotspot 2
NN2_4h_depth_hsp2 <- subset(NN2_4h_depth, start > 5577311 & end < 5578609)
NN2_4h_gtf_hsp2 <- subset(import_gtf_NN2, start > 5577311 & end < 5578609)
NN2_4h_gtf_hsp2
NN2_4h_gtf_hsp2$gene_name4 <- c("tRNA-Lys-Anticodon-TTT-pseudo", "soj")
NN2_4h_depth_hsp2 <- NN2_4h_depth_hsp2[order(NN2_4h_depth_hsp2$start),]
rownames(NN2_4h_depth_hsp2) <- NULL
NN2_4h_depth_hsp2$num <- rownames(NN2_4h_depth_hsp2)
NN2_4h_depth_hsp2$num <- as.numeric(as.character(NN2_4h_depth_hsp2$num))

nn2_hsp2_plot <-
  ggplot() +
  geom_segment(data=NN2_4h_depth_hsp2, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "darkorange2", "blue")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_4h_gtf_hsp2, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-6, yend=-6, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-12, x = 5577799), 
             label="Lys-Anticodon-TTT-pseudo", size=2, color="blue") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-12, x = 5578200), 
             label="soj", size=2, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(5577311, 5578609))


# Antisense hotspot 3
NN2_4h_depth_hsp3 <- subset(NN2_4h_depth, start > 5600220 & end < 5606993)
NN2_4h_gtf_hsp3 <- subset(import_gtf_NN2, start > 5600220 & end < 5606993)
NN2_4h_depth_hsp3$length <- NN2_4h_depth_hsp3$end - NN2_4h_depth_hsp3$start

NN2_4h_gtf_hsp3
NN2_4h_gtf_hsp3$gene_name4 <- c("CDS (undefined)","rapA", 
                                "CDS (undefined)", "CDS (undefined)", "CDS (undefined)")
NN2_4h_depth_hsp3 <- NN2_4h_depth_hsp3[order(NN2_4h_depth_hsp3$start),]
rownames(NN2_4h_depth_hsp3) <- NULL
NN2_4h_depth_hsp3$num <- rownames(NN2_4h_depth_hsp3)
NN2_4h_depth_hsp3$num <- as.numeric(as.character(NN2_4h_depth_hsp3$num))

nn2_hsp3_plot <-
  ggplot() +
  geom_segment(data=NN2_4h_depth_hsp3, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "forestgreen", "darkorange2", 
                              "forestgreen", "forestgreen", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_4h_gtf_hsp3, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-18, yend=-18, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=1) +
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-35, x = 5601520), 
             label="CDS (undefined)", size=2, color="forestgreen") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-35, x = 5603420), 
             label="rapA", size=2, color="darkorange2") + 
  geom_label(data=NN2_4h_gtf_hsp1, aes(y=-35, x = 5605620), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(5600220, 5606993))


nn2_antisense_hsp_4h <-
  ggarrange(nn2_hsp1_plot, nn2_hsp2_plot, nn2_hsp3_plot, nrow = 1,
          labels=c("A", "B", "C"))


# NN2 ####
NN2_8h_BR1_depth_fa <- read_delim("NN2_CS_8h_BR1.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR1_depth_fa <- data.frame(NN2_8h_BR1_depth_fa)
NN2_8h_BR1_depth_fa$side <- "forward_antisense"
NN2_8h_BR1_depth_ra <- read_delim("NN2_CS_8h_BR1.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR1_depth_ra <- data.frame(NN2_8h_BR1_depth_ra)
NN2_8h_BR1_depth_ra$side <- "reverse_antisense"
NN2_8h_BR1_depth <- data.frame(rbind(NN2_8h_BR1_depth_fa, NN2_8h_BR1_depth_ra))
NN2_8h_BR1_depth$X6 <- NULL
NN2_8h_BR1_depth$replicate <- "BR1"

NN2_8h_BR2_depth_fa <- read_delim("NN2_CS_8h_BR2.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR2_depth_fa <- data.frame(NN2_8h_BR2_depth_fa)
NN2_8h_BR2_depth_fa$side <- "forward_antisense"
NN2_8h_BR2_depth_ra <- read_delim("NN2_CS_8h_BR2.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR2_depth_ra <- data.frame(NN2_8h_BR2_depth_ra)
NN2_8h_BR2_depth_ra$side <- "reverse_antisense"
NN2_8h_BR2_depth <- data.frame(rbind(NN2_8h_BR2_depth_fa, NN2_8h_BR2_depth_ra))
NN2_8h_BR2_depth$X6 <- NULL
NN2_8h_BR2_depth$replicate <- "BR2"


NN2_8h_BR3_depth_fa <- read_delim("NN2_CS_8h_BR3.sorted.forward_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR3_depth_fa <- data.frame(NN2_8h_BR3_depth_fa)
NN2_8h_BR3_depth_fa$side <- "forward_antisense"
NN2_8h_BR3_depth_ra <- read_delim("NN2_CS_8h_BR3.sorted.reverse_sense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NN2_8h_BR3_depth_ra <- data.frame(NN2_8h_BR3_depth_ra)
NN2_8h_BR3_depth_ra$side <- "reverse_antisense"
NN2_8h_BR3_depth <- data.frame(rbind(NN2_8h_BR3_depth_fa, NN2_8h_BR3_depth_ra))
NN2_8h_BR3_depth$X6 <- NULL
NN2_8h_BR3_depth$replicate <- "BR3"

NN2_8h_depth <- data.frame(rbind(NN2_8h_BR1_depth, NN2_8h_BR2_depth, NN2_8h_BR3_depth))
colnames(NN2_8h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")


# Antisense hotspot 1
NN2_8h_depth_hsp1 <- subset(NN2_8h_depth, start > 766300 & end < 767903)
NN2_8h_gtf_hsp1 <- subset(import_gtf_NN2, start > 766300 & end < 767903)
NN2_8h_gtf_hsp1$gene_name4 <- c("tRNA-Trp-Anticodon-CCA", "CDS (undefined)", "nusG")
NN2_8h_depth_hsp1 <- NN2_8h_depth_hsp1[order(NN2_8h_depth_hsp1$start),]
rownames(NN2_8h_depth_hsp1) <- NULL
NN2_8h_depth_hsp1$num <- rownames(NN2_8h_depth_hsp1)
NN2_8h_depth_hsp1$num <- as.numeric(as.character(NN2_8h_depth_hsp1$num))

nn2_hsp1_plot_8h <-
  ggplot() +
  geom_segment(data=NN2_8h_depth_hsp1, 
               aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3", 
                              "forestgreen", "darkorange2", "blue")) +
  xlab("\n\n") + ylab("Read IDs (NN2-8h antisense transcripts)") + 
  geom_segment(data=NN2_8h_gtf_hsp1, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-5, yend=-5, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-10, x = 766488), 
             label="Trp-Anticodon-CCA", size=2, color="blue") + 
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-10, x = 766910), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-10, x = 767260), 
             label="nusG", size=2, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(766300, 767903))

# Antisense hotspot 2
NN2_8h_depth_hsp2 <- subset(NN2_8h_depth, start > 5577311 & end < 5578609)
NN2_8h_gtf_hsp2 <- subset(import_gtf_NN2, start > 5577311 & end < 5578609)

NN2_8h_gtf_hsp2$gene_name4 <- c("tRNA-Lys-Anticodon-TTT-pseudo", "soj")
NN2_8h_depth_hsp2 <- NN2_8h_depth_hsp2[order(NN2_8h_depth_hsp2$start),]
rownames(NN2_8h_depth_hsp2) <- NULL
NN2_8h_depth_hsp2$num <- rownames(NN2_8h_depth_hsp2)
NN2_8h_depth_hsp2$num <- as.numeric(as.character(NN2_8h_depth_hsp2$num))

nn2_hsp2_plot_8h <-
  ggplot() +
  geom_segment(data=NN2_8h_depth_hsp2, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3", 
                              "darkorange2", "blue")) +
  xlab("\nGenome position") + ylab(" ") + 
  geom_segment(data=NN2_8h_gtf_hsp2, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-5, yend=-5, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-16, x = 5577799), 
             label="Lys-Anticodon-TTT-pseudo", size=2, color="blue") + 
  geom_label(data=NN2_8h_gtf_hsp1, aes(y=-16, x = 5578200), 
             label="soj", size=2, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(5577311, 5578609))


# Antisense hotspot 3
NN2_8h_depth_hsp3 <- subset(NN2_8h_depth, start > 5600220 & end < 5606993)
NN2_8h_gtf_hsp3 <- subset(import_gtf_NN2, start > 5600220 & end < 5606993)

NN2_8h_gtf_hsp3$gene_name4 <- c("CDS (undefined)","rapA", 
                                "CDS (undefined)", "CDS (undefined)", "CDS (undefined)")
NN2_8h_depth_hsp3 <- NN2_8h_depth_hsp3[order(NN2_8h_depth_hsp3$start),]
rownames(NN2_8h_depth_hsp3) <- NULL
NN2_8h_depth_hsp3$num <- rownames(NN2_8h_depth_hsp3)
NN2_8h_depth_hsp3$num <- as.numeric(as.character(NN2_8h_depth_hsp3$num))

NN2_4h_depth_hsp3 <- subset(NN2_4h_depth, start > 5600220 & end < 5606993)
NN2_4h_gtf_hsp3 <- subset(import_gtf_NN2, start > 5600220 & end < 5606993)
NN2_4h_depth_hsp3$length <- NN2_4h_depth_hsp3$end - NN2_4h_depth_hsp3$start


nn2_hsp3_plot_8h <-
  ggplot() +
  geom_segment(data=NN2_8h_depth_hsp3, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "forestgreen", "darkorange2", 
                              "forestgreen", "forestgreen", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_8h_gtf_hsp3, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-18, yend=-18, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=1) +
  geom_label(data=NN2_8h_gtf_hsp3, aes(y=-30, x = 5601520), 
             label="CDS (undefined)", size=2, color="forestgreen") + 
  geom_label(data=NN2_8h_gtf_hsp3, aes(y=-30, x = 5603420), 
             label="rapA", size=2, color="darkorange2") + 
  geom_label(data=NN2_8h_gtf_hsp3, aes(y=-30, x = 5605620), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(5600220, 5606993))


nn2_antisense_hsp_8h <-
  ggarrange(nn2_hsp1_plot_8h, nn2_hsp2_plot_8h, nn2_hsp3_plot_8h, nrow = 1,
          labels=c("D", "E", "F"))

antisense_hotspots_nn2 <- ggarrange(nn2_antisense_hsp_4h, nn2_antisense_hsp_8h, nrow=2)

#ggsave(antisense_hotspots_nn2, filename="save_figures/antisense_transcription_hotspots_nn2.tif",
 #     dpi=600, device="tiff", units="cm", width=25.5, height=21.6)


# SG17M ####
SG17M_4h_BR1_depth_fa <- read_delim("SG17M_CS_4h_BR1.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR1_depth_fa <- data.frame(SG17M_4h_BR1_depth_fa)
SG17M_4h_BR1_depth_fa$side <- "forward_antisense"
SG17M_4h_BR1_depth_ra <- read_delim("SG17M_CS_4h_BR1.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR1_depth_ra <- data.frame(SG17M_4h_BR1_depth_ra)
SG17M_4h_BR1_depth_ra$side <- "reverse_antisense"
SG17M_4h_BR1_depth <- data.frame(rbind(SG17M_4h_BR1_depth_fa, SG17M_4h_BR1_depth_ra))
SG17M_4h_BR1_depth$X6 <- NULL
SG17M_4h_BR1_depth$replicate <- "BR1"

SG17M_4h_BR2_depth_fa <- read_delim("SG17M_CS_4h_BR2.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR2_depth_fa <- data.frame(SG17M_4h_BR2_depth_fa)
SG17M_4h_BR2_depth_fa$side <- "forward_antisense"
SG17M_4h_BR2_depth_ra <- read_delim("SG17M_CS_4h_BR2.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR2_depth_ra <- data.frame(SG17M_4h_BR2_depth_ra)
SG17M_4h_BR2_depth_ra$side <- "reverse_antisense"
SG17M_4h_BR2_depth <- data.frame(rbind(SG17M_4h_BR2_depth_fa, SG17M_4h_BR2_depth_ra))
SG17M_4h_BR2_depth$X6 <- NULL
SG17M_4h_BR2_depth$replicate <- "BR2"


SG17M_4h_BR3_depth_fa <- read_delim("SG17M_CS_4h_BR3.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR3_depth_fa <- data.frame(SG17M_4h_BR3_depth_fa)
SG17M_4h_BR3_depth_fa$side <- "forward_antisense"
SG17M_4h_BR3_depth_ra <- read_delim("SG17M_CS_4h_BR3.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_4h_BR3_depth_ra <- data.frame(SG17M_4h_BR3_depth_ra)
SG17M_4h_BR3_depth_ra$side <- "reverse_antisense"
SG17M_4h_BR3_depth <- data.frame(rbind(SG17M_4h_BR3_depth_fa, SG17M_4h_BR3_depth_ra))
SG17M_4h_BR3_depth$X6 <- NULL
SG17M_4h_BR3_depth$replicate <- "BR3"
SG17M_4h_BR3_depth

SG17M_4h_depth <- data.frame(rbind(SG17M_4h_BR1_depth, SG17M_4h_BR2_depth, SG17M_4h_BR3_depth))
colnames(SG17M_4h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")



SG17M_8h_BR1_depth_fa <- read_delim("SG17M_CS_8h_BR1.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR1_depth_fa <- data.frame(SG17M_8h_BR1_depth_fa)
SG17M_8h_BR1_depth_fa$side <- "forward_antisense"
SG17M_8h_BR1_depth_ra <- read_delim("SG17M_CS_8h_BR1.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR1_depth_ra <- data.frame(SG17M_8h_BR1_depth_ra)
SG17M_8h_BR1_depth_ra$side <- "reverse_antisense"
SG17M_8h_BR1_depth <- data.frame(rbind(SG17M_8h_BR1_depth_fa, SG17M_8h_BR1_depth_ra))
SG17M_8h_BR1_depth$X6 <- NULL
SG17M_8h_BR1_depth$replicate <- "BR1"

SG17M_8h_BR2_depth_fa <- read_delim("SG17M_CS_8h_BR2.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR2_depth_fa <- data.frame(SG17M_8h_BR2_depth_fa)
SG17M_8h_BR2_depth_fa$side <- "forward_antisense"
SG17M_8h_BR2_depth_ra <- read_delim("SG17M_CS_8h_BR2.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR2_depth_ra <- data.frame(SG17M_8h_BR2_depth_ra)
SG17M_8h_BR2_depth_ra$side <- "reverse_antisense"
SG17M_8h_BR2_depth <- data.frame(rbind(SG17M_8h_BR2_depth_fa, SG17M_8h_BR2_depth_ra))
SG17M_8h_BR2_depth$X6 <- NULL
SG17M_8h_BR2_depth$replicate <- "BR2"


SG17M_8h_BR3_depth_fa <- read_delim("SG17M_CS_8h_BR3.sorted.forward_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR3_depth_fa <- data.frame(SG17M_8h_BR3_depth_fa)
SG17M_8h_BR3_depth_fa$side <- "forward_antisense"
SG17M_8h_BR3_depth_ra <- read_delim("SG17M_CS_8h_BR3.sorted.reverse_antisense.bed",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SG17M_8h_BR3_depth_ra <- data.frame(SG17M_8h_BR3_depth_ra)
SG17M_8h_BR3_depth_ra$side <- "reverse_antisense"
SG17M_8h_BR3_depth <- data.frame(rbind(SG17M_8h_BR3_depth_fa, SG17M_8h_BR3_depth_ra))
SG17M_8h_BR3_depth$X6 <- NULL
SG17M_8h_BR3_depth$replicate <- "BR3"


SG17M_8h_depth <- data.frame(rbind(SG17M_8h_BR1_depth, SG17M_8h_BR2_depth, SG17M_8h_BR3_depth))
colnames(SG17M_8h_depth) <- c("isolate", "start", "end", "read_id", "qual", "side", "replicate")

# Antisense hotspot 1, 4h
SG17M_4h_depth_hsp1 <- subset(SG17M_4h_depth, start > 990853 & end < 992602)
SG17M_4h_gtf_hsp1 <- subset(import_gtf_sg17m, start > 990853 & end < 992602)

SG17M_4h_gtf_hsp1$gene_name4 <- c("tRNA-Trp-Anticodon-CCA", "CDS (undefined)", "nusG")
SG17M_4h_depth_hsp1 <- SG17M_4h_depth_hsp1[order(SG17M_4h_depth_hsp1$start),]
rownames(SG17M_4h_depth_hsp1) <- NULL
SG17M_4h_depth_hsp1$num <- rownames(SG17M_4h_depth_hsp1)
SG17M_4h_depth_hsp1$num <- as.numeric(as.character(SG17M_4h_depth_hsp1$num))

SG17M_hsp1_plot <-
  ggplot() +
  geom_segment(data=SG17M_4h_depth_hsp1, 
               aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "forestgreen", "darkorange2", "blue", "blue")) +
  xlab("\n") + ylab("Read IDs (SG17M-4h antisense transcripts)") + 
  geom_segment(data=SG17M_4h_gtf_hsp1, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-30, yend=-30, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-100, x = 991105), 
             label="Trp-Anticodon-CCA", size=2, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-100, x = 991577), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  geom_label(data=SG17M_4h_gtf_hsp1, aes(y=-100, x = 991977), 
             label="nusG", size=2, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(990853, 992602))


# Antisense hotspot 2
SG17M_4h_depth_hsp2 <- subset(SG17M_4h_depth, start > 5814083 & end < 5815379)
SG17M_4h_gtf_hsp2 <- subset(import_gtf_sg17m, start > 5814083 & end < 5815379)
SG17M_4h_gtf_hsp2$gene_name4 <- c("tRNA-Lys-Anticodon-TTT-pseudo", "soj")
SG17M_4h_depth_hsp2 <- SG17M_4h_depth_hsp2[order(SG17M_4h_depth_hsp2$start),]
rownames(SG17M_4h_depth_hsp2) <- NULL
SG17M_4h_depth_hsp2$num <- rownames(SG17M_4h_depth_hsp2)
SG17M_4h_depth_hsp2$num <- as.numeric(as.character(SG17M_4h_depth_hsp2$num))

SG17M_hsp2_plot <-
  ggplot() +
  geom_segment(data=SG17M_4h_depth_hsp2, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "darkorange2", "blue")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=SG17M_4h_gtf_hsp2, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-25, yend=-25, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-100, x = 5814343), 
             label="Lys-Anticodon-TTT-pseudo", size=2, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-100, x = 5814923), 
             label="soj", size=2, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(5814083, 5815379))


# Antisense hotspot 3
SG17M_4h_depth_hsp3 <- subset(SG17M_4h_depth, start > 5836292 & end < 5844865)
SG17M_4h_gtf_hsp3 <- subset(import_gtf_sg17m, start > 5836292 & end < 5844865)
SG17M_4h_depth_hsp3$length <- SG17M_4h_depth_hsp3$end - SG17M_4h_depth_hsp3$start

SG17M_4h_gtf_hsp3$gene_name4 <- c("CDS (undefined)","rapA", 
                                "CDS (undefined)", "CDS (undefined)", "CDS (undefined)")
SG17M_4h_depth_hsp3 <- SG17M_4h_depth_hsp3[order(SG17M_4h_depth_hsp3$start),]
rownames(SG17M_4h_depth_hsp3) <- NULL
SG17M_4h_depth_hsp3$num <- rownames(SG17M_4h_depth_hsp3)
SG17M_4h_depth_hsp3$num <- as.numeric(as.character(SG17M_4h_depth_hsp3$num))

SG17M_hsp3_plot <-
  ggplot() +
  geom_segment(data=SG17M_4h_depth_hsp3, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "forestgreen", "darkorange2", 
                              "forestgreen", "forestgreen", "forestgreen")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=SG17M_4h_gtf_hsp3, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-50, yend=-50, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=1) +
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-180, x = 5838292), 
             label="CDS (undefined)", size=2, color="forestgreen") + 
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-180, x = 5840292), 
             label="rapA", size=2, color="darkorange2") + 
  geom_label(data=SG17M_4h_gtf_hsp3, aes(y=-180, x = 5842392), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(label=scales::comma, limits=c(5836292, 5844865))


# Antisense hotspot 1, 8h
SG17M_8h_depth_hsp1 <- subset(SG17M_8h_depth, start > 990853 & end < 992602)
SG17M_8h_gtf_hsp1 <- subset(import_gtf_sg17m, start > 990853 & end < 992602)


SG17M_8h_gtf_hsp1$gene_name4 <- c("tRNA-Trp-Anticodon-CCA", "CDS (undefined)", "nusG")
SG17M_8h_depth_hsp1 <- SG17M_8h_depth_hsp1[order(SG17M_8h_depth_hsp1$start),]
rownames(SG17M_8h_depth_hsp1) <- NULL
SG17M_8h_depth_hsp1$num <- rownames(SG17M_8h_depth_hsp1)
SG17M_8h_depth_hsp1$num <- as.numeric(as.character(SG17M_8h_depth_hsp1$num))

SG17M_hsp1_plot_8h <-
  ggplot() +
  geom_segment(data=SG17M_8h_depth_hsp1, 
               aes(x=start, xend=end, y=num, yend=num, color=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "forestgreen", "darkorange2", "blue", "blue")) +
  xlab("\n  ") + ylab("Read IDs (SG17M-8h antisense transcripts)") + 
  geom_segment(data=SG17M_8h_gtf_hsp1, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-8, yend=-8, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=SG17M_8h_gtf_hsp1, aes(y=-20, x = 991105), 
             label="Trp-Anticodon-CCA", size=2, color="blue") + 
  geom_label(data=SG17M_8h_gtf_hsp1, aes(y=-20, x = 991577), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  geom_label(data=SG17M_8h_gtf_hsp1, aes(y=-20, x = 991977), 
             label="nusG", size=2, color="darkorange2") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(990853, 992602))


# Antisense hotspot 2
SG17M_8h_depth_hsp2 <- subset(SG17M_8h_depth, start > 5814008 & end < 5815379)
SG17M_8h_gtf_hsp2 <- subset(import_gtf_sg17m, start > 5814008 & end < 5815379)
SG17M_8h_gtf_hsp2$gene_name4 <- c("tRNA-Lys-Anticodon-TTT-pseudo", "soj")
SG17M_8h_depth_hsp2 <- SG17M_8h_depth_hsp2[order(SG17M_8h_depth_hsp2$start),]
rownames(SG17M_8h_depth_hsp2) <- NULL
SG17M_8h_depth_hsp2$num <- rownames(SG17M_8h_depth_hsp2)
SG17M_8h_depth_hsp2$num <- as.numeric(as.character(SG17M_8h_depth_hsp2$num))

SG17M_hsp2_plot_8h <-
  ggplot() +
  geom_segment(data=SG17M_8h_depth_hsp2, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "darkorange2", "blue")) +
  xlab("Genome position") + ylab(" ") + 
  geom_segment(data=SG17M_4h_gtf_hsp2, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-10, yend=-10, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2) +
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-18, x = 5814343), 
             label="Lys-Anticodon-TTT-pseudo", size=2, color="blue") + 
  geom_label(data=SG17M_4h_gtf_hsp2, aes(y=-18, x = 5814923), 
             label="soj", size=2, color="darkorange2") +
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

SG17M_8h_gtf_hsp3$gene_name4 <- c("CDS (undefined)","rapA", 
                                "CDS (undefined)", "CDS (undefined)", "CDS (undefined)")
SG17M_8h_depth_hsp3 <- SG17M_8h_depth_hsp3[order(SG17M_8h_depth_hsp3$start),]
rownames(SG17M_8h_depth_hsp3) <- NULL
SG17M_8h_depth_hsp3$num <- rownames(SG17M_8h_depth_hsp3)
SG17M_8h_depth_hsp3$num <- as.numeric(as.character(SG17M_8h_depth_hsp3$num))

SG17M_hsp3_plot_8h <-
  ggplot() +
  geom_segment(data=SG17M_8h_depth_hsp3, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "forestgreen", "darkorange2", 
                              "forestgreen", "forestgreen", "forestgreen")) +
  xlab(" ") + ylab(" ") + 
  geom_segment(data=SG17M_8h_gtf_hsp3, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-20, yend=-20, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=1) +
  geom_label(data=SG17M_8h_gtf_hsp3, aes(y=-45, x = 5838292), 
             label="CDS (undefined)", size=2, color="forestgreen") + 
  geom_label(data=SG17M_8h_gtf_hsp3, aes(y=-45, x = 5840292), 
             label="rapA", size=2, color="darkorange2") + 
  geom_label(data=SG17M_8h_gtf_hsp3, aes(y=-45, x = 5842392), 
             label="CDS (undefined)", size=2, color="forestgreen") +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(5836292, 5844865))

sg17m_4h_antisense_hotspots <-
  ggarrange(SG17M_hsp1_plot, SG17M_hsp2_plot, SG17M_hsp3_plot, nrow= 1,
          labels=c("A", "B", "C"))
sg17m_8h_antisense_hotspots <-
  ggarrange(SG17M_hsp1_plot_8h, SG17M_hsp2_plot_8h, SG17M_hsp3_plot_8h, nrow= 1,
          labels=c("D", "E", "F"))

sg17m_antisense_hsp <-
  ggarrange(sg17m_4h_antisense_hotspots, sg17m_8h_antisense_hotspots, nrow=2)


#ggsave(sg17m_antisense_hsp, filename="save_figures/antisense_transcription_hotspots_sg17m.tif",
  #    dpi=600, device="tiff", units="cm", width=25.5, height=21.6)


# differentially expressed
NN2_4h_depth_lasR <- subset(NN2_4h_depth, start > 4223560 & end < 4226100)
NN2_4h_gtf_lasR <- subset(import_gtf_NN2, start > 4223560 & end < 4226100)
NN2_4h_gtf_lasR
NN2_4h_gtf_lasR$gene_name4 <- c("lasI","exon","CDS","lasR")
NN2_4h_depth_lasR <- NN2_4h_depth_lasR[order(NN2_4h_depth_lasR$start),]
rownames(NN2_4h_depth_lasR) <- NULL
NN2_4h_depth_lasR$num <- rownames(NN2_4h_depth_lasR)
NN2_4h_depth_lasR$num <- as.numeric(as.character(NN2_4h_depth_lasR$num))


nn2_lasR_plot_4h <-
  ggplot() +
  geom_segment(data=NN2_4h_depth_lasR, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "darkorange2", "darkblue", 
                              "darkorange2", "darkorange2", "darkorange2")) +
  xlab("\n") + ylab("Antisense transcripts (NN2)") + 
  geom_segment(data=NN2_4h_gtf_lasR, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-2, yend=-2, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=0.7) +
  geom_label(data=NN2_4h_gtf_lasR, aes(y=-6, x = 4224021), 
             label="lasI", size=3, color="darkorange2") + 
  geom_label(data=NN2_4h_gtf_lasR, aes(y=-6, x = 4224350), 
             label="ncRNA", size=3, color="darkblue") + 
  geom_label(data=NN2_4h_gtf_lasR, aes(y=-6, x = 4224591), 
             label="rsaL", size=3, color="darkorange2") +
  geom_label(data=NN2_4h_gtf_lasR, aes(y=-6, x = 4224980), 
             label="lasR", size=3, color="darkorange2") +
  geom_label(data=NN2_4h_gtf_lasR, aes(x=4223821, y=54), label="NN2-4h", size=3) +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(4223650, 4225700))


# differentially expressed
NN2_8h_depth_lasR <- subset(NN2_8h_depth, start > 4223560 & end < 4226100)
NN2_8h_gtf_lasR <- subset(import_gtf_NN2, start > 4223560 & end < 4226100)
NN2_8h_gtf_lasR
NN2_8h_gtf_lasR$gene_name4 <- c("lasI","exon","CDS","lasR")
NN2_8h_depth_lasR <- NN2_8h_depth_lasR[order(NN2_8h_depth_lasR$start),]
rownames(NN2_8h_depth_lasR) <- NULL
NN2_8h_depth_lasR$num <- rownames(NN2_8h_depth_lasR)
NN2_8h_depth_lasR$num <- as.numeric(as.character(NN2_8h_depth_lasR$num))


nn2_lasR_plot_8h <-
  ggplot() +
  geom_segment(data=NN2_8h_depth_lasR, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "grey", "lightcyan3",
                              "darkorange2", "darkblue", 
                              "darkorange2", "darkorange2", "darkorange2")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=NN2_8h_gtf_lasR, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=-2, yend=-2, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=0.7) +
  geom_label(data=NN2_8h_gtf_lasR, aes(y=-8, x = 4224021), 
             label="lasI", size=3, color="darkorange2") + 
  geom_label(data=NN2_8h_gtf_lasR, aes(y=-8, x = 4224350), 
             label="ncRNA", size=3, color="darkblue") + 
  geom_label(data=NN2_8h_gtf_lasR, aes(y=-8, x = 4224591), 
             label="rsaL", size=3, color="darkorange2") +
  geom_label(data=NN2_8h_gtf_lasR, aes(y=-8, x = 4224980), 
             label="lasR", size=3, color="darkorange2") +
  geom_label(data=NN2_8h_gtf_lasR, aes(x=4223821,y=90), label="NN2-8h", size=3) +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(4223650, 4225700))

SG17M_4h_depth_lasR <- subset(SG17M_4h_depth, start > 4691366 & end < 4693885)
SG17M_4h_gtf_lasR <- subset(import_gtf_sg17m, start > 4691366 & end < 4693885)
SG17M_4h_gtf_lasR$gene_name4 <- c("lasI","exon","CDS","lasR")
SG17M_4h_depth_lasR <- SG17M_4h_depth_lasR[order(SG17M_4h_depth_lasR$start),]
rownames(SG17M_4h_depth_lasR) <- NULL
SG17M_4h_depth_lasR$num <- rownames(SG17M_4h_depth_lasR)
SG17M_4h_depth_lasR$num <- as.numeric(as.character(SG17M_4h_depth_lasR$num))


SG17M_lasR_plot_4h <-
  ggplot() +
  geom_segment(data=SG17M_4h_depth_lasR, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", 
                              "darkorange2", "darkblue", "darkorange2", 
                              "darkorange2")) +
  xlab("\n") + ylab("Antisense transcripts (SG17M)") + 
  geom_segment(data=SG17M_4h_gtf_lasR, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=0.42, yend=0.42, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=0.7) +
  geom_label(data=SG17M_4h_gtf_lasR, aes(y=0, x = 4691766), 
             label="lasI", size=3, color="darkorange2") + 
  geom_label(data=SG17M_4h_gtf_lasR, aes(y=0, x = 4692150), 
             label="ncRNA", size=3, color="darkblue") + 
  geom_label(data=SG17M_4h_gtf_lasR, aes(y=0, x = 4692350), 
             label="rsaL", size=3, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_lasR, aes(y=0, x = 4692815), 
             label="lasR", size=3, color="darkorange2") +
  geom_label(data=SG17M_4h_gtf_lasR, aes(x=4691600, y=6.5), label="SG17M-4h", size=3) +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(4691466, 4693385))


SG17M_8h_depth_lasR <- subset(SG17M_8h_depth, start > 4691366 & end < 4693885)
SG17M_8h_gtf_lasR <- subset(import_gtf_sg17m, start > 4691366 & end < 4693885)

SG17M_8h_gtf_lasR$gene_name4 <- c("lasI","exon","CDS","lasR")
SG17M_8h_depth_lasR <- SG17M_8h_depth_lasR[order(SG17M_8h_depth_lasR$start),]
rownames(SG17M_8h_depth_lasR) <- NULL
SG17M_8h_depth_lasR$num <- rownames(SG17M_8h_depth_lasR)
SG17M_8h_depth_lasR$num <- as.numeric(as.character(SG17M_8h_depth_lasR$num))

SG17M_lasR_plot_8h <-
  ggplot() +
  geom_segment(data=SG17M_8h_depth_lasR, 
               aes(x=start, xend=end, y=num, yend=num, colour=replicate), size=1, alpha=1) + 
  scale_color_manual(values=c("gray22", "lightcyan3",
                              "darkorange2", "darkblue", "darkorange2", 
                              "darkorange2")) +
  xlab("\n") + ylab(" ") + 
  geom_segment(data=SG17M_8h_gtf_lasR, 
               aes(x=ifelse(strandType == "+", start, ifelse(strandType == "\\.", start, end)), 
                   xend=ifelse(strandType == "+", end, ifelse(strandType == "\\.", end, start)),
                   y=0.32, yend=0.32, color=gene_name4), 
               arrow = arrow(length = unit(0.06, "inches")), size = 2, alpha=0.7) +
  geom_label(data=SG17M_8h_gtf_lasR, aes(y=-1.8, x = 4691766), 
             label="lasI", size=3, color="darkorange2") + 
  geom_label(data=SG17M_8h_gtf_lasR, aes(y=-1.8, x = 4692150), 
             label="ncRNA", size=3, color="darkblue") + 
  geom_label(data=SG17M_8h_gtf_lasR, aes(y=-1.8, x = 4692350), 
             label="rsaL", size=3, color="darkorange2") +
  geom_label(data=SG17M_8h_gtf_lasR, aes(y=-1.8, x = 4692815), 
             label="lasR", size=3, color="darkorange2") +
  geom_label(data=SG17M_8h_gtf_lasR, aes(x=4691600, y=30), label="SG17M-8h", size=3) +
  theme_pubr(border=TRUE, base_size=8, legend="none") + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_x_continuous(label=scales::comma, limits=c(4691466, 4693385))

lasR_compare <-
  ggarrange(nn2_lasR_plot_4h, nn2_lasR_plot_8h,
          SG17M_lasR_plot_4h, SG17M_lasR_plot_8h, labels = c("A", "B", "C", "D"),
          nrow=2, ncol=2)

#ggsave(lasR_compare, filename="save_figures/antisense_transcription_lasR.tif",
 #     dpi=600, device="tiff", units="cm", width=25.5, height=21.6)


table(SG17M_4h_depth_lasR$replicate) # 7, 0, 0 = 2.3
table(SG17M_8h_depth_lasR$replicate) # 12, 20, 0 = 10.7
table(NN2_4h_depth_lasR$replicate) # 17, 12, 29 = 19.3
table(NN2_8h_depth_lasR$replicate) # 13 #42 #42 = 32.3
