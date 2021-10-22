# title: "Nanopore-Illumina comparison"
# author: "Marie Pust"
# date: "26 9 2021"


############################################################################################################
# clean global environment
rm(list=ls())

# set working directory
setwd("/R")

############################################################################################################

# define global variables
# set global seed
set.seed(111)

# function for installing/importing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

# bootstrap the standard error of the median
bootstrap_median <- function(data, num) {
  resamples <- lapply(1:num, function(i) sample(data, replace=T))
  r.median <- sapply(resamples, median)
  std.err <- sqrt(var(r.median))
  list(std.err=std.err, resamples=resamples, medians=r.median)}

# function for correlation test
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
  }

# store required packages
packages <- c('rlist','ggplot2','readr','readxl','Rsubread','edgeR','stringr','tidyr','ggpubr','circlize','dplyr','Rsamtools',
              'GenomicRanges','GenomicAlignments','vegan','rcompanion', 'Hmisc', 'corrplot', 'pheatmap', 'TmCalculator', 'seqinr', 
              'plyr', 'DESeq2', 'RColorBrewer', 'ggplotify', 'ashr')

# store paths to files
ont_nn2_bam <- 'NN2_run/d_raw_bam'
ont_sg17m_bam <- 'SG17M_run/d_raw_bam'
illumina_nn2_bam <- 'Illumina_data/bam_files'
UserDefinedAnnotationRef_NN2 <- 'Illumina_data/NN2_ENO_curated.gtf'
UserDefinedAnnotationRef_SG17M <- 'Illumina_data/SG17M_ENO_curated.gtf'

############################################################################################################
# load required R packages
ipak(packages)
############################################################################################################

# import and clean annotation file, NN2
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
nn2_rRNA_id <- subset(import_gtf_NN2, gene_name3 == "rRNA")
nn2_rRNA_id_vector <- nn2_rRNA_id$transcript_id


# import and clean annotation file, SG17M
import_gtf_SG17M <- read.table(UserDefinedAnnotationRef_SG17M, sep = '\t', header = FALSE)
colnames(import_gtf_SG17M) <- c("Isolate", "database", "featureType", "start", "end", "V6", "strandType", "V8", "V9")
import_gtf_SG17M$V6 <- NULL
import_gtf_SG17M$V8 <- NULL
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "transcript_id ", "")
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "gene_name ", "")
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "tRNA_type ", "")
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "rRNA_type ", "")
import_gtf_SG17M <- separate(data = import_gtf_SG17M, col = V9, into = c("transcript_id", "gene_name2"), sep = "\\;")

import_gtf_SG17M$length <- import_gtf_SG17M$end - import_gtf_SG17M$start
import_gtf_SG17M$gene_name3 <-
  with(import_gtf_SG17M, 
       ifelse(database == "barrnap:0.9", "rRNA",
              ifelse(database == "Infernal:001001", "ncRNA",
                     ifelse(database == "Aragorn:001002", "tRNA",
                            ifelse(database == "IslandViewer4", "GI",
                                   ifelse(database == "Prodigal:002006", "CDS",database))))))

sg17m_accessory <- subset(import_gtf_SG17M, gene_name3 == "GI")
sg17m_accessory_bed <- select(sg17m_accessory, c(start, end))
sg17m_tRNA <- subset(import_gtf_SG17M, gene_name3 == "tRNA")
sg17m_tRNA_bed <- select(sg17m_accessory, c(start, end))
sg17m_rRNA_id <- subset(import_gtf_SG17M, gene_name3 == "rRNA")
sg17m_rRNA_id_vector <- sg17m_rRNA_id$transcript_id


############################################################################################################
# import illumina datasets
ill_input_bam_NN2 <- list.files(path = illumina_nn2_bam, pattern = ".bam", full.names = TRUE)

# obtain count data for illumina antisense reads
ill_featureCounts_antisenseNN2 <- 
  featureCounts(ill_input_bam_NN2, annot.ext = UserDefinedAnnotationRef_NN2, 
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 1, isLongRead = FALSE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

# extract counts and convert to data frame
ill_featureCounts_antisenseNN2 <- ill_featureCounts_antisenseNN2$counts
ill_featureCounts_antisenseNN2 <- data.frame(ill_featureCounts_antisenseNN2)
# rename columns
colnames(ill_featureCounts_antisenseNN2) <- c("TEX_NN2_8h_1", "TEX_NN2_8h_2", "TEX_NN2_4h_1", "TEX_NN2_4h_2", 
                                              "0TEX_NN2_8h_1", "0TEX_NN2_8h_2", "0TEX_NN2_4h_1", "0TEX_NN2_4h_2")

# obtain count data for illumina sense reads
ill_featureCounts_senseNN2 <- 
  featureCounts(ill_input_bam_NN2, annot.ext = UserDefinedAnnotationRef_NN2,
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 2, isLongRead = FALSE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

# extract counts and convert to data frame
ill_featureCounts_senseNN2 <- ill_featureCounts_senseNN2$counts
ill_featureCounts_senseNN2 <- data.frame(ill_featureCounts_senseNN2)

# rename columns
colnames(ill_featureCounts_senseNN2) <- c("TEX_NN2_8h_1", "TEX_NN2_8h_2", "TEX_NN2_4h_1", "TEX_NN2_4h_2", 
                                          "0TEX_NN2_8h_1", "0TEX_NN2_8h_2", "0TEX_NN2_4h_1", "0TEX_NN2_4h_2")


############################################################################################################
# import ont datasets, NN2
ont_input_bam_NN2 <- list.files(path = ont_nn2_bam, pattern = ".bam", full.names = TRUE)

# obtain count data for  antisense reads
ont_featureCounts_antisenseNN2 <- 
  featureCounts(ont_input_bam_NN2, annot.ext = UserDefinedAnnotationRef_NN2, 
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 2, isLongRead = TRUE,  
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

# extract counts and convert to data frame
ont_featureCounts_antisenseNN2 <- ont_featureCounts_antisenseNN2$counts
ont_featureCounts_antisenseNN2 <- data.frame(ont_featureCounts_antisenseNN2)
# rename columns
colnames(ont_featureCounts_antisenseNN2) <- c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3", 
                                              "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3")

# obtain count data for ont sense reads
ont_featureCounts_senseNN2 <- 
  featureCounts(ont_input_bam_NN2, annot.ext = UserDefinedAnnotationRef_NN2,
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 1, isLongRead = FALSE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

# extract counts and convert to data frame
ont_featureCounts_senseNN2 <- ont_featureCounts_senseNN2$counts
ont_featureCounts_senseNN2 <- data.frame(ont_featureCounts_senseNN2)

# rename columns
colnames(ont_featureCounts_senseNN2) <- c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3", 
                                          "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3")

############################################################################################################

# import ont datasets, SG17M
ont_input_bam_SG17M <- list.files(path = ont_sg17m_bam, pattern = ".bam", full.names = TRUE)

# obtain count data for illumina antisense reads
ont_featureCounts_antisenseSG17M <- 
  featureCounts(ont_input_bam_SG17M, annot.ext = UserDefinedAnnotationRef_SG17M, 
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 2, isLongRead = TRUE,  
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

# extract counts and convert to data frame
ont_featureCounts_antisenseSG17M <- ont_featureCounts_antisenseSG17M$counts
ont_featureCounts_antisenseSG17M <- data.frame(ont_featureCounts_antisenseSG17M)

# rename columns
colnames(ont_featureCounts_antisenseSG17M) <- c("SG17M_4h_BR1", "SG17M_4h_BR2", "SG17M_4h_BR3", 
                                                "SG17M_8h_BR1", "SG17M_8h_BR2", "SG17M_8h_BR3")

# obtain count data for ont sense reads
ont_featureCounts_senseSG17M <- 
  featureCounts(ont_input_bam_SG17M, annot.ext = UserDefinedAnnotationRef_SG17M,
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("CDS","exon"), strandSpecific = 1, isLongRead = FALSE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

# extract counts and convert to data frame
ont_featureCounts_senseSG17M <- ont_featureCounts_senseSG17M$counts
ont_featureCounts_senseSG17M <- data.frame(ont_featureCounts_senseSG17M)

# rename columns
colnames(ont_featureCounts_senseSG17M) <- c("SG17M_4h_BR1", "SG17M_4h_BR2", "SG17M_4h_BR3", 
                                            "SG17M_8h_BR1", "SG17M_8h_BR2", "SG17M_8h_BR3")

############################################################################################################

# merge illumina and ont output, antisense
featureCounts_antisenseNN2 <- data.frame(cbind(ont_featureCounts_antisenseNN2, ill_featureCounts_antisenseNN2))
featureCounts_antisenseNN2_sorted_0 <- select(featureCounts_antisenseNN2, 
                                              c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3", "TEX_NN2_4h_1", "TEX_NN2_4h_2", 
                                                "X0TEX_NN2_4h_1", "X0TEX_NN2_4h_2", "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3",
                                                "TEX_NN2_8h_1", "TEX_NN2_8h_2", "X0TEX_NN2_8h_1", "X0TEX_NN2_8h_2"))

featureCounts_antisenseNN2_sorted_1 <- featureCounts_antisenseNN2_sorted_0[rownames(featureCounts_antisenseNN2_sorted_0) %in% import_gtf_NN2$transcript_id,]
library_sizes = colSums(featureCounts_antisenseNN2_sorted_1)
featureCounts_antisenseNN2_sorted_2 <- featureCounts_antisenseNN2_sorted_1 / library_sizes

add_groups <- c("4h", "4h", "4h", "4", "4h", "4h", "4h", "8h", "8h", "8h", "8h", "8h_", "8h", "8h")
featureCounts_antisenseNN2_sorted_3 <- DGEList(counts=featureCounts_antisenseNN2_sorted_2, lib.size = library_sizes, group = add_groups, remove.zeros = TRUE)
featureCounts_antisenseNN2_sorted_4 <- calcNormFactors(featureCounts_antisenseNN2_sorted_3, method="TMM")
featureCounts_antisenseNN2_sorted_5 <- cpm(featureCounts_antisenseNN2_sorted_4, log=TRUE)
featureCounts_antisenseNN2_sorted_6 <- as.data.frame(featureCounts_antisenseNN2_sorted_5)

ONT_4h <- select(featureCounts_antisenseNN2_sorted_6, c(NN2_4h_BR1, NN2_4h_BR2, NN2_4h_BR3))
ONT_4h$ONT_4h <- apply(ONT_4h, 1, median, na.rm=T)

ONT_8h <- select(featureCounts_antisenseNN2_sorted_6, c(NN2_8h_BR1, NN2_8h_BR2, NN2_8h_BR3))
ONT_8h$ONT_8h <- apply(ONT_8h, 1, median, na.rm=T)

Ill_TEX_4h <- select(featureCounts_antisenseNN2_sorted_6, c(TEX_NN2_4h_1, TEX_NN2_4h_2))
Ill_TEX_4h$Ill_TEX_4h <- apply(Ill_TEX_4h, 1, median, na.rm=T)

Ill_TEX_8h <- select(featureCounts_antisenseNN2_sorted_6, c(TEX_NN2_8h_1, TEX_NN2_8h_2))
Ill_TEX_8h$Ill_TEX_8h <- apply(Ill_TEX_8h, 1, median, na.rm=T)

Ill_0TEX_4h <- select(featureCounts_antisenseNN2_sorted_6, c(X0TEX_NN2_4h_1, X0TEX_NN2_4h_2))
Ill_0TEX_4h$Ill_0TEX_4h <- apply(Ill_0TEX_4h, 1, median, na.rm=T)

Ill_0TEX_8h <- select(featureCounts_antisenseNN2_sorted_6, c(X0TEX_NN2_8h_1, X0TEX_NN2_8h_2))
Ill_0TEX_8h$Ill_0TEX_8h <- apply(Ill_0TEX_8h, 1, median, na.rm=T)


featureCounts_antisenseNN2_sorted_7 <- data.frame(cbind(ONT_4h$ONT_4h, ONT_8h$ONT_8h, Ill_TEX_4h$Ill_TEX_4h, 
                                                        Ill_TEX_8h$Ill_TEX_8h, Ill_0TEX_4h$Ill_0TEX_4h, 
                                                        Ill_0TEX_8h$Ill_0TEX_8h))


###########################################################################################################

# merge illumina and ont output, sense
featureCounts_senseNN2 <- data.frame(cbind(ont_featureCounts_senseNN2, ill_featureCounts_senseNN2))
featureCounts_senseNN2_sorted_0 <- select(featureCounts_senseNN2, c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3",
                                                                    "TEX_NN2_4h_1", "TEX_NN2_4h_2", 
                                                                    "X0TEX_NN2_4h_1", "X0TEX_NN2_4h_2",
                                                                    "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3",
                                                                    "TEX_NN2_8h_1", "TEX_NN2_8h_2", 
                                                                    "X0TEX_NN2_8h_1", "X0TEX_NN2_8h_2"))
# remove rRNA information
featureCounts_senseNN2_sorted_1 <- featureCounts_senseNN2_sorted_0[rownames(featureCounts_senseNN2_sorted_0) %in% import_gtf_NN2$transcript_id,]
library_sizes = colSums(featureCounts_senseNN2_sorted_1)
featureCounts_senseNN2_sorted_2 <- featureCounts_senseNN2_sorted_1 / library_sizes
featureCounts_senseNN2_sorted_3 <- DGEList(counts=featureCounts_senseNN2_sorted_2, lib.size = library_sizes, group = add_groups, remove.zeros = TRUE)
featureCounts_senseNN2_sorted_4 <- calcNormFactors(featureCounts_senseNN2_sorted_3, method="TMM")
featureCounts_senseNN2_sorted_5 <- cpm(featureCounts_senseNN2_sorted_4, log=TRUE)
featureCounts_senseNN2_sorted_6 <- as.data.frame(featureCounts_senseNN2_sorted_5)

ONT_4h_s <- select(featureCounts_senseNN2_sorted_6, c(NN2_4h_BR1, NN2_4h_BR2, NN2_4h_BR3))
ONT_4h_s$ONT_4h_s <- apply(ONT_4h_s, 1, median, na.rm=T)

ONT_8h_s <- select(featureCounts_senseNN2_sorted_6, c(NN2_8h_BR1, NN2_8h_BR2, NN2_8h_BR3))
ONT_8h_s$ONT_8h_s_s <- apply(ONT_8h_s, 1, median, na.rm=T)

Ill_TEX_4h_s <- select(featureCounts_senseNN2_sorted_6, c(TEX_NN2_4h_1, TEX_NN2_4h_2))
Ill_TEX_4h_s$Ill_TEX_4h_s <- apply(Ill_TEX_4h_s, 1, median, na.rm=T)

Ill_TEX_8h_s <- select(featureCounts_senseNN2_sorted_6, c(TEX_NN2_8h_1, TEX_NN2_8h_2))
Ill_TEX_8h_s$Ill_TEX_8h_s <- apply(Ill_TEX_8h_s, 1, median, na.rm=T)


Ill_0TEX_4h_s <- select(featureCounts_senseNN2_sorted_6, c(X0TEX_NN2_4h_1, X0TEX_NN2_4h_2))
Ill_0TEX_4h_s$Ill_0TEX_4h_s <- apply(Ill_0TEX_4h_s, 1, median, na.rm=T)


Ill_0TEX_8h_s <- select(featureCounts_senseNN2_sorted_6, c(X0TEX_NN2_8h_1, X0TEX_NN2_8h_2))
Ill_0TEX_8h_s$Ill_0TEX_8h_s <- apply(Ill_0TEX_8h_s, 1, median, na.rm=T)


featureCounts_senseNN2_sorted_7 <- data.frame(cbind(ONT_4h_s$ONT_4h_s, ONT_8h_s$ONT_8h_s, Ill_TEX_4h_s$Ill_TEX_4h_s, 
                                                    Ill_TEX_8h_s$Ill_TEX_8h_s, Ill_0TEX_4h_s$Ill_0TEX_4h_s, Ill_0TEX_8h_s$Ill_0TEX_8h_s))

rownames(featureCounts_senseNN2_sorted_7) <- rownames(ONT_4h_s)
colnames(featureCounts_senseNN2_sorted_7) <- c('ONT_4h','ONT_8h','Ill_TEX_4h', 'Ill_TEX_8h', 'Ill_0TEX_4h', 'Ill_0TEX_8h')


# merge tRNAs, NN2
import_gtf_NN2$gene_name4 <- ifelse(import_gtf_NN2$gene_name2 == "", import_gtf_NN2$featureType, import_gtf_NN2$gene_name2)
my_list_2_NN2 <- import_gtf_NN2$gene_name4
names(my_list_2_NN2) <- import_gtf_NN2$transcript_id

mylist_antisense_x_NN2 <- c()
for(i in rownames(featureCounts_antisenseNN2)){
    newelem <- ifelse(i %in% names(my_list_2_NN2), my_list_2_NN2[[i]], i)
    mylist_antisense_x_NN2 <- append(mylist_antisense_x_NN2, newelem)}
featureCounts_antisenseNN2$name2 <- mylist_antisense_x_NN2
tRNA_antisense_NN2 <- featureCounts_antisenseNN2[grepl( "tRNA" , featureCounts_antisenseNN2$name2 ), ]
tRNA_antisense_NN2 <- tRNA_antisense_NN2 %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_antisense_NN2 <- data.frame(tRNA_antisense_NN2)
rownames(tRNA_antisense_NN2) <- tRNA_antisense_NN2$name2
tRNA_antisense_NN2$name2 <- NULL

featureCounts_antisenseNN2 <- featureCounts_antisenseNN2[!grepl("tRNA", featureCounts_antisenseNN2$name2),]
featureCounts_antisenseNN2$name2 <- NULL
featureCounts_antisenseNN2 <- data.frame(rbind(featureCounts_antisenseNN2, tRNA_antisense_NN2))

mylist_sense_x_NN2 <- c()
for(i in rownames(featureCounts_senseNN2)){
    newelem <- ifelse(i %in% names(my_list_2_NN2), my_list_2_NN2[[i]], i)
    mylist_sense_x_NN2 <- append(mylist_sense_x_NN2, newelem)}
featureCounts_senseNN2$name2 <- mylist_sense_x_NN2
tRNA_sense_NN2 <- featureCounts_senseNN2[grepl( "tRNA" , featureCounts_senseNN2$name2 ), ]
tRNA_sense_NN2 <- tRNA_sense_NN2 %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_sense_NN2 <- data.frame(tRNA_sense_NN2)
rownames(tRNA_sense_NN2) <- tRNA_sense_NN2$name2
tRNA_sense_NN2$name2 <- NULL

featureCounts_senseNN2 <- featureCounts_senseNN2[!grepl("tRNA", featureCounts_senseNN2$name2),]
featureCounts_senseNN2$name2 <- NULL
featureCounts_senseNN2 <- data.frame(rbind(featureCounts_senseNN2, tRNA_sense_NN2))


import_gtf_SG17M$gene_name4 <- ifelse(import_gtf_SG17M$gene_name2 == "", import_gtf_SG17M$featureType, import_gtf_SG17M$gene_name2)
my_list_2_SG17M <- import_gtf_SG17M$gene_name4
names(my_list_2_SG17M) <- import_gtf_SG17M$transcript_id

mylist_antisense_x_SG17M <- c()
for(i in rownames(ont_featureCounts_antisenseSG17M)){
    newelem <- ifelse(i %in% names(my_list_2_SG17M), my_list_2_SG17M[[i]], i)
    mylist_antisense_x_SG17M <- append(mylist_antisense_x_SG17M, newelem)}
ont_featureCounts_antisenseSG17M$name2 <- mylist_antisense_x_SG17M
ont_featureCounts_antisenseSG17M
tRNA_antisense_SG17M <- ont_featureCounts_antisenseSG17M[grepl( "tRNA" , ont_featureCounts_antisenseSG17M$name2 ), ]
tRNA_antisense_SG17M <- tRNA_antisense_SG17M %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_antisense_SG17M <- data.frame(tRNA_antisense_SG17M)
rownames(tRNA_antisense_SG17M) <- tRNA_antisense_SG17M$name2
tRNA_antisense_SG17M$name2 <- NULL

ont_featureCounts_antisenseSG17M <- ont_featureCounts_antisenseSG17M[!grepl("tRNA", ont_featureCounts_antisenseSG17M$name2),]
ont_featureCounts_antisenseSG17M$name2 <- NULL
ont_featureCounts_antisenseSG17M <- data.frame(rbind(ont_featureCounts_antisenseSG17M, tRNA_antisense_SG17M))


mylist_sense_x_SG17M <- c()
for(i in rownames(ont_featureCounts_senseSG17M)){
    newelem <- ifelse(i %in% names(my_list_2_SG17M), my_list_2_SG17M[[i]], i)
    mylist_sense_x_SG17M <- append(mylist_sense_x_SG17M, newelem)}
ont_featureCounts_senseSG17M$name2 <- mylist_sense_x_SG17M
tRNA_sense_SG17M <- ont_featureCounts_senseSG17M[grepl( "tRNA" , ont_featureCounts_senseSG17M$name2 ), ]
tRNA_sense_SG17M <- tRNA_sense_SG17M %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_sense_SG17M <- data.frame(tRNA_sense_SG17M)
rownames(tRNA_sense_SG17M) <- tRNA_sense_SG17M$name2
tRNA_sense_SG17M$name2 <- NULL

ont_featureCounts_senseSG17M <- ont_featureCounts_senseSG17M[!grepl("tRNA", ont_featureCounts_senseSG17M$name2),]
ont_featureCounts_senseSG17M$name2 <- NULL
ont_featureCounts_senseSG17M <- data.frame(rbind(ont_featureCounts_senseSG17M, tRNA_sense_SG17M))


# Convert count data to a matrix of appropriate form that DEseq2 can read
row_sub = apply(featureCounts_senseNN2, 1, function(row) all(row > -2 ))
featureCounts_senseNN20 <- featureCounts_senseNN2[row_sub,]
geneID_sense <- rownames(featureCounts_senseNN20)
sampleIndex_sense <- colnames(featureCounts_senseNN20)
rawCounts_sense <- as.matrix(featureCounts_senseNN20)
rownames(rawCounts_sense) <- geneID_sense

# make coldata
sample_sense <- colnames(featureCounts_senseNN2)
time_sense <- c("4h", "4h", "4h", "8h", "8h", "8h", "8h", "8h", "4h", "4h", "8h", "8h", "4h", "4h")
platform_sense <- c("ONT", "ONT", "ONT", "ONT", "ONT", "ONT", "TEX", "TEX", "TEX", "TEX", "0TEX", "0TEX", "0TEX", "0TEX")
colData_ont_s <- data.frame(cbind(sample_sense, time_sense, platform_sense))
rownames(colData_ont_s) <- colData_ont_s$sample
colData_ont_s$platform_sense <- factor(colData_ont_s$platform, levels=c("ONT", "TEX", "0TEX"))
colData_ont_s$time_sense <- factor(colData_ont_s$time, levels=c("4h", "8h"))

# Create the DEseq2 DataSet object
de_ont_s_deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts_sense, colData=colData_ont_s, design= ~ time_sense+platform_sense)

# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
de_ont_s_deseq2Data <- DESeq(de_ont_s_deseq2Data)
deseq2Results_sense__tex_0tex <- lfcShrink(de_ont_s_deseq2Data, contrast = c("platform_sense", "TEX", "0TEX"), 
                                           type="ashr", format = "DataFrame")
deseq2Results_sense__tex_0tex$padj_value <- ifelse(deseq2Results_sense__tex_0tex$padj < 0.01, "different", "not-different")
deseq2Results_sense__tex_0tex <- na.omit(deseq2Results_sense__tex_0tex)
deseq2Results_sense__tex_0tex$compare <- "TEX (reference) - 0TEX"
deseq2Results_sense__tex_0tex_diff <- subset(deseq2Results_sense__tex_0tex, padj_value == "different")
deseq2Results_sense__tex_0tex_diff_vector <- rownames(deseq2Results_sense__tex_0tex_diff)
deseq2Results_sense__tex_0tex_ndiff <- subset(deseq2Results_sense__tex_0tex, padj_value == "not-different")
deseq2Results_sense__tex_0tex_ndiff_vector <- rownames(deseq2Results_sense__tex_0tex_ndiff)


deseq2Results_sense__tex_ont <- lfcShrink(de_ont_s_deseq2Data, contrast = c("platform_sense", "ONT", "TEX"), 
                                          type="ashr", format = "DataFrame")
deseq2Results_sense__tex_ont$padj_value <- ifelse(deseq2Results_sense__tex_ont$padj < 0.01, "different", "not-different")
deseq2Results_sense__tex_ont <- na.omit(deseq2Results_sense__tex_ont)
deseq2Results_sense__tex_ont$compare <- "ONT (reference) - TEX"
deseq2Results_sense__tex_ont_diff <- subset(deseq2Results_sense__tex_ont, padj_value == "different")
deseq2Results_sense__tex_ont_diff_vector <- rownames(deseq2Results_sense__tex_ont_diff)
deseq2Results_sense__tex_ont_ndiff <- subset(deseq2Results_sense__tex_ont, padj_value == "not-different")
deseq2Results_sense__tex_ont_ndiff_vector <- rownames(deseq2Results_sense__tex_ont_ndiff)


deseq2Results_sense__0tex_ont <- lfcShrink(de_ont_s_deseq2Data, contrast = c("platform_sense", "ONT", "0TEX"), 
                                           type="ashr", format = "DataFrame")
deseq2Results_sense__0tex_ont$padj_value <- ifelse(deseq2Results_sense__0tex_ont$padj < 0.01, "different", "not-different")
deseq2Results_sense__0tex_ont <- na.omit(deseq2Results_sense__0tex_ont)
deseq2Results_sense__0tex_ont$compare <- "ONT (reference) - 0TEX"
deseq2Results_sense__0tex_ont_diff <- subset(deseq2Results_sense__0tex_ont, padj_value == "different")
deseq2Results_sense__0tex_ont_diff_vector <- rownames(deseq2Results_sense__0tex_ont_diff)
deseq2Results_sense__0tex_ont_ndiff <- subset(deseq2Results_sense__0tex_ont, padj_value == "not-different")
deseq2Results_sense__0tex_ont_ndiff_vector <- rownames(deseq2Results_sense__0tex_ont_ndiff)

boxplot(log10(assays(de_ont_s_deseq2Data)[["cooks"]]), range=0, las=2)

# merge it
deseq2Results_sense__ont <- data.frame(rbind(deseq2Results_sense__tex_0tex, 
                                             deseq2Results_sense__tex_ont, 
                                             deseq2Results_sense__0tex_ont))
deseq2Results_sense__ont$col <- with(deseq2Results_sense__ont, 
                              ifelse(padj_value == "different" & log2FoldChange > 0, "up-regulated",
                                     ifelse(padj_value == "different" & log2FoldChange < 0, "down-regulated", "not different")))

deseq_sense_table <- data.frame(table(deseq2Results_sense__ont$col, deseq2Results_sense__ont$compare))
deseq_sense_table$Perc <- deseq_sense_table$Freq / 6545
deseq_sense_table$Type <- "Sense"


################################################################################################
# Convert count data to a matrix of appropriate form that DEseq2 can read

row_sub = apply(featureCounts_antisenseNN2, 1, function(row) all(row > -2 ))
featureCounts_antisenseNN20 <- featureCounts_antisenseNN2[row_sub,]
geneID_antisense <- rownames(featureCounts_antisenseNN20)
sampleIndex_antisense <- colnames(featureCounts_antisenseNN20)
rawCounts_antisense <- as.matrix(featureCounts_antisenseNN20)
rownames(rawCounts_antisense) <- geneID_antisense


sample_antisense <- colnames(featureCounts_antisenseNN2)
time_antisense <- c("4h", "4h", "4h", "8h", "8h", "8h", "8h", "8h", "4h", "4h", "8h", "8h", "4h", "4h")
platform_antisense <- c("ONT", "ONT", "ONT", "ONT", "ONT", "ONT", "TEX", "TEX", "TEX", "TEX", "0TEX", "0TEX", "0TEX", "0TEX")
colData_ont_as <- data.frame(cbind(sample_antisense, time_antisense, platform_antisense))
rownames(colData_ont_as) <- colData_ont_as$sample
colData_ont_as$platform_antisense <- factor(colData_ont_s$platform, levels=c("ONT", "TEX", "0TEX"))
colData_ont_as$time_antisense <- factor(colData_ont_as$time, levels=c("4h", "8h"))

# Create the DEseq2 DataSet object
de_ont_as_deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts_antisense, colData=colData_ont_as, 
                                               design= ~ time_antisense+platform_antisense)


rld_antisense <- rlog(de_ont_as_deseq2Data, blind=T)
rld_sense <- rlog(de_ont_s_deseq2Data, blind=T)



######################################################################
# PCA plots
pca_antisense <- plotPCA(rld_antisense, intgroup="platform_antisense") + theme_classic() + ylim(-50,50) + xlim(-50,50) +
  theme(axis.title = element_text(size=9), legend.title = element_blank())
pca_sense <- plotPCA(rld_sense, intgroup="platform_sense") + theme_classic() + ylim(-50,50) + xlim(-50,50) +
  theme(axis.title = element_text(size=9), legend.title = element_blank())

pca_plots <- ggarrange(pca_sense, pca_antisense, labels=c("A", "B"), common.legend = TRUE,
                       font.label = list(size = 12, color = "black"))
ggsave(pca_plots, filename="figures_all/pca_plots.tif", device = "tiff", dpi=300, width=6, height=6, bg="white")

######################################################################

# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
de_ont_as_deseq2Data <- DESeq(de_ont_as_deseq2Data)

deseq2Results_antisense__tex_0tex <- lfcShrink(de_ont_as_deseq2Data, contrast = c("platform_antisense", "TEX", "0TEX"), 
                                               type="ashr", format = "DataFrame")
deseq2Results_antisense__tex_0tex$padj_value <- ifelse(deseq2Results_antisense__tex_0tex$padj < 0.01, "different", "not-different")
deseq2Results_antisense__tex_0tex <- na.omit(deseq2Results_antisense__tex_0tex)
deseq2Results_antisense__tex_0tex$compare <- "TEX (reference) - 0TEX"
nrow_tex_0tex <- nrow(deseq2Results_antisense__tex_0tex)

deseq2Results_antisense__tex_ont <- lfcShrink(de_ont_as_deseq2Data, contrast = c("platform_antisense", "ONT", "TEX"), 
                                              type="ashr", format = "DataFrame")
deseq2Results_antisense__tex_ont$padj_value <- ifelse(deseq2Results_antisense__tex_ont$padj < 0.01, "different", "not-different")
deseq2Results_antisense__tex_ont <- na.omit(deseq2Results_antisense__tex_ont)
deseq2Results_antisense__tex_ont$compare <- "ONT (reference) - TEX"
nrow_tex_ont <- nrow(deseq2Results_antisense__tex_ont)


deseq2Results_antisense__0tex_ont <- lfcShrink(de_ont_as_deseq2Data, contrast = c("platform_antisense", "ONT", "0TEX"), 
                                               type="ashr", format = "DataFrame")

deseq2Results_antisense__0tex_ont$padj_value <- ifelse(deseq2Results_antisense__0tex_ont$padj < 0.01, "different", "not-different")
deseq2Results_antisense__0tex_ont$compare <- "ONT (reference) - 0TEX"
nrow_0tex_ont <- nrow(deseq2Results_antisense__0tex_ont)

############################################################################################################################################
# merge it
deseq2Results_antisense__ont <- data.frame(rbind(deseq2Results_antisense__tex_0tex, 
                                                 deseq2Results_antisense__tex_ont,
                                                 deseq2Results_antisense__0tex_ont))
deseq2Results_antisense__ont$col <- with(deseq2Results_antisense__ont, 
                              ifelse(padj_value == "different" & log2FoldChange > 0, "up-regulated",
                                     ifelse(padj_value == "different" & log2FoldChange < 0, "down-regulated", "not different")))

deseq_antisense_table <- data.frame(table(deseq2Results_antisense__ont$col, deseq2Results_antisense__ont$compare))
deseq_antisense_table$Total_sum <- c(nrow_0tex_ont, nrow_0tex_ont, nrow_0tex_ont,
                                     nrow_tex_ont, nrow_tex_ont, nrow_tex_ont,
                                     nrow_tex_0tex, nrow_tex_0tex, nrow_tex_0tex)
deseq_antisense_table$Perc <- deseq_antisense_table$Freq / deseq_antisense_table$Total_sum
deseq_antisense_table$Type <- "Antisense"


sense_log2FoldChange <- 
  ggplot(deseq2Results_sense__ont) +
  geom_point(aes(x=log2(baseMean), y=log2FoldChange, colour=col), size=0.5, alpha=1) +
  facet_wrap(~compare, nrow=3) + 
  xlim(0,15) + 
  theme_bw() +
  scale_y_continuous(breaks = c(-3, 0, 3),
                     labels = c(-3, 0, 3),
                     limits = c(-6,6)) +
  geom_hline(yintercept = 0, colour="black", size=0.4, linetype="longdash") + 
  labs(x=expression("Log"[2]*" base mean"), y=expression("Log"[2]*" fold change")) + 
  scale_colour_manual(values = c("steelblue4","snow4","lightsalmon")) +
  geom_hline(yintercept = 1, colour="black", size=0.2, linetype="dotted") + 
  geom_hline(yintercept = -1, colour="black", size=0.2, linetype="dotted") + 
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(), 
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size= 7),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

antisense_log2FoldChange <- 
  ggplot(deseq2Results_antisense__ont) +
  geom_point(aes(x=log2(baseMean), y=log2FoldChange, colour=col), size=0.5, alpha=1) +
  facet_wrap(~compare, nrow=3) + 
  xlim(0,15) + 
  theme_bw() +
  scale_y_continuous(breaks = c(-3, 0, 3),
                     labels = c(-3, 0, 3),
                     limits = c(-6,6)) +
  geom_hline(yintercept = 0, colour="black", size=0.4, linetype="longdash") + 
  geom_hline(yintercept = 1, colour="black", size=0.2, linetype="dotted") + 
  geom_hline(yintercept = -1, colour="black", size=0.2, linetype="dotted") + 
  labs(x=expression("Log"[2]*" base mean"), y=" ") + 
  scale_colour_manual(values = c("steelblue4","snow4","lightsalmon")) +
  guides(colour = guide_legend(override.aes = list(size =2))) +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(), 
        strip.background =element_rect(fill="beige"),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))


sense_barplot <- 
  ggplot(deseq_sense_table) +
  geom_col(aes(x=Var2, y=Perc, fill=reorder(Var1, Perc)), position = "stack", width=0.8) +
  scale_fill_manual(values = c("steelblue4", "lightsalmon", "snow4")) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1, size=7), 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  xlab(" ") + 
  ylab("\nRelative abundance\n") + facet_wrap(~Type)


antisense_barplot <- 
  ggplot(deseq_antisense_table) +
  geom_col(aes(x=Var2, y=Perc, fill=reorder(Var1, Perc)), position = "stack", width=0.8) +
  scale_fill_manual(values = c("steelblue4", "lightsalmon", "snow4")) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1, size=7),
        strip.background = element_rect(fill="beige"), 
        strip.text = element_text(size = 7),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text.y = element_text(size = 8)) + 
  facet_wrap(~Type) +
  xlab(" ") + 
  ylab("\n\n ") 

################################################################################################
# Differential expression plots
merge_barplots <-
  ggarrange(sense_barplot, antisense_barplot, labels=c("C", "D"), 
          nrow=1, ncol=2, widths = c(1,1), common.legend = TRUE, legend = "none",
          font.label = list(size = 9))

merge_de_plots <-
  ggarrange(sense_log2FoldChange, antisense_log2FoldChange, labels=c("A", "B"), 
          nrow=1, ncol=2, widths = c(1,1), common.legend = TRUE, legend = "bottom",
          font.label = list(size = 9))

merged_de <- ggarrange(merge_de_plots, merge_barplots, common.legend = TRUE, widths = c(1,0.6))

tiff(filename="figures_all/differential_expression.tif", res=300, width=20, height=12, units = "cm")
merged_de
dev.off()
################################################################################################


ftsH_genes_NN2 <- import_gtf_NN2[grepl("ftsA|ftsQ|murC|murG|ftsW|murD|ftsZ", import_gtf_NN2$gene_name4),]
ftsH_genes_NN2 <- ftsH_genes_NN2[order(ftsH_genes_NN2$transcript_id),]

ftsH_genes_SG17M <- import_gtf_SG17M[grepl("ftsA|ftsQ|murC|murG|ftsW|murD|ftsZ", import_gtf_SG17M$gene_name4),]
ftsH_genes_SG17M <- ftsH_genes_SG17M[order(ftsH_genes_SG17M$transcript_id),]

featureCounts_senseNN2_ftsH <- featureCounts_senseNN2[rownames(featureCounts_senseNN2)%in%ftsH_genes_NN2$transcript_id,]
featureCounts_senseNN2_ftsH <- featureCounts_senseNN2_ftsH[order(rownames(featureCounts_senseNN2_ftsH)),]

featureCounts_senseNN2_ftsH$ftsH_name <- ftsH_genes_NN2$gene_name4
featureCounts_senseNN2_ftsH
featureCounts_senseNN2_ftsH$ftsH_name <- str_replace(featureCounts_senseNN2_ftsH$ftsH_name, "_1", "")
featureCounts_senseNN2_ftsH$ftsH_name <- str_replace(featureCounts_senseNN2_ftsH$ftsH_name, "_2", "")
featureCounts_senseNN2_ftsH$ftsH_name <- str_replace(featureCounts_senseNN2_ftsH$ftsH_name, "_3", "")
featureCounts_senseNN2_ftsH$ftsH_name <- str_replace(featureCounts_senseNN2_ftsH$ftsH_name, "_4", "")
featureCounts_senseNN2_ftsH <- plyr::ddply(featureCounts_senseNN2_ftsH, "ftsH_name", numcolwise(sum))
rownames(featureCounts_senseNN2_ftsH) <- featureCounts_senseNN2_ftsH$ftsH_name
featureCounts_senseNN2_ftsH$ftsH_name <- NULL

featureCounts_senseSG17M_ftsH <- ont_featureCounts_senseSG17M[rownames(ont_featureCounts_senseSG17M)%in%ftsH_genes_SG17M$transcript_id,]
featureCounts_senseSG17M_ftsH <- featureCounts_senseSG17M_ftsH[order(rownames(featureCounts_senseSG17M_ftsH)),]
featureCounts_senseSG17M_ftsH$ftsH_name <- ftsH_genes_SG17M$gene_name4
featureCounts_senseSG17M_ftsH
featureCounts_senseSG17M_ftsH$ftsH_name <- str_replace(featureCounts_senseSG17M_ftsH$ftsH_name, "_1", "")
featureCounts_senseSG17M_ftsH$ftsH_name <- str_replace(featureCounts_senseSG17M_ftsH$ftsH_name, "_2", "")
featureCounts_senseSG17M_ftsH$ftsH_name <- str_replace(featureCounts_senseSG17M_ftsH$ftsH_name, "_3", "")
featureCounts_senseSG17M_ftsH$ftsH_name <- str_replace(featureCounts_senseSG17M_ftsH$ftsH_name, "_4", "")
featureCounts_senseSG17M_ftsH <- plyr::ddply(featureCounts_senseSG17M_ftsH, "ftsH_name", numcolwise(sum))
rownames(featureCounts_senseSG17M_ftsH) <- featureCounts_senseSG17M_ftsH$ftsH_name
featureCounts_senseSG17M_ftsH$ftsH_name <- NULL
#featureCounts_senseSG17M_ftsH$SG17M_4h_BR1 <- featureCounts_senseSG17M_ftsH$SG17M_4h_BR1 * 100
#featureCounts_senseSG17M_ftsH$SG17M_4h_BR2 <- featureCounts_senseSG17M_ftsH$SG17M_4h_BR1 * 100
#featureCounts_senseSG17M_ftsH$SG17M_8h_BR1 <- featureCounts_senseSG17M_ftsH$SG17M_8h_BR1 * 100
#featureCounts_senseSG17M_ftsH$SG17M_8h_BR2 <- featureCounts_senseSG17M_ftsH$SG17M_8h_BR1 * 100
#featureCounts_senseSG17M_ftsH$SG17M_8h_BR3 <- featureCounts_senseSG17M_ftsH$SG17M_8h_BR3 * 100

featureCounts_ftsH <- cbind(featureCounts_senseNN2_ftsH, featureCounts_senseSG17M_ftsH)

# Convert count data to a matrix of appropriate form that DEseq2 can read
row_sub = apply(featureCounts_ftsH, 1, function(row) all(row > -2 ))
featureCounts_ftsH0 <- featureCounts_ftsH[row_sub,]
geneID_sense_ftsH <- rownames(featureCounts_ftsH0)
sampleIndex_sense_ftsH <- colnames(featureCounts_ftsH0)
rawCounts_sense_ftsH <- as.matrix(featureCounts_ftsH0)
rownames(rawCounts_sense_ftsH) <- geneID_sense_ftsH

# make coldata
sample_sense_ftsH <- colnames(featureCounts_ftsH)
time_sense_ftsH <- c("4h", "4h", "4h", "8h", "8h", "8h", 
                     "8h", "8h", "4h", "4h", 
                     "8h", "8h", "4h", "4h", 
                     "4h", "4h", "4h", "8h", "8h", "8h")
platform_sense_ftsH <- c("ONT", "ONT", "ONT", "ONT", "ONT", "ONT", 
                         "TEX", "TEX", "TEX", "TEX", 
                         "0TEX", "0TEX", "0TEX", "0TEX", 
                         "SG17M", "SG17M", "SG17M", "SG17M", "SG17M", "SG17M")
colData_ont_s_ftsH <- data.frame(cbind(sample_sense_ftsH, time_sense_ftsH, platform_sense_ftsH))
rownames(colData_ont_s_ftsH) <- colData_ont_s_ftsH$sample
colData_ont_s_ftsH$platform_sense_ftsH <- factor(colData_ont_s_ftsH$platform, levels=c("ONT", "TEX", "0TEX", "SG17M"))
colData_ont_s_ftsH$time_sense_ftsH <- factor(colData_ont_s_ftsH$time, levels=c("4h", "8h"))

# Create the DEseq2 DataSet object
de_ont_s_deseq2Data_ftsH <- DESeqDataSetFromMatrix(countData=rawCounts_sense_ftsH, colData=colData_ont_s_ftsH, design= ~ time_sense_ftsH+platform_sense_ftsH)

# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
de_ont_s_deseq2Data_ftsH <- DESeq(de_ont_s_deseq2Data_ftsH)
deseq2Results_sense__tex_0tex_ftSH <- lfcShrink(de_ont_s_deseq2Data_ftsH, contrast = c("platform_sense_ftsH", "TEX", "0TEX"),  type="ashr", format = "DataFrame")
deseq2Results_sense__tex_0tex_ftSH$padj_value <- ifelse(deseq2Results_sense__tex_0tex_ftSH$padj < 0.01, "different", "not-different")
deseq2Results_sense__tex_0tex_ftSH <- na.omit(deseq2Results_sense__tex_0tex_ftSH)
deseq2Results_sense__tex_0tex_ftSH$compare <- "TEX (reference) - 0TEX"

deseq2Results_sense__tex_ont_ftsH <- lfcShrink(de_ont_s_deseq2Data_ftsH, contrast = c("platform_sense_ftsH", "ONT", "TEX"), type="ashr", format = "DataFrame")
deseq2Results_sense__tex_ont_ftsH$padj_value <- ifelse(deseq2Results_sense__tex_ont_ftsH$padj < 0.01, "different", "not-different")
deseq2Results_sense__tex_ont_ftsH <- na.omit(deseq2Results_sense__tex_ont_ftsH)
deseq2Results_sense__tex_ont_ftsH$compare <- "ONT (reference) - TEX"

deseq2Results_sense__0tex_ont_ftsH <- lfcShrink(de_ont_s_deseq2Data_ftsH, contrast = c("platform_sense_ftsH", "ONT", "0TEX"), type="ashr", format = "DataFrame")
deseq2Results_sense__0tex_ont_ftsH$padj_value <- ifelse(deseq2Results_sense__0tex_ont_ftsH$padj < 0.01, "different", "not-different")
deseq2Results_sense__0tex_ont_ftsH <- na.omit(deseq2Results_sense__0tex_ont_ftsH)
deseq2Results_sense__0tex_ont_ftsH$compare <- "ONT (reference) - 0TEX"


deseq2Results_sense__ont_sg17m_ftSH <- lfcShrink(de_ont_s_deseq2Data_ftsH, contrast = c("platform_sense_ftsH", "ONT", "SG17M"), type="ashr", format = "DataFrame")
deseq2Results_sense__ont_sg17m_ftSH$padj_value <- ifelse(deseq2Results_sense__ont_sg17m_ftSH$padj < 0.01, "different", "not-different")
deseq2Results_sense__ont_sg17m_ftSH <- na.omit(deseq2Results_sense__ont_sg17m_ftSH)
deseq2Results_sense__ont_sg17m_ftSH$compare <- "ONT (reference) - SG17M"

deseq2Results_sense__tex_SG17M_ftsH <- lfcShrink(de_ont_s_deseq2Data_ftsH, contrast = c("platform_sense_ftsH", "SG17M", "TEX"), type="ashr", format = "DataFrame")
deseq2Results_sense__tex_SG17M_ftsH$padj_value <- ifelse(deseq2Results_sense__tex_SG17M_ftsH$padj < 0.01, "different", "not-different")
deseq2Results_sense__tex_SG17M_ftsH <- na.omit(deseq2Results_sense__tex_SG17M_ftsH)
deseq2Results_sense__tex_SG17M_ftsH$compare <- "SG17M (reference) - TEX"


deseq2Results_sense__0tex_SG17M_ftsH <- lfcShrink(de_ont_s_deseq2Data_ftsH, contrast = c("platform_sense_ftsH", "SG17M", "0TEX"), type="ashr", format = "DataFrame")
deseq2Results_sense__0tex_SG17M_ftsH$padj_value <- ifelse(deseq2Results_sense__0tex_SG17M_ftsH$padj < 0.01, "different", "not-different")
deseq2Results_sense__0tex_SG17M_ftsH <- na.omit(deseq2Results_sense__0tex_SG17M_ftsH)
deseq2Results_sense__0tex_SG17M_ftsH$compare <- "SG17M (reference) - 0TEX"


deseq2Results_sense__ont_ftsH <- data.frame(rbind(deseq2Results_sense__0tex_ont_ftsH, deseq2Results_sense__tex_ont_ftsH, deseq2Results_sense__tex_0tex_ftSH,
                                             deseq2Results_sense__ont_sg17m_ftSH, deseq2Results_sense__tex_SG17M_ftsH, deseq2Results_sense__0tex_SG17M_ftsH))

# make a boxplot of the Cook’s distances to see if one sample is 
# consistently higher than others (here this is not the case):
boxplot(log10(assays(de_ont_s_deseq2Data)[["cooks"]]), range=0, las=2)

# merge it
deseq2Results_sense__ont_ftsH$col <- with(deseq2Results_sense__ont_ftsH, 
                              ifelse(padj_value == "different" & log2FoldChange > 0, "up-regulated",
                                     ifelse(padj_value == "different" & log2FoldChange < 0, "down-regulated", 
                                            "not different")))


nn2_tRNA_no0 <- nn2_tRNA[grep("tRNA-", nn2_tRNA$gene_name2),]
sg17m_tRNA_no0 <- sg17m_tRNA[grep("tRNA-", sg17m_tRNA$gene_name2),]

group <- c("SG17M_ONT_4h", "SG17M_4h_ONT", "SG17M_4h_ONT", "SG17M_8h_ONT", "SG17M_8h_ONT", "SG17M_8h_ONT",
           "NN2_ONT_4h", "NN2_ONT_4h", "NN2_ONT_4h", "NN2_ONT_8h", "NN2_ONT_8h", "NN2_ONT_8h",
           "NN2_TEX_8h", "NN2_TEX_8h", "NN2_TEX_4h", "NN2_TEX_4h", "NN2_0TEX_8h", "NN2_0TEX_8h", "NN2_0TEX_4h", "NN2_0TEX_4h")

featureCounts_senseSG17M <- ont_featureCounts_senseSG17M
featureCounts_antisenseSG17M <- ont_featureCounts_antisenseSG17M

featureCounts_senseNN2_tRNA <- featureCounts_senseNN2[rownames(featureCounts_senseNN2) %in% nn2_tRNA_no0$gene_name2,]
featureCounts_antisenseNN2_tRNA <- featureCounts_antisenseNN2[rownames(featureCounts_antisenseNN2) %in% nn2_tRNA_no0$gene_name2,]

featureCounts_senseSG17M_tRNA <- featureCounts_senseSG17M[rownames(featureCounts_senseSG17M) %in% sg17m_tRNA_no0$gene_name2,]
featureCounts_antisenseSG17M_tRNA <- featureCounts_antisenseSG17M[rownames(featureCounts_antisenseSG17M) %in% sg17m_tRNA_no0$gene_name2,]

# merged sense and antisense
rownames(featureCounts_senseNN2_tRNA) <- paste("(+) ", rownames(featureCounts_senseNN2_tRNA), sep="")
rownames(featureCounts_antisenseNN2_tRNA) <- paste("(-) ", rownames(featureCounts_antisenseNN2_tRNA), sep="")

rownames(featureCounts_senseSG17M_tRNA) <- paste("(+) ", rownames(featureCounts_senseSG17M_tRNA), sep="")
rownames(featureCounts_antisenseSG17M_tRNA) <- paste("(-) ", rownames(featureCounts_antisenseSG17M_tRNA), sep="")

featureCounts_tRNA <- data.frame(rbind(featureCounts_senseNN2_tRNA, featureCounts_antisenseNN2_tRNA))
featureCounts_tRNA <- featureCounts_tRNA[order(rownames(featureCounts_tRNA)),]
featureCounts_tRNA_SG17M <- data.frame(rbind(featureCounts_senseSG17M_tRNA, featureCounts_antisenseSG17M_tRNA))

featureCounts_tRNA_SG17M <- featureCounts_tRNA_SG17M[order(rownames(featureCounts_tRNA_SG17M)),]
merged_featureCounts_tRNA <- data.frame(cbind(featureCounts_tRNA_SG17M, featureCounts_tRNA))


y_both <- DGEList(counts=merged_featureCounts_tRNA, group=group)
# We filter out lowly expressed genes using the following commands:
keep_both <- filterByExpr(y_both, min.count=10)
y_both <- y_both[keep_both, keep.lib.sizes=TRUE]
all_counts <- y_both$counts
all_counts <- data.frame(all_counts)

write.csv(all_counts, file="all_counts.csv", sep=";")
y_both <- calcNormFactors(y_both, method="RLE")
logcpm_both <- cpm(y_both, log=TRUE)
logcpm_both_df <- data.frame(logcpm_both)

strReverse <- function(x){sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")}
NN2_codon_usage <- read_delim("NN2_codon_usage.csv", ";", escape_double = FALSE, trim_ws = TRUE)
NN2_codon_usage <- NN2_codon_usage[,1:5]
NN2_codon_usage$Anticodon1 <- strReverse(NN2_codon_usage$Codon)
nn2_comp <- c()
for (items in NN2_codon_usage$Anticodon1) { a = c2s(comp(s2c(items), forceToLower = FALSE))
   nn2_comp <- append(nn2_comp, a)}
NN2_codon_usage$Anticodon2 <- nn2_comp
NN2_codon_usage$Anticodon3 <- paste("tRNA-",NN2_codon_usage$AmAcid,NN2_codon_usage$Anticodon2)
NN2_codon_usage$Anticodon3 <- str_replace_all(NN2_codon_usage$Anticodon3, " ", "")
NN2_codon_usage <- NN2_codon_usage[1:61,]
NN2_codon_usage1 <- NN2_codon_usage[order(NN2_codon_usage$Anticodon3),]
NN2_codon_usage1 <- data.frame(NN2_codon_usage1)

NN2_codon_usage1_tmRNA <- c("Undet", "NNN", NA, NA, NA, "NNN", "NNN", "tRNA-UndetNNN")
NN2_codon_usage1_LysTTT_pseudo <- c("Lys-pseudo", "TTT", NA, NA, NA, "AAA", "TTT", "tRNA-LysTTT-pseudo")
NN2_codon_usage1_SerCGA_pseudo <- c("Ser-pseudo", "CGA", NA, NA, NA, "GCT", "CGA", "tRNA-SerCGA-pseudo")

NN2_codon_usage2 <- data.frame(rbind(NN2_codon_usage1,
                                     NN2_codon_usage1_tmRNA,
                                     NN2_codon_usage1_LysTTT_pseudo,
                                     NN2_codon_usage1_SerCGA_pseudo))
NN2_codon_usage3 <- NN2_codon_usage2[order(NN2_codon_usage2$Anticodon3),]


missing_genes_tRNA <- c("tRNA-AlaAGC", 
                        "tRNA-AlaCGC", 
                        "tRNA-ArgGCG", 
                        "tRNA-ArgTCG", 
                        "tRNA-AsnATT",
                        "tRNA-AspATC", 
                        "tRNA-CysACA", 
                        "tRNA-GlnCTG", 
                        "tRNA-GluCTC",
                        "tRNA-GlyACC",
                        "tRNA-HisATG",
                        "tRNA-IleAAT",
                        "tRNA-IleTAT", 
                        "tRNA-LeuAAG", 
                        "tRNA-LeuTAA",
                        "tRNA-LysCTT", 
                        "tRNA-PheAAA",  
                        "tRNA-ProAGG",
                        "tRNA-SerACT", 
                        "tRNA-SerAGA", 
                        "tRNA-ThrAGT", 
                        "tRNA-TyrATA",
                        "tRNA-ValAAC", 
                        "tRNA-ValCAC")


y_both_unfiltered_df <- logcpm_both_df
y_both_unfiltered_df$Anticodon3 <- rownames(y_both_unfiltered_df)
y_both_unfiltered_df$Anticodon3 <- str_replace_all(y_both_unfiltered_df$Anticodon3, " ", "")

rownames(y_both_unfiltered_df) <- NULL
y_both_unfiltered_df$Anticodon3 <- str_replace_all(y_both_unfiltered_df$Anticodon3, "tRNA-fMetCAT", "tRNA-MetCAT")
y_both_unfiltered_df$Anticodon3 <- str_replace_all(y_both_unfiltered_df$Anticodon3, "tRNA-SeCTCA", "tRNA-SerTGA")
y_both_unfiltered_df$Anticodon3 <- str_replace_all(y_both_unfiltered_df$Anticodon3, "tRNA-Ile2CAT", "tRNA-IleGAT")

y_both_unfiltered_df_plus <- subset(y_both_unfiltered_df, grepl("\\(\\+\\)", y_both_unfiltered_df$Anticodon3 ))
y_both_unfiltered_df_plus$Anticodon3 <- str_replace_all(y_both_unfiltered_df_plus$Anticodon3, "\\(\\+\\)", "")
y_both_unfiltered_df_plus <- plyr::ddply(y_both_unfiltered_df_plus, "Anticodon3", numcolwise(sum))

y_both_unfiltered_df_minus <- subset(y_both_unfiltered_df, grepl("\\(\\-\\)", y_both_unfiltered_df$Anticodon3 ))
y_both_unfiltered_df_minus$Anticodon3 <- str_replace_all(y_both_unfiltered_df_minus$Anticodon3, "\\(\\-\\)", "")
y_both_unfiltered_df_minus <- plyr::ddply(y_both_unfiltered_df_minus, "Anticodon3", numcolwise(sum))

# how many tRNAs are detected in samples that is not in genome sequence? n = 0
missing_ab1_minus <- setdiff(y_both_unfiltered_df_minus$Anticodon3, NN2_codon_usage3$Anticodon3)

# how many tRNAs are detected in genome sequence but not in samples? n = 24
missing_df1_minus <- setdiff(NN2_codon_usage3$Anticodon3, y_both_unfiltered_df_minus$Anticodon3)
lengthn <- length(missing_df1_minus)
# convert missing to data frame with 0's set
missing_df2_minus <- rep(0,lengthn) #22
missing_df3_minus <- data.frame(cbind(missing_df1_minus, missing_df2_minus, missing_df2_minus, missing_df2_minus, 
                                      missing_df2_minus, missing_df2_minus, missing_df2_minus, 
                                      missing_df2_minus, missing_df2_minus, missing_df2_minus, 
                                      missing_df2_minus, missing_df2_minus, missing_df2_minus,
                                      missing_df2_minus, missing_df2_minus, missing_df2_minus, 
                                      missing_df2_minus, missing_df2_minus, missing_df2_minus, 
                                      missing_df2_minus, missing_df2_minus))

colnames(missing_df3_minus) <- colnames(y_both_unfiltered_df_minus)

#y_both_unfiltered_df_minus <- plyr::ddply(y_both_unfiltered_df_minus, "Anticodon3", numcolwise(sum))
rownames(y_both_unfiltered_df_minus) <- y_both_unfiltered_df_minus$Anticodon3
y_both_unfiltered_df_minus$Anticodon3 <- NULL

#missing_df3_minus <- plyr::ddply(missing_df3_minus, "Anticodon3", numcolwise(sum))
rownames(missing_df3_minus) <- missing_df3_minus$Anticodon3
missing_df3_minus$Anticodon3 <- NULL
citation("DESeq2")
y_both_minus <- data.frame(rbind(y_both_unfiltered_df_minus, missing_df3_minus))
y_both_minus_1 <- y_both_minus[order(rownames(y_both_minus)),]
y_both_minus_2 <- data.frame(y_both_minus_1)
y_both_minus_2$ref_NN2 <- (NN2_codon_usage3$div1000)

y_both_minus_2$NN2_4h_BR1 <- as.numeric(as.character(y_both_minus_2$NN2_4h_BR1))
y_both_minus_2$NN2_4h_BR2 <- as.numeric(as.character(y_both_minus_2$NN2_4h_BR2))
y_both_minus_2$NN2_4h_BR3 <- as.numeric(as.character(y_both_minus_2$NN2_4h_BR3))
y_both_minus_2$NN2_8h_BR1 <- as.numeric(as.character(y_both_minus_2$NN2_8h_BR1))
y_both_minus_2$NN2_8h_BR2 <- as.numeric(as.character(y_both_minus_2$NN2_8h_BR1))
y_both_minus_2$NN2_8h_BR3 <- as.numeric(as.character(y_both_minus_2$NN2_8h_BR1))

y_both_minus_2$SG17M_4h_BR1 <- as.numeric(as.character(y_both_minus_2$SG17M_4h_BR1))
y_both_minus_2$SG17M_4h_BR2 <- as.numeric(as.character(y_both_minus_2$SG17M_4h_BR2))
y_both_minus_2$SG17M_4h_BR3 <- as.numeric(as.character(y_both_minus_2$SG17M_4h_BR3))
y_both_minus_2$SG17M_8h_BR1 <- as.numeric(as.character(y_both_minus_2$SG17M_8h_BR1))
y_both_minus_2$SG17M_8h_BR2 <- as.numeric(as.character(y_both_minus_2$SG17M_8h_BR1))
y_both_minus_2$SG17M_8h_BR3 <- as.numeric(as.character(y_both_minus_2$SG17M_8h_BR1))

y_both_minus_2$TEX_NN2_8h_1 <- as.numeric(as.character(y_both_minus_2$TEX_NN2_8h_1))
y_both_minus_2$TEX_NN2_8h_2 <- as.numeric(as.character(y_both_minus_2$TEX_NN2_8h_2))
y_both_minus_2$TEX_NN2_4h_1 <- as.numeric(as.character(y_both_minus_2$TEX_NN2_4h_1))
y_both_minus_2$TEX_NN2_4h_2 <- as.numeric(as.character(y_both_minus_2$TEX_NN2_4h_2))
y_both_minus_2$X0TEX_NN2_8h_1 <- as.numeric(as.character(y_both_minus_2$X0TEX_NN2_8h_1))
y_both_minus_2$X0TEX_NN2_8h_2 <- as.numeric(as.character(y_both_minus_2$X0TEX_NN2_8h_2))
y_both_minus_2$X0TEX_NN2_4h_1 <- as.numeric(as.character(y_both_minus_2$X0TEX_NN2_4h_1))
y_both_minus_2$X0TEX_NN2_4h_2 <- as.numeric(as.character(y_both_minus_2$X0TEX_NN2_4h_2))
y_both_minus_2$ref_NN2 <- as.numeric(as.character(y_both_minus_2$ref_NN2))

# how many tRNAs are detected in samples that is not in genome sequence? n = 0
missing_ab1_plus <- setdiff(y_both_unfiltered_df_plus$Anticodon3, NN2_codon_usage3$Anticodon3)

# how many tRNAs are detected in genome sequence but not in samples? n = 24
missing_df1_plus <- setdiff(NN2_codon_usage3$Anticodon3, y_both_unfiltered_df_plus$Anticodon3)
nlength_plus <- length(missing_df1_plus)

# convert missing to data frame with 0's set
missing_df2_plus <- rep(0,nlength_plus)
missing_df3_plus <- data.frame(cbind(missing_df1_plus, 
                                     missing_df2_plus, missing_df2_plus, missing_df2_plus, 
                                      missing_df2_plus, missing_df2_plus, missing_df2_plus, 
                                      missing_df2_plus, missing_df2_plus, missing_df2_plus, 
                                      missing_df2_plus, missing_df2_plus, missing_df2_plus,
                                      missing_df2_plus, missing_df2_plus, missing_df2_plus, 
                                      missing_df2_plus, missing_df2_plus, missing_df2_plus, 
                                      missing_df2_plus, missing_df2_plus)) 
rownames(missing_df3_plus) <- missing_df1_plus
colnames(missing_df3_plus) <- colnames(y_both_unfiltered_df_plus)
missing_df3_plus$Anticodon3 <- NULL

y_both_unfiltered_df_plus <- plyr::ddply(y_both_unfiltered_df_plus, "Anticodon3", numcolwise(sum))
rownames(y_both_unfiltered_df_plus) <- y_both_unfiltered_df_plus$Anticodon3
y_both_unfiltered_df_plus$Anticodon3 <- NULL

y_both_plus <- data.frame(rbind(y_both_unfiltered_df_plus, missing_df3_plus))
y_both_plus_1 <- y_both_plus[order(rownames(y_both_plus)),]
y_both_plus_2 <- data.frame(y_both_plus_1)

y_both_plus_2$ref_NN2 <- (NN2_codon_usage3$div1000) 

y_both_plus_2$NN2_4h_BR1 <- as.numeric(as.character(y_both_plus_2$NN2_4h_BR1))
y_both_plus_2$NN2_4h_BR2 <- as.numeric(as.character(y_both_plus_2$NN2_4h_BR2))
y_both_plus_2$NN2_4h_BR3 <- as.numeric(as.character(y_both_plus_2$NN2_4h_BR3))
y_both_plus_2$NN2_8h_BR1 <- as.numeric(as.character(y_both_plus_2$NN2_8h_BR1))
y_both_plus_2$NN2_8h_BR2 <- as.numeric(as.character(y_both_plus_2$NN2_8h_BR1))
y_both_plus_2$NN2_8h_BR3 <- as.numeric(as.character(y_both_plus_2$NN2_8h_BR1))

y_both_plus_2$SG17M_4h_BR1 <- as.numeric(as.character(y_both_plus_2$SG17M_4h_BR1))
y_both_plus_2$SG17M_4h_BR2 <- as.numeric(as.character(y_both_plus_2$SG17M_4h_BR2))
y_both_plus_2$SG17M_4h_BR3 <- as.numeric(as.character(y_both_plus_2$SG17M_4h_BR3))
y_both_plus_2$SG17M_8h_BR1 <- as.numeric(as.character(y_both_plus_2$SG17M_8h_BR1))
y_both_plus_2$SG17M_8h_BR2 <- as.numeric(as.character(y_both_plus_2$SG17M_8h_BR1))
y_both_plus_2$SG17M_8h_BR3 <- as.numeric(as.character(y_both_plus_2$SG17M_8h_BR1))

y_both_plus_2$TEX_NN2_8h_1 <- as.numeric(as.character(y_both_plus_2$TEX_NN2_8h_1))
y_both_plus_2$TEX_NN2_8h_2 <- as.numeric(as.character(y_both_plus_2$TEX_NN2_8h_2))
y_both_plus_2$TEX_NN2_4h_1 <- as.numeric(as.character(y_both_plus_2$TEX_NN2_4h_1))
y_both_plus_2$TEX_NN2_4h_2 <- as.numeric(as.character(y_both_plus_2$TEX_NN2_4h_2))
y_both_plus_2$X0TEX_NN2_8h_1 <- as.numeric(as.character(y_both_plus_2$X0TEX_NN2_8h_1))
y_both_plus_2$X0TEX_NN2_8h_2 <- as.numeric(as.character(y_both_plus_2$X0TEX_NN2_8h_2))
y_both_plus_2$X0TEX_NN2_4h_1 <- as.numeric(as.character(y_both_plus_2$X0TEX_NN2_4h_1))
y_both_plus_2$X0TEX_NN2_4h_2 <- as.numeric(as.character(y_both_plus_2$X0TEX_NN2_4h_2))
y_both_plus_2$ref_NN2 <- as.numeric(as.character(y_both_plus_2$ref_NN2))


ties_method = "min"
na_method = "keep"
y_both_plus_2  <- y_both_plus_2[!rownames(y_both_plus_2)%in%missing_genes_tRNA,]
y_both_plus_2$NN2_4h_BR1_rank <- rank(y_both_plus_2$NN2_4h_BR1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$NN2_4h_BR2_rank <- rank(y_both_plus_2$NN2_4h_BR2, ties.method = ties_method, na.last = na_method)
y_both_plus_2$NN2_4h_BR3_rank <- rank(y_both_plus_2$NN2_4h_BR3, ties.method = ties_method, na.last = na_method)
y_both_plus_2$NN2_8h_BR1_rank <- rank(y_both_plus_2$NN2_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$NN2_8h_BR2_rank <- rank(y_both_plus_2$NN2_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$NN2_8h_BR3_rank <- rank(y_both_plus_2$NN2_8h_BR1, ties.method = ties_method, na.last = na_method)

y_both_plus_2$SG17M_4h_BR1_rank <- rank(y_both_plus_2$SG17M_4h_BR1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$SG17M_4h_BR2_rank <- rank(y_both_plus_2$SG17M_4h_BR2, ties.method = ties_method, na.last = na_method)
y_both_plus_2$SG17M_4h_BR3_rank <- rank(y_both_plus_2$SG17M_4h_BR3, ties.method = ties_method, na.last = na_method)
y_both_plus_2$SG17M_8h_BR1_rank <- rank(y_both_plus_2$SG17M_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$SG17M_8h_BR2_rank <- rank(y_both_plus_2$SG17M_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$SG17M_8h_BR3_rank <- rank(y_both_plus_2$SG17M_8h_BR1, ties.method = ties_method, na.last = na_method)

y_both_plus_2$TEX_NN2_8h_1_rank <- rank(y_both_plus_2$TEX_NN2_8h_1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$TEX_NN2_8h_2_rank <- rank(y_both_plus_2$TEX_NN2_8h_2, ties.method = ties_method, na.last = na_method)
y_both_plus_2$TEX_NN2_4h_1_rank <- rank(y_both_plus_2$TEX_NN2_4h_1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$TEX_NN2_4h_2_rank <- rank(y_both_plus_2$TEX_NN2_4h_2, ties.method = ties_method, na.last = na_method)
y_both_plus_2$X0TEX_NN2_8h_1_rank <- rank(y_both_plus_2$X0TEX_NN2_8h_1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$X0TEX_NN2_8h_2_rank <- rank(y_both_plus_2$X0TEX_NN2_8h_2, ties.method = ties_method, na.last = na_method)
y_both_plus_2$X0TEX_NN2_4h_1_rank <- rank(y_both_plus_2$X0TEX_NN2_4h_1, ties.method = ties_method, na.last = na_method)
y_both_plus_2$X0TEX_NN2_4h_2_rank <- rank(y_both_plus_2$X0TEX_NN2_4h_2, ties.method = ties_method, na.last = na_method)

y_both_plus_2$ref_NN2_rank <- rank(y_both_plus_2$ref_NN2, ties.method = "average", na.last = na_method)
y_both_plus_2$ref_SG17M_rank <- rank(y_both_plus_2$ref_NN2, ties.method = "average", na.last = na_method)

y_both_plus_3 <- y_both_plus_2[,-c(1:20)]
y_both_plus_3$NN2_4h_BR3_rank <- NULL
y_both_plus_3$NN2_8h_BR3_rank <- NULL
y_both_plus_3$SG17M_4h_BR3_rank <- NULL
y_both_plus_3$SG17M_8h_BR3_rank <- NULL

y_both_minus_2  <- y_both_minus_2[!rownames(y_both_minus_2)%in%missing_genes_tRNA,]
y_both_minus_2$NN2_4h_BR1_rank_minus <- rank(y_both_minus_2$NN2_4h_BR1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$NN2_4h_BR2_rank_minus <- rank(y_both_minus_2$NN2_4h_BR2, ties.method = ties_method, na.last = na_method)
y_both_minus_2$NN2_4h_BR3_rank_minus <- rank(y_both_minus_2$NN2_4h_BR3, ties.method = ties_method, na.last = na_method)
y_both_minus_2$NN2_8h_BR1_rank_minus <- rank(y_both_minus_2$NN2_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$NN2_8h_BR2_rank_minus <- rank(y_both_minus_2$NN2_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$NN2_8h_BR3_rank_minus <- rank(y_both_minus_2$NN2_8h_BR1, ties.method = ties_method, na.last = na_method)

y_both_minus_2$SG17M_4h_BR1_rank_minus <- rank(y_both_minus_2$SG17M_4h_BR1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$SG17M_4h_BR2_rank_minus <- rank(y_both_minus_2$SG17M_4h_BR2, ties.method = ties_method, na.last = na_method)
y_both_minus_2$SG17M_4h_BR3_rank_minus <- rank(y_both_minus_2$SG17M_4h_BR3, ties.method = ties_method, na.last = na_method)
y_both_minus_2$SG17M_8h_BR1_rank_minus <- rank(y_both_minus_2$SG17M_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$SG17M_8h_BR2_rank_minus <- rank(y_both_minus_2$SG17M_8h_BR1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$SG17M_8h_BR3_rank_minus <- rank(y_both_minus_2$SG17M_8h_BR1, ties.method = ties_method, na.last = na_method)

y_both_minus_2$TEX_NN2_8h_1_rank_minus <- rank(y_both_minus_2$TEX_NN2_8h_1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$TEX_NN2_8h_2_rank_minus <- rank(y_both_minus_2$TEX_NN2_8h_2, ties.method = ties_method, na.last = na_method)
y_both_minus_2$TEX_NN2_4h_1_rank_minus <- rank(y_both_minus_2$TEX_NN2_4h_1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$TEX_NN2_4h_2_rank_minus <- rank(y_both_minus_2$TEX_NN2_4h_2, ties.method = ties_method, na.last = na_method)
y_both_minus_2$X0TEX_NN2_8h_1_rank_minus <- rank(y_both_minus_2$X0TEX_NN2_8h_1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$X0TEX_NN2_8h_2_rank_minus <- rank(y_both_minus_2$X0TEX_NN2_8h_2, ties.method = ties_method, na.last = na_method)
y_both_minus_2$X0TEX_NN2_4h_1_rank_minus <- rank(y_both_minus_2$X0TEX_NN2_4h_1, ties.method = ties_method, na.last = na_method)
y_both_minus_2$X0TEX_NN2_4h_2_rank_minus <- rank(y_both_minus_2$X0TEX_NN2_4h_2, ties.method = ties_method, na.last = na_method)

y_both_minus_3 <- y_both_minus_2[,-c(1:21)]
y_both_minus_3$NN2_4h_BR1_rank_minus <- NULL
y_both_minus_3$NN2_8h_BR1_rank_minus <- NULL
y_both_minus_3$SG17M_4h_BR1_rank_minus <- NULL
y_both_minus_3$SG17M_8h_BR1_rank_minus <- NULL

y_both_plus_minus <- data.frame(cbind(y_both_minus_3, y_both_plus_3))
y_both_plus_minus$ref_NN2 <-NULL
colnames(y_both_plus_minus)
cluster_names_rows <- c("Ala-Anticodon-GGC", 
                        "Ala-Anticodon-TGC",
                        "Arg-Anticodon-ACG", 
                        "Arg-Anticodon-CCG", 
                        "Arg-Anticodon-CCT", 
                        "Arg-Anticodon-TCT", 
                        "Asn-Anticodon-GTT",
                        "Asp-Anticodon-GTC",
                        "Cys-Anticodon-GCA",
                        "Gln-Anticodon-TTG",
                        "Glu-Anticodon-TTC",
                        "Gly-Anticodon-CCC",
                        "Gly-Anticodon-GCC",
                        "Gly-Anticodon-TCC",
                        "His-Anticodon-GTG",
                        "Ile-Anticodon-GAT",
                        "Leu-Anticodon-CAA",
                        "Leu-Anticodon-CAG",
                        "Leu-Anticodon-GAG",
                        "Leu-Anticodon-TAG",    
                        "Lys-Anticodon-TTT",
                        "Lys-Anticodon-TTT (pKLC102)", 
                        "Met-Anticodon-CAT",
                        "Phe-Anticodon-GAA",
                        "Pro-Anticodon-CGG",
                        "Pro-Anticodon-GGG",
                        "Pro-Anticodon-TGG",       
                        "Ser-Anticodon-CGA",
                        "Ser-Anticodon-CGA (pseudo)",
                        "Ser-Anticodon-GCT",
                        "Ser-Anticodon-GGA",
                        "Ser-Anticodon-TGA",
                        "Thr-Anticodon-CGT", 
                        "Thr-Anticodon-GGT",
                        "Thr-Anticodon-TGT",
                        "Trp-Anticodon-CCA", 
                        "Tyr-Anticodon-GTA",
                        "tmRNA-SsrA",
                        "Val-Anticodon-GAC", 
                        "Val-Anticodon-TAC")

cluster_names_columns <- c("NN2-4h-ONT (antisense)", 
                           "NN2-4h-ONT (antisense)",
                           "NN2-8h-ONT (antisense)", 
                           "NN2-8h-ONT (antisense)",
                           "SG17M-4h-ONT (antisense)", 
                           "SG17M-4h-ONT (antisense)",
                           "SG17M-8h-ONT (antisense)", 
                           "SG17M-8h-ONT (antisense)",
                           "NN2-8h-TEX (antisense)", 
                           "NN2-8h-TEX (antisense)",  
                           "NN2-4h-TEX (antisense)", 
                           "NN2-4h-TEX (antisense)", 
                           "NN2-8h-0TEX (antisense)", 
                           "NN2-8h-0TEX (antisense)",
                           "NN2-4h-0TEX (antisense)", 
                           "NN2-4h-0TEX (antisense)",
                           "NN2-4h-ONT (sense)", 
                           "NN2-4h-ONT (sense)", 
                           "NN2-8h-ONT (sense)", 
                           "NN2-8h-ONT (sense)", 
                           "SG17M-4h-ONT (sense)", 
                           "SG17M-4h-ONT (sense)",
                           "SG17M-8h-ONT (sense)", 
                           "SG17M-8h-ONT (sense)", 
                           "NN2-8h-TEX (sense)", 
                           "NN2-8h-TEX (sense)", 
                           "NN2-4h-TEX (sense)", 
                           "NN2-4h-TEX (sense)", 
                           "NN2-8h-0TEX (sense)", 
                           "NN2-8h-0TEX (sense)", 
                           "NN2-4h-0TEX (sense)", 
                           "NN2-4h-0TEX (sense)", 
                           "NN2 (codon usage)", 
                           "SG17M (codon usage)")

codonUsage_plot <- pheatmap(y_both_plus_minus, scale = "none", treeheight_row = 30,
                            cluster_rows = TRUE, cutree_rows = 5,
                            cutree_cols = 4, cluster_cols = TRUE, treeheight_col = 30,
                            clustering_distance_cols = "euclidean", 
                            color = colorRampPalette((brewer.pal(n = 7, name ="OrRd")))(100),
                            cellheight = 8, cellwidth = 11, 
                            border_color = "white", na_col = "white",
                            labels_row = cluster_names_rows,
                            labels_col = cluster_names_columns, 
                            angle_col = 90, fontsize = 7)

codonUsage_tRNA_heatmap_gg <- as.ggplot(codonUsage_plot, scale = 1.2)

ggsave(codonUsage_tRNA_heatmap_gg, filename="figures_all/Figure_03.tif", device="tiff",dpi=300,units="cm",width=23,height=17,bg="white") 

NN2_codon_usage3_sub <- NN2_codon_usage3[order(NN2_codon_usage3$Anticodon3),]
tRNA_sense_corr <- y_both_plus
tRNA_sense_corr <- tRNA_sense_corr[order(rownames(tRNA_sense_corr)),]
tRNA_sense_corr$ref_NN2 <- NN2_codon_usage3_sub$div1000

tRNA_sense_corr <- subset(tRNA_sense_corr, rownames(tRNA_sense_corr) !="tRNA-UndetNNN")
tRNA_sense_corr <- subset(tRNA_sense_corr, rownames(tRNA_sense_corr) !="tRNA-SerCGA-pseudo")
tRNA_sense_corr[is.na(tRNA_sense_corr)] <- 0
tRNA_sense_corr$SG17M_4h_BR1 <- as.numeric(as.character(tRNA_sense_corr$SG17M_4h_BR1))
tRNA_sense_corr$SG17M_4h_BR2 <- as.numeric(as.character(tRNA_sense_corr$SG17M_4h_BR2))
tRNA_sense_corr$SG17M_4h_BR3 <- as.numeric(as.character(tRNA_sense_corr$SG17M_4h_BR3))
tRNA_sense_corr$SG17M_8h_BR1 <- as.numeric(as.character(tRNA_sense_corr$SG17M_8h_BR1))
tRNA_sense_corr$SG17M_8h_BR2 <- as.numeric(as.character(tRNA_sense_corr$SG17M_8h_BR2))
tRNA_sense_corr$SG17M_8h_BR3 <- as.numeric(as.character(tRNA_sense_corr$SG17M_8h_BR3))
tRNA_sense_corr$NN2_4h_BR1 <- as.numeric(as.character(tRNA_sense_corr$NN2_4h_BR1))
tRNA_sense_corr$NN2_4h_BR2 <- as.numeric(as.character(tRNA_sense_corr$NN2_4h_BR2))
tRNA_sense_corr$NN2_4h_BR3 <- as.numeric(as.character(tRNA_sense_corr$NN2_4h_BR3))
tRNA_sense_corr$NN2_8h_BR1 <- as.numeric(as.character(tRNA_sense_corr$NN2_8h_BR1))
tRNA_sense_corr$NN2_8h_BR2 <- as.numeric(as.character(tRNA_sense_corr$NN2_8h_BR2))
tRNA_sense_corr$NN2_8h_BR3 <- as.numeric(as.character(tRNA_sense_corr$NN2_8h_BR3))
tRNA_sense_corr$TEX_NN2_8h_1 <- as.numeric(as.character(tRNA_sense_corr$TEX_NN2_8h_1))
tRNA_sense_corr$TEX_NN2_8h_2 <- as.numeric(as.character(tRNA_sense_corr$TEX_NN2_8h_2))
tRNA_sense_corr$TEX_NN2_4h_1 <- as.numeric(as.character(tRNA_sense_corr$TEX_NN2_4h_1))
tRNA_sense_corr$TEX_NN2_4h_2 <- as.numeric(as.character(tRNA_sense_corr$TEX_NN2_4h_2))
tRNA_sense_corr$X0TEX_NN2_8h_1 <- as.numeric(as.character(tRNA_sense_corr$X0TEX_NN2_8h_1))
tRNA_sense_corr$X0TEX_NN2_8h_2 <- as.numeric(as.character(tRNA_sense_corr$X0TEX_NN2_8h_2))
tRNA_sense_corr$X0TEX_NN2_4h_1 <- as.numeric(as.character(tRNA_sense_corr$X0TEX_NN2_4h_1))
tRNA_sense_corr$X0TEX_NN2_4h_2 <- as.numeric(as.character(tRNA_sense_corr$X0TEX_NN2_4h_2))
tRNA_sense_corr$ref_SG17M <- as.numeric(as.character(tRNA_sense_corr$ref_NN2))
tRNA_sense_corr$ref_NN2<- as.numeric(as.character(tRNA_sense_corr$ref_NN2))
tRNA_sense_corr <- round(tRNA_sense_corr,2)
data_matrix_1 <- cor(tRNA_sense_corr, method = "spearman")

testRes_Sense <- corrplot::cor.mtest(data_matrix_1, conf.level = 0.95)
colnames(data_matrix_1) <- c("ONT_SG17M_4h","ONT_SG17M_4h","ONT_SG17M_4h",
                             "ONT_SG17M_8h","ONT_SG17M_8h","ONT_SG17M_8h",
                             "ONT_NN2_4h",  "ONT_NN2_4h",  "ONT_NN2_4h", 
                             "ONT_NN2_8h",  "ONT_NN2_8h",  "ONT_NN2_8h",
                             "TEX_NN2_8h","TEX_NN2_8h","TEX_NN2_4h","TEX_NN2_4h",
                             "0TEX_NN2_8h", "0TEX_NN2_8h", "0TEX_NN2_4h", "0TEX_NN2_4h",
                             "NN2 (codon usage)", "SG17M (codon usage)")

rownames(data_matrix_1) <- c("ONT_SG17M_4h","ONT_SG17M_4h","ONT_SG17M_4h",
                             "ONT_SG17M_8h","ONT_SG17M_8h","ONT_SG17M_8h",
                             "ONT_NN2_4h",  "ONT_NN2_4h",  "ONT_NN2_4h", 
                             "ONT_NN2_8h",  "ONT_NN2_8h",  "ONT_NN2_8h",
                             "TEX_NN2_8h","TEX_NN2_8h","TEX_NN2_4h","TEX_NN2_4h",
                             "0TEX_NN2_8h", "0TEX_NN2_8h", "0TEX_NN2_4h", "0TEX_NN2_4h",
                             "NN2 (codon usage)", "SG17M (codon usage)")

tRNA_antisense_corr <- y_both_minus
tRNA_antisense_corr <- tRNA_antisense_corr[order(rownames(tRNA_antisense_corr)),]
tRNA_antisense_corr$ref_NN2 <- NN2_codon_usage3_sub$div1000

tRNA_antisense_corr <- subset(tRNA_antisense_corr, rownames(tRNA_antisense_corr) !="tRNA-UndetNNN")
tRNA_antisense_corr <- subset(tRNA_antisense_corr, rownames(tRNA_antisense_corr) !="tRNA-SerCGA-pseudo")
tRNA_antisense_corr[is.na(tRNA_antisense_corr)] <- 0
tRNA_antisense_corr$SG17M_4h_BR1 <- as.numeric(as.character(tRNA_antisense_corr$SG17M_4h_BR1))
tRNA_antisense_corr$SG17M_4h_BR2 <- as.numeric(as.character(tRNA_antisense_corr$SG17M_4h_BR2))
tRNA_antisense_corr$SG17M_4h_BR3 <- as.numeric(as.character(tRNA_antisense_corr$SG17M_4h_BR3))
tRNA_antisense_corr$SG17M_8h_BR1 <- as.numeric(as.character(tRNA_antisense_corr$SG17M_8h_BR1))
tRNA_antisense_corr$SG17M_8h_BR2 <- as.numeric(as.character(tRNA_antisense_corr$SG17M_8h_BR2))
tRNA_antisense_corr$SG17M_8h_BR3 <- as.numeric(as.character(tRNA_antisense_corr$SG17M_8h_BR3))
tRNA_antisense_corr$NN2_4h_BR1 <- as.numeric(as.character(tRNA_antisense_corr$NN2_4h_BR1))
tRNA_antisense_corr$NN2_4h_BR2 <- as.numeric(as.character(tRNA_antisense_corr$NN2_4h_BR2))
tRNA_antisense_corr$NN2_4h_BR3 <- as.numeric(as.character(tRNA_antisense_corr$NN2_4h_BR3))
tRNA_antisense_corr$NN2_8h_BR1 <- as.numeric(as.character(tRNA_antisense_corr$NN2_8h_BR1))
tRNA_antisense_corr$NN2_8h_BR2 <- as.numeric(as.character(tRNA_antisense_corr$NN2_8h_BR2))
tRNA_antisense_corr$NN2_8h_BR3 <- as.numeric(as.character(tRNA_antisense_corr$NN2_8h_BR3))
tRNA_antisense_corr$TEX_NN2_8h_1 <- as.numeric(as.character(tRNA_antisense_corr$TEX_NN2_8h_1))
tRNA_antisense_corr$TEX_NN2_8h_2 <- as.numeric(as.character(tRNA_antisense_corr$TEX_NN2_8h_2))
tRNA_antisense_corr$TEX_NN2_4h_1 <- as.numeric(as.character(tRNA_antisense_corr$TEX_NN2_4h_1))
tRNA_antisense_corr$TEX_NN2_4h_2 <- as.numeric(as.character(tRNA_antisense_corr$TEX_NN2_4h_2))
tRNA_antisense_corr$X0TEX_NN2_8h_1 <- as.numeric(as.character(tRNA_antisense_corr$X0TEX_NN2_8h_1))
tRNA_antisense_corr$X0TEX_NN2_8h_2 <- as.numeric(as.character(tRNA_antisense_corr$X0TEX_NN2_8h_2))
tRNA_antisense_corr$X0TEX_NN2_4h_1 <- as.numeric(as.character(tRNA_antisense_corr$X0TEX_NN2_4h_1))
tRNA_antisense_corr$X0TEX_NN2_4h_2 <- as.numeric(as.character(tRNA_antisense_corr$X0TEX_NN2_4h_2))
tRNA_antisense_corr$ref_SG17M <- as.numeric(as.character(tRNA_antisense_corr$ref_NN2))
tRNA_antisense_corr$ref_NN2<- as.numeric(as.character(tRNA_antisense_corr$ref_NN2))
tRNA_antisense_corr <- round(tRNA_antisense_corr,2)
data_matrix_1_antisense <- cor(tRNA_antisense_corr, method = "spearman")

testRes <- corrplot::cor.mtest(data_matrix_1_antisense, conf.level = 0.99)


colnames(data_matrix_1_antisense) <- c("ONT_SG17M_4h","ONT_SG17M_4h","ONT_SG17M_4h",
                                       "ONT_SG17M_8h","ONT_SG17M_8h","ONT_SG17M_8h",
                                       "ONT_NN2_4h",  "ONT_NN2_4h",  "ONT_NN2_4h", 
                                       "ONT_NN2_8h",  "ONT_NN2_8h",  "ONT_NN2_8h",
                                       "TEX_NN2_8h","TEX_NN2_8h","TEX_NN2_4h","TEX_NN2_4h",
                                       "0TEX_NN2_8h", "0TEX_NN2_8h", "0TEX_NN2_4h", "0TEX_NN2_4h",
                                       "NN2 (codon usage)", "SG17M (codon usage)")

rownames(data_matrix_1_antisense) <- c("ONT_SG17M_4h","ONT_SG17M_4h","ONT_SG17M_4h",
                                       "ONT_SG17M_8h","ONT_SG17M_8h","ONT_SG17M_8h",
                                       "ONT_NN2_4h",  "ONT_NN2_4h",  "ONT_NN2_4h", 
                                       "ONT_NN2_8h",  "ONT_NN2_8h",  "ONT_NN2_8h",
                                       "TEX_NN2_8h","TEX_NN2_8h","TEX_NN2_4h","TEX_NN2_4h",
                                       "0TEX_NN2_8h", "0TEX_NN2_8h", "0TEX_NN2_4h", "0TEX_NN2_4h",
                                       "NN2 (codon usage)", "SG17M (codon usage)")

col2 <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
tiff(filename="figures_all/correlation_tRNA_sense.tif", res = 300, height = 17.5, width = 17, units = "cm")
corrplot(data_matrix_1, method="shade", order="hclust", type="lower",number.cex = 0.5, number.digits = 1, 
         addCoef.col = "black", tl.col="black", tl.cex = 0.6, tl.srt = 45, diag=TRUE, addrect = 4,
          rect.col = 'white', rect.lwd = 1,  col=brewer.pal(n=8, name="Spectral"))
dev.off()



library(readr)
trna_quantity_2 <- read_delim("trna_quantity_2.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
trna_quantity_2 <- data.frame(trna_quantity_2)

trna_quantity_2$Expression <- factor(trna_quantity_2$Expression, levels = c("absent", "low", "high"))

ggplot(trna_quantity_2, aes(x=transcript_type, y=sqrt(MFE*-1))) +
  geom_boxplot(outlier.size = 0, width=0.4) +
  geom_jitter(width=0.3) +
  facet_wrap(~Expression) + 
  stat_compare_means()



trna_quantity_2_high <- subset(trna_quantity_2, Expression == "high")
trna_high_R <- wilcoxonR(trna_quantity_2_high$MFE, g=trna_quantity_2_high$transcript_type, ci=TRUE)
trna_high_p <- wilcox.test(trna_quantity_2_high$MFE ~ trna_quantity_2_high$transcript_type)


trna_quantity_2_low <- subset(trna_quantity_2, Expression == "low")
trna_low_R <- wilcoxonR(trna_quantity_2_low$MFE, g=trna_quantity_2_low$transcript_type, ci=TRUE)
trna_low_p <- wilcox.test(trna_quantity_2_low$MFE ~ trna_quantity_2_low$transcript_type)


trna_quantity_2_absent <- subset(trna_quantity_2, Expression == "absent")
trna_absent_R <- wilcoxonR(trna_quantity_2_absent$MFE, g=trna_quantity_2_absent$transcript_type, ci=TRUE)
trna_absent_p <- wilcox.test(trna_quantity_2_absent$MFE ~ trna_quantity_2_absent$transcript_type)

trna_quantity_2_absent$tRNA.anticodon <- factor(trna_quantity_2_absent$tRNA.anticodon)


trna_quantity_3 <- trna_quantity_2[!trna_quantity_2$tRNA.anticodon%in%remove_ids,]

tRNA_paired <-
  ggpaired(trna_quantity_2, x="transcript_type", y="MFE", 
         color="transcript_type", id="tRNA.anticodon",palette="lancet", line.size = 0.3) +
  stat_compare_means(paired=TRUE, label.x.npc = "left", label.y = -15, size=3) + 
  facet_wrap(~Expression) + theme_pubr(border=TRUE, base_size = 10) +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.title = element_text(size=10),
        strip.text = element_text(size=10),
        axis.text = element_text(size=10)) + xlab(expression(bold("RNA transcript"))) + 
  ylab(expression(bold("Free energy (kcal/mol)"))) + ylim(-50,-10)


tiff(filename="paired_test.tif", units="cm", pointsize=10, width=19, height=11, type="cairo", res=300)
tRNA_paired
dev.off()
