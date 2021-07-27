
# title: "tRNA analysis"
# author: "Marie Pust"
# date: "08 07 2021"

# Clean environment
rm(list = ls())

# Packages ####
library('Rsubread')
library('ggplot2')
library('edgeR')
library('pheatmap')
library('plyr')
library('dplyr')
library('tidyr')
library('ggpubr')
library('vegan')
library('lemon')
library('factoextra')
library('stringr')
library('ggplotify')
library('tidyverse')
library('RColorBrewer')
library('purrr')
library('readr')
library('ggdendro')
library('ggrepel')
library('ggcorrplot')
library('Hmisc')
library('GenomicRanges')
library('GenomicAlignments')
library('Rsamtools')
library('seqinr')
library('corrplot')

# set working directory
setwd("/MinionExperiments202009/d_Ranalysis")

# load input
input_bam_NN2<- list.files(path = "/NN2_run/d_raw_bam", pattern = ".bam", full.names = TRUE)
input_bam_SG17M <- list.files(path = "/SG17M_run/d_raw_bam", pattern = ".bam", full.names = TRUE)

UserDefinedAnnotationRef <- "NN2_ENO_SG17M_curated.gtf"
UserDefinedAnnotationRef_NN2 <- "NN2_ENO_curated.gtf"
UserDefinedAnnotationRef_SG17M <- "SG17M_ENO_curated.gtf"

# prepare annotation files, NN2
import_gtf_NN2 <- read.table(UserDefinedAnnotationRef_NN2, sep = '\t', header = FALSE)
colnames(import_gtf_NN2) <- c("Isolate", "database", "featureType", "start", "end", "V6", "strandType", "V8", "V9")
import_gtf_NN2$V6 <- NULL
import_gtf_NN2$V8 <- NULL
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "transcript_id ", "")
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "gene_name ", "")
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "tRNA_type ", "")
import_gtf_NN2$V9 <- str_replace_all(import_gtf_NN2$V9, "rRNA_type ", "")
import_gtf_NN2 <- separate(data = import_gtf_NN2, col = V9, into = c("transcript_id", "gene_name2"), sep = "\\;")
import_gtf_NN2$gene_name3 <- ifelse(import_gtf_NN2$gene_name2 == "", import_gtf_NN2$featureType, import_gtf_NN2$gene_name2)
import_gtf_NN2$gene_name <- ifelse(import_gtf_NN2$gene_name3 == "CDS", "CDS-undefined",
                               ifelse(import_gtf_NN2$gene_name3 == "exon", "ncRNA-undefined", 
                                      ifelse(grepl("tRNA", import_gtf_NN2$gene_name3), "tRNA",
                                             ifelse(grepl("rRNA", import_gtf_NN2$gene_name3), "rRNA", "CDS"))))
import_gtf_NN2$gene_name3 <- NULL
import_gtf_NN2$length <- import_gtf_NN2$end - import_gtf_NN2$start
import_gtf_NN2_island <- import_gtf_NN2 
import_gtf_NN2 <- subset(import_gtf_NN2, database != "IslandViewer4")

# prepare annotation files, SG17M
import_gtf_SG17M <- read.table(UserDefinedAnnotationRef_SG17M, sep = '\t', header = FALSE)
colnames(import_gtf_SG17M) <- c("Isolate", "database", "featureType", "start", "end", "V6", "strandType", "V8", "V9")
import_gtf_SG17M$V6 <- NULL
import_gtf_SG17M$V8 <- NULL
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "transcript_id ", "")
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "gene_name ", "")
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "tRNA_type ", "")
import_gtf_SG17M$V9 <- str_replace_all(import_gtf_SG17M$V9, "rRNA_type ", "")
import_gtf_SG17M <- separate(data = import_gtf_SG17M, col = V9, into = c("transcript_id", "gene_name2"), sep = "\\;")
import_gtf_SG17M$gene_name3 <- ifelse(import_gtf_SG17M$gene_name2 == "", import_gtf_SG17M$featureType, import_gtf_SG17M$gene_name2)
import_gtf_SG17M$gene_name <- ifelse(import_gtf_SG17M$gene_name3 == "CDS", "CDS-undefined",
                               ifelse(import_gtf_SG17M$gene_name3 == "exon", "ncRNA-undefined", 
                                      ifelse(grepl("tRNA", import_gtf_SG17M$gene_name3), "tRNA",
                                             ifelse(grepl("rRNA", import_gtf_SG17M$gene_name3), "rRNA", "CDS"))))
import_gtf_SG17M$gene_name3 <- NULL
import_gtf_SG17M$length <- import_gtf_SG17M$end - import_gtf_SG17M$start
import_gtf_SG17M_island <- import_gtf_SG17M 
import_gtf_SG17M <- subset(import_gtf_SG17M, database != "IslandViewer4")

# Obtain count data
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

featureCounts_PA_antisenseSG17M <- 
  featureCounts(input_bam_SG17M, annot.ext = UserDefinedAnnotationRef_SG17M, 
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 2, isLongRead = TRUE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

featureCounts_PA_senseSG17M <- 
  featureCounts(input_bam_SG17M, annot.ext = UserDefinedAnnotationRef_SG17M,
                isGTFAnnotationFile = TRUE, GTF.attrType="transcript_id", 
                GTF.featureType = c("exon", "CDS"), strandSpecific = 1, isLongRead = TRUE, 
                useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                countMultiMappingReads = TRUE)

featureCounts_PA_antisense_df_NN2 <- featureCounts_PA_antisenseNN2$counts
featureCounts_PA_antisense_df_NN2 <- data.frame(featureCounts_PA_antisense_df_NN2)
featureCounts_PA_sense_df_NN2 <- featureCounts_PA_senseNN2$counts
featureCounts_PA_sense_df_NN2 <- data.frame(featureCounts_PA_sense_df_NN2)
colnames(featureCounts_PA_antisense_df_NN2) <- c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3", "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3")
colnames(featureCounts_PA_sense_df_NN2) <- c("NN2_4h_BR1", "NN2_4h_BR2", "NN2_4h_BR3", "NN2_8h_BR1", "NN2_8h_BR2", "NN2_8h_BR3")
tRNA_antisense_NN2 <- featureCounts_PA_antisense_df_NN2
mRNA_antisense_NN2 <- featureCounts_PA_antisense_df_NN2
tRNA_sense_NN2 <- featureCounts_PA_sense_df_NN2
mRNA_sense_NN2 <- featureCounts_PA_sense_df_NN2

featureCounts_PA_antisense_df_SG17M <- featureCounts_PA_antisenseSG17M$counts
featureCounts_PA_antisense_df_SG17M <- data.frame(featureCounts_PA_antisense_df_SG17M)
featureCounts_PA_sense_df_SG17M <- featureCounts_PA_senseSG17M$counts
featureCounts_PA_sense_df_SG17M <- data.frame(featureCounts_PA_sense_df_SG17M)

colnames(featureCounts_PA_antisense_df_SG17M) <-c("SG17M_4h_BR1", "SG17M_4h_BR2", "SG17M_4h_BR3", "SG17M_8h_BR1","SG17M_8h_BR2","SG17M_8h_BR3")
colnames(featureCounts_PA_sense_df_SG17M) <-c("SG17M_4h_BR1", "SG17M_4h_BR2", "SG17M_4h_BR3", "SG17M_8h_BR1","SG17M_8h_BR2","SG17M_8h_BR3")
tRNA_antisense_SG17M <- featureCounts_PA_antisense_df_SG17M
mRNA_antisense_SG17M <- featureCounts_PA_antisense_df_SG17M
tRNA_sense_SG17M <- featureCounts_PA_sense_df_SG17M
mRNA_sense_SG17M <- featureCounts_PA_sense_df_SG17M

gene_name_list_2_NN2 <- import_gtf_NN2$gene_name
names(gene_name_list_2_NN2) <- import_gtf_NN2$transcript_id

mylist_antisense_NN2 <- c()
for(i in rownames(featureCounts_PA_antisense_df_NN2)){
    newelem <- ifelse(i %in% names(gene_name_list_2_NN2), gene_name_list_2_NN2[[i]], i)
    mylist_antisense_NN2 <- append(mylist_antisense_NN2, newelem)}
featureCounts_PA_antisense_df_NN2$name2 <- mylist_antisense_NN2
featureCounts_PA_antisense_df_NN2$name <- map(strsplit(featureCounts_PA_antisense_df_NN2$name2, "_"),1)
featureCounts_PA_antisense_df_NN2$name2 <- NULL

rownames(featureCounts_PA_antisense_df_NN2) <- NULL
featureCounts_PA_antisense_df_NN2 <- featureCounts_PA_antisense_df_NN2 %>% group_by(name) %>% summarise_all(funs(sum))
featureCounts_PA_antisense_df_NN2 <- data.frame(featureCounts_PA_antisense_df_NN2)
rownames(featureCounts_PA_antisense_df_NN2) <- featureCounts_PA_antisense_df_NN2$name
featureCounts_PA_antisense_df_NN2$name <- NULL

featureCounts_PA_antisense_df_t_NN2 <- data.frame(t(featureCounts_PA_antisense_df_NN2))
featureCounts_PA_antisense_df_t_NN2$rowSums <- rowSums(featureCounts_PA_antisense_df_t_NN2)
featureCounts_PA_antisense_df_t_rel_NN2 <- featureCounts_PA_antisense_df_t_NN2[,1:5] / featureCounts_PA_antisense_df_t_NN2$rowSums
featureCounts_PA_antisense_df_t_rel_NN2 <- round(featureCounts_PA_antisense_df_t_rel_NN2,3)
featureCounts_PA_antisense_df_t_rel_NN2$featureType <- rownames(featureCounts_PA_antisense_df_t_rel_NN2)
featureCounts_PA_antisense_df_t_rel_L_NN2 <- gather(featureCounts_PA_antisense_df_t_rel_NN2, key="type", value="abundance", -featureType)
featureCounts_PA_antisense_df_t_rel_L_NN2$strand_type = "Antisense RNA"
featureCounts_PA_antisense_df_t_rel_L_NN2$isolate = "NN2"

mylist_sense_NN2 <- c()
for(i in rownames(featureCounts_PA_sense_df_NN2)){
    newelem <- ifelse(i %in% names(gene_name_list_2_NN2), gene_name_list_2_NN2[[i]], i)
    mylist_sense_NN2 <- append(mylist_sense_NN2, newelem)}
featureCounts_PA_sense_df_NN2$name2 <- mylist_sense_NN2
featureCounts_PA_sense_df_NN2$name <- map(strsplit(featureCounts_PA_sense_df_NN2$name2, "_"),1)
featureCounts_PA_sense_df_NN2$name2 <- NULL

rownames(featureCounts_PA_sense_df_NN2) <- NULL
featureCounts_PA_sense_df_NN2 <- featureCounts_PA_sense_df_NN2 %>% group_by(name) %>% summarise_all(funs(sum))
featureCounts_PA_sense_df_NN2 <- data.frame(featureCounts_PA_sense_df_NN2)
rownames(featureCounts_PA_sense_df_NN2) <- featureCounts_PA_sense_df_NN2$name
featureCounts_PA_sense_df_NN2$name <- NULL

featureCounts_PA_sense_df_t_NN2 <- data.frame(t(featureCounts_PA_sense_df_NN2))
featureCounts_PA_sense_df_t_NN2$rowSums <- rowSums(featureCounts_PA_sense_df_t_NN2)
featureCounts_PA_sense_df_t_rel_NN2 <- featureCounts_PA_sense_df_t_NN2[,1:5] / featureCounts_PA_sense_df_t_NN2$rowSums
featureCounts_PA_sense_df_t_rel_NN2 <- round(featureCounts_PA_sense_df_t_rel_NN2,3)
featureCounts_PA_sense_df_t_rel_NN2$featureType <- rownames(featureCounts_PA_sense_df_t_rel_NN2)
featureCounts_PA_sense_df_t_rel_L_NN2 <- gather(featureCounts_PA_sense_df_t_rel_NN2, key="type", value="abundance", -featureType)
featureCounts_PA_sense_df_t_rel_L_NN2$strand_type = "Sense RNA"
featureCounts_PA_sense_df_t_rel_L_NN2$isolate = "NN2"

gene_name_list_2_SG17M <- import_gtf_SG17M$gene_name
names(gene_name_list_2_SG17M) <- import_gtf_SG17M$transcript_id

mylist_antisense_SG17M <- c()
for(i in rownames(featureCounts_PA_antisense_df_SG17M)){
    newelem <- ifelse(i %in% names(gene_name_list_2_SG17M), gene_name_list_2_SG17M[[i]], i)
    mylist_antisense_SG17M <- append(mylist_antisense_SG17M, newelem)}
featureCounts_PA_antisense_df_SG17M$name2 <- mylist_antisense_SG17M
featureCounts_PA_antisense_df_SG17M$name <- map(strsplit(featureCounts_PA_antisense_df_SG17M$name2, "_"),1)
featureCounts_PA_antisense_df_SG17M$name2 <- NULL

rownames(featureCounts_PA_antisense_df_SG17M) <- NULL
featureCounts_PA_antisense_df_SG17M <- featureCounts_PA_antisense_df_SG17M %>% group_by(name) %>% summarise_all(funs(sum))
featureCounts_PA_antisense_df_SG17M <- data.frame(featureCounts_PA_antisense_df_SG17M)
rownames(featureCounts_PA_antisense_df_SG17M) <- featureCounts_PA_antisense_df_SG17M$name
featureCounts_PA_antisense_df_SG17M$name <- NULL

featureCounts_PA_antisense_df_t_SG17M <- data.frame(t(featureCounts_PA_antisense_df_SG17M))
featureCounts_PA_antisense_df_t_SG17M$rowSums <- rowSums(featureCounts_PA_antisense_df_t_SG17M)
featureCounts_PA_antisense_df_t_rel_SG17M <- featureCounts_PA_antisense_df_t_SG17M[,1:5] / featureCounts_PA_antisense_df_t_SG17M$rowSums
featureCounts_PA_antisense_df_t_rel_SG17M <- round(featureCounts_PA_antisense_df_t_rel_SG17M,3)
featureCounts_PA_antisense_df_t_rel_SG17M$featureType <- rownames(featureCounts_PA_antisense_df_t_rel_SG17M)
featureCounts_PA_antisense_df_t_rel_L_SG17M <- gather(featureCounts_PA_antisense_df_t_rel_SG17M, key="type", value="abundance", -featureType)
featureCounts_PA_antisense_df_t_rel_L_SG17M$strand_type = "Antisense RNA"
featureCounts_PA_antisense_df_t_rel_L_SG17M$isolate = "SG17M"

mylist_sense_SG17M <- c()
for(i in rownames(featureCounts_PA_sense_df_SG17M)){
    newelem <- ifelse(i %in% names(gene_name_list_2_SG17M), gene_name_list_2_SG17M[[i]], i)
    mylist_sense_SG17M <- append(mylist_sense_SG17M, newelem)}
featureCounts_PA_sense_df_SG17M$name2 <- mylist_sense_SG17M
featureCounts_PA_sense_df_SG17M$name <- map(strsplit(featureCounts_PA_sense_df_SG17M$name2, "_"),1)
featureCounts_PA_sense_df_SG17M$name2 <- NULL

rownames(featureCounts_PA_sense_df_SG17M) <- NULL
featureCounts_PA_sense_df_SG17M <- featureCounts_PA_sense_df_SG17M %>% group_by(name) %>% summarise_all(funs(sum))
featureCounts_PA_sense_df_SG17M <- data.frame(featureCounts_PA_sense_df_SG17M)
rownames(featureCounts_PA_sense_df_SG17M) <- featureCounts_PA_sense_df_SG17M$name
featureCounts_PA_sense_df_SG17M$name <- NULL

featureCounts_PA_sense_df_t_SG17M <- data.frame(t(featureCounts_PA_sense_df_SG17M))
featureCounts_PA_sense_df_t_SG17M$rowSums <- rowSums(featureCounts_PA_sense_df_t_SG17M)
featureCounts_PA_sense_df_t_rel_SG17M <- featureCounts_PA_sense_df_t_SG17M[,1:5] / featureCounts_PA_sense_df_t_SG17M$rowSums
featureCounts_PA_sense_df_t_rel_SG17M <- round(featureCounts_PA_sense_df_t_rel_SG17M,3)
featureCounts_PA_sense_df_t_rel_SG17M$featureType <- rownames(featureCounts_PA_sense_df_t_rel_SG17M)
featureCounts_PA_sense_df_t_rel_L_SG17M <- gather(featureCounts_PA_sense_df_t_rel_SG17M, key="type", value="abundance", -featureType)
featureCounts_PA_sense_df_t_rel_L_SG17M$strand_type = "Sense RNA"
featureCounts_PA_sense_df_t_rel_L_SG17M$isolate = "SG17M"

merge_sense_df_NN2_SG17M <- data.frame(rbind(featureCounts_PA_sense_df_t_rel_L_NN2, 
                                             featureCounts_PA_sense_df_t_rel_L_SG17M))
merge_sense_df_NN2_SG17M$type1 <- paste(merge_sense_df_NN2_SG17M$type,"_sense")

merge_antisense_df_NN2_SG17M <- data.frame(rbind(featureCounts_PA_antisense_df_t_rel_L_NN2, 
                                             featureCounts_PA_antisense_df_t_rel_L_SG17M))
merge_antisense_df_NN2_SG17M$type1 <- paste(merge_antisense_df_NN2_SG17M$type,"_antisense")

featureCounts_PA_sense_antisense_df_t_rel_L <- data.frame(rbind(merge_sense_df_NN2_SG17M, merge_antisense_df_NN2_SG17M))

featureCounts_PA_sense_antisense_df_t_rel_L$type1b <- with(
  featureCounts_PA_sense_antisense_df_t_rel_L,
  ifelse(type1 == "CDS _antisense", "antisense to CDS",
         ifelse(type1 == "rRNA _antisense", "antisense to rRNA",
                ifelse(type1 == "CDS.undefined _antisense", "antisense to CDS (undefined)",
                       ifelse(type1 == "tRNA _antisense", "antisense to tRNA",
                              ifelse(type1 == "ncRNA.undefined _antisense", 
                                     "antisense to ncRNA (undefined)", type1))))))

table(featureCounts_PA_sense_antisense_df_t_rel_L$type1b)
featureCounts_PA_sense_antisense_df_t_rel_L$type2 <- with(
  featureCounts_PA_sense_antisense_df_t_rel_L,
  ifelse(type1b == "CDS _sense", " CDS",
         ifelse(type1b == "CDS.undefined _sense", " CDS (undefined)",
                ifelse(type1b == "ncRNA.undefined _sense", " ncRNA (undefined)",
                       ifelse(type1b == "rRNA _sense", " rRNA",
                              ifelse(type1b == "tRNA _sense", " tRNA", type1b))))))


featureCounts_PA_sense_antisense_df_t_rel_L$type3 <- factor(
  featureCounts_PA_sense_antisense_df_t_rel_L$type2,
  levels = c(" rRNA", " tRNA", " ncRNA (undefined)", " CDS (undefined)", " CDS",
             "antisense to rRNA", "antisense to tRNA", "antisense to ncRNA (undefined)",
             "antisense to ncRNA", "antisense to CDS (undefined)", "antisense to CDS"))
                                                              
featureCounts_PA_sense_antisense_df_t_rel_L$strand_type2 <- ifelse(
  featureCounts_PA_sense_antisense_df_t_rel_L$strand_type == "Sense RNA", " Sense RNA ", "Antisense RNA") 

amount_rna <- ggplot(featureCounts_PA_sense_antisense_df_t_rel_L) +
  geom_col(aes(x=featureType, y=abundance, fill=reorder(type3))) +
  facet_wrap(~strand_type2) +
  theme_pubr(border=TRUE, base_size=11) + 
  theme(legend.title = element_blank(), legend.position = "bottom", legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"), axis.text.y = element_text(size=10)) + 
  xlab(" ") + ylab("Relative abundance") +
  scale_fill_manual(values=c("skyblue4", "darkolivegreen4","gray5",  "tan2", "brown3", 
                             "skyblue3", "darkolivegreen3", "gray60", "tan1", "brown1")) + 
  coord_flip() +
  guides(fill=guide_legend(ncol=2)) 

amount_rna2 <- ggarrange(amount_rna, labels=c("A"))

featureCounts_PA_antisense_df_t_2_NN2 <- featureCounts_PA_antisense_df_t_NN2
featureCounts_PA_antisense_df_t_2_NN2$rowSums <- NULL
featureCounts_PA_antisense_df_t_2_NN2$coding <- 
  featureCounts_PA_antisense_df_t_2_NN2$CDS + featureCounts_PA_antisense_df_t_2_NN2$CDS.undefined
featureCounts_PA_antisense_df_t_2_NN2$ncRNA <- 
  featureCounts_PA_antisense_df_t_2_NN2$ncRNA.undefined + featureCounts_PA_antisense_df_t_2_NN2$tRNA + featureCounts_PA_antisense_df_t_2_NN2$rRNA
featureCounts_PA_antisense_df_t_2_NN2$CDS <- NULL
featureCounts_PA_antisense_df_t_2_NN2$CDS.undefined <- NULL
featureCounts_PA_antisense_df_t_2_NN2$ncRNA.undefined <- NULL
featureCounts_PA_antisense_df_t_2_NN2$tRNA <- NULL
featureCounts_PA_antisense_df_t_2_NN2$rRNA <- NULL
featureCounts_PA_antisense_df_t_2_NN2$sample <- rownames(featureCounts_PA_antisense_df_t_2_NN2)

featureCounts_PA_antisense_df_t_2_SG17M <- featureCounts_PA_antisense_df_t_SG17M
featureCounts_PA_antisense_df_t_2_SG17M$rowSums <- NULL
featureCounts_PA_antisense_df_t_2_SG17M$coding <- 
  featureCounts_PA_antisense_df_t_2_SG17M$CDS + featureCounts_PA_antisense_df_t_2_SG17M$CDS.undefined
featureCounts_PA_antisense_df_t_2_SG17M$ncRNA <- 
  featureCounts_PA_antisense_df_t_2_SG17M$ncRNA.undefined + featureCounts_PA_antisense_df_t_2_SG17M$tRNA + featureCounts_PA_antisense_df_t_2_SG17M$rRNA
featureCounts_PA_antisense_df_t_2_SG17M$CDS <- NULL
featureCounts_PA_antisense_df_t_2_SG17M$CDS.undefined <- NULL
featureCounts_PA_antisense_df_t_2_SG17M$ncRNA.undefined <- NULL
featureCounts_PA_antisense_df_t_2_SG17M$tRNA <- NULL
featureCounts_PA_antisense_df_t_2_SG17M$rRNA <- NULL
featureCounts_PA_antisense_df_t_2_SG17M$sample <- rownames(featureCounts_PA_antisense_df_t_2_SG17M)

featureCounts_PA_antisense_df_t_2_NN2$strain <- "NN2"
featureCounts_PA_antisense_df_t_2_SG17M$strain <- "SG17M"
featureCounts_PA_antisense_df_t_2_NN2$time <- c("4h", "4h", "4h", "8h", "8h", "8h")
featureCounts_PA_antisense_df_t_2_SG17M$time <- c("4h", "4h", "4h", "8h", "8h", "8h")

featureCounts_PA_antisense_df_t_2_both <- data.frame(rbind(featureCounts_PA_antisense_df_t_2_NN2, featureCounts_PA_antisense_df_t_2_SG17M))
featureCounts_PA_antisense_df_t_2_L <- gather(featureCounts_PA_antisense_df_t_2_both, key="RNA_type", "abundance", -c(sample,strain,time))

featureCounts_PA_sense_df_t_2_NN2 <- featureCounts_PA_sense_df_t_NN2
featureCounts_PA_sense_df_t_2_NN2$rowSums <- NULL
featureCounts_PA_sense_df_t_2_NN2$coding <- 
  featureCounts_PA_sense_df_t_2_NN2$CDS + featureCounts_PA_sense_df_t_2_NN2$CDS.undefined
featureCounts_PA_sense_df_t_2_NN2$ncRNA <- 
  featureCounts_PA_sense_df_t_2_NN2$ncRNA.undefined + featureCounts_PA_sense_df_t_2_NN2$tRNA + featureCounts_PA_sense_df_t_2_NN2$rRNA
featureCounts_PA_sense_df_t_2_NN2$CDS <- NULL
featureCounts_PA_sense_df_t_2_NN2$CDS.undefined <- NULL
featureCounts_PA_sense_df_t_2_NN2$ncRNA.undefined <- NULL
featureCounts_PA_sense_df_t_2_NN2$tRNA <- NULL
featureCounts_PA_sense_df_t_2_NN2$rRNA <- NULL
featureCounts_PA_sense_df_t_2_NN2$sample <- rownames(featureCounts_PA_sense_df_t_2_NN2)

featureCounts_PA_sense_df_t_2_SG17M <- featureCounts_PA_sense_df_t_SG17M
featureCounts_PA_sense_df_t_2_SG17M$rowSums <- NULL
featureCounts_PA_sense_df_t_2_SG17M$coding <- 
  featureCounts_PA_sense_df_t_2_SG17M$CDS + featureCounts_PA_sense_df_t_2_SG17M$CDS.undefined
featureCounts_PA_sense_df_t_2_SG17M$ncRNA <- 
  featureCounts_PA_sense_df_t_2_SG17M$ncRNA.undefined + featureCounts_PA_sense_df_t_2_SG17M$tRNA + featureCounts_PA_sense_df_t_2_SG17M$rRNA
featureCounts_PA_sense_df_t_2_SG17M$CDS <- NULL
featureCounts_PA_sense_df_t_2_SG17M$CDS.undefined <- NULL
featureCounts_PA_sense_df_t_2_SG17M$ncRNA.undefined <- NULL
featureCounts_PA_sense_df_t_2_SG17M$tRNA <- NULL
featureCounts_PA_sense_df_t_2_SG17M$rRNA <- NULL
featureCounts_PA_sense_df_t_2_SG17M$sample <- rownames(featureCounts_PA_sense_df_t_2_SG17M)

featureCounts_PA_sense_df_t_2_NN2$strain <- "NN2"
featureCounts_PA_sense_df_t_2_SG17M$strain <- "SG17M"
featureCounts_PA_sense_df_t_2_NN2$time <- c("4h", "4h", "4h", "8h", "8h", "8h")
featureCounts_PA_sense_df_t_2_SG17M$time <- c("4h", "4h", "4h", "8h", "8h", "8h")

featureCounts_PA_sense_df_t_2_both <- data.frame(rbind(featureCounts_PA_sense_df_t_2_NN2, featureCounts_PA_sense_df_t_2_SG17M))
featureCounts_PA_sense_df_t_2_L <- gather(featureCounts_PA_sense_df_t_2_both, key="RNA_type", "abundance", -c(sample,strain,time))
featureCounts_PA_antisense_df_t_2_L$RNA_type2 <- with(featureCounts_PA_antisense_df_t_2_L, ifelse(RNA_type == "ncRNA", "ncRNA", "CDS"))
featureCounts_PA_sense_df_t_2_L$RNA_type2 <- with(featureCounts_PA_sense_df_t_2_L, ifelse(RNA_type == "ncRNA", "ncRNA", "CDS"))

scientific_10 <- function(x) {
    xout <- gsub("1e", "10^{", format(x),fixed=TRUE)
    xout <- gsub("{-0", "{-", xout,fixed=TRUE)
    xout <- gsub("{+", "{", xout,fixed=TRUE)
    xout <- gsub("{0", "{", xout,fixed=TRUE)
    xout <- paste(xout,"}",sep="")
    return(parse(text=xout))}

scale_y_log10nice <- function(name=NULL,...) {
    breaks10 <- c(100, 10000, 1000000)
    scale_y_log10(name, breaks=breaks10, 
                  labels=scientific_10(breaks10), limits=c(10, 1900000),...)}

quantify_antisense_RNA <- ggplot(featureCounts_PA_antisense_df_t_2_L, aes(x=RNA_type2, y=abundance)) +
  geom_violin() + geom_jitter(aes(colour=time), width = 0.05, alpha=0.4) + 
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange", color = "red") +
  stat_compare_means(label = "p.signif", size=6, label.x = 1.44, label.y = 5.9) +  
  theme_pubr(border=TRUE, base_size=11)  + xlab(" ") + 
  facet_wrap(~strain) + scale_y_log10nice() +
  theme(legend.position = "none", axis.text.x = element_text(size=10)) +
  scale_color_manual(values=c("black","blue"))

quantify_sense_RNA <- ggplot(featureCounts_PA_sense_df_t_2_L, aes(x=RNA_type2, y=abundance)) +
  geom_violin() + geom_jitter(aes(colour=time), width = 0.05, alpha=0.4) + 
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange", color = "red") +
  stat_compare_means(label = "p.signif", size=6, label.x = 1.44, label.y = 5.9) + 
  theme_pubr(border=TRUE, base_size=11) + xlab(" ") + 
  facet_wrap(~strain) + 
  theme(legend.position = "none", axis.text.x = element_text(size=10)) +
  scale_y_log10nice("Read count") + scale_color_manual(values=c("black","blue"))

merge_quantify <- ggarrange(quantify_sense_RNA, quantify_antisense_RNA, labels=c("B", "C"), widths = c(1,0.9))
quantify_sense_antisense <- ggarrange(amount_rna2, merge_quantify, nrow=2, heights = c(1,0.5))
# ggsave(quantify_sense_antisense, filename = "save_figures/Figure_02.jpeg", device="jpeg", dpi=800, units="cm", width = 21.1, height = 21.6)

featureCounts_PA_sense_df_t_2_L_NN2 <- subset(featureCounts_PA_sense_df_t_2_L, strain == "NN2")
featureCounts_PA_sense_df_t_2_L_SG17M <- subset(featureCounts_PA_sense_df_t_2_L, strain == "SG17M")
wilcox.test(featureCounts_PA_sense_df_t_2_L_NN2$abundance ~ featureCounts_PA_sense_df_t_2_L_NN2$RNA_type2)
rcompanion::wilcoxonR(x=featureCounts_PA_sense_df_t_2_L_NN2$abundance, g=featureCounts_PA_sense_df_t_2_L_NN2$RNA_type2, ci=TRUE)

wilcox.test(featureCounts_PA_sense_df_t_2_L_SG17M$abundance ~ featureCounts_PA_sense_df_t_2_L_SG17M$RNA_type2)
rcompanion::wilcoxonR(x=featureCounts_PA_sense_df_t_2_L_SG17M$abundance, g=featureCounts_PA_sense_df_t_2_L_SG17M$RNA_type2, ci=TRUE)

featureCounts_PA_antisense_df_t_2_L_NN2 <- subset(featureCounts_PA_antisense_df_t_2_L, strain == "NN2")
featureCounts_PA_antisense_df_t_2_L_SG17M <- subset(featureCounts_PA_antisense_df_t_2_L, strain == "SG17M")
wilcox.test(featureCounts_PA_antisense_df_t_2_L_NN2$abundance ~ featureCounts_PA_antisense_df_t_2_L_NN2$RNA_type2)
rcompanion::wilcoxonR(x=featureCounts_PA_antisense_df_t_2_L_NN2$abundance, g=featureCounts_PA_antisense_df_t_2_L_NN2$RNA_type2, ci=TRUE)

wilcox.test(featureCounts_PA_antisense_df_t_2_L_SG17M$abundance ~ featureCounts_PA_antisense_df_t_2_L_SG17M$RNA_type2)
rcompanion::wilcoxonR(x=featureCounts_PA_antisense_df_t_2_L_SG17M$abundance, g=featureCounts_PA_antisense_df_t_2_L_SG17M$RNA_type2, ci=TRUE)

# Within and between group variance
featureCounts_PA_antisense_df_t_SG17M$mixed <- c("SG17M_4h", "SG17M_4h", "SG17M_4h", "SG17M_8h", "SG17M_8h", "SG17M_8h")
featureCounts_PA_antisense_df_t_SG17M$strandType <- "Antisense RNA"

featureCounts_PA_sense_df_t_SG17M$mixed <- c("SG17M_4h", "SG17M_4h", "SG17M_4h", "SG17M_8h", "SG17M_8h", "SG17M_8h")
featureCounts_PA_sense_df_t_SG17M$strandType <- "Sense RNA"

featureCounts_PA_antisense_df_t_NN2$mixed <- c("NN2_4h", "NN2_4h", "NN2_4h", "NN2_8h", "NN2_8h", "NN2_8h")
featureCounts_PA_antisense_df_t_NN2$strandType <- "Antisense RNA"

featureCounts_PA_sense_df_t_NN2$mixed <- c("NN2_4h", "NN2_4h", "NN2_4h", "NN2_8h", "NN2_8h", "NN2_8h")
featureCounts_PA_sense_df_t_NN2$strandType <- "Sense RNA"

featureCounts_PA_sense_variance <- data.frame(
  rbind(featureCounts_PA_sense_df_t_NN2, featureCounts_PA_sense_df_t_SG17M))
featureCounts_PA_sense_variance$logCDS <- log2(featureCounts_PA_sense_variance$CDS)
featureCounts_PA_sense_variance$logCDS.undefined <- log2(featureCounts_PA_sense_variance$CDS.undefined)
featureCounts_PA_sense_variance$ncRNA.undefined <- log2(featureCounts_PA_sense_variance$ncRNA.undefined)
featureCounts_PA_sense_variance$tRNA <- log2(featureCounts_PA_sense_variance$tRNA)
featureCounts_PA_sense_variance$rRNA <- log2(featureCounts_PA_sense_variance$rRNA)

featureCounts_PA_sense_variance_log <- select(featureCounts_PA_sense_variance, c(logCDS, logCDS.undefined, 
                                                ncRNA.undefined, tRNA, rRNA, mixed, strandType))
hc_sense <- hclust(dist(featureCounts_PA_sense_variance_log))
ddata_sense <- dendro_data(hc_sense, type = "rectangle")
featureCounts_PA_antisense_variance <- data.frame(rbind(featureCounts_PA_antisense_df_t_NN2, featureCounts_PA_antisense_df_t_SG17M))

featureCounts_PA_antisense_variance$logCDS <- log2(featureCounts_PA_antisense_variance$CDS)
featureCounts_PA_antisense_variance$logCDS.undefined <- log2(featureCounts_PA_antisense_variance$CDS.undefined)
featureCounts_PA_antisense_variance$log10ncRNA.undefined <- log2(featureCounts_PA_antisense_variance$ncRNA.undefined)
featureCounts_PA_antisense_variance$log10tRNA <- log2(featureCounts_PA_antisense_variance$tRNA)

featureCounts_PA_antisense_variance_log <- select(featureCounts_PA_antisense_variance, 
                                              c(logCDS, logCDS.undefined, log10ncRNA.undefined,
                                                log10tRNA, mixed, strandType))

pca_antisense_1 <- prcomp(featureCounts_PA_antisense_variance_log[1:4], center = TRUE, scale. = FALSE)
eig_var_antisense <- fviz_eig(pca_antisense_1)
eig_var_antisense$data
vis <- featureCounts_PA_antisense_variance_log$mixed

pcaInd_antisense_1_ind <- fviz_pca_ind(pca_antisense_1)
pcaInd_antisense_1_ind$data$vis <- vis

# Results for Variables (sense)
res.var.antisense <- get_pca_var(pca_antisense_1)
res.var.antisense <- data.frame(res.var.antisense$contrib)  
res.var.antisense <- res.var.antisense[order(-res.var.antisense$Dim.1),]
res.var.antisense <- data.frame(res.var.antisense)

res.var.antisense.sub <- select(res.var.antisense, c("Dim.1", "Dim.2"))
res.var.antisense.sub$type <- rownames(res.var.antisense.sub)
res.var.antisense.sub$type <- str_replace(res.var.antisense.sub$type, "log", "")
res.var.antisense.sub$type <- str_replace(res.var.antisense.sub$type, "10", "")
res.var.antisense.sub$strand <- "Antisense"
res.var.antisense.sub$type2 <- with(
  res.var.antisense.sub, 
  ifelse(type == "ncRNA.undefined", "antisense to ncRNA (undefined)",
         ifelse(type == "CDS.undefined", "antisense to CDS (undefined)",
                ifelse(type == "CDS", "antisense to CDS",
                       ifelse(type == "tRNA", "antisense to tRNA", NA)))))
rRNA_manual <- c(0,0,"rRNA", "Antisense", "antisense to rRNA")
res.var.antisense.sub <- data.frame(rbind(
  res.var.antisense.sub, rRNA_manual))

pca_sense_1 <- prcomp(featureCounts_PA_sense_variance_log[1:5], center = TRUE, scale. = FALSE)
eig_var_sense <- fviz_eig(pca_sense_1)
eig_var_sense$data
vis <- featureCounts_PA_sense_variance_log$mixed

pcaInd_sense_1_ind <- fviz_pca_ind(pca_sense_1)
pcaInd_sense_1_ind$data$vis <- vis

# Results for Variables (sense)
res.var.sense <- get_pca_var(pca_sense_1)
res.var.sense <- data.frame(res.var.sense$contrib)        
res.var.sense <- res.var.sense[order(-res.var.sense$Dim.1),]
res.var.sense <- data.frame(res.var.sense)

res.var.sense.sub <- select(res.var.sense, c("Dim.1", "Dim.2"))
res.var.sense.sub$type <- rownames(res.var.sense.sub)
res.var.sense.sub$type <- str_replace(res.var.sense.sub$type, "log", "")
res.var.sense.sub$type <- str_replace(res.var.sense.sub$type, "10", "")
res.var.sense.sub$strand <- " Sense "
res.var.sense.sub$type2 <- with(
  res.var.sense.sub, 
  ifelse(type == "ncRNA.undefined", " ncRNA (undefined)",
         ifelse(type == "CDS.undefined", " CDS (undefined)",
                ifelse(type == "CDS", " CDS",
                       ifelse(type == "tRNA", " tRNA", " rRNA")))))

res.var.both <- data.frame(rbind(res.var.antisense.sub, res.var.sense.sub))
rownames(res.var.both) <- NULL

res.var.both.L <- gather(res.var.both, key="Dim", value="PCA", -c("type", "strand", "type2"))
variance_antisense_pca_plot <- ggplot(pcaInd_antisense_1_ind$data) + 
  geom_text_repel(aes(x=x, y=y, label=name, color=vis), size= 3, angle=0) + 
  xlab("Dim1 (92.4 %)") + ylab("Dim2 (6.7 %)") +
  theme_pubr(border=TRUE, legend="bottom", base_size=11) + ylim(-2.3, 2.3) + xlim(-7, 8) +
  scale_color_manual(values=c("deeppink4", "violet", "navyblue", "lightseagreen")) +
  theme(legend.position = "none") + 
  geom_label(aes(x=6, y = 2, label = "Antisense"), fill="white")

variance_sense_pca_plot <- ggplot(pcaInd_sense_1_ind$data) + 
  geom_text_repel(aes(x=x, y=y, label=name, color=vis), size= 3, angle=0) +
  xlab("Dim1 (77.3 %)") + ylab("Dim2 (13.7 %)") +
  theme_pubr(border=TRUE, legend="bottom", base_size=11) +
  xlim(-7, 8) +  ylim(-2.3, 2.3) +
  scale_color_manual(values=c("deeppink4", "violet", "navyblue", "lightseagreen")) +
  theme(legend.position = "none") + 
  geom_label(aes(x=6, y = 2, label = "Sense"), fill="white")

res.var.both.L$PCA <- as.numeric(as.character(res.var.both.L$PCA))

variance_both <- ggplot(res.var.both.L) + geom_col(aes(x=Dim, y=PCA, fill=type2)) + 
  facet_grid(~strand) + theme_pubr(base_size=11, border=TRUE) + 
  theme(legend.title = element_blank(), legend.position = "top", axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
        legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
        axis.text.y = element_text(size=10)) + 
  xlab(" ") + ylab("Relative abundance") + coord_flip() + guides(fill=guide_legend(ncol=2)) +
  scale_fill_manual(values=c("skyblue4", "darkolivegreen4","gray5",  "tan2", "brown3", 
                             "skyblue3", "darkolivegreen3", "gray60", "tan1", "brown1"))


ab <- ggarrange(variance_sense_pca_plot, variance_antisense_pca_plot, labels=c("A", "B"))
ab2 <- ggarrange(ab, variance_both, nrow = 2, labels=c("A", "C"), heights=c(1,0.5))

# ggsave(ab2, filename = "save_figures/Supplementary_Figure_05.jpeg", device="jpeg", dpi=600, units="cm", width=21.1, height=20.3)

import_gtf_NN2$gene_name3 <- ifelse(import_gtf_NN2$gene_name2 == "", import_gtf_NN2$featureType, import_gtf_NN2$gene_name2)
my_list_2_NN2 <- import_gtf_NN2$gene_name3
names(my_list_2_NN2) <- import_gtf_NN2$transcript_id

mylist_antisense_x_NN2 <- c()
for(i in rownames(tRNA_antisense_NN2)){
    newelem <- ifelse(i %in% names(my_list_2_NN2), my_list_2_NN2[[i]], i)
    mylist_antisense_x_NN2 <- append(mylist_antisense_x_NN2, newelem)}
tRNA_antisense_NN2$name2 <- mylist_antisense_x_NN2
tRNA_antisense_NN2 <- tRNA_antisense_NN2[grepl( "tRNA" , tRNA_antisense_NN2$name2 ), ]
tRNA_antisense_NN2 <- tRNA_antisense_NN2 %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_antisense_NN2 <- data.frame(tRNA_antisense_NN2)
rownames(tRNA_antisense_NN2) <- tRNA_antisense_NN2$name2
tRNA_antisense_NN2$name2 <- NULL
remove_empty <- c(" tRNA-SerGGA")
tRNA_antisense_NN2 <- tRNA_antisense_NN2[(!rownames(tRNA_antisense_NN2) %in% remove_empty),]

mylist_sense_x_NN2 <- c()
for(i in rownames(tRNA_sense_NN2)){
    newelem <- ifelse(i %in% names(my_list_2_NN2), my_list_2_NN2[[i]], i)
    mylist_sense_x_NN2 <- append(mylist_sense_x_NN2, newelem)}
tRNA_sense_NN2$name2 <- mylist_sense_x_NN2
tRNA_sense_NN2 <- tRNA_sense_NN2[grepl( "tRNA" , tRNA_sense_NN2$name2 ), ]
tRNA_sense_NN2 <- tRNA_sense_NN2 %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_sense_NN2 <- data.frame(tRNA_sense_NN2)
rownames(tRNA_sense_NN2) <- tRNA_sense_NN2$name2
tRNA_sense_NN2$name2 <- NULL
tRNA_sense_NN2 <- tRNA_sense_NN2[(!rownames(tRNA_sense_NN2) %in% remove_empty),]

import_gtf_SG17M$gene_name3 <- ifelse(import_gtf_SG17M$gene_name2 == "", import_gtf_SG17M$featureType, import_gtf_SG17M$gene_name2)
my_list_2_SG17M <- import_gtf_SG17M$gene_name3
names(my_list_2_SG17M) <- import_gtf_SG17M$transcript_id

mylist_antisense_x_SG17M <- c()
for(i in rownames(tRNA_antisense_SG17M)){
    newelem <- ifelse(i %in% names(my_list_2_SG17M), my_list_2_SG17M[[i]], i)
    mylist_antisense_x_SG17M <- append(mylist_antisense_x_SG17M, newelem)}
tRNA_antisense_SG17M$name2 <- mylist_antisense_x_SG17M
tRNA_antisense_SG17M <- tRNA_antisense_SG17M[grepl( "tRNA" , tRNA_antisense_SG17M$name2 ), ]
tRNA_antisense_SG17M <- tRNA_antisense_SG17M %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_antisense_SG17M <- data.frame(tRNA_antisense_SG17M)
rownames(tRNA_antisense_SG17M) <- tRNA_antisense_SG17M$name2
tRNA_antisense_SG17M$name2 <- NULL

mylist_sense_x_SG17M <- c()
for(i in rownames(tRNA_sense_SG17M)){
    newelem <- ifelse(i %in% names(my_list_2_SG17M), my_list_2_SG17M[[i]], i)
    mylist_sense_x_SG17M <- append(mylist_sense_x_SG17M, newelem)}
tRNA_sense_SG17M$name2 <- mylist_sense_x_SG17M
tRNA_sense_SG17M <- tRNA_sense_SG17M[grepl( "tRNA" , tRNA_sense_SG17M$name2 ), ]
tRNA_sense_SG17M <- tRNA_sense_SG17M %>% group_by(name2) %>% summarise_all(funs(sum))
tRNA_sense_SG17M <- data.frame(tRNA_sense_SG17M)
rownames(tRNA_sense_SG17M) <- tRNA_sense_SG17M$name2
tRNA_sense_SG17M$name2 <- NULL

tRNA_sense <- data.frame(cbind(tRNA_sense_NN2, tRNA_sense_SG17M))
tRNA_antisense <- data.frame(cbind(tRNA_antisense_NN2, tRNA_antisense_SG17M))

#A grouping factor can be added at the same time:
group <- c("NN2_4h", "NN2_4h", "NN2_4h", "NN2_8h", "NN2_8h", "NN2_8h",
           "SG17M_4h", "SG17M_4h", "SG17M_4h", "SG17M_8h", "SG17M_8h", "SG17M_8h")
y_antisense <- DGEList(counts=tRNA_antisense, group=group)
y_sense <- DGEList(counts=tRNA_sense, group=group)

# We filter out lowly expressed genes using the following commands:
keep_antisense <- filterByExpr(y_antisense)
keep_sense <- filterByExpr(y_sense)
y_antisense_unfiltered <- y_antisense
y_sense_unfiltered <- y_sense
y_antisense <- y_antisense[keep_antisense, keep.lib.sizes=FALSE]
y_sense <- y_sense[keep_sense, keep.lib.sizes=FALSE]

# TMM normalisation
y_antisense <- calcNormFactors(y_antisense, method="TMM")
logcpm_antisense <- cpm(y_antisense, log=TRUE)
logcpm_antisense_df <- data.frame(logcpm_antisense)
y_sense <- calcNormFactors(y_sense, method="TMM")
logcpm_sense <- cpm(y_sense, log=TRUE)
logcpm_sense_df <- data.frame(logcpm_sense)

range <- max(abs(logcpm_sense_df))
rownames(logcpm_antisense_df) <- paste("(-) ", rownames(logcpm_antisense_df), sep="")
rownames(logcpm_sense_df) <- paste("(+) ", rownames(logcpm_sense_df), sep="")

logcpm_heatmap_df <- data.frame(rbind(logcpm_antisense_df, logcpm_sense_df))
rownames(logcpm_heatmap_df) <- c("(-) Lys-anticodon TTT", 
                                 "(-) Lys-anticodon TTT-pseudo", 
                                 "(+) Ala-anticodon GGC",
                                 "(+) Ala-anticodon TGC",
                                 "(+) Arg-anticodon ACG", 
                                 "(+) Arg-anticodon TCT", 
                                 "(+) Asn-anticodon GTT", 
                                 "(+) Asp-anticodon GTC", 
                                 "(+) fMet-anticodon CAT", 
                                 "(+) Gln-anticodon TTG", 
                                 "(+) Glu-anticodon TTC", 
                                 "(+) Gly-anticodon GCC", 
                                 "(+) Gly-anticodon TCC", 
                                 "(+) His-anticodon GTG",
                                 "(+) Ile-anticodon GAT", 
                                 "(+) Leu-anticodon CAG", 
                                 "(+) Leu-anticodon GAG", 
                                 "(+) Lys-anticodon TTT",
                                 "(+) Pro-anticodon GGG", 
                                 "(+) Pro-anticodon TGG", 
                                 "(+) Ser-anticodon CGA-pseudo", 
                                 "(+) Ser-anticodon GCT",
                                 "(+) Thr-anticodon CGT", 
                                 "(+) Thr-anticodon GGT", 
                                 "(+) Trp-anticodon CCA", 
                                 "(+) Tyr-anticodon GTA", 
                                 "(+) tmRNA-SsrA", 
                                 "(+) Val-anticodon TAC")

tRNA_heatmap <- pheatmap(logcpm_heatmap_df, cutree_rows  = 5, cluster_cols = FALSE, scale = "none", 
                         color = colorRampPalette(c("beige", "navy", "firebrick3"))(10),
                         clustering_distance_cols = "euclidean", angle_col = 45, fontsize = 10)

tRNA_heatmap_gg <- as.ggplot(tRNA_heatmap, scale=1)
# ggsave(tRNA_heatmap_gg, filename = "save_figures/Figure_03.jpeg", dpi=800, units="cm", width=18.2, height=18.9, bg="white")

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

SG17M_codon_usage <- read_delim("SG17M_codon_usage.csv", ";", escape_double = FALSE, trim_ws = TRUE)
SG17M_codon_usage <- SG17M_codon_usage[,1:5]
SG17M_codon_usage$Anticodon1 <- strReverse(SG17M_codon_usage$Codon)

sg17m_comp <- c()
for (items in SG17M_codon_usage$Anticodon1) { a = c2s(comp(s2c(items), forceToLower = FALSE))
   sg17m_comp <- append(sg17m_comp, a)}

SG17M_codon_usage$Anticodon2 <- sg17m_comp
SG17M_codon_usage$Anticodon3 <- paste("tRNA-",SG17M_codon_usage$AmAcid,SG17M_codon_usage$Anticodon2)
SG17M_codon_usage$Anticodon3 <- str_replace_all(SG17M_codon_usage$Anticodon3, " ", "")
SG17M_codon_usage <- SG17M_codon_usage[1:61,]

# TMM normalisation
y_sense_unfiltered_df <- data.frame(y_sense_unfiltered$counts)
y_sense_unfiltered_df$Anticodon3 <- rownames(y_sense_unfiltered_df)
y_sense_unfiltered_df$Anticodon3 <- str_replace_all(y_sense_unfiltered_df$Anticodon3, "-pseudo", "")
y_sense_unfiltered_df$Anticodon3 <- str_replace_all(y_sense_unfiltered_df$Anticodon3, " ", "")
y_antisense_unfiltered_df <- data.frame(y_antisense_unfiltered$counts)
y_antisense_unfiltered_df$Anticodon3 <- rownames(y_antisense_unfiltered_df)
y_antisense_unfiltered_df$Anticodon3 <- str_replace_all(y_antisense_unfiltered_df$Anticodon3, "-pseudo", "")
y_antisense_unfiltered_df$Anticodon3 <- str_replace_all(y_antisense_unfiltered_df$Anticodon3, " ", "")

missing_ab1 <- setdiff(y_sense_unfiltered_df$Anticodon3, NN2_codon_usage$Anticodon3)
missing_ab2 <- setdiff(y_sense_unfiltered_df$Anticodon3, SG17M_codon_usage$Anticodon3)
missing_df1 <- setdiff(NN2_codon_usage$Anticodon3, y_sense_unfiltered_df$Anticodon3)
missing_df2 <- rep(0,25)
missing_df3 <- data.frame(cbind(missing_df2, missing_df2, missing_df2, missing_df2, 
                                missing_df2, missing_df2, missing_df2, missing_df2, 
                                missing_df2, missing_df2, missing_df2, missing_df2, missing_df1)) 
colnames(missing_df3) <- colnames(y_sense_unfiltered_df)
y_sense_unfiltered_df2 <- data.frame(rbind(y_sense_unfiltered_df, missing_df3))
y_sense_unfiltered_df2 <- y_sense_unfiltered_df2[!y_sense_unfiltered_df2$Anticodon3 %in% missing_ab1,]
y_sense_unfiltered_df2$NN2_4h_BR1 <- as.numeric(as.character(y_sense_unfiltered_df2$NN2_4h_BR1)) 
y_sense_unfiltered_df2$NN2_4h_BR2 <- as.numeric(as.character(y_sense_unfiltered_df2$NN2_4h_BR2))
y_sense_unfiltered_df2$NN2_4h_BR3 <- as.numeric(as.character(y_sense_unfiltered_df2$NN2_4h_BR3))
y_sense_unfiltered_df2$NN2_8h_BR1 <- as.numeric(as.character(y_sense_unfiltered_df2$NN2_8h_BR1)) 
y_sense_unfiltered_df2$NN2_8h_BR2 <- as.numeric(as.character(y_sense_unfiltered_df2$NN2_8h_BR2))
y_sense_unfiltered_df2$NN2_8h_BR3 <- as.numeric(as.character(y_sense_unfiltered_df2$NN2_8h_BR3))

y_sense_unfiltered_df2$SG17M_4h_BR1 <- as.numeric(as.character(y_sense_unfiltered_df2$SG17M_4h_BR1)) 
y_sense_unfiltered_df2$SG17M_4h_BR2 <- as.numeric(as.character(y_sense_unfiltered_df2$SG17M_4h_BR2))
y_sense_unfiltered_df2$SG17M_4h_BR3 <- as.numeric(as.character(y_sense_unfiltered_df2$SG17M_4h_BR3))
y_sense_unfiltered_df2$SG17M_8h_BR1 <- as.numeric(as.character(y_sense_unfiltered_df2$SG17M_8h_BR1)) 
y_sense_unfiltered_df2$SG17M_8h_BR2 <- as.numeric(as.character(y_sense_unfiltered_df2$SG17M_8h_BR2))
y_sense_unfiltered_df2$SG17M_8h_BR3 <- as.numeric(as.character(y_sense_unfiltered_df2$SG17M_8h_BR3))
y_sense_unfiltered_df2
y_sense_unfiltered_df2_NN2 <- y_sense_unfiltered_df2
y_sense_unfiltered_df2_NN2_1 <- plyr::ddply(y_sense_unfiltered_df2_NN2, "Anticodon3", numcolwise(sum))
rownames(y_sense_unfiltered_df2_NN2_1) <- y_sense_unfiltered_df2_NN2_1$Anticodon3
y_sense_unfiltered_df2_NN2_1$Anticodon3 <- NULL
y_sense_unfiltered_df2_NN2_1 <- y_sense_unfiltered_df2_NN2_1[order(rownames(y_sense_unfiltered_df2_NN2_1)),]

NN2_codon_usage1 <- NN2_codon_usage[order(NN2_codon_usage$Anticodon3),]
NN2_codon_usage1 <- data.frame(NN2_codon_usage1)
SG17M_codon_usage1 <- SG17M_codon_usage[order(SG17M_codon_usage$Anticodon3),]
SG17M_codon_usage1 <- data.frame(SG17M_codon_usage1)

y_sense_unfiltered_df2_NN2_2 <- prop.table(as.matrix(y_sense_unfiltered_df2_NN2_1), margin = 2)
y_sense_unfiltered_df2_NN2_2 <- data.frame(y_sense_unfiltered_df2_NN2_2)
y_sense_unfiltered_df2_NN2_2$ref_NN2 <- (NN2_codon_usage1$div1000) / 10
y_sense_unfiltered_df2_NN2_2$ref_SG17M <- (SG17M_codon_usage1$div1000) / 10

y_sense_unfiltered_df2_NN2_1$ref_NN2 <- NN2_codon_usage1$div1000
y_sense_unfiltered_df2_NN2_1$ref_SG17M <- SG17M_codon_usage1$div1000

y_sense_unfiltered_df2_NN2_2$NN2_4h_BR1_rank <- rank(y_sense_unfiltered_df2_NN2_2$NN2_4h_BR1)
y_sense_unfiltered_df2_NN2_2$NN2_4h_BR2_rank <- rank(y_sense_unfiltered_df2_NN2_2$NN2_4h_BR2)
y_sense_unfiltered_df2_NN2_2$NN2_4h_BR3_rank <- rank(y_sense_unfiltered_df2_NN2_2$NN2_4h_BR3)

y_sense_unfiltered_df2_NN2_2$NN2_8h_BR1_rank <- rank(y_sense_unfiltered_df2_NN2_2$NN2_8h_BR1)
y_sense_unfiltered_df2_NN2_2$NN2_8h_BR2_rank <- rank(y_sense_unfiltered_df2_NN2_2$NN2_8h_BR2)
y_sense_unfiltered_df2_NN2_2$NN2_8h_BR3_rank <- rank(y_sense_unfiltered_df2_NN2_2$NN2_8h_BR3)

y_sense_unfiltered_df2_NN2_2$SG17M_4h_BR1_rank <- rank(y_sense_unfiltered_df2_NN2_2$SG17M_4h_BR1)
y_sense_unfiltered_df2_NN2_2$SG17M_4h_BR2_rank <- rank(y_sense_unfiltered_df2_NN2_2$SG17M_4h_BR2)
y_sense_unfiltered_df2_NN2_2$SG17M_4h_BR3_rank <- rank(y_sense_unfiltered_df2_NN2_2$SG17M_4h_BR3)

y_sense_unfiltered_df2_NN2_2$SG17M_8h_BR1_rank <- rank(y_sense_unfiltered_df2_NN2_2$SG17M_8h_BR1)
y_sense_unfiltered_df2_NN2_2$SG17M_8h_BR2_rank <- rank(y_sense_unfiltered_df2_NN2_2$SG17M_8h_BR2)
y_sense_unfiltered_df2_NN2_2$SG17M_8h_BR3_rank <- rank(y_sense_unfiltered_df2_NN2_2$SG17M_8h_BR3)
y_sense_unfiltered_df2_NN2_2$ref_NN2_rank <- rank(y_sense_unfiltered_df2_NN2_2$ref_NN2)
y_sense_unfiltered_df2_NN2_2$ref_SG17M_rank <- rank(y_sense_unfiltered_df2_NN2_2$ref_SG17M)
y_sense_unfiltered_df2_NN2_2_rank <- y_sense_unfiltered_df2_NN2_2[,c(15:ncol(y_sense_unfiltered_df2_NN2_2))]
colnames(y_sense_unfiltered_df2_NN2_2_rank) <- str_replace_all(colnames(y_sense_unfiltered_df2_NN2_2_rank), "_rank", "")

rownames(y_sense_unfiltered_df2_NN2_2_rank) <- c("Ala-anticodon AGC", 
                                                 "Ala-anticodon CGC", 
                                                 "Ala-anticodon GGC", "Ala-anticodon TGC", 
                                                 "Arg-anticodon ACG",
                                                 "Arg-anticodon CCG", "Arg-anticodon CCT", 
                                                 "Arg-anticodon GCG", "Arg-anticodon TCG",
                                                 "Arg-anticodon TCT", "Asn-anticodon ATT", "Asn-anticodon GTT", 
                                                 "Asp-anticodon ATC", "Asp-anticodon GTC", 
                                                 "Cys-anticodon ACA", "Cys-anticodon GCA", 
                                                 "Gln-anticodon CTG", "Gln-anticodon TTG",
                                                 "Glu-anticodon CTC", "Glu-anticodon TTC", "Gly-anticodon ACC", 
                                                 "Gly-anticodon CCC", "Gly-anticodon GCC", 
                                                 "Gly-anticodon TCC", "His-anticodon ATG", 
                                                 "His-anticodon GTG", "Ile-anticodon AAT",
                                                 "Ile-anticodon GAT", "Ile-anticodon TAT", "Leu-anticodon AAG", 
                                                 "Leu-anticodon CAA", "Leu-anticodon CAG", 
                                                 "Leu-anticodon GAG", "Leu-anticodon TAA", 
                                                 "Leu-anticodon TAG", "Lys-anticodon CTT",
                                                 "Lys-anticodon TTT", "Met-anticodon CAT", 
                                                 "Phe-anticodon AAA", "Phe-anticodon GAA",
                                                 "Pro-anticodon AGG", "Pro-anticodon CGG", "Pro-anticodon GGG", 
                                                 "Pro-anticodon TGG", "Ser-anticodon ACT",
                                                 "Ser-anticodon AGA", "Ser-anticodon CGA", "Ser-anticodon GCT", 
                                                 "Ser-anticodon GGA", "Ser-anticodon TGA", 
                                                 "Thr-anticodon AGT", "Thr-anticodon CGT", 
                                                 "Thr-anticodon GGT", "Thr-anticodon TGT",
                                                 "Trp-anticodon CCA", "Tyr-anticodon ATA", 
                                                 "Tyr-anticodon GTA", "Val-anticodon AAC", 
                                                 "Val-anticodon CAC", "Val-anticodon GAC", "Val-anticodon TAC")

codonUsage_trna_heatmap <- pheatmap(y_sense_unfiltered_df2_NN2_2_rank, scale = "none", cutree_cols = 2,
           color = colorRampPalette(c("grey", "beige",  "firebrick3"))(10), 
         cluster_rows = TRUE, cutree_rows = 8 ,cellheight = 9, cellwidth = 16, 
         clustering_distance_cols = "euclidean", treeheight_col = 3, treeheight_row = 0,
         angle_col = 45, fontsize = 8)
codonUsage_tRNA_heatmap_gg <- as.ggplot(codonUsage_trna_heatmap, scale=1)

codonUsage_trna_heatmap_alphabet <- pheatmap(y_sense_unfiltered_df2_NN2_2_rank, scale = "none", cutree_cols = 2,
                                             color = colorRampPalette(c("grey",  "beige", "firebrick3"))(10), 
         cluster_rows = FALSE, cutree_rows = 8 ,cellheight = 9, cellwidth = 16, 
         clustering_distance_cols = "euclidean", treeheight_col = 3, treeheight_row = 0, angle_col = 45, fontsize = 8)

codonUsage_tRNA_heatmap_gg <- as.ggplot(codonUsage_trna_heatmap, scale=1)
codonUsage_tRNA_heatmap_alphabet_gg <- as.ggplot(codonUsage_trna_heatmap_alphabet, scale=1)

#ggsave(codonUsage_tRNA_heatmap_gg, filename="save_figures/Figure_04.jpeg", 
  #     device="jpeg", dpi=800, units="cm", width=22, height=24, bg="white") 

#ggsave(codonUsage_tRNA_heatmap_alphabet_gg,
    #   filename="save_figures/Supplementary_Figure_07.jpeg", 
     #  device="jpeg", dpi=800, units="cm", width=22, height=24, bg="white") 

corPlot1 <- cor(y_sense_unfiltered_df2_NN2_1, method = "spearman")
corPlot1pval <- rcorr(as.matrix(y_sense_unfiltered_df2_NN2_1), type = "spearman")
corPlot1pval_df <- data.frame(corPlot1pval$P)
corPlot1pval_df$a <- rownames(corPlot1pval_df)
corPlot1pval_df_L <- gather(corPlot1pval_df, key="Sample", value="scale", -a)
corPlot1pval_df_L_noNA <- na.omit(corPlot1pval_df_L)

corPlot1_df <- data.frame(corPlot1)
corPlot1_df$a <- rownames(corPlot1_df)
corPlot1_df_L <- gather(corPlot1_df, key="Sample", value="scale", -a)
corPlot1_df_L <- subset(corPlot1_df_L, scale < 1.0000000)

pval_cor <- ggplot(corPlot1pval_df_L_noNA) +
  geom_jitter(aes(x=a, y=scale, colour=Sample), size=3, alpha=0.8, width = 0.2) +
  scale_colour_manual(values = c("yellow", "yellow", "yellow", 
                                 "orange", "orange", "orange", "red", "pink",
                                 "blue", "blue", "blue", 
                                 "black", "black", "black")) + 
  theme_pubr(border=TRUE, base_size=8) +
  coord_flip() + ylab("Spearman's, p-value") + xlab("") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme(legend.title=element_blank())

coff_cor <-ggplot(corPlot1_df_L) +
  geom_jitter(aes(x=a, y=scale, colour=Sample), size=3, alpha=0.8, width = 0.2) +
  scale_colour_manual(values = c("yellow", "yellow", "yellow", 
                                 "orange", "orange", "orange", "red", "pink",
                                 "blue", "blue", "blue", 
                                 "black", "black", "black")) + 
  theme_pubr(border=TRUE, base_size=8, legend = "bottom") +
  coord_flip() + ylab("Spearman's correlation coefficient") + xlab("") +
  geom_hline(yintercept = 0.6, linetype = "dashed") +
  geom_hline(yintercept = 0.2, linetype = "dashed") + ylim(0,1) +
  geom_label(aes(x=14, y=0.05), label="low", size=3.2) +
  geom_label(aes(x=14, y=0.4), label="medium", size=3.2) +
  geom_label(aes(x=14, y=0.8), label="high", size=3.2)+
  theme(legend.title=element_blank(), axis.text.y = element_blank())

data_matrix_1 <- cor(y_sense_unfiltered_df2_NN2_1, method = "spearman")

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
  p.mat}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(data_matrix_1)
#jpeg(filename = "save_figures/Supplementary_figure_08.jpeg", units="cm", res=800, width=24, height=24)
corrplot(data_matrix_1, method="square", type="lower", order="hclust", number.cex = 0.7, 
         addCoef.col = "black", tl.col="black", tl.cex = 0.7, tl.srt = 40, diag=FALSE)
# dev.off()
