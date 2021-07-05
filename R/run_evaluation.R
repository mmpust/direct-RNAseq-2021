# title run_evaluation
# author: "Marie Pust"
# date: "05 07 2021"

# clean environment
rm(list = ls())

# set working directory
setwd("/MinionExperiments202009/d_Ranalysis")

# load libraries
library(Rsamtools)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyr)
library(readr)
library(rcompanion)

# load input files
input_bam <- list.files(
  path = "b_minimap2/",
  pattern = ".bam", full.names = TRUE)

UserDefinedAnnotationRef <- "NN2_ENO_SG17M_curated.gtf"
Ref_genome <- "NN2_SG17M_CS.fna"

ont_NN2_4h_BR1 <-  scanBam("NN2_CS_4h_BR1.bam")
ont_NN2_4h_BR1 <- data.frame(ont_NN2_4h_BR1)
ont_NN2_4h_BR1$isolate <- "NN2"
ont_NN2_4h_BR1$time <- "4h"
ont_NN2_4h_BR1$mix <- "NN2_4h"
ont_NN2_4h_BR1$replicate <- "BR1"
ont_NN2_4h_BR1$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_NN2_4h_BR1 <- subset(ont_NN2_4h_BR1, flag == 0 | flag == 16)
 

ont_NN2_4h_BR2 <-  scanBam("NN2_CS_4h_BR2.bam")
ont_NN2_4h_BR2 <- data.frame(ont_NN2_4h_BR2)
ont_NN2_4h_BR2$isolate <- "NN2"
ont_NN2_4h_BR2$time <- "4h"
ont_NN2_4h_BR2$mix <- "NN2_4h"
ont_NN2_4h_BR2$replicate <- "BR2"
ont_NN2_4h_BR2$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_NN2_4h_BR2 <- subset(ont_NN2_4h_BR2, flag == 0 | flag == 16)


ont_NN2_4h_BR3 <-  scanBam("NN2_CS_4h_BR3.bam")
ont_NN2_4h_BR3 <- data.frame(ont_NN2_4h_BR3)
ont_NN2_4h_BR3$isolate <- "NN2"
ont_NN2_4h_BR3$time <- "4h"
ont_NN2_4h_BR3$mix <- "NN2_4h"
ont_NN2_4h_BR3$replicate <- "BR3"
ont_NN2_4h_BR3$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_NN2_4h_BR3 <- subset(ont_NN2_4h_BR3, flag == 0 | flag == 16)


ont_NN2_8h_BR1 <-  scanBam("NN2_CS_8h_BR1.bam")
ont_NN2_8h_BR1 <- data.frame(ont_NN2_8h_BR1)
ont_NN2_8h_BR1$isolate <- "NN2"
ont_NN2_8h_BR1$time <- "8h"
ont_NN2_8h_BR1$mix <- "NN2_8h"
ont_NN2_8h_BR1$replicate <- "BR1"
ont_NN2_8h_BR1$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_NN2_8h_BR1 <- subset(ont_NN2_8h_BR1, flag == 0 | flag == 16)


ont_NN2_8h_BR2 <-  scanBam("NN2_CS_8h_BR2.bam")
ont_NN2_8h_BR2 <- data.frame(ont_NN2_8h_BR2)
ont_NN2_8h_BR2$isolate <- "NN2"
ont_NN2_8h_BR2$time <- "8h"
ont_NN2_8h_BR2$mix <- "NN2_8h"
ont_NN2_8h_BR2$replicate <- "BR2"
ont_NN2_8h_BR2$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_NN2_8h_BR2 <- subset(ont_NN2_8h_BR2, flag == 0 | flag == 16)


ont_NN2_8h_BR3 <-  scanBam("NN2_CS_8h_BR3.bam")
ont_NN2_8h_BR3 <- data.frame(ont_NN2_8h_BR3)
ont_NN2_8h_BR3$isolate <- "NN2"
ont_NN2_8h_BR3$time <- "8h"
ont_NN2_8h_BR3$mix <- "NN2_8h"
ont_NN2_8h_BR3$replicate <- "BR3"
ont_NN2_8h_BR3$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_NN2_8h_BR3 <- subset(ont_NN2_8h_BR3, flag == 0 | flag == 16)


ont_SG17M_4h_BR1 <-  scanBam("SG17M_CS_4h_BR1.bam")
ont_SG17M_4h_BR1 <- data.frame(ont_SG17M_4h_BR1)
ont_SG17M_4h_BR1$isolate <- "SG17M"
ont_SG17M_4h_BR1$time <- "4h"
ont_SG17M_4h_BR1$mix <- "SG17M_4h"
ont_SG17M_4h_BR1$replicate <- "BR1"
ont_SG17M_4h_BR1$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_SG17M_4h_BR1 <- subset(ont_SG17M_4h_BR1, flag == 0 | flag == 16)


ont_SG17M_4h_BR2 <-  scanBam("SG17M_CS_4h_BR2.bam")
ont_SG17M_4h_BR2 <- data.frame(ont_SG17M_4h_BR2)
ont_SG17M_4h_BR2$isolate <- "SG17M"
ont_SG17M_4h_BR2$time <- "4h"
ont_SG17M_4h_BR2$mix <- "SG17M_4h"
ont_SG17M_4h_BR2$replicate <- "BR2"
ont_SG17M_4h_BR2$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_SG17M_4h_BR2 <- subset(ont_SG17M_4h_BR2, flag == 0 | flag == 16)


ont_SG17M_4h_BR3 <-  scanBam("SG17M_CS_4h_BR3.bam")
ont_SG17M_4h_BR3 <- data.frame(ont_SG17M_4h_BR3)
ont_SG17M_4h_BR3$isolate <- "SG17M"
ont_SG17M_4h_BR3$time <- "4h"
ont_SG17M_4h_BR3$mix <- "SG17M_4h"
ont_SG17M_4h_BR3$replicate <- "BR3"
ont_SG17M_4h_BR3$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_SG17M_4h_BR3 <- subset(ont_SG17M_4h_BR3, flag == 0 | flag == 16)


ont_SG17M_8h_BR1 <-  scanBam("SG17M_CS_8h_BR1.bam")
ont_SG17M_8h_BR1 <- data.frame(ont_SG17M_8h_BR1)
ont_SG17M_8h_BR1$isolate <- "SG17M"
ont_SG17M_8h_BR1$time <- "8h"
ont_SG17M_8h_BR1$mix <- "SG17M_8h"
ont_SG17M_8h_BR1$replicate <- "BR1"
ont_SG17M_8h_BR1$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_SG17M_8h_BR1 <- subset(ont_SG17M_8h_BR1, flag == 0 | flag == 16)


ont_SG17M_8h_BR2 <-  scanBam("SG17M_CS_8h_BR2.bam")
ont_SG17M_8h_BR2 <- data.frame(ont_SG17M_8h_BR2)
ont_SG17M_8h_BR2$isolate <- "SG17M"
ont_SG17M_8h_BR2$time <- "8h"
ont_SG17M_8h_BR2$mix <- "SG17M_8h"
ont_SG17M_8h_BR2$replicate <- "BR2"
ont_SG17M_8h_BR2$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_SG17M_8h_BR2 <- subset(ont_SG17M_8h_BR2, flag == 0 | flag == 16)


ont_SG17M_8h_BR3 <-  scanBam("SG17M_CS_8h_BR3.bam")
ont_SG17M_8h_BR3 <- data.frame(ont_SG17M_8h_BR3)
ont_SG17M_8h_BR3$isolate <- "SG17M"
ont_SG17M_8h_BR3$time <- "8h"
ont_SG17M_8h_BR3$mix <- "SG17M_8h"
ont_SG17M_8h_BR3$replicate <- "BR3"
ont_SG17M_8h_BR3$platform <- "Nanopore"
# filter for primary mapping reads only, flags 0 and 16
#ont_SG17M_8h_BR3 <- subset(ont_SG17M_8h_BR3, flag == 0 | flag == 16)


# start analysis
ont_NN2_4h <- data.frame(rbind(ont_NN2_4h_BR1, ont_NN2_4h_BR2, ont_NN2_4h_BR3))
ont_NN2_4h_noNA <- ont_NN2_4h[!is.na(ont_NN2_4h$qwidth), ]
ont_NN2_4h_noNA$rname2 <- ifelse(ont_NN2_4h_noNA$rname == "RNA_CS_ENO2", "ENO2","PA")
ont_NN2_4h_noNA$rname2 <- as.factor(as.character(ont_NN2_4h_noNA$rname2))

ont_SG17M_4h <- data.frame(rbind(ont_SG17M_4h_BR1, ont_SG17M_4h_BR2, ont_SG17M_4h_BR3))
ont_SG17M_4h_noNA <- ont_SG17M_4h[!is.na(ont_SG17M_4h$qwidth), ]
ont_SG17M_4h_noNA$rname2 <- ifelse(ont_SG17M_4h_noNA$rname == "RNA_CS_ENO2", "ENO2", "PA")
ont_SG17M_4h_noNA$rname2 <- as.factor(as.character(ont_SG17M_4h_noNA$rname2))

ont_SG17M_NN2_4h_noNA <- data.frame(rbind(ont_NN2_4h_noNA, ont_SG17M_4h_noNA))

ont_NN2_8h <- data.frame(rbind(ont_NN2_8h_BR1, ont_NN2_8h_BR2, ont_NN2_8h_BR3))
ont_NN2_8h_noNA <- ont_NN2_8h[!is.na(ont_NN2_8h$qwidth), ]
ont_NN2_8h_noNA$rname2 <- ifelse(ont_NN2_8h_noNA$rname == "RNA_CS_ENO2", "ENO2", "PA")
ont_NN2_8h_noNA$rname2 <- as.factor(as.character(ont_NN2_8h_noNA$rname2))

ont_SG17M_8h <- data.frame(rbind(ont_SG17M_8h_BR1, ont_SG17M_8h_BR2, ont_SG17M_8h_BR3))
ont_SG17M_8h_noNA <- ont_SG17M_8h[!is.na(ont_SG17M_8h$qwidth), ]
ont_SG17M_8h_noNA$rname2 <- ifelse(ont_SG17M_8h_noNA$rname == "RNA_CS_ENO2", "ENO2", "PA")
ont_SG17M_8h_noNA$rname2 <- as.factor(as.character(ont_SG17M_8h_noNA$rname2))

ont_SG17M_NN2_8h_noNA <- data.frame(rbind(ont_NN2_8h_noNA, ont_SG17M_8h_noNA))

nreads_PA_sense <- c((table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[1] + table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[4]),
                     (table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[1] + table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[4]),
                     (table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[1] + table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[4]),
                     (table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[1] + table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[4]),
                     (table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[1] + table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[4]),
                     (table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[1] + table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[4]),
                     (table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[1] + table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[4]),
                     (table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[1] + table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[4]),
                     (table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[1] + table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[4]),
                     (table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[1] + table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[4]),
                     (table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[1] + table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[4]),
                     (table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[1] + table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[4]))

nreads_ENO_sense  <- c(table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[7],
                     table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[7],
                     table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[7],
                     table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[7],
                     table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[7],
                     table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[7],
                     table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[7],
                     table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[7],
                     table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[7],
                     table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[7],
                     table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[7],
                     table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[7])

nreads_PA_antisense <- c((table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[2] + table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[5]),
                     (table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[2] + table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[5]),
                     (table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[2] + table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[5]),
                     (table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[2] + table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[5]),
                     (table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[2] + table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[5]),
                     (table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[2] + table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[5]),
                     (table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[2] + table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[5]),
                     (table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[2] + table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[5]),
                     (table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[2] + table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[5]),
                     (table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[2] + table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[5]),
                     (table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[2] + table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[5]),
                     (table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[2] + table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[5]))

nreads_ENO_antisense  <- c(table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[8],
                     table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[8],
                     table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[8],
                     table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[8],
                     table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[8],
                     table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[8],
                     table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[8],
                     table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[8],
                     table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[8],
                     table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[8],
                     table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[8],
                     table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[8])

# percentage of spurious antisense RNA
sum(nreads_ENO_antisense) / (sum(nreads_ENO_sense + sum(nreads_ENO_antisense)))
n_reads_sense_antisense <- data.frame(cbind(nreads_PA_sense, nreads_PA_antisense, nreads_ENO_sense, nreads_ENO_antisense))
rownames(n_reads_sense_antisense) <- c("ONT_NN2_4h_BR1", "ONT_NN2_4h_BR2", "ONT_NN2_4h_BR3", "ONT_NN2_8h_BR1", "ONT_NN2_8h_BR2", "ONT_NN2_8h_BR3",
                                       "ONT_SG17M_4h_BR1", "ONT_SG17M_4h_BR2", "ONT_SG17M_4h_BR3","ONT_SG17M_8h_BR1", "ONT_SG17M_8h_BR2", "ONT_SG17M_8h_BR3")
n_reads_sense_antisense$sum <- rowSums(n_reads_sense_antisense)

n_reads_sense_antisense$sense_PA_perc <- (n_reads_sense_antisense$nreads_PA_sense / n_reads_sense_antisense$sum) * 100
n_reads_sense_antisense$antisense_PA_perc <- (n_reads_sense_antisense$nreads_PA_antisense / n_reads_sense_antisense$sum) * 100
n_reads_sense_antisense$sense_ENO_perc <- (n_reads_sense_antisense$nreads_ENO_sense / n_reads_sense_antisense$sum) * 100
n_reads_sense_antisense$antisense_ENO_perc <- (n_reads_sense_antisense$nreads_ENO_antisense / n_reads_sense_antisense$sum) * 100
n_reads_sense_antisense$sense_PA_perc <- round(n_reads_sense_antisense$sense_PA_perc,2)
n_reads_sense_antisense$antisense_PA_perc <- round(n_reads_sense_antisense$antisense_PA_perc,2)
n_reads_sense_antisense$sense_ENO_perc <- round(n_reads_sense_antisense$sense_ENO_perc,2)
n_reads_sense_antisense$antisense_ENO_perc <- round(n_reads_sense_antisense$antisense_ENO_perc,2)
n_reads_sense_antisense$id <- rownames(n_reads_sense_antisense)
n_reads_sense_antisense_L <- gather(n_reads_sense_antisense, key="isolate", value="value", -c("id"))
n_reads_sense_antisense_L_total <- subset(n_reads_sense_antisense_L, 
                                         isolate == "nreads_PA_sense" | isolate == "nreads_PA_antisense" | 
                                           isolate == "nreads_ENO_sense" | isolate == "nreads_ENO_antisense")
n_reads_sense_antisense_L_total$isolate2 <- with(n_reads_sense_antisense_L_total,
                                                ifelse(isolate == "nreads_PA_sense", "PA, forward strand",
                                                       ifelse(isolate == "nreads_PA_antisense", "PA, reverse strand",
                                                              ifelse(isolate == "nreads_ENO_sense", "ENO, forward strand", 
                                                                     ifelse(isolate == "nreads_ENO_antisense", "ENO, reverse strand", NA)))))



libsize_plot <- ggplot(n_reads_sense_antisense_L_total, aes(x=id, y=sqrt(value), fill=isolate2)) +
  geom_col() + theme_pubr(border=TRUE, base_size=11) + coord_flip() + ylab("Read count (square-root scale)\n\n") + xlab("") +
  theme(legend.title = element_blank(), legend.position = c(0.7,0.2)) + 
  scale_fill_manual(values=c("orange", "red", "cyan", "darkblue")) +
  scale_y_continuous(labels = scales::comma, breaks = c(0, 1000, 2000))

summary_names <- c("ONT_NN2_4h_BR1", "ONT_NN2_4h_BR2", "ONT_NN2_4h_BR3", 
                   "ONT_NN2_8h_BR1", "ONT_NN2_8h_BR2", "ONT_NN2_8h_BR3", 
                   "ONT_SG17M_4h_BR1", "ONT_SG17M_4h_BR2", "ONT_SG17M_4h_BR3", 
                   "ONT_SG17M_8h_BR1", "ONT_SG17M_8h_BR2", "ONT_SG17M_8h_BR3")
nraw_reads_total <- c(164300, 124269, 254652, 142493, 91097, 
                      456000, 323314, 267750, 1401492, 448244, 
                      411871, 334857)
nbases <- c(91645543, 67762786, 103114020, 65434441, 31702988, 
            260692826, 150006926, 160857070, 547180145, 252533658, 225927528, 98126611)

df_library_size_antisense <- data.frame(cbind(summary_names, nraw_reads_total, nbases))
df_library_size_antisense$sense_ENO <- n_reads_sense_antisense$nreads_ENO_sense
df_library_size_antisense_ont <- df_library_size_antisense[1:12,]
df_library_size_antisense_ont$nraw_reads_total <- as.numeric(as.character(df_library_size_antisense_ont$nraw_reads_total))
df_library_size_antisense_ont$nbases <- as.numeric(as.character(df_library_size_antisense_ont$nbases))

dataset2 <- read_table2("input_data/pore_capacity_2.csv")
dataset2 <- data.frame(dataset2)

centroids <- aggregate(cbind(pore_capacity,cumsum)~Run,dataset2,median)
cor_test <- cor.test(dataset2$pore_capacity, dataset2$cumsum)

pore_capacity_plot <-
  ggplot(dataset2) +
  geom_jitter(aes(x=pore_capacity, y=cumsum, color=Run), size=2) +
  geom_point(data=centroids,aes(x=pore_capacity, y=cumsum, color=Run), size=6) +
  geom_line(aes(x=pore_capacity, y=cumsum, color=Run)) +
  theme_pubr(border=TRUE,legend = "none",base_size=11) +
  scale_colour_manual(values=c("lightgreen", "yellow", "darkgreen", "gold", "black", "blue", "orange", 
                               "cyan", "grey", "firebrick", "violet")) +
  xlab("Number of available pores\n\n") + ylab("Number of reads") + scale_y_continuous(labels = scales::comma)

size_capacity_plot <- ggarrange(libsize_plot, pore_capacity_plot, nrow=1, labels = c("A", "B"))

ont_noNA <- rbind(ont_NN2_4h_noNA, ont_NN2_8h_noNA, ont_SG17M_4h_noNA, ont_SG17M_8h_noNA)
ont_noNA$loq_qwith <- sqrt(ont_noNA$qwidth)
ont_noNA$mix2 <- paste(ont_noNA$mix,"_",ont_noNA$replicate)
ont_noNA_4h <- subset(ont_noNA, time == "4h")
ont_noNA_8h <- subset(ont_noNA, time == "8h")

length_4h_plot <-
  ggplot(ont_noNA_4h) +
  geom_density(aes(x=loq_qwith, color=strand)) + facet_grid(~isolate) + 
  theme_pubr(border=TRUE, base_size=11) + xlab("Read length (square-root scale)\n") +
  scale_color_manual(values=c("blue", "cyan")) + ylab("Density") + xlim(1,100) + 
  theme(legend.title = element_blank(), legend.position = c(0.1,0.75)) 

length_8h_plot <-
  ggplot(ont_noNA_8h) +
  geom_density(aes(x=loq_qwith, color=strand)) + facet_grid(~isolate) + 
  theme_pubr(border=TRUE, base_size=11, legend = "none") + xlab("Read length (square-root scale)\n") +
  scale_color_manual(values=c("blue", "cyan")) + ylab(" ") + xlim(1,100) + 
  theme(axis.ticks.y = element_blank())

length_sense_antisense <- ggarrange(length_4h_plot, length_8h_plot, nrow=1, labels=c("C", "D"))

merge_plots2 <- 
  ggarrange(size_capacity_plot, length_sense_antisense, nrow=2, heights = c(1,0.5))

#ggsave(merge_plots2, filename="save_figures/run_evaluation.tif", 
 #      device="tiff", dpi=600, units="cm", width=25.2, height=21.6)


## NN2, 4h
nn2_4h_kw <- kruskal.test(log10(ont_NN2_4h_noNA$qwidth), g=ont_NN2_4h_noNA$replicate)
nn2_4h_e2 <- epsilonSquared(log10(ont_NN2_4h_noNA$qwidth), g=ont_NN2_4h_noNA$replicate, ci=TRUE)

## NN2, 8h
nn2_8h_kw <- kruskal.test(log10(ont_NN2_8h_noNA$qwidth), g=ont_NN2_8h_noNA$replicate)
nn2_8h_e2 <- epsilonSquared(log10(ont_NN2_8h_noNA$qwidth), g=ont_NN2_8h_noNA$replicate, ci=TRUE)

## NN2, 4h-8h
ont_NN2_4h_8h_noNA <- data.frame(rbind(ont_NN2_4h_noNA, ont_NN2_8h_noNA))
ont_NN2_4h_8h_noNA$merged <- paste(ont_NN2_4h_8h_noNA$mix, ont_NN2_4h_8h_noNA$replicate, sep="_")
NN2_4h_8h_kw <- kruskal.test(log10(ont_NN2_4h_8h_noNA$qwidth), g=ont_NN2_4h_8h_noNA$merged)
NN2_4h_8h_e2 <- epsilonSquared(log10(ont_NN2_4h_8h_noNA$qwidth), g=ont_NN2_4h_8h_noNA$merged, ci=TRUE)

## SG17M, 4h
sg17m_4h_kw <- kruskal.test(log10(ont_SG17M_4h_noNA$qwidth), g=ont_SG17M_4h_noNA$replicate)
sg17m_4h_e2 <- epsilonSquared(log10(ont_SG17M_4h_noNA$qwidth), g=ont_SG17M_4h_noNA$replicate, ci=TRUE)

## SG17m, 8h
sg17m_8h_kw <- kruskal.test(log10(ont_SG17M_8h_noNA$qwidth), g=ont_SG17M_8h_noNA$replicate)
sg17m_8h_e2 <- epsilonSquared(log10(ont_SG17M_8h_noNA$qwidth), g=ont_SG17M_8h_noNA$replicate, ci=TRUE)

## SG17M, 4h-8h
ont_SG17M_4h_8h_noNA <- data.frame(rbind(ont_SG17M_4h_noNA, ont_SG17M_8h_noNA))
ont_SG17M_4h_8h_noNA$merged <- paste(ont_SG17M_4h_8h_noNA$mix, ont_SG17M_4h_8h_noNA$replicate, sep="_")
SG17M_4h_8h_kw <- kruskal.test(log10(ont_SG17M_4h_8h_noNA$qwidth), g=ont_SG17M_4h_8h_noNA$merged)
SG17M_4h_8h_e2 <- epsilonSquared(log10(ont_SG17M_4h_8h_noNA$qwidth), g=ont_SG17M_4h_8h_noNA$merged, ci=TRUE)

## NN2 vs SG17M, 4h
nn2_sg17m_4h <- data.frame(rbind(ont_NN2_4h_noNA, ont_SG17M_NN2_4h_noNA))
nn2_sg17m_4h$merged <- paste(nn2_sg17m_4h$mix, nn2_sg17m_4h$replicate, sep="_")
nn2_sg17m_4h_kw <- kruskal.test(log10(nn2_sg17m_4h$qwidth), g=nn2_sg17m_4h$merged)
nn2_sg17m_4h_e2 <- epsilonSquared(log10(nn2_sg17m_4h$qwidth), g=nn2_sg17m_4h$merged, ci=TRUE)

## NN2 vs SG17M, 8h 
nn2_sg17m_8h <- data.frame(rbind(ont_NN2_8h_noNA, ont_SG17M_8h_noNA))
nn2_sg17m_8h$merged <- paste(nn2_sg17m_8h$mix, nn2_sg17m_8h$replicate, sep="_")
nn2_sg17m_8h_kw <- kruskal.test(log10(nn2_sg17m_8h$qwidth), g=nn2_sg17m_8h$merged)
nn2_sg17m_8h_e2 <- epsilonSquared(log10(nn2_sg17m_8h$qwidth), g=nn2_sg17m_8h$merged, ci=TRUE)

# plots
nn2_4h_plot_length <- ggplot(ont_NN2_4h_noNA, aes(x=replicate, y=log10(qwidth))) + 
  geom_jitter(width = 0.1) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_pubr(border = TRUE, base_size=11) +  xlab("NN2-4h replicates") + ylab("Read length (log10 scale)")

nn2_8h_plot_length <- ggplot(ont_NN2_8h_noNA, aes(x=replicate, y=log10(qwidth))) +  
  geom_jitter(width = 0.1) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_pubr(border = TRUE, base_size=11) +  xlab("NN2-8h replicates") + ylab(" ") + 
  theme(axis.text.y = element_blank())

sg17m_4h_plot_length <-ggplot(ont_SG17M_4h_noNA, aes(x=replicate, y=log10(qwidth))) +  
  geom_jitter(width = 0.1) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_pubr(border = TRUE, base_size=11) + xlab("SG17M-4h replicates") + ylab(" ") + 
  theme(axis.text.y = element_blank())

sg18m_8h_plot_length <- ggplot(ont_SG17M_8h_noNA, aes(x=replicate, y=log10(qwidth))) +  
  geom_jitter(width = 0.1) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_pubr(border = TRUE, base_size=11) + xlab("SG17M-8h replicates") + ylab(" ") + 
  theme(axis.text.y = element_blank())

#               kruskal_statistic kruskal_pvalue e2_values e2_lower_ci e2_upper_ci
#NN2_4h                 20112.455              0   0.06640     0.06470     0.06840
#NN2_8h                 26434.382              0   0.05500     0.05380     0.05620
#SG17M_4h              153333.714              0   0.10500     0.10400     0.10600
#SG17M_8h                5118.860              0   0.00571     0.00539     0.00601
#NN2-SG17M (4h)         66040.651              0   0.03200     0.03150     0.03240
#NN2-SG17M (8h)          4352.079              0   0.00316     0.03150     0.00334

sample_names <- c("NN2-4h", "NN2-8h", "NN2-4h-8h", 
                  "SG17M-4h", "SG17M-8h", "SG17M-4h-8h", 
                  "NN2-SG17M (4h)", "NN2-SG17M (8h)")

kruskal_statistic <- c(nn2_4h_kw$statistic, nn2_8h_kw$statistic, NN2_4h_8h_kw$statistic, 
                       sg17m_4h_kw$statistic, sg17m_8h_kw$statistic, SG17M_4h_8h_kw$statistic,
                       nn2_sg17m_4h_kw$statistic, nn2_sg17m_8h_kw$statistic)

kruskal_pvalue <- c(nn2_4h_kw$p.value, nn2_8h_kw$p.value, NN2_4h_8h_kw$p.value,
                    sg17m_4h_kw$p.value, sg17m_8h_kw$p.value, SG17M_4h_8h_kw$p.value,
                    nn2_sg17m_4h_kw$p.value, nn2_sg17m_8h_kw$p.value)

e2_values <- c(nn2_4h_e2$epsilon.squared, nn2_8h_e2$epsilon.squared, NN2_4h_8h_e2$epsilon.squared,
               sg17m_4h_e2$epsilon.squared, sg17m_8h_e2$epsilon.squared, SG17M_4h_8h_e2$epsilon.squared,
               nn2_sg17m_4h_e2$epsilon.squared, nn2_sg17m_8h_e2$epsilon.squared)

e2_lower_ci <- c(nn2_4h_e2$lower.ci, nn2_8h_e2$lower.ci, NN2_4h_8h_e2$lower.ci,
                 sg17m_4h_e2$lower.ci, sg17m_8h_e2$lower.ci, SG17M_4h_8h_e2$lower.ci,
                 nn2_sg17m_4h_e2$lower.ci, nn2_sg17m_4h_e2$lower.ci)

e2_upper_ci <- c(nn2_4h_e2$upper.ci, nn2_8h_e2$upper.ci, NN2_4h_8h_e2$upper.ci,
                 sg17m_4h_e2$upper.ci, sg17m_8h_e2$upper.ci, SG17M_4h_8h_e2$upper.ci,
                 nn2_sg17m_4h_e2$upper.ci, nn2_sg17m_8h_e2$upper.ci)

length_stats_df <- data.frame(cbind(kruskal_statistic, kruskal_pvalue, e2_values, e2_lower_ci, e2_upper_ci))
rownames(length_stats_df) <- sample_names
length_stats_df$e2_values <- round(length_stats_df$e2_values,2)
length_stats_df$e2_lower_ci <- round(length_stats_df$e2_lower_ci, 2)
length_stats_df$e2_upper_ci <- round(length_stats_df$e2_upper_ci, 2)

stats_plot <- 
  ggplot(length_stats_df) +
  geom_point(aes(x=e2_values, y=rownames(length_stats_df)), size=3, colour="red") +
  geom_point(aes(x=e2_lower_ci, y=rownames(length_stats_df)), size=1, colour="black") +
  geom_point(aes(x=e2_upper_ci, y=rownames(length_stats_df)), size=1, colour="black") +
  theme_pubr(border=TRUE, legend = "none", base_size=11) + xlim(-1,1) +
  geom_vline(xintercept=-0.2) + geom_vline(xintercept=0.2) +
  xlab("Epsilon-squared effect size") + ylab("")

plot_length1 <-
  ggarrange(nn2_4h_plot_length, nn2_8h_plot_length, 
          sg17m_4h_plot_length, sg18m_8h_plot_length, 
          labels = c("A", "B", "C", "D"), nrow = 1)

merged_plot <- ggarrange(plot_length1, stats_plot, labels=c("A","E"), nrow=2, heights = c(1,0.5))

#ggsave(merged_plot, file="save_figures/plot_length2.tif", device="tiff", dpi=300, units="cm") 
       # width=9.77, height=7.64)
