
# title: "Nanopore Run Evaluation"
# author: "Marie Pust"
# date: "04 05 2021"


# Clean environment
rm(list = ls())

setwd("/R")

library(Rsamtools)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyr)
library(readr)
library(rcompanion)

ont_NN2_4h_BR1 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/d_raw_bam/NN2_CS_4h_BR1.bam")
ont_NN2_4h_BR1 <- data.frame(ont_NN2_4h_BR1)
ont_NN2_4h_BR1$isolate <- "NN2"
ont_NN2_4h_BR1$time <- "4h"
ont_NN2_4h_BR1$mix <- "NN2_4h"
ont_NN2_4h_BR1$replicate <- "BR1"
ont_NN2_4h_BR1$platform <- "Nanopore"

ont_NN2_4h_BR2 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/d_raw_bam/NN2_CS_4h_BR2.bam")
ont_NN2_4h_BR2 <- data.frame(ont_NN2_4h_BR2)
ont_NN2_4h_BR2$isolate <- "NN2"
ont_NN2_4h_BR2$time <- "4h"
ont_NN2_4h_BR2$mix <- "NN2_4h"
ont_NN2_4h_BR2$replicate <- "BR2"
ont_NN2_4h_BR2$platform <- "Nanopore"

ont_NN2_4h_BR3 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/d_raw_bam/NN2_CS_4h_BR3.bam")
ont_NN2_4h_BR3 <- data.frame(ont_NN2_4h_BR3)
ont_NN2_4h_BR3$isolate <- "NN2"
ont_NN2_4h_BR3$time <- "4h"
ont_NN2_4h_BR3$mix <- "NN2_4h"
ont_NN2_4h_BR3$replicate <- "BR3"
ont_NN2_4h_BR3$platform <- "Nanopore"

ont_NN2_8h_BR1 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/d_raw_bam/NN2_CS_8h_BR1.bam")
ont_NN2_8h_BR1 <- data.frame(ont_NN2_8h_BR1)
ont_NN2_8h_BR1$isolate <- "NN2"
ont_NN2_8h_BR1$time <- "8h"
ont_NN2_8h_BR1$mix <- "NN2_8h"
ont_NN2_8h_BR1$replicate <- "BR1"
ont_NN2_8h_BR1$platform <- "Nanopore"

ont_NN2_8h_BR2 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/d_raw_bam/NN2_CS_8h_BR2.bam")
ont_NN2_8h_BR2 <- data.frame(ont_NN2_8h_BR2)
ont_NN2_8h_BR2$isolate <- "NN2"
ont_NN2_8h_BR2$time <- "8h"
ont_NN2_8h_BR2$mix <- "NN2_8h"
ont_NN2_8h_BR2$replicate <- "BR2"
ont_NN2_8h_BR2$platform <- "Nanopore"

ont_NN2_8h_BR3 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/d_raw_bam/NN2_CS_8h_BR3.bam")
ont_NN2_8h_BR3 <- data.frame(ont_NN2_8h_BR3)
ont_NN2_8h_BR3$isolate <- "NN2"
ont_NN2_8h_BR3$time <- "8h"
ont_NN2_8h_BR3$mix <- "NN2_8h"
ont_NN2_8h_BR3$replicate <- "BR3"
ont_NN2_8h_BR3$platform <- "Nanopore"

ont_SG17M_4h_BR1 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/d_raw_bam/SG17M_CS_4h_BR1.bam")
ont_SG17M_4h_BR1 <- data.frame(ont_SG17M_4h_BR1)
ont_SG17M_4h_BR1$isolate <- "SG17M"
ont_SG17M_4h_BR1$time <- "4h"
ont_SG17M_4h_BR1$mix <- "SG17M_4h"
ont_SG17M_4h_BR1$replicate <- "BR1"
ont_SG17M_4h_BR1$platform <- "Nanopore"

ont_SG17M_4h_BR2 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/d_raw_bam/SG17M_CS_4h_BR2.bam")
ont_SG17M_4h_BR2 <- data.frame(ont_SG17M_4h_BR2)
ont_SG17M_4h_BR2$isolate <- "SG17M"
ont_SG17M_4h_BR2$time <- "4h"
ont_SG17M_4h_BR2$mix <- "SG17M_4h"
ont_SG17M_4h_BR2$replicate <- "BR2"
ont_SG17M_4h_BR2$platform <- "Nanopore"

ont_SG17M_4h_BR3 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/d_raw_bam/SG17M_CS_4h_BR3.bam")
ont_SG17M_4h_BR3 <- data.frame(ont_SG17M_4h_BR3)
ont_SG17M_4h_BR3$isolate <- "SG17M"
ont_SG17M_4h_BR3$time <- "4h"
ont_SG17M_4h_BR3$mix <- "SG17M_4h"
ont_SG17M_4h_BR3$replicate <- "BR3"
ont_SG17M_4h_BR3$platform <- "Nanopore"

ont_SG17M_8h_BR1 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/d_raw_bam/SG17M_CS_8h_BR1.bam")
ont_SG17M_8h_BR1 <- data.frame(ont_SG17M_8h_BR1)
ont_SG17M_8h_BR1$isolate <- "SG17M"
ont_SG17M_8h_BR1$time <- "8h"
ont_SG17M_8h_BR1$mix <- "SG17M_8h"
ont_SG17M_8h_BR1$replicate <- "BR1"
ont_SG17M_8h_BR1$platform <- "Nanopore"

ont_SG17M_8h_BR2 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/d_raw_bam/SG17M_CS_8h_BR2.bam")
ont_SG17M_8h_BR2 <- data.frame(ont_SG17M_8h_BR2)
ont_SG17M_8h_BR2$isolate <- "SG17M"
ont_SG17M_8h_BR2$time <- "8h"
ont_SG17M_8h_BR2$mix <- "SG17M_8h"
ont_SG17M_8h_BR2$replicate <- "BR2"
ont_SG17M_8h_BR2$platform <- "Nanopore"

ont_SG17M_8h_BR3 <-  scanBam("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/d_raw_bam/SG17M_CS_8h_BR3.bam")
ont_SG17M_8h_BR3 <- data.frame(ont_SG17M_8h_BR3)
ont_SG17M_8h_BR3$isolate <- "SG17M"
ont_SG17M_8h_BR3$time <- "8h"
ont_SG17M_8h_BR3$mix <- "SG17M_8h"
ont_SG17M_8h_BR3$replicate <- "BR3"
ont_SG17M_8h_BR3$platform <- "Nanopore"

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
nreads_PA_sense <- c((table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[1]),
                     (table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[1]),
                     (table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[1]),
                     (table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[1]),
                     (table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[1]),
                     (table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[1]),
                     (table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[1]),
                     (table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[1]),
                     (table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[1]),
                     (table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[1]),
                     (table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[1]),
                     (table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[1]))

nreads_ENO_sense  <- c(table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[4],
                     table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[4],
                     table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[4],
                     table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[4],
                     table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[4],
                     table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[4],
                     table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[4],
                     table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[4],
                     table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[4],
                     table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[4],
                     table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[4],
                     table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[4])

nreads_PA_antisense <- c((table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[2]),
                     (table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[2]),
                     (table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[2]),
                     (table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[2]),
                     (table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[2]),
                     (table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[2]),
                     (table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[2]),
                     (table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[2]),
                     (table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[2]),
                     (table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[2]),
                     (table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[2]),
                     (table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[2]))

nreads_ENO_antisense  <- c(table(ont_NN2_4h_BR1$strand, ont_NN2_4h_BR1$rname)[5],
                     table(ont_NN2_4h_BR2$strand, ont_NN2_4h_BR2$rname)[5],
                     table(ont_NN2_4h_BR3$strand, ont_NN2_4h_BR3$rname)[5],
                     table(ont_NN2_8h_BR1$strand, ont_NN2_8h_BR1$rname)[5],
                     table(ont_NN2_8h_BR2$strand, ont_NN2_8h_BR2$rname)[5],
                     table(ont_NN2_8h_BR3$strand, ont_NN2_8h_BR3$rname)[5],
                     table(ont_SG17M_4h_BR1$strand, ont_SG17M_4h_BR1$rname)[5],
                     table(ont_SG17M_4h_BR2$strand, ont_SG17M_4h_BR2$rname)[5],
                     table(ont_SG17M_4h_BR3$strand, ont_SG17M_4h_BR3$rname)[5],
                     table(ont_SG17M_8h_BR1$strand, ont_SG17M_8h_BR1$rname)[5],
                     table(ont_SG17M_8h_BR2$strand, ont_SG17M_8h_BR2$rname)[5],
                     table(ont_SG17M_8h_BR3$strand, ont_SG17M_8h_BR3$rname)[5])


# percentage of spurious antisense RNA
sum(nreads_ENO_antisense)
sum(nreads_ENO_sense) + sum(nreads_ENO_antisense)
sum(nreads_ENO_antisense) / (sum(nreads_ENO_sense + sum(nreads_ENO_antisense)))

n_reads_sense_antisense <- data.frame(cbind((nreads_PA_sense), 
                                            (nreads_PA_antisense), 
                                            (nreads_ENO_sense), 
                                            (nreads_ENO_antisense)))
n_reads_sense_antisense <- round(n_reads_sense_antisense,1)

colnames(n_reads_sense_antisense) <- c("nreads_PA_sense", "nreads_PA_antisense", 
                                       "nreads_ENO_sense", "nreads_ENO_antisense")
rownames(n_reads_sense_antisense) <- c("ONT_NN2_4h_BR1", "ONT_NN2_4h_BR2", "ONT_NN2_4h_BR3", 
                                       "ONT_NN2_8h_BR1", "ONT_NN2_8h_BR2", "ONT_NN2_8h_BR3",
                                       "ONT_SG17M_4h_BR1", "ONT_SG17M_4h_BR2", "ONT_SG17M_4h_BR3",
                                       "ONT_SG17M_8h_BR1", "ONT_SG17M_8h_BR2", "ONT_SG17M_8h_BR3")
n_reads_sense_antisense$sum <- rowSums(n_reads_sense_antisense)
n_reads_sense_antisense$id <- rownames(n_reads_sense_antisense)
n_reads_sense_antisense_L <- gather(n_reads_sense_antisense, key="isolate", value="value", -c("id"))
n_reads_sense_antisense_L_total <- subset(n_reads_sense_antisense_L, 
                                         isolate == "nreads_PA_sense" | isolate == "nreads_PA_antisense" | 
                                           isolate == "nreads_ENO_sense" | isolate == "nreads_ENO_antisense")
n_reads_sense_antisense_L_total$isolate2 <- with(n_reads_sense_antisense_L_total,
                                                ifelse(isolate == "nreads_PA_sense", "PA, forward strand",
                                                       ifelse(isolate == "nreads_PA_antisense", "PA, reverse strand",
                                                              ifelse(isolate == "nreads_ENO_sense", "ENO, forward strand", 
                                                                     ifelse(isolate == "nreads_ENO_antisense", 
                                                                            "ENO, reverse strand", NA)))))


summary_names <- c("ONT_NN2_4h_BR1","ONT_NN2_4h_BR2","ONT_NN2_4h_BR3",
                   "ONT_NN2_8h_BR1","ONT_NN2_8h_BR2","ONT_NN2_8h_BR3", 
                   "ONT_SG17M_4h_BR1","ONT_SG17M_4h_BR2","ONT_SG17M_4h_BR3", 
                   "ONT_SG17M_8h_BR1","ONT_SG17M_8h_BR2","ONT_SG17M_8h_BR3")

nraw_reads_total <- c(164300,124269,254652, 
                      142493,91097,456000, 
                      323314,267750,1401492, 
                      448244,411871,334857)

nbases <- c(91645543,67762786,103114020, 
            65434441,31702988,260692826, 
            150006926,160857070,547180145, 
            252533658,225927528,98126611)

df_library_size_antisense <- data.frame(cbind(summary_names, nraw_reads_total, nbases))
df_library_size_antisense$sense_ENO <- n_reads_sense_antisense$nreads_ENO_sense
df_library_size_antisense_ont <- df_library_size_antisense[1:12,]
df_library_size_antisense_ont$nraw_reads_total <- as.numeric(as.character(df_library_size_antisense_ont$nraw_reads_total))
df_library_size_antisense_ont$nbases <- as.numeric(as.character(df_library_size_antisense_ont$nbases))

dataset2 <- read_table2("input_data/pore_capacity_2.csv")
dataset2 <- data.frame(dataset2)
dataset2 <- subset(dataset2, minute < 361)
aggregate(dataset2$cumsum, by=list(Category=dataset2$Run), FUN=sum)
centroids <- aggregate(cbind(pore_capacity,cumsum)~Run,dataset2,median)
cor_test <- cor.test(dataset2$pore_capacity, dataset2$cumsum)

libsize_plot <- ggplot(n_reads_sense_antisense_L_total, aes(x=id, color=isolate2)) +
  geom_jitter(aes(x=id, y=sqrt(value), fill=isolate2), shape = 15, size=2.5, width=0.05) +
  theme_pubr(border=TRUE, base_size=11) + coord_flip() + 
  ylab("Number of reads") + 
  xlab("") +
  theme(legend.title = element_blank(), 
        legend.position = c(0.7,0.2), 
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10)) + 
  scale_color_manual(values=c("orange", "red", "cyan", "darkblue")) + 
  scale_y_continuous(breaks = c(0, 500, 1000), 
                     labels=c("0", expression("500"^2), expression("1,000"^2)))

pore_capacity_plot <- ggplot(dataset2) + 
  geom_jitter(aes(x=pore_capacity, y=cumsum, color=Run), size=2) +
  geom_point(data=centroids, aes(x=pore_capacity, y=cumsum, color=Run), size=6) + 
  geom_line(aes(x=pore_capacity, y=cumsum, color=Run)) + 
  theme_pubr(border=TRUE, legend = "none", base_size=11) +
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10)) +
  scale_colour_manual(values=c("lightgreen", "aquamarine4", 
                               "darkgreen", "yellow", 
                               "black", "navajowhite2", 
                               "purple", "pink", 
                               "grey", "firebrick", "violet")) +
  xlab("Number of available pores") + ylab("Number of reads") + 
  scale_y_continuous(labels = scales::comma)


size_capacity_plot <- ggarrange(libsize_plot, pore_capacity_plot, nrow=1, labels = c("A", "B"),
                                font.label = list(size = 14, color = "black"), widths = c(1,0.8))

ggsave(size_capacity_plot, filename="save_figures/Figure_01.tif", device="tiff", dpi=300, units="cm", width=23, height=15)

# Load input to evaluate read length####
ont_NN2_4h_BR1_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_4h_BR1.sorted.bam.antisense_summary.bed", 
                             "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_4h_BR1_antisense <- data.frame(ont_NN2_4h_BR1_antisense)
ont_NN2_4h_BR1_antisense$isolate <- "NN2"
ont_NN2_4h_BR1_antisense$time <- "4h"
ont_NN2_4h_BR1_antisense$mix <- "NN2_4h"
ont_NN2_4h_BR1_antisense$replicate <- "BR1"
ont_NN2_4h_BR1_antisense$platform <- "Nanopore"
ont_NN2_4h_BR1_antisense$read_length <- nchar(ont_NN2_4h_BR1_antisense$X4)
ont_NN2_4h_BR1_antisense$strand_transcript <- "Antisense"

ont_NN2_4h_BR2_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_4h_BR2.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_4h_BR2_antisense <- data.frame(ont_NN2_4h_BR2_antisense)
ont_NN2_4h_BR2_antisense$isolate <- "NN2"
ont_NN2_4h_BR2_antisense$time <- "4h"
ont_NN2_4h_BR2_antisense$mix <- "NN2_4h"
ont_NN2_4h_BR2_antisense$replicate <- "BR2"
ont_NN2_4h_BR2_antisense$platform <- "Nanopore"
ont_NN2_4h_BR2_antisense$read_length <- nchar(ont_NN2_4h_BR2_antisense$X4)
ont_NN2_4h_BR2_antisense$strand_transcript <- "Antisense"

ont_NN2_4h_BR3_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_4h_BR3.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_4h_BR3_antisense <- data.frame(ont_NN2_4h_BR3_antisense)
ont_NN2_4h_BR3_antisense$isolate <- "NN2"
ont_NN2_4h_BR3_antisense$time <- "4h"
ont_NN2_4h_BR3_antisense$mix <- "NN2_4h"
ont_NN2_4h_BR3_antisense$replicate <- "BR3"
ont_NN2_4h_BR3_antisense$platform <- "Nanopore"
ont_NN2_4h_BR3_antisense$read_length <- nchar(ont_NN2_4h_BR3_antisense$X4)
ont_NN2_4h_BR3_antisense$strand_transcript <- "Antisense"

ont_NN2_4h_BR1_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_4h_BR1.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_4h_BR1_sense <- data.frame(ont_NN2_4h_BR1_sense)
ont_NN2_4h_BR1_sense$isolate <- "NN2"
ont_NN2_4h_BR1_sense$time <- "4h"
ont_NN2_4h_BR1_sense$mix <- "NN2_4h"
ont_NN2_4h_BR1_sense$replicate <- "BR1"
ont_NN2_4h_BR1_sense$platform <- "Nanopore"
ont_NN2_4h_BR1_sense$read_length <- nchar(ont_NN2_4h_BR1_sense$X4)
ont_NN2_4h_BR1_sense$strand_transcript <- "Sense"

ont_NN2_4h_BR2_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_4h_BR2.sorted.bam.sense_summary.bed","\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_4h_BR2_sense <- data.frame(ont_NN2_4h_BR2_sense)
ont_NN2_4h_BR2_sense$isolate <- "NN2"
ont_NN2_4h_BR2_sense$time <- "4h"
ont_NN2_4h_BR2_sense$mix <- "NN2_4h"
ont_NN2_4h_BR2_sense$replicate <- "BR2"
ont_NN2_4h_BR2_sense$platform <- "Nanopore"
ont_NN2_4h_BR2_sense$read_length <- nchar(ont_NN2_4h_BR2_sense$X4)
ont_NN2_4h_BR2_sense$strand_transcript <- "Sense"

ont_NN2_4h_BR3_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_4h_BR3.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_4h_BR3_sense <- data.frame(ont_NN2_4h_BR3_sense)
ont_NN2_4h_BR3_sense$isolate <- "NN2"
ont_NN2_4h_BR3_sense$time <- "4h"
ont_NN2_4h_BR3_sense$mix <- "NN2_4h"
ont_NN2_4h_BR3_sense$replicate <- "BR3"
ont_NN2_4h_BR3_sense$platform <- "Nanopore"
ont_NN2_4h_BR3_sense$read_length <- nchar(ont_NN2_4h_BR3_sense$X4)
ont_NN2_4h_BR3_sense$strand_transcript <- "Sense"

ont_NN2_8h_BR1_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_8h_BR1.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_8h_BR1_antisense <- data.frame(ont_NN2_8h_BR1_antisense)
ont_NN2_8h_BR1_antisense$isolate <- "NN2"
ont_NN2_8h_BR1_antisense$time <- "8h"
ont_NN2_8h_BR1_antisense$mix <- "NN2_8h"
ont_NN2_8h_BR1_antisense$replicate <- "BR1"
ont_NN2_8h_BR1_antisense$platform <- "Nanopore"
ont_NN2_8h_BR1_antisense$read_length <- nchar(ont_NN2_8h_BR1_antisense$X4)
ont_NN2_8h_BR1_antisense$strand_transcript <- "Antisense"

ont_NN2_8h_BR2_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_8h_BR2.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_8h_BR2_antisense <- data.frame(ont_NN2_8h_BR2_antisense)
ont_NN2_8h_BR2_antisense$isolate <- "NN2"
ont_NN2_8h_BR2_antisense$time <- "8h"
ont_NN2_8h_BR2_antisense$mix <- "NN2_8h"
ont_NN2_8h_BR2_antisense$replicate <- "BR2"
ont_NN2_8h_BR2_antisense$platform <- "Nanopore"
ont_NN2_8h_BR2_antisense$read_length <- nchar(ont_NN2_8h_BR2_antisense$X4)
ont_NN2_8h_BR2_antisense$strand_transcript <- "Antisense"

ont_NN2_8h_BR3_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_8h_BR3.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_8h_BR3_antisense <- data.frame(ont_NN2_8h_BR3_antisense)
ont_NN2_8h_BR3_antisense$isolate <- "NN2"
ont_NN2_8h_BR3_antisense$time <- "8h"
ont_NN2_8h_BR3_antisense$mix <- "NN2_8h"
ont_NN2_8h_BR3_antisense$replicate <- "BR3"
ont_NN2_8h_BR3_antisense$platform <- "Nanopore"
ont_NN2_8h_BR3_antisense$read_length <- nchar(ont_NN2_8h_BR3_antisense$X4)
ont_NN2_8h_BR3_antisense$strand_transcript <- "Antisense"

ont_NN2_8h_BR1_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_8h_BR1.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_8h_BR1_sense <- data.frame(ont_NN2_8h_BR1_sense)
ont_NN2_8h_BR1_sense$isolate <- "NN2"
ont_NN2_8h_BR1_sense$time <- "8h"
ont_NN2_8h_BR1_sense$mix <- "NN2_8h"
ont_NN2_8h_BR1_sense$replicate <- "BR1"
ont_NN2_8h_BR1_sense$platform <- "Nanopore"
ont_NN2_8h_BR1_sense$read_length <- nchar(ont_NN2_8h_BR1_sense$X4)
ont_NN2_8h_BR1_sense$strand_transcript <- "Sense"

ont_NN2_8h_BR2_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_8h_BR2.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_8h_BR2_sense <- data.frame(ont_NN2_8h_BR2_sense)
ont_NN2_8h_BR2_sense$isolate <- "NN2"
ont_NN2_8h_BR2_sense$time <- "8h"
ont_NN2_8h_BR2_sense$mix <- "NN2_8h"
ont_NN2_8h_BR2_sense$replicate <- "BR2"
ont_NN2_8h_BR2_sense$platform <- "Nanopore"
ont_NN2_8h_BR2_sense$read_length <- nchar(ont_NN2_8h_BR2_sense$X4)
ont_NN2_8h_BR2_sense$strand_transcript <- "Sense"

ont_NN2_8h_BR3_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/NN2_run/c_raw_bam_sam/antisense_transcription/NN2_CS_8h_BR3.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_NN2_8h_BR3_sense <- data.frame(ont_NN2_8h_BR3_sense)
ont_NN2_8h_BR3_sense$isolate <- "NN2"
ont_NN2_8h_BR3_sense$time <- "8h"
ont_NN2_8h_BR3_sense$mix <- "NN2_8h"
ont_NN2_8h_BR3_sense$replicate <- "BR3"
ont_NN2_8h_BR3_sense$platform <- "Nanopore"
ont_NN2_8h_BR3_sense$read_length <- nchar(ont_NN2_8h_BR3_sense$X4)
ont_NN2_8h_BR3_sense$strand_transcript <- "Sense"

ont_SG17M_4h_BR1_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_4h_BR1.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_4h_BR1_antisense <- data.frame(ont_SG17M_4h_BR1_antisense)
ont_SG17M_4h_BR1_antisense$isolate <- "SG17M"
ont_SG17M_4h_BR1_antisense$time <- "4h"
ont_SG17M_4h_BR1_antisense$mix <- "SG17M_4h"
ont_SG17M_4h_BR1_antisense$replicate <- "BR1"
ont_SG17M_4h_BR1_antisense$platform <- "Nanopore"
ont_SG17M_4h_BR1_antisense$read_length <- nchar(ont_SG17M_4h_BR1_antisense$X4)
ont_SG17M_4h_BR1_antisense$strand_transcript <- "Antisense"

ont_SG17M_4h_BR2_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_4h_BR2.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_4h_BR2_antisense <- data.frame(ont_SG17M_4h_BR2_antisense)
ont_SG17M_4h_BR2_antisense$isolate <- "SG17M"
ont_SG17M_4h_BR2_antisense$time <- "4h"
ont_SG17M_4h_BR2_antisense$mix <- "SG17M_4h"
ont_SG17M_4h_BR2_antisense$replicate <- "BR2"
ont_SG17M_4h_BR2_antisense$platform <- "Nanopore"
ont_SG17M_4h_BR2_antisense$read_length <- nchar(ont_SG17M_4h_BR2_antisense$X4)
ont_SG17M_4h_BR2_antisense$strand_transcript <- "Antisense"

ont_SG17M_4h_BR3_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_4h_BR3.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_4h_BR3_antisense <- data.frame(ont_SG17M_4h_BR3_antisense)
ont_SG17M_4h_BR3_antisense$isolate <- "SG17M"
ont_SG17M_4h_BR3_antisense$time <- "4h"
ont_SG17M_4h_BR3_antisense$mix <- "SG17M_4h"
ont_SG17M_4h_BR3_antisense$replicate <- "BR3"
ont_SG17M_4h_BR3_antisense$platform <- "Nanopore"
ont_SG17M_4h_BR3_antisense$read_length <- nchar(ont_SG17M_4h_BR3_antisense$X4)
ont_SG17M_4h_BR3_antisense$strand_transcript <- "Antisense"

ont_SG17M_4h_BR1_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_4h_BR1.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_4h_BR1_sense <- data.frame(ont_SG17M_4h_BR1_sense)
ont_SG17M_4h_BR1_sense$isolate <- "SG17M"
ont_SG17M_4h_BR1_sense$time <- "4h"
ont_SG17M_4h_BR1_sense$mix <- "SG17M_4h"
ont_SG17M_4h_BR1_sense$replicate <- "BR1"
ont_SG17M_4h_BR1_sense$platform <- "Nanopore"
ont_SG17M_4h_BR1_sense$read_length <- nchar(ont_SG17M_4h_BR1_sense$X4)
ont_SG17M_4h_BR1_sense$strand_transcript <- "Sense"

ont_SG17M_4h_BR2_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_4h_BR2.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_4h_BR2_sense <- data.frame(ont_SG17M_4h_BR2_sense)
ont_SG17M_4h_BR2_sense$isolate <- "SG17M"
ont_SG17M_4h_BR2_sense$time <- "4h"
ont_SG17M_4h_BR2_sense$mix <- "SG17M_4h"
ont_SG17M_4h_BR2_sense$replicate <- "BR2"
ont_SG17M_4h_BR2_sense$platform <- "Nanopore"
ont_SG17M_4h_BR2_sense$read_length <- nchar(ont_SG17M_4h_BR2_sense$X4)
ont_SG17M_4h_BR2_sense$strand_transcript <- "Sense"

ont_SG17M_4h_BR3_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_4h_BR3.sorted.bam.sense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_4h_BR3_sense <- data.frame(ont_SG17M_4h_BR3_sense)
ont_SG17M_4h_BR3_sense$isolate <- "SG17M"
ont_SG17M_4h_BR3_sense$time <- "4h"
ont_SG17M_4h_BR3_sense$mix <- "SG17M_4h"
ont_SG17M_4h_BR3_sense$replicate <- "BR3"
ont_SG17M_4h_BR3_sense$platform <- "Nanopore"
ont_SG17M_4h_BR3_sense$read_length <- nchar(ont_SG17M_4h_BR3_sense$X4)
ont_SG17M_4h_BR3_sense$strand_transcript <- "Sense"

ont_SG17M_8h_BR1_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_8h_BR1.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_8h_BR1_antisense <- data.frame(ont_SG17M_8h_BR1_antisense)
ont_SG17M_8h_BR1_antisense$isolate <- "SG17M"
ont_SG17M_8h_BR1_antisense$time <- "8h"
ont_SG17M_8h_BR1_antisense$mix <- "SG17M_8h"
ont_SG17M_8h_BR1_antisense$replicate <- "BR1"
ont_SG17M_8h_BR1_antisense$platform <- "Nanopore"
ont_SG17M_8h_BR1_antisense$read_length <- nchar(ont_SG17M_8h_BR1_antisense$X4)
ont_SG17M_8h_BR1_antisense$strand_transcript <- "Antisense"

ont_SG17M_8h_BR2_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_8h_BR2.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_8h_BR2_antisense <- data.frame(ont_SG17M_8h_BR2_antisense)
ont_SG17M_8h_BR2_antisense$isolate <- "SG17M"
ont_SG17M_8h_BR2_antisense$time <- "8h"
ont_SG17M_8h_BR2_antisense$mix <- "SG17M_8h"
ont_SG17M_8h_BR2_antisense$replicate <- "BR2"
ont_SG17M_8h_BR2_antisense$platform <- "Nanopore"
ont_SG17M_8h_BR2_antisense$read_length <- nchar(ont_SG17M_8h_BR2_antisense$X4)
ont_SG17M_8h_BR2_antisense$strand_transcript <- "Antisense"

ont_SG17M_8h_BR3_antisense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_8h_BR3.sorted.bam.antisense_summary.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_8h_BR3_antisense <- data.frame(ont_SG17M_8h_BR3_antisense)
ont_SG17M_8h_BR3_antisense$isolate <- "SG17M"
ont_SG17M_8h_BR3_antisense$time <- "8h"
ont_SG17M_8h_BR3_antisense$mix <- "SG17M_8h"
ont_SG17M_8h_BR3_antisense$replicate <- "BR3"
ont_SG17M_8h_BR3_antisense$platform <- "Nanopore"
ont_SG17M_8h_BR3_antisense$read_length <- nchar(ont_SG17M_8h_BR3_antisense$X4)
ont_SG17M_8h_BR3_antisense$strand_transcript <- "Antisense"

ont_SG17M_8h_BR1_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_8h_BR1.sorted.bam.sense_summary.bed", 
                                     "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_8h_BR1_sense <- data.frame(ont_SG17M_8h_BR1_sense)
ont_SG17M_8h_BR1_sense$isolate <- "SG17M"
ont_SG17M_8h_BR1_sense$time <- "8h"
ont_SG17M_8h_BR1_sense$mix <- "SG17M_8h"
ont_SG17M_8h_BR1_sense$replicate <- "BR1"
ont_SG17M_8h_BR1_sense$platform <- "Nanopore"
ont_SG17M_8h_BR1_sense$read_length <- nchar(ont_SG17M_8h_BR1_sense$X4)
ont_SG17M_8h_BR1_sense$strand_transcript <- "Sense"

ont_SG17M_8h_BR2_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_8h_BR2.sorted.bam.sense_summary.bed", 
                             "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_8h_BR2_sense <- data.frame(ont_SG17M_8h_BR2_sense)
ont_SG17M_8h_BR2_sense$isolate <- "SG17M"
ont_SG17M_8h_BR2_sense$time <- "8h"
ont_SG17M_8h_BR2_sense$mix <- "SG17M_8h"
ont_SG17M_8h_BR2_sense$replicate <- "BR2"
ont_SG17M_8h_BR2_sense$platform <- "Nanopore"
ont_SG17M_8h_BR2_sense$read_length <- nchar(ont_SG17M_8h_BR2_sense$X4)
ont_SG17M_8h_BR2_sense$strand_transcript <- "Sense"

ont_SG17M_8h_BR3_sense <- read_delim("/ngsssd1/tuem_mp/projects_MariePust_2018/rna_ont_nn2_sg17m_2021/MinionExperiments202009/b_minimap2/repeat_mapping/SG17M_run/c_raw_bam_sam/antisense_transcription/SG17M_CS_8h_BR3.sorted.bam.sense_summary.bed", 
                             "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ont_SG17M_8h_BR3_sense <- data.frame(ont_SG17M_8h_BR3_sense)
ont_SG17M_8h_BR3_sense$isolate <- "SG17M"
ont_SG17M_8h_BR3_sense$time <- "8h"
ont_SG17M_8h_BR3_sense$mix <- "SG17M_8h"
ont_SG17M_8h_BR3_sense$replicate <- "BR3"
ont_SG17M_8h_BR3_sense$platform <- "Nanopore"
ont_SG17M_8h_BR3_sense$read_length <- nchar(ont_SG17M_8h_BR3_sense$X4)
ont_SG17M_8h_BR3_sense$strand_transcript <- "Sense"

ont_NN2_4h <- data.frame(rbind(ont_NN2_4h_BR1_antisense, 
                               ont_NN2_4h_BR2_antisense,
                               ont_NN2_4h_BR3_antisense,
                               ont_NN2_4h_BR1_sense, 
                               ont_NN2_4h_BR2_sense, 
                               ont_NN2_4h_BR3_sense))
ont_NN2_4h_noNA <- ont_NN2_4h[!is.na(ont_NN2_4h$read_length), ]
ont_NN2_4h_noNA$rname2 <- ifelse(ont_NN2_4h_noNA$X1 == "RNA_CS_ENO2", "ENO2","PA")
ont_NN2_4h_noNA$rname2 <- as.factor(as.character(ont_NN2_4h_noNA$rname2))

ont_NN2_8h <- data.frame(rbind(ont_NN2_8h_BR1_antisense, 
                               ont_NN2_8h_BR2_antisense, 
                               ont_NN2_8h_BR3_antisense,
                               ont_NN2_8h_BR1_sense, 
                               ont_NN2_8h_BR2_sense, 
                               ont_NN2_8h_BR3_sense))
ont_NN2_8h_noNA <- ont_NN2_8h[!is.na(ont_NN2_8h$read_length), ]
ont_NN2_8h_noNA$rname2 <- ifelse(ont_NN2_8h_noNA$X1 == "RNA_CS_ENO2", "ENO2","PA")
ont_NN2_8h_noNA$rname2 <- as.factor(as.character(ont_NN2_8h_noNA$rname2))

ont_SG17M_4h <- data.frame(rbind(ont_SG17M_4h_BR1_antisense, 
                                 ont_SG17M_4h_BR2_antisense, 
                                 ont_SG17M_4h_BR3_antisense,
                                 ont_SG17M_4h_BR1_sense, 
                                 ont_SG17M_4h_BR2_sense, 
                                 ont_SG17M_4h_BR3_sense))
ont_SG17M_4h_noNA <- ont_SG17M_4h[!is.na(ont_SG17M_4h$read_length), ]
ont_SG17M_4h_noNA$rname2 <- ifelse(ont_SG17M_4h_noNA$X1 == "RNA_CS_ENO2", "ENO2","PA")
ont_SG17M_4h_noNA$rname2 <- as.factor(as.character(ont_SG17M_4h_noNA$rname2))

ont_SG17M_8h <- data.frame(rbind(ont_SG17M_8h_BR1_antisense, 
                                 ont_SG17M_8h_BR2_antisense, 
                                 ont_SG17M_8h_BR3_antisense,
                                 ont_SG17M_8h_BR1_sense, 
                                 ont_SG17M_8h_BR2_sense, 
                                 ont_SG17M_8h_BR3_sense))
ont_SG17M_8h_noNA <- ont_SG17M_8h[!is.na(ont_SG17M_8h$read_length), ]
ont_SG17M_8h_noNA$rname2 <- ifelse(ont_SG17M_8h_noNA$X1 == "RNA_CS_ENO2", "ENO2","PA")
ont_SG17M_8h_noNA$rname2 <- as.factor(as.character(ont_SG17M_8h_noNA$rname2))

ont_noNA <- rbind(ont_NN2_4h_noNA, 
                  ont_NN2_8h_noNA, 
                  ont_SG17M_4h_noNA, 
                  ont_SG17M_8h_noNA)
ont_noNA <- subset(ont_noNA, read_length > 50)

ont_noNA$sqrt_qwith <- sqrt(ont_noNA$read_length)
ont_noNA$mix2 <- paste(ont_noNA$mix,"_",ont_noNA$replicate)
ont_noNA_4h <- subset(ont_noNA, time == "4h")
ont_noNA_8h <- subset(ont_noNA, time == "8h")

length_4h_plot <- ggplot(ont_noNA_4h) + 
  geom_density(aes(x=sqrt_qwith, colour=strand_transcript)) + 
  facet_grid(~isolate) + 
  theme_pubr(border=TRUE, base_size=11) + 
  xlab("RNA read length") + 
  scale_color_manual(values=c("cyan", "blue")) + ylab("Density") + 
  scale_x_continuous(breaks=c(30,60,90), 
                     labels=c(expression("30"^2),
                              expression("60"^2),
                              expression("90"^2)),
                     limits=c(10,90)) +
  theme(legend.title = element_blank(), 
        legend.position = c(0.2,0.75)) +
  ylim(0,0.3)

ont_noNA_8h <- subset(ont_noNA_8h, mix2 != "NN2_8h _ BR1" & mix2 != "NN2_8h _ BR3")

length_8h_plot <- ggplot(ont_noNA_8h) + 
  geom_density(aes(x=sqrt_qwith, colour=strand_transcript)) + 
  facet_grid(~isolate) + 
  theme_pubr(border=TRUE, base_size=11, legend = "none") + 
  xlab("RNA read length") +
  scale_color_manual(values=c("cyan", "blue")) + ylab(" ") +
  scale_x_continuous(breaks=c(30,60,90), 
                     labels=c(expression("30"^2),
                              expression("60"^2),
                              expression("90"^2)),
                     limits=c(10,90)) +
  theme(axis.ticks.y = element_blank())+
  ylim(0,0.3)

length_sense_antisense <- ggarrange(length_4h_plot, length_8h_plot, nrow=1, labels=c("A", "B"))
length_sense_antisense
ggsave(length_sense_antisense, file="save_figures/length_sense_antisense.tif", device="tiff", dpi=300, units="cm", width=25, height=13.4)


ont_noNA_4h_SG17M <- subset(ont_noNA_4h, isolate == "SG17M")
set.seed(111)
ont_noNA_4h_SG17M_n1000 <- ont_noNA_4h_SG17M[sample(nrow(ont_noNA_4h_SG17M), 10000), ]
wilcoxonR(ont_noNA_4h_SG17M_n1000$sqrt_qwith, g=ont_noNA_4h_SG17M_n1000$strand_transcript, ci=TRUE)
# -0.0182   -0.085    0.043

ont_noNA_4h_NN2 <- subset(ont_noNA_4h, isolate == "NN2")                       
ont_noNA_4h_NN2_n1000 <- ont_noNA_4h_NN2[sample(nrow(ont_noNA_4h_NN2), 10000), ]
wilcoxonR(ont_noNA_4h_NN2_n1000$sqrt_qwith, g=ont_noNA_4h_NN2_n1000$strand_transcript, ci=TRUE)
# 0.0383  -0.0239    0.099

ont_noNA_8h_SG17M <- subset(ont_noNA_8h, isolate == "SG17M")
ont_noNA_8h_SG17M_n1000 <- ont_noNA_8h_SG17M[sample(nrow(ont_noNA_8h_SG17M), 10000), ]
wilcoxonR(ont_noNA_8h_SG17M_n1000$sqrt_qwith, g=ont_noNA_8h_SG17M_n1000$strand_transcript, ci=TRUE)
# -0.0522   -0.112  0.00947

ont_noNA_8h_NN2 <- subset(ont_noNA_8h, isolate == "NN2")    
ont_noNA_8h_NN2_n1000 <- ont_noNA_8h_NN2[sample(nrow(ont_noNA_8h_NN2), 10000), ]
wilcoxonR(ont_noNA_8h_NN2_n1000$sqrt_qwith, g=ont_noNA_8h_NN2_n1000$strand_transcript, ci=TRUE)
#     r lower.ci upper.ci
# 1 -0.0693   -0.136 -0.00132

## NN2, 4h
ont_NN2_4h_noNA <- subset(ont_NN2_4h_noNA, read_length > 50)
ont_NN2_8h_noNA <- subset(ont_NN2_8h_noNA, read_length > 50)
ont_SG17M_4h_noNA <- subset(ont_SG17M_4h_noNA, read_length > 50)
ont_SG17M_8h_noNA <- subset(ont_SG17M_8h_noNA, read_length > 50)

nn2_4h_kw <- kruskal.test(log10(ont_NN2_4h_noNA$read_length), g=ont_NN2_4h_noNA$replicate)
nn2_4h_e2 <- epsilonSquared(log10(ont_NN2_4h_noNA$read_length), g=ont_NN2_4h_noNA$replicate, ci=FALSE)

## NN2, 8h
nn2_8h_kw <- kruskal.test(log10(ont_NN2_8h_noNA$read_length), g=ont_NN2_8h_noNA$replicate)
nn2_8h_e2 <- epsilonSquared(log10(ont_NN2_8h_noNA$read_length), g=ont_NN2_8h_noNA$replicate, ci=FALSE)

## NN2, 4h-8h
ont_NN2_4h_8h_noNA <- data.frame(rbind(ont_NN2_4h_noNA, ont_NN2_8h_noNA))
ont_NN2_4h_8h_noNA$merged <- paste(ont_NN2_4h_8h_noNA$mix, ont_NN2_4h_8h_noNA$replicate, sep="_")
NN2_4h_8h_kw <- kruskal.test(log10(ont_NN2_4h_8h_noNA$read_length), g=ont_NN2_4h_8h_noNA$merged)
NN2_4h_8h_e2 <- epsilonSquared(log10(ont_NN2_4h_8h_noNA$read_length), g=ont_NN2_4h_8h_noNA$merged, ci=FALSE)

## SG17M, 4h
sg17m_4h_kw <- kruskal.test(log10(ont_SG17M_4h_noNA$read_length), g=ont_SG17M_4h_noNA$replicate)
sg17m_4h_e2 <- epsilonSquared(log10(ont_SG17M_4h_noNA$read_length), g=ont_SG17M_4h_noNA$replicate, ci=FALSE)

## SG17m, 8h
sg17m_8h_kw <- kruskal.test(log10(ont_SG17M_8h_noNA$read_length), g=ont_SG17M_8h_noNA$replicate)
sg17m_8h_e2 <- epsilonSquared(log10(ont_SG17M_8h_noNA$read_length), g=ont_SG17M_8h_noNA$replicate, ci=FALSE)

## SG17M, 4h-8h
ont_SG17M_4h_8h_noNA <- data.frame(rbind(ont_SG17M_4h_noNA, ont_SG17M_8h_noNA))
ont_SG17M_4h_8h_noNA$merged <- paste(ont_SG17M_4h_8h_noNA$mix, ont_SG17M_4h_8h_noNA$replicate, sep="_")
SG17M_4h_8h_kw <- kruskal.test(log10(ont_SG17M_4h_8h_noNA$read_length), g=ont_SG17M_4h_8h_noNA$merged)
SG17M_4h_8h_e2 <- epsilonSquared(log10(ont_SG17M_4h_8h_noNA$read_length), g=ont_SG17M_4h_8h_noNA$merged, ci=FALSE)

## NN2 vs SG17M, 4h
nn2_sg17m_4h <- data.frame(rbind(ont_NN2_4h_noNA, ont_SG17M_NN2_4h_noNA))
nn2_sg17m_4h$merged <- paste(nn2_sg17m_4h$mix, nn2_sg17m_4h$replicate, sep="_")
nn2_sg17m_4h_kw <- kruskal.test(log10(nn2_sg17m_4h$read_length), g=nn2_sg17m_4h$merged)
nn2_sg17m_4h_e2 <- epsilonSquared(log10(nn2_sg17m_4h$read_length), g=nn2_sg17m_4h$merged, ci=FALSE)

## NN2 vs SG17M, 8h 
nn2_sg17m_8h <- data.frame(rbind(ont_NN2_8h_noNA, ont_SG17M_8h_noNA))
nn2_sg17m_8h$merged <- paste(nn2_sg17m_8h$mix, nn2_sg17m_8h$replicate, sep="_")
nn2_sg17m_8h_kw <- kruskal.test(log10(nn2_sg17m_8h$read_length), g=nn2_sg17m_8h$merged)
nn2_sg17m_8h_e2 <- epsilonSquared(log10(nn2_sg17m_8h$read_length), g=nn2_sg17m_8h$merged, ci=FALSE)

# plots
nn2_4h_plot_length <- ggplot(ont_NN2_4h_noNA, aes(x=replicate, y=log10(read_length))) + geom_jitter(width = 0.1) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_pubr(border = TRUE, base_size=11) +  
  xlab("NN2-4h replicates") + ylab("Read length (log10 scale)")

nn2_8h_plot_length <- ggplot(ont_NN2_8h_noNA, aes(x=replicate, y=log10(read_length))) +  geom_jitter(width = 0.1) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_pubr(border = TRUE, base_size=11) +  
  xlab("NN2-8h replicates") + ylab(" ") + theme(axis.text.y = element_blank())

sg17m_4h_plot_length <-ggplot(ont_SG17M_4h_noNA, aes(x=replicate, y=log10(read_length))) +  geom_jitter(width = 0.1) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_pubr(border = TRUE, base_size=11) + 
  xlab("SG17M-4h replicates") + ylab(" ") + theme(axis.text.y = element_blank())

sg18m_8h_plot_length <- ggplot(ont_SG17M_8h_noNA, aes(x=replicate, y=log10(read_length))) +  geom_jitter(width = 0.1) + geom_violin() +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_pubr(border = TRUE, base_size=11) + 
  xlab("SG17M-8h replicates") + ylab(" ") + theme(axis.text.y = element_blank())

plot_length_merge <- ggarrange(nn2_4h_plot_length, nn2_8h_plot_length, 
                               sg17m_4h_plot_length, sg18m_8h_plot_length, labels=c("A", "B", "C", "D"), nrow=1)

sample_names <- c("NN2-4h", "NN2-8h", "NN2-4h-8h", "SG17M-4h", "SG17M-8h", "SG17M-4h-8h", "NN2-SG17M (4h)", "NN2-SG17M (8h)")

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

length_stats_df <- data.frame(cbind(kruskal_statistic, 
                                    kruskal_pvalue, e2_values, e2_lower_ci, e2_upper_ci))
rownames(length_stats_df) <- sample_names
length_stats_df$e2_values <- round(length_stats_df$e2_values,2)
length_stats_df$e2_lower_ci <- round(length_stats_df$e2_lower_ci, 2)
length_stats_df$e2_upper_ci <- round(length_stats_df$e2_upper_ci, 2)

stats_plot <- ggplot(length_stats_df) + geom_point(aes(x=e2_values, y=rownames(length_stats_df)), size=3, colour="red") +
  geom_point(aes(x=e2_lower_ci, y=rownames(length_stats_df)), size=1, colour="black") +
  geom_point(aes(x=e2_upper_ci, y=rownames(length_stats_df)), size=1, colour="black") +
  theme_pubr(border=TRUE, legend = "none", base_size=11) + xlim(-1,1) + geom_vline(xintercept=-0.2) + geom_vline(xintercept=0.2) +
  xlab("Epsilon-squared effect size") + ylab("")

plot_length1 <- ggarrange(nn2_4h_plot_length, nn2_8h_plot_length, 
                          sg17m_4h_plot_length, sg18m_8h_plot_length, labels = c("A", "B", "C", "D"), nrow = 1)

merged_plot <- ggarrange(plot_length1, stats_plot, labels=c("A","E"), nrow=2, heights = c(1,0.5))
merged_plot
#ggsave(merged_plot, file="save_figures/plot_length2.tif", device="tiff", dpi=300, units="cm", width=9.77, height=7.64))


