# Analysis of Weigel scRNAseq data from 78 A375 parental single cells 
## Date: 07.19.2018
## Author: Michael Chimenti
## Organism: hg38
## Aligners: kallisto
## Design: 78 single-cell libraries from fluidigm C1
## Contrasts:  none
## Reps: none
## Updated based on: https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html

##########
## Imports
##########

#source("http://bioconductor.org/biocLite.R")
#biocLite("COMBINE-lab/wasabi")       #install wasabi the first time

#pseudoalignment analysis

library(sleuth)
library(cowplot)

#annotation
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library(tidyverse)

#pathway
library(pathview)
library(gage)
library(gageData)

setwd("~/iihg/sc_rna_seq/weigel/C1_A375_scRNA/")

#########################
## kallisto > sleuth
#########################

base_dir <- "~/iihg/sc_rna_seq/weigel/C1_A375_scRNA"

kal_files <- list.files(path = '.', pattern = '*output')
kal_dirs <- file.path(base_dir, kal_files)

sample <- unlist(lapply(kal_files, strsplit, split='\\.'))
sample <- sort(sample, decreasing = TRUE)[1:78]

## create a R matrix containing sample names and conditions from a text file
s2c <- as.data.frame(sample)

## add a column called "kal_dirs" containing paths to the data
s2c <- dplyr::mutate(s2c, path = sort(kal_dirs, decreasing = TRUE))

## add a column called cell number 
s2c <- dplyr::mutate(s2c, cell = as.factor(paste0("cell", 1:78)))

## Get common gene names for transcripts

## this section queries Ensemble online database for gene names associated with transcript IDs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
t2g <- getBM(attributes = c("ensembl_transcript_id", "transcript_version",
                            "ensembl_gene_id", "external_gene_name", "description",
                            "transcript_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)



## testing 
s2c_small <- head(s2c)
so_small <- sleuth_prep(s2c_small, target_mapping = t2g, extra_bootstrap_summary = TRUE)


plot_pca(so_small, color_by = 'cell', text_labels = FALSE)
filter(kal_tab_small, target_id == 'ENST00000245479.2')   ## there is no SOX9 expression, so it is filtered out

mygenes <- c("SOX9", "GPX3","IGFBP5","AIM1","SEMA3B","EPHB3","S100A2","FOXD3","NRN1","PRDM7","YPEL4","RASGRP1")
my_t2g <- filter(t2g, ext_gene %in% mygenes)
my_t2g <- mutate(my_t2g, target_id_vers = paste0(target_id,'.', transcript_version))

png(filename = 'test.png', width = 1200, height = 1500, res = 150)
plot_transcript_heatmap(so_small, my_t2g$target_id_vers, 
                        units = 'tpm', cluster_transcripts = TRUE, trans = 'log', offset = 1,
                        annotation_cols = NULL, 
                        color_high = 'yellow', color_mid = 'magenta', color_low = 'black')
dev.off()





## big 'so' with all 78 samples 
so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE)

png(filename = 'A375_single_cell_exp_heatmap.png', width = 1200, height = 1500, res = 150)
p <- plot_transcript_heatmap(so, my_t2g$target_id_vers, 
                        units = 'tpm', 
                        cluster_transcripts = TRUE, 
                        trans = 'log', 
                        offset = 1,
                        annotation_cols = NULL,
                        color_high = 'yellow', color_mid = 'magenta', color_low = 'black',
                        main = "Gene expression in 78 A375 parental single cells by Fluidigm C1 (log TPM)",
                        fontsize = 8)
dev.off()

kal_tab <- kallisto_table(so, include_covariates = FALSE)
kal_tab <- select(kal_tab, c("target_id", "sample", "tpm"))

kal_tab_my_genes <- filter(kal_tab, target_id %in% my_t2g$target_id_vers)
kal_tab_my_genes <- mutate(kal_tab_my_genes, ext_gene = t2g$ext_gene)

png(filename = 'A375_parental_scRNA_gene_boxplot.png', width = 1200, height = 1200, res = 150)
p <- ggplot(kal_tab_my_genes, aes(x = as.factor(target_id), y = as.numeric(tpm))) + 
       geom_boxplot(outlier.colour = 'red', outlier.size = 1)
p <- p + coord_flip(ylim = c(0,200))
p <- p + theme(axis.text.y = element_text(size = 8, face = 'bold'))
p <- p + ylab("transcripts per million (tpm)") + xlab("transcript isoform") + ggtitle("Boxplot of TPM expression in 78 A375 single cells: Fluidigm C1 data")
print(p)
dev.off()

write.csv(file = "transcript_to_gene_tab.csv", x = my_t2g)
