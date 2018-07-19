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
library(dplyr)

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

samples <- unlist(lapply(kal_files, strsplit, split='\\.'))
samples <- sort(samples, decreasing = TRUE)[1:78]

## create a R matrix containing sample names and conditions from a text file
s2c <- read.table(file.path(base_dir, "design.txt"), header = TRUE, stringsAsFactors=FALSE)

## add a column called "kal_dirs" containing paths to the data
s2c <- dplyr::mutate(s2c, path = kal_dirs)

## Get common gene names for transcripts

## this section queries Ensemble online database for gene names associated with transcript IDs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
t2g <- getBM(attributes = c("ensembl_transcript_id", "transcript_version",
                            "ensembl_gene_id", "external_gene_name", "description",
                            "transcript_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE)
plot_pca(so, color_by = 'type', text_labels = TRUE)
plot_pca(so, color_by = 'id', text_labels = FALSE)

plot_loadings(so, pc_input=2)
plot_bootstrap(so, 'ENSMUST00000082402.1', color_by = 'type')
plot_bootstrap(so, 'ENSMUST00000042235.14', color_by = 'type')
