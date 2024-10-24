library(stringr)
library(RFLPtools)
library(dplyr)
library(GenomicFeatures)
library(tximport)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)
library(pheatmap)
library(pvclust)
library(RMariaDB)
library(edgeR)
source("scripts/Usefull_functions.R")
Strata_pos <- data.frame(strata = c("Strata1", "Strata2", "Strata3"), 
                         start = c(8*10e5,40.5*10e5, 71.2*10e5), end = c(19.55*10e5, 44.35*10e5, 72.95*10e5))


############ DEKETE 
txdb_coords <- makeTxDbFromGFF("indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_genomic.gff.gz")


k <- keys(txdb_coords, keytype = "GENEID")
chrs <- read.table('indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_assembly_report.txt', sep = '\t')[,c(3,7)]
chrs <- chrs[which(chrs$V3 != "na" & chrs$V3 != "MT"),]
txdf <- AnnotationDbi::select(txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
coords_genes <- as.data.frame(genes(txdb_coords, c("TXCHROM", "GENEID")))
coords_genesZ <- filter(coords_genes, TXCHROM == chrs$V7[chrs$V3 == "Z"])

#### SALMON DATA READING IN ----
sg_temp <- filter(coords_genesZ, TXCHROM == chrs$V7[chrs$V3 == "Z"]])


inv_temp <- filter(coords_genes, TXCHROM == chrs$V7[chrs$V3 == "Z"]) %>% 
  filter((start < 6.5e6 | end < 6.5e6) | (start > 7.01e7| end > 7.01e7))


fixAin <- 181
fixBin <- 182
fixAout <- 0
fixBout <- 0

OI_genes <- unique(inv_temp$GENEID)
Igenes <- coords_genesZ$GENEID[!(coords_genesZ$GENEID %in% OI_genes)] %>% unlist()
genes_in <- 885
genes_out <- 199
fixed_diff_in <- 363
fixed_diff_out <- 0
total <- genes_in + genes_out
exp_outA <- (199/885)*181
exp_outB <- (199/885)*182

Odata <- matrix(c(fixAin, fixBin, fixAout, fixBout), nrow = 2, ncol = 2)
colnames(Odata) <- c("In", "Out")
rownames(Odata) <- c("A", "B")

Edata <- matrix(c(fixAin, fixBin, exp_outA, exp_outB), nrow = 2, ncol = 2)
colnames(Edata) <- c("In", "Out")
rownames(Edata) <- c("A", "B")


q <- (((Odata - Edata)^2 )/Edata) %>% 
  sum()

pvalue <- pchisq(q, df = 1, lower.tail = FALSE)

