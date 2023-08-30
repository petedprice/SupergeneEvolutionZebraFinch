library(ggpubr)
library(ggplot2)
library(tidyverse)
library(edgeR)
library(GenomicFeatures)

load("outdata/DEG_analyses/ABvsAAet.RData")
DNDS <- read.table("indata/work/ZChromGenesMKT.csv", sep = ",", 
                   header = T, row.names = 1)
txdb_coords <- makeTxDbFromGFF("indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_genomic.gff.gz")

##### GETTING EXON LENGTH PER GENE as better proxy for CDS than gene length -----
exons.list.per.gene <- exonsBy(txdb_coords,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons.list.per.gene))))
exonic.gene.sizes$gene <- rownames(exonic.gene.sizes)
colnames(exonic.gene.sizes) <- c('exon_length', 'Gene_Name')

DNDS <- merge(DNDS, exonic.gene.sizes, all.y = T)
DNDS$length <- abs(DNDS$FeatureStart - DNDS$FeatureEnd)
DNDS$Dtotal <- DNDS$D_nonsyn + DNDS$D_syn
rownames(DNDS) <- DNDS$Gene_Name
DNDS[is.na(DNDS$Dtotal) == T,c("D_nonsyn", "D_syn", "Dtotal")] <- 0

DNDS$Dn_cont <- DNDS$D_nonsyn/DNDS$exon_length
DNDS$Ds_cont <- DNDS$D_syn/DNDS$exon_length
DNDS$Dtotal_cont <- DNDS$Dtotal/DNDS$exon_length
DNDS$dif_or_not <- "NOT"
DNDS$dif_or_not[DNDS$Dtotal > 0] <- "DIF"

DNDS$dnds <- (DNDS$D_nonsyn)/(DNDS$D_syn + 1)


#### DEG #########
DEG_DNDS <- DNDS%>% 
  merge(ABAet$table[ABAet$table$Chromosome == "Z",], by=0, all = T) %>% 
  filter(is.na(fdr) == F) 
DEG_DNDS$logFC[DEG_DNDS$logFC == 0 ]  <- min(DEG_DNDS$logFC * 0.1)
write.table(DEG_DNDS, "outdata/FINAL_STATS/DEG_DNDS.txt")


######### ASE ###########
ASE <- read.table("outdata/ASE_DGE_data_males.csv")
ASE <- ASE[duplicated(ASE) == F,] %>% 
  filter(contig == "NC_044241.2" & stop > 6.5e6 & start < 70.1e6)

ASE_DNDS <- merge(ASE, DNDS, by.x = 'name', by.y = 'Gene_Name')
ASE_DNDS[ASE_DNDS$DTotal == 0,c('Dn_cont', 'Ds_cont', 'Dtotal_cont', 'dnds')] <- 0.0001
ASE_DNDS$log2_aFC[ASE_DNDS$log2_aFC == 0 ]  <- min(ASE_DNDS$log2_aFC * 0.1)

write.table(ASE_DNDS, "outdata/FINAL_STATS/ASE_DNDS.txt")

