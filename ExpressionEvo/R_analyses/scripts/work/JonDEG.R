library(GenomicFeatures)
library(edgeR)
library(tidyverse)

load("outdata/DEG_analyses/ABvsAAet.RData")
DNDS <- read.table("indata/work/ZChromGenesMKT.csv", sep = ",", 
                   header = T, row.names = 1)
txdb_coords <- makeTxDbFromGFF("indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_genomic.gff.gz")
sample_info <- read.table("outdata/project_data_gen.csv", sep = ",", header = T)
##### GETTING EXON LENGTH PER GENE as better proxy for CDS than gene length -----
exons.list.per.gene <- exonsBy(txdb_coords,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons.list.per.gene))))
exonic.gene.sizes$gene <- rownames(exonic.gene.sizes)
colnames(exonic.gene.sizes) <- c('exon_length', 'Gene_Name')

DNDS <- merge(DNDS, exonic.gene.sizes)
DNDS$length <- abs(DNDS$FeatureStart - DNDS$FeatureEnd)
DNDS$Dtotal <- DNDS$D_nonsyn + DNDS$D_syn
rownames(DNDS) <- DNDS$Gene_Name

DNDS$Dn_cont <- DNDS$D_nonsyn/DNDS$exon_length
DNDS$Ds_cont <- DNDS$D_syn/DNDS$exon_length
DNDS$Dtotal_cont <- DNDS$Dtotal/DNDS$exon_length
DNDS$dif_or_not <- "DIF"
DNDS$dnds <- (DNDS$D_nonsyn)/(DNDS$D_syn + 1)

tidy_data <- ABArpkm %>% mutate(gene = rownames(ABArpkm)) %>% 
  filter(gene %in% c(sex_genes, tip_genes)) %>% 
  pivot_longer(colnames(ABArpkm), names_to = 'sample', values_to = 'reads')
tidy_data$genotype <- sapply(tidy_data$sample, function(x)(return(
  sample_info$genotype[sample_info$sample == x])))

tidy_data <- merge(tidy_data, DNDS, by.x = 'gene', by.y = 'Gene_Name', all.x = T) 
tidy_data[is.na(tidy_data$Dtotal) == T,c('D_nonsyn', 'D_syn', 'Dn_cont', 'Ds_cont', 'Dtotal_cont', 'dnds')] <- 0
tidy_data$dif_or_not[is.na(tidy_data$dif_or_not) == T] <- "NOT"
tidy_data$genotype <- as.factor(tidy_data$genotype)

tidy_data$position = "tip"
tidy_data$position[tidy_data$gene %in% sex_genes] <- "inv"
#tidy_data$reads <- round(tidy_data$reads)
write.table(tidy_data, "outdata/work/JON_DEGDATA.txt")

ggplot(tidy_data, aes(x = Ds_cont, y = reads)) + geom_point() + 
  geom_smooth(method = 'lm') + stat_cor()
ggplot(tidy_data, aes(x = Dn_cont, y = reads)) + geom_point() + 
  geom_smooth(method = 'lm') + stat_cor()



