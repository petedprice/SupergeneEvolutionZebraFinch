install.packages("glmmTMB")
library(glmmTMB)
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

DNDS <- merge(DNDS, exonic.gene.sizes)
DNDS$length <- abs(DNDS$FeatureStart - DNDS$FeatureEnd)
DNDS$Dtotal <- DNDS$D_nonsyn + DNDS$D_syn
rownames(DNDS) <- DNDS$Gene_Name

DNDS$Dn_cont <- DNDS$D_nonsyn/DNDS$exon_length
DNDS$Ds_cont <- DNDS$D_syn/DNDS$exon_length
DNDS$Dtotal_cont <- DNDS$Dtotal/DNDS$exon_length
DNDS$dif_or_not <- "DIF"
DNDS$dnds <- (DNDS$D_nonsyn)/(DNDS$D_syn + 1)


#### DEG #########
DEG_DNDS <- DNDS%>% 
  merge(ABAet$table[ABAet$table$Chromosome == "Z",], by=0, all = T) %>% 
  filter(is.na(fdr) == F) 
DEG_DNDS$dif_or_not[is.na(DEG_DNDS$dif_or_not) == T] <- "NOT"
DEG_DNDS[is.na(DEG_DNDS$DTotal) == T,c('Dn_cont', 'Ds_cont', 'Dtotal_cont', 'dnds')] <- 0

write.table(DEG_DNDS, "outdata/DEG_DNDS.txt")

DEG_plots <- list()
DEG_plots[[1]] <- DEG_DNDS  %>% pivot_longer(
  c("Dn_cont", "Ds_cont", "Dtotal_cont"), names_to ='dif_type', values_to = "difs") %>% 
  #filter(abs(logFC) < 1 & difs > 0) %>% 
  ggplot(aes(y = log(abs(logFC)), x =  difs, colour = dif_type)) + 
  geom_point() + stat_cor() + 
  geom_smooth(method = "lm", se = T)

DEG_plots[[1]]

DEG_plots[[2]] <- DEG_DNDS %>% 
  ggplot(aes(x = abs(logFC), y = (D_nonsyn)/(D_syn))) +
  geom_point() + stat_cor() + 
  geom_smooth(method = "lm", se = T)

DEG_ploted <- ggarrange(plotlist = DEG_plots, nrow = 2)
ggsave("plots/DEGDNDS.pdf", DEG_ploted)


######### ASE ###########
ASE <- read.table("outdata/ASE_DGE_data_males.csv")
ASE <- ASE[duplicated(ASE) == F,] %>% 
  filter(contig == "NC_044241.2")

ASE_DNDS <- merge(ASE, DNDS, by.x = 'name', by.y = 'Gene_Name', all.x = T)
ASE_DNDS[is.na(ASE_DNDS$DTotal) == T,c('Dn_cont', 'Ds_cont', 'Dtotal_cont', 'dnds')] <- 0
ASE_DNDS$dif_or_not[is.na(ASE_DNDS$dif_or_not) == T] <- "NOT"

ASE_plots <- list()
ASE_plots[[1]] <- ASE_DNDS  %>% pivot_longer(
  c("Dn_cont", "Ds_cont", "Dtotal_cont"), names_to ='dif_type', values_to = "difs") %>% 
  filter(dif_type == 'D_nonsyn') %>% 
  ggplot(aes(x = (abs(log2_aFC)), y =  difs/exon_length, shape = sample)) + 
  geom_point() + stat_cor() + 
  geom_smooth(method = "lm", se = T)
  
ASE_plots[[2]] <- ASE_DNDS %>% 
  ggplot(aes(x = abs(log2_aFC), y = (D_nonsyn)/(D_syn), colour = sample)) +
  geom_point() + stat_cor() + 
  geom_smooth(method = "lm", se = T)  

ASE_plotted <- ggarrange(plotlist = ASE_plots, nrow = 2)
ggsave("plots/ASEDNDS.pdf", ASE_plotted)


### MODELS ----
ASE_DNDS %>% 
  glm(abs(log2_aFC) ~ Dtotal_cont + dif_or_not + genotype + sample, ., family = gaussian(link = 'log')) %>% 
  summary()

glm(abs(log2_aFC) ~ Dtotal_cont, ASE_DNDS, family = gaussian(link = "log"), start = 1)

ASE_DNDS %>% 
  lm(log(abs(log2_aFC) + 0.0001) ~ Dtotal_cont + dif_or_not + genotype + sample,.) %>% summary()

ASE_DNDS %>% 
  glm(abs(log2_aFC) ~ dnds + dif_or_not+ genotype + sample, ., family = gaussian(link = 'log', start = 1)) %>% 
  summary()

ASE_DNDS %>% 
  lm(log(abs(log2_aFC) + 0.0001) ~ dnds + dif_or_not+ genotype + 
       sample,.) %>% summary()






