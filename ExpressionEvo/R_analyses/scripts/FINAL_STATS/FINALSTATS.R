library(ggplot2)
library(ggpubr)
library(tidyverse)

############################# FINAL STATS --------------------------------------

## DIFFERENTIAL GENE EXPRESSION ANALYSES 
# PER GENE BASIS USING abs LF
#PLOTS WITH RANK CORRELATIONS DO NOT FILTER FOR DN_CONT < 0
DEG_data <- read.table("outdata/FINAL_STATS/DEG_DNDS.txt")


### SPEARMAN RANK ANALYSIS ----
## DN
cor.test(abs(DEG_data$logFC), DEG_data$Dn_cont, method = 'spearman', exact = FALSE )
## DS
cor.test(abs(DEG_data$logFC), DEG_data$Ds_cont, method = 'spearman', exact = FALSE )

#GLM: SAME DIRECTION, DIFFERENT END RESULTS 
#model1 <- glm(abs(logFC) ~ Dn_cont, DEG_data, family = gaussian(link = 'log')) %>% 
# summary()

#model2 <- glm(abs(logFC) ~ Ds_cont, DEG_data, family = gaussian(link = 'log')) %>% 
#  summary()


DEGplots <- list()
DEGplots[[1]] <- ggplot(DEG_data, aes(x = Dn_cont, y = abs(logFC))) + geom_point() + 
  stat_cor(method = 'spearman') +
  labs(x = '', y = '')
DEGplots[[2]] <- ggplot(DEG_data, aes(x = Ds_cont, y = abs(logFC))) + geom_point() + 
  stat_cor(method = 'spearman') +
  labs(x = '', y = '')


DEGplots_arranged <- ggarrange(plotlist = DEGplots, nrow = 1, labels = c("A", "B"))

ggsave("plots/DEGplots_SPEARMAN_arranged.pdf", DEGplots_arranged, height = 5, width = 10)
## ALLELE SPECIFIC EXPRESSION EXPRESSION ANALYSIS ----
#PER GENE BASIS USING abs LF
#PLOTS WITH RANK CORRELATIONS DO NOT FILTER FOR DN_CONT < 0
# REPORT MEASURES FOR HOM BIRDS and HET BIRDS
ASE_DATA <- read.table("outdata/FINAL_STATS/ASE_DNDS.txt")
#### SPEARMANS RANK FOR HETEROZYGOTES 
hetASE <- ASE_DATA[ASE_DATA$genotype == "AB",]
cor.test(abs(hetASE$log2_aFC), hetASE$Ds_cont, method = 'spearman', exact = FALSE)
cor.test(abs(hetASE$log2_aFC), hetASE$Dn_cont, method = 'spearman', exact = FALSE)

### SPEARMNS RANK FOR HOMOZYGOTES 
homASE <- ASE_DATA[ASE_DATA$genotype == "AA",]
cor.test(abs(homASE$log2_aFC), homASE$Ds_cont, method = 'spearman', exact = FALSE)
cor.test(abs(homASE$log2_aFC), homASE$Dn_cont, method = 'spearman', exact = FALSE)

ASEplots <- list()
ASEplots[[1]] <- ggplot(hetASE, aes(x = (Ds_cont), y = (abs(log2_aFC)))) + 
  geom_point(alpha = 0.1) + 
  stat_cor(method = 'spearman') + labs(x = "", y = "") +
  ylim(0,11) + theme_bw()
ASEplots[[2]] <- ggplot(hetASE, aes(x = Dn_cont, y = abs(log2_aFC))) + 
  geom_point(alpha = 0.1) + 
  stat_cor(method = 'spearman') + labs(x = "", y = "")+ylim(0,11)+ theme_bw()

ASEplots[[3]] <- ggplot(homASE, aes(x = (Ds_cont), y = (abs(log2_aFC)))) + 
  geom_point(alpha = 0.1) + 
  stat_cor(method = 'spearman') + 
  labs(x = bquote(~D[s]~ '/ gene length'), y = "") +ylim(0,11)+ theme_bw()

ASEplots[[4]] <- ggplot(homASE, aes(x = Dn_cont, y = abs(log2_aFC))) + 
  geom_point(alpha = 0.1) + 
  stat_cor(method = 'spearman') + 
  labs(x = bquote(~D[n]~ '/ gene length'), y = "") +ylim(0,11)+ theme_bw()


ASEplots_arranged <- ggarrange(plotlist = ASEplots, labels = c("a", "b", "c", "d"))
ASEplots_arranged <- annotate_figure(ASEplots_arranged, left = 
                                       text_grob(bquote('|'~log[2]~ '(fold change) |'), rot = 90))
ggsave("plots/ASEplots_SPEARMAN_arranged.pdf", ASEplots_arranged, height = 6, width = 7)
system("open plots/ASEplots_SPEARMAN_arranged.pdf")


## Check candidates of Viitaniemi paper of ASE 
candidates <- c("DMGDH", "FREM1", "ENSTGUG00000021693")
DEG_data[DEG_data$Gene_Name %in% candidates,] %>% View()
ASE_DATA[ASE_DATA$name %in% candidates,] %>% View()


VitENb <- read.table("indata/Viitaniemi_candidates/btg1.4_pro_blastout.tsv")
VitEnE <- read.table("indata/Viitaniemi_candidates/ENSTGUP00000034343.1_blastout.tsv")
colnames(VitENb) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                     "qend", "sstart", "send", "evalue", "bitscore")
colnames(VitEnE) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                      "qend", "sstart", "send", "evalue", "bitscore")

# filter both for top 5 e values 
VitENb <- VitENb[order(VitENb$evalue),]
VitEnE <- VitEnE[order(VitEnE$evalue),]

VitENb <- VitENb[1:5,]
VitEnE <- VitEnE[1:5,]

# get the gene names
intersect(VitENb$sseqid, VitEnE$qseqid)
View(VitENb)
View(VitEnE)
