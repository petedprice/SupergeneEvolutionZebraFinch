### libraries and functions ----
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

sample_info <- read.table("outdata/project_data_gen.csv", sep = ',', header = T)

#READING IN GENE COORDINATES
txdb_coords <- makeTxDbFromGFF("indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_genomic.gff.gz")
k <- keys(txdb_coords, keytype = "GENEID")
chrs <- read.table('indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_assembly_report.txt', sep = '\t')[,c(3,7)]
chrs <- chrs[which(chrs$V3 != "na" & chrs$V3 != "MT"),]
txdf <- AnnotationDbi::select(txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
coords_genes <- as.data.frame(genes(txdb_coords, c("TXCHROM", "GENEID")))
coords_genesZ <- filter(coords_genes, TXCHROM == chrs$V7[chrs$V3 == "Z"])

#### SALMON DATA READING IN ----
files <- Sys.glob(file.path(paste("indata/salmon_quant/bTG1.4/", sample_info$code, "*/quant.sf", sep = "")))
names(files) <- sample_info$sample
#### RUNNINGE TRANSCRIPT TO GENE 
tx2gene <- txdf[,c(2,1)]
tx2gene$TXNAME <- str_split(paste("rna-", tx2gene$TXNAME, sep = ""), "[.]", simplify = TRUE)[,1]
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE, 
                ignoreAfterBar = TRUE,
                countsFromAbundance = "no")

sg_temp <- unique(txdf$GENEID[txdf$TXCHROM == chrs$V7[chrs$V3 == "Z"]])
ag_temp <- unique(txdf$GENEID[txdf$TXCHROM != chrs$V7[chrs$V3 == "Z"]])
inv_temp <- filter(coords_genes, TXCHROM == chrs$V7[chrs$V3 == "Z"]) %>% 
  filter((start < 6.5e5 | end < 6.5e5) | (start > 7.01e7| end > 7.01e7))

autosome_genes <- intersect(ag_temp, rownames(txi$counts))
tip_genes <- intersect(inv_temp$GENEID, rownames(txi$counts)) %>% unlist()
sgs <- intersect(sg_temp, rownames(txi$counts))
sex_genes <- sgs[sgs %in% tip_genes == FALSE]
salmon_DE_data <- data.frame("MF_lfc" = rep(NA, length(sex_genes)), "MF_fdrp" = rep(NA, length(sex_genes)), "MF_bias" = rep("UNBIAS", length(sex_genes)),
                             "M_AB_A_lfc" = rep(NA, length(sex_genes)), "M_AB_A_fdrp" = rep(NA, length(sex_genes)), "M_AB_A_bias" = rep("UNBIAS", length(sex_genes)),
                             "F_A_B_lfc" = rep(NA, length(sex_genes)), "F_A_B_fdrp" = rep(NA, length(sex_genes)), "F_A_B_bias" = rep("UNBIAS", length(sex_genes)))
rownames(salmon_DE_data) <- sex_genes


#Male vs Female                            
group = sample_info$sex
design <- model.matrix(~group)
inds = sample_info$sample
y <- DGEList(counts = txi$counts[,inds], genes = txi$length[,inds], group = group)
y <- calcNormFactors(y)
temp_rpkm <- rpkm(y, gene.length = y$genes, log = T)
Mrs <- rowSums(temp_rpkm[,sample_info$sample[sample_info$sex == "M"]] >= 2)
Frs <- rowSums(temp_rpkm[,sample_info$sample[sample_info$sex == "F"]] >= 2)

Mrm <- Mrs > (0.5 * length(which(sample_info$sex == "M")))
Frm <- Frs > (0.5 * length(which(sample_info$sex == "M")))
keep <- Mrm == TRUE | Frm == TRUE
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
MFrpkm <- rpkm(y, gene.length = y$genes, log = T)

et <- exactTest(y, pair = c("M", "F")) #Positive values are Male bias (first in pair stated)
et$table$fdr <- p.adjust(et$table$PValue, "fdr")
et$table$Chromosome <- "NA"
et$table[sex_genes,]$Chromosome <- "Z"
et$table[autosome_genes,]$Chromosome <- "aut"
et$table$Bias <- "UNBIAS"
et$table$Bias[which(et$table$fdr < 0.05 & et$table$logFC < -1)] <- "Male_bias"
et$table$Bias[which(et$table$fdr < 0.05 & et$table$logFC > 1)] <- "Female_bias"

et$table %>% 
  count(Chromosome, Bias)

set <-(et$table[sex_genes,])
salmon_DE_data[,c(1:2)] <- set[,c(1,4)]
MB <- rownames(set)[which(set$fdr < 0.05 & set$logFC < -1)]
FB <- rownames(set)[which(set$fdr < 0.05 & set$logFC > 1)]
salmon_DE_data[MB,3] <- "Male_Bias" 
salmon_DE_data[FB,3] <- "Female_Bias" 





#AB vs AA
group = sample_info$genotype[sample_info$sex == "M"]
design <- model.matrix(~group)
inds = sample_info$sample[sample_info$sex == "M"]
y <- DGEList(counts = txi$counts[,inds], genes = txi$length[,inds], group = group)
y <- calcNormFactors(y)
temp_rpkm <- rpkm(y, gene.length = y$genes, log = T)
ABrs <- rowSums(temp_rpkm[,sample_info$sample[sample_info$genotype == "AB"]] >= 2)
AArs <- rowSums(temp_rpkm[,sample_info$sample[sample_info$genotype == "AA"]] >= 2)
ABrm <- ABrs > 0.5 * length(which(sample_info$genotype == "AB"))
AArm <- AArs > 0.5 * length(which(sample_info$genotype == "AA"))
keep <- ABrm == TRUE | AArm == TRUE

y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
ABArpkm <- rpkm(y, gene.length = y$genes, log = T)
ABAet <- exactTest(y, pair = c("AB", "AA")) 
ABAet$table$fdr <- p.adjust(ABAet$table$PValue, "fdr")
ABAet$table$Chromosome <- NA
ABAet$table[sex_genes,]$Chromosome <- "Z"
ABAet$table[autosome_genes,]$Chromosome <- "aut"
ABAet$table[tip_genes,]$Chromosome <- "Ztip"

ABAet$table$Bias <- "Unbias"
ABAet$table$Bias[which(ABAet$table$fdr < 0.05 & ABAet$table$logFC < -1)] <- "AB_bias"
ABAet$table$Bias[which(ABAet$table$fdr < 0.05 & ABAet$table$logFC > 1)] <- "AA_bias"

AB_bias <- rownames(ABAet$table[which(ABAet$table$Bias == "AB_bias"),])



ABAA_summary <- ABAet$table %>% 
  filter(is.na(logFC) == F) %>% 
  count(Chromosome, Bias)
ABAA_summary$total <- 0
ABAA_summary$total[ABAA_summary$Chromosome == "Z"] <- 
  sum(ABAA_summary$n[ABAA_summary$Chromosome == "Z"])
ABAA_summary$total[ABAA_summary$Chromosome == "Ztip"] <- 
  sum(ABAA_summary$n[ABAA_summary$Chromosome == "Ztip"])
ABAA_summary$total[ABAA_summary$Chromosome == "aut"] <- 
  sum(ABAA_summary$n[ABAA_summary$Chromosome == "aut"])
ABAA_summary$pntg <- (ABAA_summary$n/ABAA_summary$total) * 100
ABAA_summary$comp <- "ABvsAA"

ABAset <-ABAet$table[sex_genes,]
ABAaet <-ABAet$table[autosome_genes,]

chisq_table <- ABAA_summary %>% 
  group_by(Chromosome) %>% 
  summarise(DEG = sum(n[Bias != "Unbias"]), 
            UNBIAS = sum(n[Bias == "Unbias"]),
            total = total[1], 
            pntg = sum(pntg[Bias != "Unbias"]))

chisq_result <- chisq.test(chisq_table[c(1,3),c(2:3)])
chisq_result$expected
chisq_result$observed

salmon_DE_data[,c(4,5)] <- ABAset[,c(1,4)]
AB_B <- rownames(ABAset)[which(ABAset$fdr < 0.05 & ABAset$logFC < -1)]
AA_B <- rownames(ABAset)[which(ABAset$fdr < 0.05 & ABAset$logFC > 1)]
salmon_DE_data[AB_B,6] <- "AB_Bias" 
salmon_DE_data[AA_B,6] <- "A_Bias" 

ABa_B <- rownames(ABAaet)[which(ABAaet$fdr < 0.05 & ABAaet$logFC < -1)]
AAa_B <- rownames(ABAaet)[which(ABAaet$fdr < 0.05 & ABAaet$logFC > 1)]


write.table(intersect(rownames(ABArpkm), autosome_genes), "outdata/DEG_analyses/DGE_AB_AA_Expressed_autosomal_genes.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(intersect(rownames(ABArpkm), sex_genes), "outdata/DEG_analyses/DGE_AB_AA_Expressed_sex_genes.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

for (f in c(0,0.5, 1)){
  write.table(intersect(rownames(filter(ABAet$table, fdr < 0.05 & abs(logFC) >= f)), autosome_genes), paste("outdata/DEG_analyses/DGE_AB_AA_Bias_lfc", f, "_autosomal_genes.txt", sep = ""), 
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(intersect(rownames(filter(ABAet$table, fdr < 0.05 & abs(logFC) >= f)), sex_genes), paste("outdata/DEG_analyses/DGE_AB_AA_Bias_lfc", f, "_sex_genes.txt", sep = ""), 
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}



#A vs B (females)
group = sample_info$genotype[sample_info$sex == "F"]
design <- model.matrix(~group)
inds = sample_info$sample[sample_info$sex == "F"]
y <- DGEList(counts = txi$counts[,inds], genes = txi$length[,inds], group = group)
y <- calcNormFactors(y)
temp_rpkm <- rpkm(y, gene.length = y$genes, log = T)
Ars <- rowSums(temp_rpkm[,sample_info$sample[sample_info$genotype == "A"]] >= 2)
Brs <- rowSums(temp_rpkm[,sample_info$sample[sample_info$genotype == "B"]] >= 2)
Arm <- Ars > 0.5 * length(which(sample_info$genotype == "A"))
Brm <- Brs > 0.5 * length(which(sample_info$genotype == "B"))
keep <- Arm == TRUE | Brm == TRUE
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
ABrpkm <- rpkm(y, gene.length = y$genes, log = T)

ABet <- exactTest(y, pair = c("A", "B")) 
ABet$table$fdr <- p.adjust(ABet$table$PValue, "fdr")

ABet$table$Chromosome <- NA
ABet$table[sex_genes,]$Chromosome <- "Z"
ABet$table[autosome_genes,]$Chromosome <- "aut"
ABet$table[tip_genes,]$Chromosome <- "Ztip"

ABet$table$Bias <- "Unbias"
ABet$table$Bias[which(ABet$table$fdr < 0.05 & ABet$table$logFC < -1)] <- "A_bias"
ABet$table$Bias[which(ABet$table$fdr < 0.05 & ABet$table$logFC > 1)] <- "B_bias"

AB_summary <- ABet$table %>% 
  filter(is.na(logFC) == F) %>% 
  count(Chromosome, Bias)
AB_summary$total <- 0
AB_summary$total[AB_summary$Chromosome == "Z"] <- 
  sum(AB_summary$n[AB_summary$Chromosome == "Z"])
AB_summary$total[AB_summary$Chromosome == "Ztip"] <- 
  sum(AB_summary$n[AB_summary$Chromosome == "Ztip"])
AB_summary$total[AB_summary$Chromosome == "aut"] <- 
  sum(AB_summary$n[AB_summary$Chromosome == "aut"])
AB_summary$pntg <- (AB_summary$n/AB_summary$total) * 100
AB_summary$comp <- "AvsB"
chisq_table <- AB_summary %>% 
  group_by(Chromosome) %>% 
  summarise(DEG = sum(n[Bias != "Unbias"]), 
            UNBIAS = sum(n[Bias == "Unbias"]),
            total = total[1], 
            pntg = sum(pntg[Bias != "Unbias"]))

chisq_result <- chisq.test(chisq_table[c(1,3),c(2:3)])
chisq_result$expected
chisq_result$observed

ABset <-ABet$table[sex_genes,]

salmon_DE_data[,c(7,8)] <- ABset[,c(1,4)]
A_B <- rownames(ABset)[which(ABset$fdr < 0.05 & ABset$logFC < -1)]
B_B <- rownames(ABset)[which(ABset$fdr < 0.05 & ABset$logFC > 1)]

salmon_DE_data[A_B,9] <- "A_Bias" 
salmon_DE_data[B_B,9] <- "B_Bias" 

salmon_DE_data$GENE <- rownames(salmon_DE_data)


DEG_manuscript_summary <- rbind(AB_summary, ABAA_summary) %>% 
  filter(Chromosome != "Ztip")
write.table(DEG_manuscript_summary, "outdata/DEG.csv")
save(file = "outdata/DEG_analyses/ABvsAAet.RData", ABAet)
