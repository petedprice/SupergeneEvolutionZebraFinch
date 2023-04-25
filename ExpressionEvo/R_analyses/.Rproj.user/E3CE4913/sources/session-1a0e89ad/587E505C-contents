library(gridExtra)
library(grid)
library(ggpubr)
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
library(lme4)


data_path = "indata/phaser/"
siglevel <- data.frame("stars" = c("ns", "*", "**", "***", "****"), "siglevel" = c(1, 0.05, 0.01, 0.001, 0.0001), "other" = 
                         c("non-sig", "p<0.05", "p<0.01", "p<0.001", "p<0.0001"))
sigcolors <- setNames(c("red", "blue", "green", "purple", "orange"), siglevel$other)

sample_info <- read.table("outdata/project_data_gen.csv", header = T, sep = ",")

males <- sample_info$sample[sample_info$sex == "M"]
names(males) <- str_split(sample_info$code[sample_info$sex == "M"], 
                          "_", simplify = T)[,3]
plots <- list()

# GENE LIST TO LOOK AT 
si_ss <- sample_info[sample_info$sex =="M",]
files <- Sys.glob(file.path(paste("indata/salmon_quant/bTG1.4/", si_ss$code, "*/quant.sf", sep = "")))
names(files) <- si_ss$sample


#### RUNNING TRANSCRIPT TO GENE 
txdb_coords <- makeTxDbFromGFF("indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_genomic.gff.gz")
k <- keys(txdb_coords, keytype = "GENEID")
chrs <- read.table('indata/genome_files/GCF_003957565.2_bTaeGut1.4.pri_assembly_report.txt', sep = '\t')[,c(3,7)]
chrs <- chrs[which(chrs$V3 != "na" & chrs$V3 != "MT"),]
txdf <- AnnotationDbi::select(txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
tx2gene <- txdf[,c(2,1)]
tx2gene$TXNAME <- str_split(paste("rna-", tx2gene$TXNAME, sep = ""), "[.]", simplify = TRUE)[,1]
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE, 
                ignoreAfterBar = TRUE,
                countsFromAbundance = "no")
group = si_ss$genotype
design <- model.matrix(~group)
inds = si_ss$sample
y <- DGEList(counts = txi$counts[,inds], genes = txi$length[,inds], group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

rpkm <- rpkm(y, gene.length = y$genes, log = T)
include_genes <- rpkm > 2
et <- exactTest(y, pair = c("AA", "AB")) #Positive values are Male bias (first in pair stated)
et$table$fdr <- p.adjust(et$table$PValue, "fdr")
et$table$bias <- "UNBIAS" 
AB_bias <- rownames(et$table)[which(et$table$fdr < 0.05 & et$table$logFC < -1)]
AA_bias <- rownames(et$table)[which(et$table$fdr < 0.05 & et$table$logFC > 1)]
et$table[AB_bias,"bias"] <- "AB_BIAS"
et$table[AA_bias, "bias"] <- "AA_BIAS"
et$table$name = rownames(et$table)
colnames(et$table)[1:5] <- paste("DGE_", colnames(et$table)[1:5], sep = "")
######


contigs <- paste("NC_", str_split(list.files(
  path = data_path), "_", simplify = TRUE)[,2]
  , sep = "") %>% unique()
#contigs = "NC_044241.2"
chisq_table <- data.frame(Contig = NA, Focal = NA, Others = NA, lfc = NA, test = NA, pvalue = NA, 
                          sample = NA, rZ = c(NA, NA), rA = c(NA, NA), sig = c(NA, NA))

o=0.01
focal_genes <- list()
other <- list()
whole_data = list.files(path = data_path, full.names = TRUE) %>% 
  lapply(., read.delim) %>% 
  do.call("rbind", .) %>% 
  left_join(et$table, ("name" = "name"))
whole_data$sample <- str_split(whole_data$bam, "_", simplify = TRUE)[,5]



for (ct in contigs){
  print(ct)
  for (s in names(males)){
    
    
    df = whole_data %>% 
      filter(sample == s)
    # select only genes with sufficient coverage
    #cov8 = subset(df, df$totalCount >= 8)
    code <- males[s]
    Zinc <- intersect(df$name, rownames(include_genes)[include_genes[,code]])
    cov8 <- filter(df, name %in% Zinc)
    cov8 = subset(cov8, cov8$totalCount >= 8)
    cov8$log2_aFC <- log(cov8$aCount + 1, 2) - log(cov8$bCount + 1, 2)
    cov8$binom_p = apply(cov8[,c("aCount","bCount")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)$p.value)
    cov8$binom_q = p.adjust(cov8$binom_p, method = "fdr")
    cov8$ASE_bias <- "UNBIAS"
    cov8$ASE_bias[cov8$binom_q < 0.05 & abs(cov8$log2_aFC) >= 1.22] <- "BIAS"
    if (s == names(males)[1] & ct == contigs[1]){
      save_ASE <- cov8
    } else {
      save_ASE <- rbind(save_ASE, cov8)
    }
    focal_genes[[ct]][[s]] <- filter(cov8, binom_q < 0.05 & abs(log2_aFC) >= 1.22 & contig == ct)$name
    other[[ct]][[s]] <- filter(cov8, binom_q >= 0.05 & contig == ct)$name
    
    
    
    
    lf_table <- data.frame("low" = c(0,Inf), "medium" = c(0.5, Inf), 
                           "high" = c(1.22, Inf), "very high" = c(1.5, Inf), "Super_high" = c(2, Inf))
    
    
    for (lf in 1:ncol(lf_table)){
      cst <- data.frame(Contig = c(ct, ct), Focal = c(NA,NA), Others = c(NA,NA), lfc = c(NA,NA), 
                        pvalue = c(NA, NA), rZ = c(NA, NA), rA = c(NA, NA), sig = c(NA, NA))
      cst[1,2] <- nrow(filter(cov8, binom_q <= 0.05 & abs(log2_aFC) >= lf_table[1,lf] &
                                abs(log2_aFC) <= lf_table[2,lf] &
                                contig == ct))
      cst[2,2] <- nrow(filter(cov8, contig == ct)) - cst[1,2]
      cst[1,3] <- nrow(filter(cov8, binom_q < 0.05 & abs(log2_aFC) >= lf_table[1,lf] &
                                abs(log2_aFC) <= lf_table[2,lf] & 
                                contig != ct))
      cst[2,3] <- nrow(filter(cov8, contig != ct))- cst[1,3]
      
      cst$lfc = colnames(lf_table)[lf]
      cst$test <- c("ASE", "NON")
      chisq <- chisq.test(cst[1:2,2:3], simulate.p.value = TRUE, B = 1000000)
      cst$pvalue <- chisq$p.value
      cst$sample <- males[s]
      cst$rZ <- chisq$residuals[1]
      cst$rA <- chisq$residuals[2]
      cst$sig<- unlist(lapply(cst$pvalue, function(x)(return(siglevel$other[tail(which(siglevel$siglevel >=x), 1)]))))
      chisq_table <- rbind(chisq_table, cst)
      
      
    }
  }
  
  #chisq_table <- chisq_table[-which(is.na(chisq_table$Contig)),]
  chisq_table$sig <- factor(chisq_table$sig, levels = siglevel$other)
  
  pdf(paste("plots/phaser/males_PHASER_residuals_binned_", ct, ".pdf", sep = ""), height = 5, width = 7)
  if (ct == "NC_044241.2"){
    g <- ggplot(filter(chisq_table, test == "ASE" & 
                         Contig == ct), aes(x = lfc, y = rZ, shape = sample, colour = sig)) + geom_point(size = 3) + 
      theme_bw() + scale_colour_manual(values = c("grey68", "hotpink", "deeppink3", "darkslateblue", "orange")) +
      labs(colour = "Significance", shape = "Sample", x = bquote(~log[2]~ 'Fold Change'), y = "Residuals") +
      scale_shape_manual(values = c(16,17,15, 18, 19)) + scale_x_discrete(limits = colnames(lf_table))
    print(g)
    dev.off()
  } else {
    g <- ggplot(filter(chisq_table, test == "ASE" & lfc != "above_0.5" & 
                         Contig == ct), aes(x = lfc, y = rZ, shape = sample, colour = sig)) + geom_point(size = 3) + 
      theme_bw() + 
      labs(colour = "Significance", shape = "Sample", x = "", y = "", title = ct) +
      scale_shape_manual(values = c(16,17,15,18, 19)) + scale_x_discrete(limits = colnames(lf_table)) +
      scale_color_manual(values = sigcolors)
    print(g)
    dev.off()
  }
  plots[[ct]] <- g
}
chisq_table$genotype <- NA
chisq_table[chisq_table$sample %in% sample_info$sample[sample_info$genotype == "AB"],]$genotype <- "AB"
chisq_table[chisq_table$sample %in% sample_info$sample[sample_info$genotype == "AA"],]$genotype <- "AA"


codes <- str_split(sample_info$code, "_", simplify = TRUE)[,3]

save_ASE$genotype <- NA
save_ASE[save_ASE$sample %in% codes[sample_info$genotype == "AB"],]$genotype <- "AB"
save_ASE[save_ASE$sample %in% codes[sample_info$genotype == "AA"],]$genotype <- "AA"
write.table(save_ASE, "outdata/ASE_DGE_data_males.csv")





fig <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow =2, ncol = 2, common.legend = TRUE, legend = "right", labels = c("A", "B", "C", "D"))
pdf("plots/males_arranged_autosomes_phaser_binned.pdf", height = 8, width = 11)
annotate_figure(fig,
                bottom = textGrob(expression(paste("log"["2"], "Fold Change"))),
                left = text_grob("Residuals", rot = 90))
dev.off()

plot_list <- list()
for (lf in colnames(lf_table)){ 
  #pdf(paste("plots/manuscript/AvsZ/residuals_aut_Z_lfcgt", lf_table[1,lf], ".pdf", sep = ""), height = 3.5, width = 5.5)
  g <- ggplot(filter(chisq_table, lfc == lf), 
              aes(x =Contig, y = rZ, shape = genotype, colour = sig)) + geom_point() +
    theme_bw() + scale_colour_manual(values = c("grey68", "hotpink", "deeppink3", "darkslateblue", "orange")) + 
    labs(colour = "Significance", shape = "Sample", x = "", y = "ASE Enrichment (Standardised Residuals)") +
    geom_text(aes(fontface = 1, label = "Z", x = 'NC_044241.2' , y = -Inf), colour = 'black', vjust = 2, show.legend = FALSE) +
    geom_text(aes(fontface = 1, label = "Autosomes", x = 6.5 , y = -Inf), colour = 'black',vjust = 2, show.legend = FALSE) +
    coord_cartesian(clip = 'off') + geom_vline(xintercept = 12.5)+ geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = c(-1.96,1.96), linetype = 'dotted', colour = 'maroon') + 
    theme(axis.text.x=element_blank())
  plot_list[[lf]] <- g + ggtitle(paste("lf threshold = ", lf_table[1,lf]))
  ggsave(paste("plots/residuals_aut_Z_lfcgt", lf_table[1,lf], ".pdf", sep = ""), g, height = 3.5, width = 5.5)
}

g2 <- ggarrange(plotlist = plot_list)
ggsave("plots/males_residuals_aut_Z_comp.pdf", g2, height = 10, width = 20)



#### MAKING CHISQUARED TEST TABLES ####
table <- filter(chisq_table, Contig == 'NC_044241.2' & lfc == 'high')
tidy_table <- table %>% 
  group_by(test) %>% 
  summarise(
    Z = mean(Focal), 
    Autosomes = mean(Others),
    n = n())

ZvsA_chisq <- matrix(nrow = 0, ncol = 8) %>% as.data.frame()
for (s in unique(table$sample)){
  table_ss <- filter(table, sample == s)
  chisq_man <- chisq.test(table_ss[,2:3], simulate.p.value = TRUE, B = 1000000)
  name = paste("outdata/chisq_tables/hets", s, "_AvsZ_chisq_p=", round(chisq_man$p.value, 10), "_", sep = "")
  
  temp_out <- matrix(nrow =2, ncol = 10, data = c(rep(s,2), 
                                                  rep(table_ss$genotype[1], 2), 
                                                  c("Aut", "Z"), 
                                                  t(chisq_man$observed)[c(2,1),], 
                                                  round(t(chisq_man$expected)[c(2,1),]), 
                                                  rep(chisq_man$p.value, 2),
                                                  rep(chisq_man$statistic, 2), 
                                                  rep(sum(chisq_man$observed), 2))) %>% as.data.frame()
  
  
  colnames(temp_out) <- c("Sample", "genotype", "Chromosome", "Obs_ASE", "Obs_nASE", "Ex_ASE", "Ex_nASE", "MCMCpvalue", "Chi", "N")
  ZvsA_chisq <- rbind(ZvsA_chisq, temp_out)
}

write.table(ZvsA_chisq, "outdata/chisq_tables/ZvsA_chisq.tsv", quote = FALSE, row.names = FALSE, sep = '\t')


### CHISQURED USING OVERLAPPING GENES 
ase_tab_hets <- data.frame(ASE = NA, NON = NA, contig = NA)
for (c in unique(names(focal_genes))){
  ASE <- intersect(focal_genes[[c]][[1]], intersect  (focal_genes[[c]][[2]], focal_genes[[c]][[3]]))
  NON <- intersect(other[[c]][[1]], intersect(other[[c]][[2]], other[[c]][[3]]))
  ase_tab_hets <- rbind(ase_tab_hets, c(length(ASE), length(NON), c))
}
overlap_ASE <- ASE
write.table(overlap_ASE, file = "outdata/Z_ASE_overlap.tsv", col.names = FALSE, 
            row.names = FALSE, quote = FALSE)
ase_tab_hets <- ase_tab_hets[-1,] 
ase_tab_hets[,1:2] <- as.numeric(unlist(ase_tab_hets[,1:2]))
write.table(ase_tab_hets, "outdata/chisq_tables/hets_ASE_overlap.tsv",
            quote = FALSE, row.names = FALSE, sep = "\t")

ase_tab_homs <- data.frame(ASE = NA, NON = NA, contig = NA)
for (c in unique(names(focal_genes))){
  ASE <- intersect(focal_genes[[c]][[4]], focal_genes[[c]][[5]])
  NON <- intersect(other[[c]][[4]], other[[c]][[5]])
  ase_tab_homs <- rbind(ase_tab_homs, c(length(ASE), length(NON), c))
}
ase_tab_homs <- ase_tab_homs[-1,] 
ase_tab_homs[,1:2] <- as.numeric(unlist(ase_tab_homs[,1:2]))
write.table(ase_tab_homs, "outdata//chisq_tables/homs_ASE_overlap.tsv",
            quote = FALSE, row.names = FALSE, sep = "\t")



ase_tab_hets_chimat <- matrix(data = c(tail(ase_tab_hets$ASE, 1), tail(ase_tab_hets$NON, 1),
                                       sum(ase_tab_hets$ASE[1:12]), sum(ase_tab_hets$NON[1:12])), 
                              nrow = 2, byrow = FALSE)
chisq_man <- chisq.test(ase_tab_hets_chimat, simulate.p.value = TRUE, B = 1000000)

############### FINISHING EDITING SOMWEHRE AROUND HERE!!!!

ZvsA_chisq_overlap <- matrix(nrow =2, ncol = 9, data = c(rep("AB", 2), 
                                                         c("Aut", "Z"), 
                                                         t(chisq_man$observed)[c(2,1),], 
                                                         round(t(chisq_man$expected)[c(2,1),]), 
                                                         rep(chisq_man$p.value, 2), 
                                                         rep(chisq_man$statistic, 2), 
                                                         rep(sum(chisq_man$observed), 2))) %>% as.data.frame()


ase_tab_homs_chimat <- matrix(data = c(tail(ase_tab_homs$ASE, 1), tail(ase_tab_homs$NON, 1),
                                       sum(ase_tab_homs$ASE[1:12]), sum(ase_tab_homs$NON[1:12])), 
                              nrow = 2, byrow = FALSE)
chisq_man <- chisq.test(ase_tab_homs_chimat, simulate.p.value = TRUE, B = 1000000)
ZvsA_chisq_overlap <- rbind(ZvsA_chisq_overlap, matrix(nrow =2, ncol = 9, data = c(rep("AA", 2), 
                                                                                   c("Aut", "Z"), 
                                                                                   t(chisq_man$observed)[c(2,1),], 
                                                                                   round(t(chisq_man$expected)[c(2,1),]), 
                                                                                   rep(chisq_man$p.value, 2),
                                                                                   rep(chisq_man$statistic, 2), 
                                                                                   rep(sum(chisq_man$observed), 2)))) %>% as.data.frame()

colnames(ZvsA_chisq_overlap) <- c("genotype", "Chromosome", "Obs_ASE", "Obs_nASE", "Ex_ASE", "Ex_nASE", "MCMCpvalue", "Chi", "N")

write.table(ZvsA_chisq_overlap, "outdata//chisq_tables/ZvsA_overlap_chisq.tsv", quote = FALSE, row.names = FALSE, sep = '\t')



Z_compare <- ZvsA_chisq_overlap[c(2,4), c(3,4)]
Z_compare[,1] <- as.numeric(Z_compare[,1])
Z_compare[,2] <- as.numeric(Z_compare[,2])
chisq_man <- chisq.test(Z_compare, simulate.p.value = TRUE, B = 1000000)
ZvsZ_chisq_overlap <- matrix(nrow =2, ncol = 9, data = c(c("AA","AB", "Z", "Z"), 
                                                         chisq_man$observed[c(2,1),], 
                                                         round(chisq_man$expected[c(2,1),]), 
                                                         rep(chisq_man$p.value, 2),
                                                         rep(chisq_man$statistic, 2),
                                                         rep(sum(chisq_man$observed), 2)))  %>% as.data.frame()


colnames(ZvsZ_chisq_overlap) <- c("genotype", "Chromosome", "Obs_ASE", "Obs_nASE", "Ex_ASE", "Ex_nASE", "MCMCpvalue", "Chi", "N")
write.table(ZvsZ_chisq_overlap, "outdata//chisq_tables/ZvsZ_overlap_chisq.tsv", quote = FALSE, row.names = FALSE, sep = '\t')



lmtable <- data.frame(BirdID = NA, ASE = NA, Genes = NA, Chromosome = NA, Genotype = NA)
lmtable <- ZvsA_chisq[,c(1:5)]
lmtable[,4:5] <- apply((lmtable[,4:5]), 2, as.numeric)
lmtable$Genes <- lmtable$Obs_ASE + lmtable$Obs_nASE
colnames(lmtable)[4:5] <- c("ASE", "NON")
lm <- glmer(cbind(ASE, Genes) ~ genotype  * Chromosome + (1|Sample), 
            data = lmtable, family = "binomial")
summary(lm)

overlap_ASE_df <- filter(save_ASE, name %in% overlap_ASE & genotype == "AB")
overlap_ASE_df <- overlap_ASE_df[which(duplicated(overlap_ASE_df) == FALSE),-c(9:11)]
sum_stats <- overlap_ASE_df %>% 
  group_by(name) %>% 
  summarise(meanlf = mean(log2_aFC), 
            #  n = n(),
            var = var(log2_aFC),
            pos = mean(start))

