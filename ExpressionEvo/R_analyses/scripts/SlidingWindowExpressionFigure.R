### libraries ----
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
library(emojifont)

Strata_pos <- data.frame(strata = c("Strata1", "Strata2", "Strata3"), 
                         start = c(8*10e5,40.5*10e5, 71.2*10e5), end = c(19.55*10e5, 44.35*10e5, 72.95*10e5))

sample_info <- read.table("data/project_data_gen.csv", sep = ',')
fd <- read.table("data/Analyses_outputs/FULL_DATA.csv")
load("data/Analyses_outputs/bTG1.4_salmon_txi.RData")


group = sample_info$genotype
design <- model.matrix(~group)
inds = sample_info$sample

y <- DGEList(counts = txi$counts[,inds], genes = txi$length[,inds], group = group)
y <- calcNormFactors(y)
allrpkm <- rpkm(y, gene.length = y$genes, log = T)

keep <- rep(FALSE, nrow(allrpkm))
for (s in unique(sample_info$genotype)){
  temprpkm <- allrpkm[sample_info$genotype == s]
  temp_keep <- rowSums(temprpkm > 2) > (0.5*ncol(temprpkm))
  keep[temp_keep] <- TRUE
  print(length(temp_keep))
}

y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
filt_allrpkm <- rpkm(y, gene.length = y$genes)

avrpkm <- matrix(ncol = length(unique(sample_info$genotype)), nrow = nrow(filt_allrpkm))
colnames(avrpkm) <- unique(sample_info$genotype)
rownames(avrpkm) <- rownames(filt_allrpkm)
for (g in colnames(avrpkm)){
  for (l in 1:nrow(avrpkm)){
    temp_exp <- filt_allrpkm[l,which(sample_info$genotype == g)]
    avrpkm[l,g] <- mean(as.numeric(temp_exp))
  }
}

gavrpkm <- avrpkm[rownames(avrpkm) %in% sex_genes,] %>% as.data.frame()
gavrpkm$GENE <- rownames(gavrpkm)
a1 <- merge(gavrpkm, fd[,c('start', 'end', "GENE")], by.x = 'GENE', by.y = 'GENE')
a2 <- filter(a1, is.na(a1$start) == FALSE)
a3 <- a2
for (g in unique(a3$GENE)){
  dups <- which(a3$GENE == g)
  if (length(dups) > 1){
    a3 <- a3[-dups[2:length(dups)],]
  }
}

ws <- 500000
#no overlap
a4 <- a3[order(a3$start),]
no_sw_average <- matrix(ncol = 7, nrow = round(max(a4$end)/ws)) %>% as.data.frame()
colnames(no_sw_average) <- c("start", "end", "B_genotypes", "AA_genotypes", "AB_genotypes", "A_genotypes", "nogenes")
no_sw_average$start <- seq(0,ws*(nrow(no_sw_average)-1), ws)
no_sw_average$end <- seq(ws+1,(ws+1)*nrow(no_sw_average), ws)

for (w in 1:nrow(no_sw_average)){
  ss <- a3[which(a3$end < no_sw_average[w,]$end & a3$start >=no_sw_average[w,]$start),]
  avs <- colMeans(ss[,c(2:5)])
  no_genes <- nrow(ss)
  no_sw_average[w,3:7] <- c(avs, no_genes)
}
no_sw_average <- no_sw_average[is.na(no_sw_average$AB_genotypes) == FALSE,]

ws2=5000000
overlap =100000
#overlap
ol_sw_average <- matrix(ncol = 7, nrow = round(max(a4$end)/overlap)) %>% as.data.frame()
colnames(ol_sw_average) <- c("start", "end", "mid", "B_genotypes", "AA_genotypes", "AB_genotypes", "A_genotypes")

ol_sw_average$mid <- seq(0,max(a4$end), overlap)
ol_sw_average$start <- ol_sw_average$mid - ws2/2
ol_sw_average$end <- ol_sw_average$mid + ws2/2

for (w in 1:nrow(ol_sw_average)){
  ss <- a3[which(a3$end < ol_sw_average[w,]$end & a3$start >=ol_sw_average[w,]$start),]
  avs <- colMeans(ss[,c(2:5)])
  no_genes <- nrow(ss)
  ol_sw_average[w,4:7] <- c(avs)
}
ol_sw_average <- ol_sw_average[is.na(ol_sw_average$AB_genotypes) == FALSE,]



tidy_windows <- ol_sw_average %>% 
  pivot_longer(c("B_genotypes", "AB_genotypes", "A_genotypes", "AA_genotypes"), names_to = "genotypes",
               values_to = "rpkm")
tidy_windows$mid <- (tidy_windows$start + (tidy_windows$end - tidy_windows$start))/1e+06
Strata_pos_Mb <- Strata_pos
Strata_pos_Mb[,2:3] <- Strata_pos_Mb[,2:3] /1e+06

ABtopDE <- filter(fd, M_AB_A_lfc < -1 & M_AB_A_fdrp < 0.05)[,c("start", "M_AB_A_bias")]
ABtopDE$rpkm = rep(2^6.5) 
colnames(ABtopDE) <- c("mid", "genotypes", "rpkm")
ABtopDE$mid <- ABtopDE$mid/1e+06
AAtopDE <- filter(fd, M_AB_A_lfc > 1 & M_AB_A_fdrp < 0.05)[,c("start", "M_AB_A_bias")]
AAtopDE$rpkm = rep(2^6.5)  
colnames(AAtopDE) <- c("mid", "genotypes", "rpkm")
AAtopDE$mid <- AAtopDE$mid/1e+06
tidy_windows[tidy_windows == "A_genotypes"] <- "A (\u2640)"
tidy_windows[tidy_windows == "B_genotypes"] <- "B (\u2640)"
tidy_windows[tidy_windows == "AA_genotypes"] <- "AA (\u2642)"
tidy_windows[tidy_windows == "AB_genotypes"] <- "AB (\u2642)"
tidy_windows$genotypes <- factor(tidy_windows$genotypes, 
                                 levels = c("AA (\u2642)", "AB (\u2642)",
                                            "A (\u2640)", "B (\u2640)"))

pdf("plots/Zexpression.pdf", height = 3, width = 7)
tidy_windows_plot <- 
  ggplot(arrange(tidy_windows, genotypes) , aes(x = mid, y = log(rpkm, base = 2), colour = genotypes)) + geom_line() +
  geom_vline(xintercept = c(74, 7), color = 'black', linetype = 'dotted') + 
  annotate(xmin = Strata_pos_Mb[1,2], xmax = Strata_pos_Mb[1,3], ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2, fill = "red") +
  annotate(xmin = Strata_pos_Mb[2,2], xmax = Strata_pos_Mb[2,3], ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2, fill = "red") +
  annotate(xmin = Strata_pos_Mb[3,2], xmax = Strata_pos_Mb[3,3], ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2, fill = "red") +
  annotate(xmin = 7, xmax = Strata_pos_Mb[1,2], ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.1) +
  annotate(xmin = Strata_pos_Mb[1,3], xmax = Strata_pos_Mb[2,2], ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.1) +
  annotate(xmin = Strata_pos_Mb[2,3], xmax = 74, ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.1) +
  theme_classic() + xlab("Position (Mb)")  +ylab(bquote(~Log[2]~ '(rpkm)')) +
  geom_point(data = ABtopDE, col = 'grey45', shape = 'star', size = 0.5) +
  geom_point(data = AAtopDE, col = 'green4', shape = 'star', size = 0.5) + 
  scale_colour_manual(values = brewer.pal(4, "Paired")[c(4,8,1,2)]) +
  labs(colour = "Genotype")

#scale_colour_manual(values = brewer.pal(4, "Set1")[c(2,4,3,1)]) +

tidy_windows_plot
dev.off()






