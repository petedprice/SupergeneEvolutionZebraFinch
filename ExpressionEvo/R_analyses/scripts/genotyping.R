#### GENOTYPING BIRDS ----
library(tidyverse)

sample_info <- read.table('indata/project_data.csv', sep = ',', header = T)
source("scripts/Usefull_functions.R")
vcf <- "indata/Z.vcf.gz"

tmp_vcf<-readLines(vcf)
fg_vcf<-read.table(vcf, stringsAsFactors = FALSE)
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
names(fg_vcf)<-vcf_names

diog_snps <- read.table('indata/bTG1.4_chi0.9_ident.txt', header = TRUE) %>% 
  filter(Scaffold == "NC_044241")


ds2 <- merge(diog_snps, fg_vcf[,c(2,4,5)], by.x = 'LocusPOS', by.y = 'POS')

fg_vcf <- fg_vcf %>% 
  filter(POS %in% diog_snps$LocusPOS)
dim(fg_vcf)

sample <-colnames(fg_vcf)[10:ncol(fg_vcf)]
genotype_func <- function(s, fg_vcf, diog_snps){
  s2 = sample_info$sample[which(sample_info$code == gsub("_tgu", "tgu", s))]
  filtered_vcf <- vcf_filter(vcf = fg_vcf, samples = s, DP = 10, HETT = 0.90, GQ = 10)
  genotype <- str_split(filtered_vcf[,s], ":", simplify = T)[,1]
  ds3 <- cbind(ds2, genotype)
  AB <- which(ds3$genotype %in% c("0|1", "0/1"))
  AA <- c(which(ds3$genotype %in% c("0|0", "0/0") & ds3$A_ident == ds3$REF), 
          which(ds3$genotype %in% c("1|1", "1/1") & ds3$B_ident == ds3$REF)) 
  BB <- c(which(ds3$genotype %in% c("1|1", "1/1") & ds3$A_ident == ds3$REF),
          which(ds3$genotype %in% c("0|0", "0/0") & ds3$B_ident == ds3$REF))
  no_snps <- length(AB) + length(BB) + length(AA)
  AAp <- length(AA)/no_snps
  ABp <- length(AB)/no_snps
  BBp <- length(BB)/no_snps
  gen <- c("AA", "AB", "BB")[which.max(c(AAp, ABp, BBp))]
  out_df <- data.frame(AA = AAp, AB = ABp, BB = BBp, no_snps = no_snps, 
                       sample = s2, genotype = gen)
  return(out_df)
}   

out <- lapply(sample, genotype_func, fg_vcf = fg_vcf, diog_snps = diog_snps) %>% 
  bind_rows()

sample_info_gen <- merge(sample_info, out, by = 'sample')
sample_info_gen$genotype[sample_info_gen$sex == "F" & 
                           sample_info_gen$genotype == "AA"] <- "A"
sample_info_gen$genotype[sample_info_gen$sex == "F" & 
                           sample_info_gen$genotype == "BB"] <- "B"
sample_info_gen[,c(7:9)] <- round(sample_info_gen[,c(7:9)], 3)
write.table(sample_info_gen, "outdata/project_data_gen.csv", sep = ',', 
            row.names = F)
