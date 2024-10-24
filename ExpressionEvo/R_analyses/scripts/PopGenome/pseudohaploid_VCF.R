## CONVERT DIPLOID TO HAPLOID ##
library(PopGenome)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(vcfR)
library(stringr)


### SAMPLE INFO ###
sample_info_ENI <- read.table("indata/genomeR_data/filereport_read_run_PRJEB10586_tsv.txt", sep = 
                                '\t', header = T) %>% 
  select(run_accession, sample_accession, scientific_name)

sample_info_jake <- read.csv("indata/genomeR_data/jake_sample_info.csv")

sample_info <- sample_info_ENI %>% 
  merge(sample_info_jake, by.x = 'sample_accession', by.y = 'Sample.Accession') %>% 
  filter(Karyotype %in% c("A", "AA", "AB", "BB", "B"))

###################

### READING IN VCF ###
vcf_whole <- read.vcfR("indata//Genome_VCFs/fixploidy/subset_nocoding.vcf")

### GETTING SEX SAMPLES ###
Females <- intersect(sample_info$run_accession[sample_info$Sex == "F"], colnames(vcf_whole@gt))
Males <- intersect(sample_info$run_accession[sample_info$Sex == "M"], colnames(vcf_whole@gt))

### GETTING GT FUNCTION ###
get_gt <- function(samp, gt){
  gs <- str_split(gt[,samp], ":", simplify = T)[,1]
  df <- data.frame(genotype = gs)
  colnames(df) <- samp
  return(df)
}

### MALE GENOTYPES ###
just_gt_Males <- lapply(Males, get_gt, gt = vcf_whole@gt) %>% 
  bind_cols() 

### FEMALE GENOTYPES ###
just_gt_Females <- lapply(Females, get_gt, gt = vcf_whole@gt) %>% 
  bind_cols()


### FUNCTION TO CONVERT DIPLOID MALES TO PSEUDOHAPLOID ###
hap_the_dip <- function(samp, gt){
  gs <- str_split(gt[,samp], "/|\\|", simplify = T) %>% 
    as.data.frame()
  
  pseudo_dip <- lapply(gs, function(x)(return(paste0(x, "|", x)))) %>% 
    bind_cols()
  
  colnames(pseudo_dip) <- paste0(samp, "_hap", 1:2)
  return(pseudo_dip)
}

### GETTING THE PSEUDOHAPLOID MALE GT MATRIX ###
pseudo_haploid_males_gts <- lapply(Males, hap_the_dip, gt = just_gt) %>% 
  bind_cols()
pseudo_haploid_males_gts[pseudo_haploid_males_gts == "0|0"] <- "0|0:10,0:10"
pseudo_haploid_males_gts[pseudo_haploid_males_gts == "1|1"] <- "1|1:0,10:10"
pseudo_haploid_males_gts[pseudo_haploid_males_gts == ".|."] <- ".|.:0,0:0"


### FUNCTION FOR GETTING PSEUDOHAPLOID FEMALES 
hap_the_hap <- function(samp, gt){
  gs <- str_split(gt[,samp], "/|\\|", simplify = T) %>% 
    as.data.frame()
  
  pseudo_dip <- lapply(gs, function(x)(return(paste0(x, "|", x)))) %>% 
    bind_cols()
  pseudo_dip[pseudo_dip[,1] != pseudo_dip[,2],] <- ".|."
  pseudo_dip <- pseudo_dip[,1]
  colnames(pseudo_dip) <- paste0(samp, "_hapF")
  return(pseudo_dip)
}

### GETTING THE PSEUDOHAPLOID FEMALE GT MATRIX ###
pseudo_haploid_female_gts <- lapply(Females, hap_the_hap, gt = just_gt_Females) %>% 
  bind_cols()
pseudo_haploid_female_gts[pseudo_haploid_female_gts == "0|0"] <- "0|0:10,0:10"
pseudo_haploid_female_gts[pseudo_haploid_female_gts == "1|1"] <- "1|1:0,10:10"
pseudo_haploid_female_gts[pseudo_haploid_female_gts == ".|."] <- ".|.:0,0:0"


### MERGING ###
pseudo_haploid_gts <- cbind(data.frame(FORMAT = rep("GT:AD:DP"))) %>% 
                              cbind(pseudo_haploid_males_gts, pseudo_haploid_female_gts) %>% 
  as.matrix()


pseudohaploid_vcf <- vcf_whole
pseudohaploid_vcf@gt <- pseudo_haploid_gts


write.vcf(pseudohaploid_vcf, "indata/Genome_VCFs/pseudohaploid_Z.vcf.gz")



#### TEST IN POPGENOME ###
genome <- readVCF("indata/Genome_VCFs/pseudohaploid_Z.vcf.gz", numcols=10000, tid="NC_044241.2", frompos=1, topos=254285)
samples <- get.individuals(genome)[[1]]


A <- c("ERR1013156","ERR1013157","ERR1013158","ERR1013159","ERR1013160","ERR1013161",
       "ERR1013165","ERR1013166","ERR1013170","ERR1013171","ERR1013172","ERR1013173",
       "ERR1013174")
B <- c("ERR1013162","ERR1013164","ERR1013167","ERR1013169","ERR1013176","ERR1013177",
       "ERR1013178")
populations <- list(A, B)
outgroup <- c("LFT")

genome <- set.populations(genome, populations, diploid = T)
genome <- set.outgroup(genome, outgroup, diploid = T)




