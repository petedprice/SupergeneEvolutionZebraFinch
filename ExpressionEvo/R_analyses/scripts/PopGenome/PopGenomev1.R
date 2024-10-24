# Install and load necessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

library(PopGenome)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(vcfR)

filter = "NO"

#### GET SAMPLE INFORMATION ---
sample_info_ENI <- read.table("indata/genomeR_data/filereport_read_run_PRJEB10586_tsv.txt", sep = 
                            '\t', header = T) %>% 
  select(run_accession, sample_accession, scientific_name)

sample_info_jake <- read.csv("indata/genomeR_data/jake_sample_info.csv")

sample_info <- sample_info_ENI %>% 
  merge(sample_info_jake, by.x = 'sample_accession', by.y = 'Sample.Accession') %>% 
  filter(Karyotype %in% c("A", "AA", "AB", "BB", "B"))

write.table(sample_info[,c(2,6)], "indata/Genome_VCFs/fixploidy/samples.txt", col.names = F, quote = F, row.names = F)

######### INITIAL Read the VCF file and filter for coding regions -----------
if (filter == "YES"){
  vcf <- read.vcfR("indata/Genome_VCFs/Singhal_bTG1_4_Z_plusLTF_filtered_snps.vcf.gz")
  gff <- import.gff("indata/genome_files/NC_044241.2.gff.gz")
  coding_regions <- gff[gff$type == "CDS"]
  
  
  filter_vcf_func <- function(snps, start, width){
    end <- start + width-1
    snps <- as.numeric(snps)
    rm_positions <- snps[which(snps > start & snps < end)]
    return(rm_positions)
  }
  
  
  snps <- vcf@fix[,2] %>% as.numeric()
  
  rm_snps <- coding_regions@ranges %>% as.data.frame() %>% 
    apply(., 1, function(x)(filter_vcf_func(snps, x[1], width = x[3]))) %>% unlist() %>% 
    unique()
  
  rms_indx <- which(snps %in% rm_snps)
  vcf@gt <- vcf@gt[-rms_indx,]
  vcf@fix <- vcf@fix[-rms_indx,]
  
  write.vcf(vcf, "indata/Genome_VCFs/Z_nocoding.vcf.gz")
}
################################################################

#conda activate del 
# gunzip Z_nocoding.vcf.gz
#bgzip -c Z_nocoding.vcf > Z_nocoding.vcf.gz
#tabix -p vcf Z_nocoding.vcf.gz

## Filter VCF for depth < 1
#conda activate vcftools
#vcftools --gzvcf Z_nocoding.vcf.gz --recode --out Z_nocoding_filt --minQ 30 --minDP 10
#vcftools --gzvcf subset_nocoding.vcf --recode --out subset_nocoding_filt --minQ 30 --minDP 10 --minGQ 30

#bgzip subset_nocoding_filt.recode.vcf > subset_nocoding_filt.recode.vcf.gz
#tabix -p vcf subset_nocoding_filt.recode.vcf.gz

#bcftools convert --hapsample --vcf-ids subset_nocoding_filt.recode.vcf.gz -o haplotypes
#bcftools convert --haplegendsample2vcf haplotypes

#bcftools convert -haplegendsample subset_nocoding_filt.recode.vcf.gz
#bcftools convert --haplegendsample2vcf aplegendsample > pseudodiploid.vcf


########### READ INTO POPGENOME FOR NEUTRALITY AND DIVERSITY STATISTICS ---------
genome <- readVCF("indata/Genome_VCFs/subset_nocoding_filt.recode.vcf.gz", numcols=1000, tid="NC_044241.2", frompos=1, topos=254285)

samples <- get.individuals(genome)[[1]]
AA_birds <- c()
BB_birds <- c()
outgroup <- c()

# Define populations
#populations <- list(c(), c() c())
A <- c("ERR1013156","ERR1013157","ERR1013158","ERR1013159","ERR1013160","ERR1013161",
       "ERR1013165","ERR1013166","ERR1013170","ERR1013171","ERR1013172","ERR1013173",
       "ERR1013174")
B <- c("ERR1013162","ERR1013164","ERR1013167","ERR1013169","ERR1013176","ERR1013177",
       "ERR1013178")
populations <- list(A, B)
outgroup <- c("LFT")

genome <- set.populations(genome, populations, diploid = T)
genome <- set.outgroup(genome, outgroup, diploid = T)

# Define sliding windows and calculate Fst
window_size <- 10000
step_size <- 10000

genome <- sliding.window.transform(genome, width=window_size, jump=step_size, type=2)
genome <- F_ST.stats(genome)

get.F_ST(genome)

# Extract Fst values
fst_values <- genome@region.stats@nucleotide.F_ST
genome@region.stats@Pop_FSTH


########## DELETE -----
vcf <- read.vcfR("indata/Genome_VCFs/subset_nocoding.vcf.gz")
vcf_females <- vcf
vcf_females@gt <- cbind(vcf@gt[,1], vcf_females@gt[,colnames(vcf@gt) %in% sample_info$run_accession[sample_info$Sex == "F"]])
colnames(vcf_females@gt)[1] <- 'FORMAT'
samples <- colnames(vcf_females@gt[,-1])
write.table(data.frame(samples = samples, sex = rep("F", length(samples))), "indata/Genome_VCFs/fixploidy/Fsamples.txt", quote = F, col.names = F, row.names = F)
write.vcf(vcf_females, 'indata/Genome_VCFs/fixploidy/Females.vcf.gz')



vcf <- read.vcfR("indata//Genome_VCFs/fixploidy/ploidy.vcf")
keep <- colnames(vcf@gt) %in% sample_info$run_accession
keep[1] <-  TRUE
vcf@gt <- vcf@gt[,keep]
write.vcf(vcf, "indata/Genome_VCFs/fixploidy/hapdip.vcf.gz")

genome <- readVCF("indata/Genome_VCFs/fixploidy/hapdip.vcf.gz", numcols=1000, tid="NC_044241.2", frompos=1, topos=254285)







############################
genome <- readVCF("indata/Genome_VCFs/subset_nocoding_filt.recode.vcf.gz", numcols=1000, tid="NC_044241.2", frompos=1, topos=254285)

# Define populations
#populations <- list(c(), c() c())
MAs <- sample_info[sample_info$Karyotype %in% c("AA") & sample_info$Sex == "M",]$run_accession
MAs2 <- paste0(MAs, ".2")
FAs <- sample_info[sample_info$Karyotype %in% c("A") & sample_info$Sex == "F",]$run_accession

MBs <- sample_info[sample_info$Karyotype %in% c("BB", "B") & sample_info$Sex == "M",]$run_accession
MBs2 <- paste0(MBs, ".2")
FBs <- sample_info[sample_info$Karyotype %in% c("BB", "B") & sample_info$Sex == "F",]$run_accession

A <- c(MAs, MAs2, FAs)
B <- c(MBs, MBs2, FBs)
A_ss <- sample(A, length(B), replace = F)

populations <- list(A_ss, B)
outgroup <- c("LFT", "LFT.2")

slide <- sliding.window.transform(genome,width=1000,1000, type=2)

#calculate pi diversity
slide <- diversity.stats(slide, pi=T)

Pi <-  as.data.frame(slide@Pi)

slide <- F_ST.stats(slide, mode="nucleotide")
pairwise.FST <- t(slide@nuc.F_ST.pairwise)


position <- as.data.frame(seq(10000, 113930000, 10000))

all.stats_Chr <- bind_cols(position, TajimasD, FayWuH, ZengE, pairwise.FST, Pi) 

write.csv(all.stats_Chr, "all_stats_Chr1.csv", row.names = F)


