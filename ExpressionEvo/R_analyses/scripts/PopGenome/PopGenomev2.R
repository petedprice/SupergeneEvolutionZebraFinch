########### READ INTO POPGENOME FOR NEUTRALITY AND DIVERSITY STATISTICS ---------
############################
genome <- readVCF("indata/Genome_VCFs/Z_nocoding.vcf.gz", numcols=10000, tid="NC_044241.2", frompos=1, topos=7539615)
genome <- readVCF("indata/Genome_VCFs/pseudohaploid_Z.vcf.gz", numcols=10000, tid="NC_044241.2", frompos=1, topos=7539615)

genome <- sliding.window.transform(genome,width=100000,10000, type=2)

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
genome <- set.populations(genome, populations, diploid = F)
genome <- set.outgroup(genome,new.outgroup=outgroup,diploid = F)

#calculate fst stats
genome <- F_ST.stats(genome, mode="nucleotide")
pairwise.FST <- data.frame(FST = t(genome@nuc.F_ST.pairwise))
position <- data.frame(position = seq(1, 7539615, (7539615/nrow(pairwise.FST))))
all.stats_Chr <- bind_cols(position, pairwise.FST) 
colnames(all.stats_Chr) <- c("Pos", "FST")

correct_haps <- all.stats_Chr %>% 
  ggplot(aes(x =Pos, y = FST)) + geom_line()
ggsave("plots/del_correct_hapsFST.pdf", correct_haps)



#############################################################################################

populations <- list(c(FAs, MAs[1:2]), c(FBs, MBs))
populations <- list(MAs, MBs[1])

outgroup <- "LFT"

genome <- set.populations(genome, populations, diploid = T)
genome <- set.outgroup(genome,new.outgroup=outgroup,diploid = T)

#calculate fst stats
genome <- F_ST.stats(genome, mode="nucleotide")
pairwise.FSTin <- t(genome@nuc.F_ST.pairwise)
position <- data.frame(position = seq(1, 7539615, (7539615/nrow(pairwise.FSTin))))
all.stats_Chrin <- bind_cols(position, pairwise.FSTin) 
colnames(all.stats_Chrin) <- c("Pos", "FST")

incorrect_haps <- all.stats_Chrin %>% 
  ggplot(aes(x =Pos, y = FST)) + geom_line()
ggsave("plots/del_incorrect_hapsFST.pdf", incorrect_haps)


