ASE<- read.table("outdata/ASE_DGE_data_males.csv")
ASE <- ASE[duplicated(ASE) == F,] %>% 
  filter(contig == "NC_044241.2")

ASE_DNDS <- merge(ASE, DNDS, by.x = 'name', by.y = 'Gene_Name', all.x = T)
ASE_DNDS[is.na(ASE_DNDS$Dtotal) == T,c('Dn_cont', 'Ds_cont', 'Dtotal_cont', 'dnds')] <- 0
ASE_DNDS$dif_or_not[is.na(ASE_DNDS$dif_or_not) == T] <- "NOT"
write.table(ASE_DNDS, "outdata/work/JON_ASEDATA.txt")

temp <- ASE_DNDS %>% filter(totalCount > 0)

lmer <- glmer(cbind(aCount, totalCount) ~ genotype * Ds_cont +
             (1|sample) + (1|name), temp,
           family = "binomial")
summary(lmer)
