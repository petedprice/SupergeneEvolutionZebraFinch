load("outdata/DEG_analyses/ABvsAAet.RData")
DNDS <- read.table("indata/work/ZChromGenesMKT.csv", sep = ",", 
                   header = T, row.names = 1)

head(DNDS)
DNDS$P1dif <- DNDS$P1_nonsyn + DNDS$P1_syn
DNDS$P2dif <- DNDS$P2_nonsyn + DNDS$P2_syn
DNDS$TotalDif <- DNDS$D_nonsyn + DNDS$D_syn
DNDS$P1_score <- DNDS$D_nonsyn/(DNDS$P1_nonsyn + DNDS$P1_syn +1)
DNDS$P2_score <- DNDS$D_nonsyn/(DNDS$P2_nonsyn + DNDS$P2_syn +1)

rownames(DNDS) <- DNDS$Gene_Name
DNDS2 <- DNDS%>% 
  merge(ABAet$table, by=0) %>% 
  filter(is.na(fdr) == F) 

DNDS2 %>% ggplot(aes(x = D_syn + D_nonsyn, y = logFC, colour = fdr)) + 
  geom_point() + stat_cor()


DNDS2 %>% pivot_longer(
  c("D_nonsyn", "D_syn", "TotalDif"), names_to ='dif_type', values_to = "difs") %>% 
  filter(Bias != "Unbias") %>% 
  ggplot(aes(x = abs(logFC), y =  difs, colour = dif_type)) + 
  geom_point() + stat_cor() + 
  geom_smooth(method = "lm", se = FALSE)

DNDS2 %>% 
  pivot_longer(c("P1_score", "P2_score"), names_to = "score_type", values_to = "score") %>% 
  ggplot(aes(x = abs(logFC), y = log(score), colour = score_type)) + 
  geom_point() + stat_cor() + 
  geom_smooth(method = "lm", se = FALSE)

hist(log(DNDS2$P1_score))
