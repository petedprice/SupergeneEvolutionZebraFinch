#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
`%notin%` <- Negate(`%in%`)

data=read.csv(args[1], sep = " ", header = FALSE)
#data=read.csv("data/nextflow/bfc126389ec32018da5e345d53fd7a/comp_read_counts.csv", sep = " ", header=F)

m_lengths = sapply(unique(data$V2), function(x)(mean(data$V1[data$V2 == x])))
names(m_lengths)=unique(data$V2)
mml=min(m_lengths)
nt=(which((m_lengths - mml) > 0.25 *mml))
median=median(m_lengths[names(m_lengths) %notin% names(nt)])
data$V4 = rep(0)
data$V5 = rep(NA)

for (s in 1:nrow(data)){
  if (data[s,2] == names(nt)){
    data[s,4] = 1
    data[s,5] = median/data[s,1]
  }
}

write.table(data, "transform_data.csv", quote = F, row.names = F)


