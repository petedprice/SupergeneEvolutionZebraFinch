#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
data=read.csv(args[1], sep = " ", header = TRUE)
#data=read.csv("transform_data.csv", sep = " ", header=F)
sample = args[2]
#species="duck"

cat(noquote(data[,4][data[,3] == sample]))
