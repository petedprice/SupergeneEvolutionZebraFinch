library(GenomicRanges)
library(tidyverse)
library(data.table)
rm(list = ls())

gff <-ape::read.gff("indata/genome_files/intron_GCF_003957565.2_bTaeGut1.4.pri_genomic.gff") %>% 
  filter(type %in% c("gene", "intron")) %>% 
  dplyr::select(seqid, type, start, end)
metadata <- read.csv("indata/genome_files/metadata_full.csv", header = F)

out <- data.frame(contig = NA, start = NA, end = NA, length = NA, Pos = NA, window = NA)


stp <- 10000
nd=0
for (contig in unique(metadata$V2)){
  print(contig)
  ct_gff <- filter(gff, seqid == contig)
  st = 2
  Pos <- 0
  i=0
  nd = 0
  ct_length <- metadata$V3[metadata$V2 == contig]
  while (nd < ct_length){
    window <- 100000
    Pos <- Pos + stp
    i = i + 1
    if(i %% 100 == 0) print(i)
    nd = st + window -1
    if (nd > ct_length) window = ct_length - st
    window_gff <- ct_gff %>% 
      filter((start >= st & start <= nd) | (end >= st & end <= nd) | 
               start <= st & end >=nd) %>% 
      unique()
    window_gff$start[window_gff$start < st] <- st
    window_gff$end[window_gff$end > nd] <- nd
    
    introns <- filter(window_gff, type == 'intron') %>% 
      GRanges(.) %>% GenomicRanges::reduce() %>% GRanges(.) %>% as.data.table()
    gene <- filter(window_gff, type %in% c('gene')) %>% 
      GRanges(.) %>% GenomicRanges::reduce() %>% GRanges(.) %>% as.data.table()
      
    
    #length <- window - (sum(gene$width) - sum(introns$width))
    
    sd <- GenomicRanges::intersect(GRanges(introns), GRanges(gene)) %>% as.data.table()
    length = window - sum(as.data.table(gene)$width) + sum(sd$width)
    
    out_tmp <- data.frame(contig = contig, start = st, end = nd, length = length, Pos = Pos, window = window)
    out <- rbind(out, out_tmp)
    if (length > window) stop("length > window")
    st = st + stp
  }
}
out <- out[-1,]
#write.csv(out, "indata/genome_files/NCXG_window_sizes.csv", row.names = F, quote = F)


rm(list = ls())
metadata <- read.csv("indata/genome_files/metadata_full.csv", header = F)


stp <- 10000
nd=0
for (contig in unique(metadata$V2)){
  print(contig)
  ct_FZD <- read.table(paste0("indata/genome_files/", contig, "_ZFD.bed"), sep = '\t') %>% 
    dplyr::select(V1,V2,V3,V5) %>% 
    unique()
  st = 2
  Pos <- 0
  i=0
  nd = 0
  ct_length <- metadata$V3[metadata$V2 == contig]
  ZFD_out <- data.frame(contig = NA, start = NA, end = NA, 
                        FD_length = NA, ZD_length = NA, 
                        Pos = NA, window = NA)
  while (nd < ct_length){
    window <- 100000
    Pos <- Pos + stp
    i = i + 1
    if(i %% 100 == 0) print(st)
    nd = st + window -1
    if (nd > ct_length) window = ct_length - st
    
    window_FZD <- ct_FZD %>% 
      filter(V3 >= st & V3 <= nd) %>% 
      unique()
    
    FD <- length(unique(filter(window_FZD, V5 == 4)$V2))
    ZD <- length(unique(filter(window_FZD, V5 == 0)$V2))
    
    ZFD_out_tmp <- data.frame(contig = contig, start = st, end = nd, 
                              FD_length = FD, ZD_length = ZD, Pos = Pos, 
                              window = window)
    ZFD_out <- rbind(ZFD_out, ZFD_out_tmp)
    if (FD > window) stop("FD > window")
    if (ZD > window) stop("ZD > window")
    st = st + stp  
  }
  ZFD_out <- ZFD_out[-1,]
  write.table(ZFD_out, paste0("indata/genome_files/", contig, "_ZFD_window_sizes.csv"), 
              row.names = F, quote = F, sep = '\t')
}
