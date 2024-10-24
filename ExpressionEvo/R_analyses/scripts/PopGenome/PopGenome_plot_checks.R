library(tidyverse)
library(gridExtra)
library(ggtext)
library(cowplot)
library(patchwork)
library(readr)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(PopGenome)
library(openxlsx)
rm(list = ls())

## READ IN DATA FUNCTION
metadata <- read.csv("indata/PopGenome/metadata_full.csv", header = F) %>% 
  unique()
colnames(metadata) <- c("VCF", "CONTIG", "LENGTH")

metadata$Chr <- metadata$VCF %>% 
  gsub("Singhal_bTG1_4_", "", .) %>% 
  gsub("_plusLTF_filtered_snps.vcf.gz", "", .)

for (vcf_type in c("correctedZ_PH")){#, "8As_8Bs_PH")){
  
  files = list.files(paste0("indata/PopGenome/", vcf_type, "/"), full.names = T, pattern = ".csv")
  file = files[1]
  
  read_popgenome_data <- function(file){
    data <- read.csv(file, header = T)
    fn = strsplit(file, "/")[[1]][5]
    parts <- strsplit(fn, "_")[[1]]
    if (parts[1] == "all"){
      parts[1] <- paste0(parts[1], parts[2])
      parts <- parts[-2]
    } else if (parts[1] %in% c("CorZPH", "8As8BsPH")) {
      if (parts[2] == "all"){
        parts[2] <- paste0(parts[2], parts[3])
        parts <- parts[-3]
      }
      parts[1] <- paste0(parts[1], parts[2])
      parts <- parts[-2]
    }
    SNP_type <- parts[1]
    contig <- paste(parts[2], parts[3], sep = "_")
    chr <- metadata$Chr[metadata$CONTIG == contig]
    data$contig <- contig
    data$chr <- chr
    data$SNP_type <- SNP_type
    return(data)
  }
  
  
  chr_order <- c("Chr1", "Chr1a", "Chr2", "Chr3", 
                 "Chr4", "Chr4a", paste0("Chr", 5:30), 
                 "Z")
  
  PG_data <- lapply(files, read_popgenome_data) %>% bind_rows()
  PG_data$chr <- factor(PG_data$chr, levels = c(intersect(chr_order, unique(PG_data$chr))))

  min_seg_sites <- c(0,4,10,20,50,100,200) 
  system(paste0("mkdir plots/", vcf_type))
  
  cols = colnames(PG_data)[c(1,5,14,15,20,21,23,24,29,30,34,35,36)]
  
  oldnames <- cols[c(3,5,6,7,9,10)]
  newnames <- c("Tajima's D (A)", "Fay and Wu's H (A)", "Zeng's E (A)", "Tajima's D (B)", "Fay and Wu's H (B)", "Zeng's E (B)")
  ZvsA_plot <- PG_data %>% 
    filter(SNP_type == 'CorZPHNCXG') %>% 
    dplyr::select(cols) %>% 
    filter(chr %in% c("Z", "Chr1", "Chr1a", "Chr2", "Chr4", "Chr5", "Chr6")) %>% 
    mutate(chr = as.character(chr)) %>% 
    mutate(chr = ifelse(chr == "Z", "Z", "Autosomes")) %>%
    mutate(chr = ifelse(chr == "Z" & Pos > 6500000 & Pos < 70100000, "Z Inversion", chr)) %>%
    pivot_longer(cols = cols[-c(1,2,4,8,11,12,13)], names_to = "stat", values_to = "value") %>%
    mutate(chr = ifelse(chr == "Z", "Outside Z Inversion", chr)) %>%
    mutate(chr = factor(chr, levels = c("Autosomes", "Outside Z Inversion", "Z Inversion"))) %>% 
    mutate(stat = factor(stat, levels = oldnames, labels = newnames)) %>%
    group_by(chr, stat) %>% 
    #summarise(mean = mean(value, na.rm = T), var = var(value, na.rm = T), sd = sd(value, na.rm = T), n = n(), se = sd/sqrt(n)) %>% 
    #ggplot(aes(x = chr, y = mean, fill = chr)) + geom_bar(stat = 'identity') + 
    #geom_errorbar(aes(ymin = mean - var, ymax = mean + var), width = 0.2) +
    ggplot(aes(x = chr, y = value, fill = chr)) +
    geom_boxplot(outlier.alpha = 0.1) + 
    facet_wrap(~stat) + scale_fill_grey(start = 0.4) + geom_hline(yintercept = 0, linetype = 2) + labs(y = "", x = "") + 
    theme_bw()+
    theme(legend.position ="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
          axis.ticks = element_line(color = "black"))
  
  
  ggsave(filename = paste0("plots/", vcf_type, "_ZvsA.png"),
         ZvsA_plot,
         width = 5, height = 5, units = "in", device='png')
  system("open plots/correctedZ_PH_ZvsA.png")
  
  Summary_tables <- list()
  make_plots = T
  for (mss in min_seg_sites[1]){
    for (snptype in unique(PG_data$SNP_type)){
      
      All_Chroms <- PG_data %>% 
        filter(SNP_type == snptype) %>% 
        filter(!is.na(A_n.segregating.sites) & !is.na(B_n.segregating.sites))
    
      low_seg_sites <- which(All_Chroms$A_n.segregating.sites < mss & All_Chroms$B_n.segregating.sites < mss)

      All_Chroms <- All_Chroms %>% 
        mutate(Pi_All_Adj = Pi_alls, 
               Position = Pos, 
               Chr = gsub("Chr", "", chr), 
               Pi_A = Ahaps_pi, 
               Pi_B = Bhaps_pi, 
               TajimasD_A = A_Tajima.D, 
               TajimasD_B = B_Tajima.D,
               FayWuH_A = A_Fay.Wu.H,
               FayWuH_B = B_Fay.Wu.H, 
               ZengE_A = A_Zeng.E, 
               ZengE_B = B_Zeng.E, 
               TajimasD_All = Tajima.D, 
               FayWuH_All = Fay.Wu.H,
               ZengE_All = Zeng.E)
      
      if (make_plots == T){
        
    
       All_Chroms[low_seg_sites,c("FST", "dxy")] <- NA
        
        All_Chroms$Chr <- factor(All_Chroms$Chr, levels = c("1","1a", "2":"4", "4a", "5":"30", "Z"))
        All_Chroms <- All_Chroms %>%
          arrange(Chr, Position) %>%
          mutate(PlotOrder = 1:n())
        
        All_Chroms$Dxy_Adj <- All_Chroms$dxy / 100000
        
        
        
        ## GETTING PLOT DATA ##
        
        MeanPi_Auto <- All_Chroms %>%
          filter(Chr != "Z") %>%
          summarise(mean=mean(Pi_All_Adj, na.rm = T), sd=sd(Pi_All_Adj, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Autosomes", 
                 stat = "Pi")
        
        MeanPi_Notinversion <- All_Chroms %>%
          filter(Chr == "Z", Position < 6500000 | Position > 70100000) %>%
          summarise(mean=mean(Pi_All_Adj, na.rm = T), sd=sd(Pi_All_Adj, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_outside", 
                 stat = "Pi")
        
        MeanPi_Inversion <- All_Chroms %>%
          filter(Chr == "Z", Position > 6500000 & Position < 70100000) %>%
          summarise(mean=mean(Pi_All_Adj, na.rm = T), sd=sd(Pi_All_Adj, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_inside", 
                 stat = "Pi")
        
        
        MeanPi_Inversion_A <- All_Chroms %>%
          filter(Chr == "Z", Position > 6500000 & Position < 70100000) %>%
          summarise(mean=mean(Pi_A/100000, na.rm = T), sd=sd(Pi_A/100000, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_insideA", 
                 stat = "Pi")
        
        MeanPi_Inversion_B <- All_Chroms %>%
          filter(Chr == "Z", Position > 6500000 & Position < 70100000) %>%
          summarise(mean=mean(Pi_B/100000, na.rm = T), sd=sd(Pi_B/100000, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_insideB", 
                 stat = "Pi")
        
        MeanDxy_Auto <- All_Chroms %>%
          filter(Chr != "Z") %>%
          summarise(mean=mean(Dxy_Adj, na.rm = T), sd=sd(Dxy_Adj, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Autosomes", 
                 stat = 'dxy')
        
        
        MeanDxy_inversion <- All_Chroms %>%
          filter(Chr == "Z", Position > 6500000 & Position < 70100000) %>%
          summarise(mean=mean(Dxy_Adj, na.rm = T), sd=sd(Dxy_Adj, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_inside", 
                 stat = 'dxy')
        
        MeanDxy_Notinversion <- All_Chroms %>%
          filter(Chr == "Z", Position < 6500000 | Position > 70100000) %>%
          summarise(mean=mean(Dxy_Adj, na.rm = T), sd=sd(Dxy_Adj, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_outside", 
                 stat = 'dxy')
        
        
        MeanFst_Notinversion <- All_Chroms %>%
          filter(Chr == "Z", Position < 6500000 | Position > 70100000) %>%
          summarise(mean=mean(FST, na.rm = T), sd=sd(FST, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_outside", 
                 stat = 'FST')
        
        MeanFst_inversion <- All_Chroms %>%
          filter(Chr == "Z", Position > 6500000 & Position < 70100000) %>%
          summarise(mean=mean(FST, na.rm = T), sd=sd(FST, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Z_inside", 
                 stat = 'FST')
        
        MeanFst_Auto <- All_Chroms %>%
          filter(Chr != "Z") %>%
          #drop_na(FST) %>%
          summarise(mean=mean(FST, na.rm = T), sd=sd(FST, na.rm = T), n= n()) %>%
          mutate(se = sd/sqrt(n)) %>% 
          mutate(Chr = "Autosomes", 
                 stat = 'FST')
        
        Summary_table <- rbind(MeanPi_Auto, MeanPi_Notinversion, MeanPi_Inversion, 
                               MeanPi_Inversion_A, MeanPi_Inversion_B, 
                               MeanDxy_Auto, MeanDxy_Notinversion, MeanDxy_inversion, 
                               MeanFst_Auto, MeanFst_Notinversion, MeanFst_inversion)
        Summary_table$snptype = snptype
        Summary_tables[[snptype]] <- Summary_table
        MedianChrom <- All_Chroms %>%
          group_by(Chr) %>%
          summarise(Middle = median(PlotOrder, na.rm=TRUE))
        
        All_Chroms_Long <- All_Chroms %>%
          select(Position, FST, Dxy_Adj, Pi_All_Adj, Chr, PlotOrder) %>%
          pivot_longer(names_to = "Test", values_to="TestStat", cols=c(FST, Dxy_Adj, Pi_All_Adj)) %>%
          mutate(Test = str_replace_all(Test,c("FST"="Fst: A vs B", "Dxy_Adj" = "Dxy: A vs B", "Pi_All_Adj" = "Pi (Nei's)")))
        
      
        Z_Chroms_Long <- All_Chroms %>%
          filter(Chr =="Z") %>%
          select(Position, FST, Dxy_Adj, Pi_All_Adj, Chr, PlotOrder) %>%
          pivot_longer(names_to = "Test", values_to="TestStat", cols=c(FST, Dxy_Adj, Pi_All_Adj)) %>%
          mutate(Test = str_replace_all(Test,c("FST"="Fst: A vs B", "Dxy_Adj" = "Dxy: A vs B",  "Pi_All_Adj" = "Pi (Nei's)")))
        
        ###########
        
        
        ## FIGURE 1 PLOTING ##
        
        mypalette <- rep(c("grey80", "darkorange1", "black", "firebrick1", "grey60", "tan4"),length.out=40) # chr color palette
        
        AllStats <- All_Chroms_Long %>% 
          filter(Test != "Dxy: A vs B") %>% 
          ggplot(.,aes(x=PlotOrder,y = TestStat, color=Chr)) +
          geom_point(size=0.1)+
          theme_bw()+
          geom_hline(yintercept=0,linetype='solid',colour="grey50")+
          scale_colour_manual(values=mypalette) +
          scale_x_continuous(label = MedianChrom$Chr[c(1:12,33)], breaks= MedianChrom$Middle[c(1:12,33)]) +
          theme(legend.position ="none") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          xlab("") +
          ylab("") +
          facet_wrap(~Test, ncol=1, scales="free_y")
        
      # Allstats_smooth_line <-  
      #   All_Chroms_Long %>% filter(Test != "Dxy: A vs B") %>% 
      #     ggplot(.,aes(x=PlotOrder,y = TestStat, color=Chr)) +
      #     #geom_line(linewidth=0.3)+
      #     geom_point(size=0.3, alpha = 0.1)+
      #     theme_bw()+
      #     geom_hline(yintercept=0,linetype='solid',colour="grey50")+
      #     scale_colour_manual(values=mypalette) +
      #     scale_x_continuous(label = MedianChrom$Chr[c(1:12,33)], breaks= MedianChrom$Middle[c(1:12,33)]) +
      #     theme(legend.position ="none") +
      #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      #     xlab("") +
      #     ylab("") +
      #     facet_wrap(~Test, ncol=1, scales="free_y") + 
      #     geom_smooth(data=subset(All_Chroms_Long, Chr %in% c('1',"1a",'2','3','4','4a','5','6','7','8','9','10') & 
      #                               Test != "Dxy: A vs B"), 
      #                 linewidth=0.4, aes(group = Chr), colour = 'blue', se = F, 
      #                method = 'loess', span = 0.3)
      # 
      #   ggsave(filename = paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_Pi_Fst_Dxy_All_smoothed_line.png"),
      #          Allstats_smooth_line,
      #          width = 8, height = 6, units = "in", device='png')
        
        Z_labeller <- as_labeller(c(`Dxy: A vs B` = expression("Dxy~*~10^-3:A~vs~B"),
                                    `Fst: A vs B` = "`Fst: A vs B`", 
                                    `Pi (Nei's)` = "`Pi (Nei's)`"), label_parsed)
        
        AllStatsZ <- ggplot(Z_Chroms_Long,aes(x=Position/1000000,y = TestStat, color=Chr)) +
          geom_point(size=0.3)+
          theme_bw()+
          geom_hline(yintercept=0,linetype='solid',colour="grey50")+
          scale_colour_manual(values="black") +
          #scale_x_continuous(label = MedianChrom$Chr[c(33)], breaks= MedianChrom$Middle[c(33)]) +
          theme(legend.position ="none") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          geom_vline(xintercept = 6.5, colour="firebrick1", linetype = "longdash") +
          geom_vline(xintercept = 70.1, colour="firebrick1", linetype = "longdash") +
          xlab("Position (Mbp)") +
          ylab("") +
          facet_wrap(~Test, ncol=3, scales="free_y")
        
        
        Pi_Fst_Dxy_All_Z <- AllStats + AllStatsZ +
          plot_layout(ncol=1, nrow=2, widths=c(2.5), heights = c(3, 1)) +
          plot_annotation(tag_levels = 'a') 
        
        
        ggsave(filename = paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_Pi_Fst_Dxy_All_Z.png"),
               Pi_Fst_Dxy_All_Z,
               width = 8, height = 6, units = "in", device='png')
        
        ggsave(filename = paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_Pi_Fst_Dxy_All_Z.pdf"),
                Pi_Fst_Dxy_All_Z,
                width = 8, height = 6, units = "in", device='pdf')
         
        ggsave(filename = paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_Pi_Fst_Dxy_All_Z.svg"),
                Pi_Fst_Dxy_All_Z,
                width = 8, height = 6, units = "in", device='svg')
      
        
        ###########################################
        
        ## FIGURE 2 Plotting ##
        
        ChrZ_Long <- All_Chroms %>%
          filter(Chr == "Z") %>% 
          dplyr::select(Position, TajimasD_A, TajimasD_B, FayWuH_A, FayWuH_B,
                 ZengE_A, ZengE_B, Pi_A, Pi_B) %>%
          mutate(Pi_A = Pi_A / 10000) %>%
          mutate(Pi_B = Pi_B / 10000) %>%
          mutate(Position = Position / 1000000) %>%
          pivot_longer(names_to = "Test", values_to="TestStat", cols = TajimasD_A:Pi_B) %>%
          mutate(Haplotype = rep(c("Haplotype A", "Haplotype B"),length.out=nrow(.))) %>%
          mutate(Test2 = str_replace_all(Test,c("TajimasD_A" = "Tajima's D",
                                                "TajimasD_B" = "Tajima's D", 
                                                "FayWuH_A" = "Fay and Wu's H",
                                                "FayWuH_B" = "Fay and Wu's H",
                                                "ZengE_A" = "Zeng's E",
                                                "ZengE_B" = "Zeng's E",
                                                "Pi_A" = "Pi (Nei's)",
                                                "Pi_B" = "Pi (Nei's)")))
        
        #plot add ons
        pas <- theme_bw()+
            theme(legend.position ="none") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(strip.background=element_rect(fill="darkorange2"))
        
        ChrZ_Neutrality <- ###Edit to ChrZ_Pi for diversity, Chr_Neutrality for others
          ggplot(subset(ChrZ_Long,Test2 != "Pi (Nei's)"),aes(x=Position,y = TestStat)) +
          geom_point(size=0.3)+
          geom_hline(yintercept=0,linetype='solid',colour="grey50")+
          geom_vline(xintercept = 6.5, colour="firebrick1", linetype = "longdash") +
          geom_vline(xintercept = 70.1, colour="firebrick1", linetype = "longdash") +
          scale_colour_manual(values="black") +
          facet_grid(cols=vars(Haplotype), rows = vars(Test2)) + 
          ylab("") + xlab("") + pas
        
        ChrZ_Pi <- 
          ggplot(subset(ChrZ_Long,Test2 == "Pi (Nei's)"),aes(x=Position,y = TestStat)) +
          geom_point(size=0.3)+
          ylim(0,0.1) + 
          geom_hline(yintercept=0,linetype='solid',colour="grey50")+
          geom_vline(xintercept = 6.5, colour="firebrick1", linetype = "longdash") +
          geom_vline(xintercept = 70.1, colour="firebrick1", linetype = "longdash") +
          scale_colour_manual(values="black") +
          facet_grid(cols=vars(Haplotype), rows = vars(Test2)) + ylab("") +
          xlab("Position (Mbp)") +
          pas
          
        
        
        ChromZPlots <- ChrZ_Neutrality + ChrZ_Pi +
          plot_layout(ncol=1, nrow=2, axis_titles='collect_x', heights=c(2.5,1)) +
          plot_annotation(tag_levels = 'a') 
        
        ChromZPlots
        
        
        ggsave(plot=ChromZPlots, filename=paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_ChromZ_Neutrality_Diversity_Longitudinal.png"), width=5,height=8, device = 'png')
        ggsave(plot=ChromZPlots, filename=paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_ChromZ_Neutrality_Diversity_Longitudinal.pdf"), width=5,height=8, device = 'pdf')
        ggsave(plot=ChromZPlots, filename=paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_ChromZ_Neutrality_Diversity_Longitudinal.svg"), width=5,height=8, device = 'svg')
        
        
        ###########################################
        
        
        
        ##############################################
        ##### Supplementary SFS plots
      
        All_Chroms_SFS <- All_Chroms %>%
          select(Position, TajimasD_All, FayWuH_All, ZengE_All, Chr, PlotOrder) %>%
          pivot_longer(names_to = "Test", values_to="TestStat", cols=c(TajimasD_All,FayWuH_All,ZengE_All)) %>%
          mutate(Test = str_replace_all(Test,c("TajimasD_All"="Tajima's D", "FayWuH_All" = "Fay & Wu's H", "ZengE_All" = "Zeng's E")))
      
        Z_Chroms_SFS <- All_Chroms %>%
          filter(Chr =="Z") %>%
          select (Position, TajimasD_A, TajimasD_B, FayWuH_A, FayWuH_B, ZengE_A, ZengE_B, Chr, PlotOrder) %>%
          pivot_longer(names_to = "Test", values_to="TestStat", cols=c(TajimasD_A, TajimasD_B, FayWuH_A, FayWuH_B, ZengE_A, ZengE_B)) %>%
          mutate(Test = str_replace_all(Test,c("TajimasD_A" = "Tajima's D: A",
                                               "TajimasD_B" = "Tajima's D: B",
                                               "FayWuH_A" = "Fay & Wu's H: A",
                                               "FayWuH_B" = "Fay & Wu's H: B",
                                               "ZengE_A" = "Zeng's E: A",
                                               "ZengE_B" = "Zeng's E: B")))
      
        All_Chroms_SFS$Test <- factor(All_Chroms_SFS$Test, levels = c("Tajima's D", "Fay & Wu's H", "Zeng's E"))
        Z_Chroms_SFS$Test <- factor(Z_Chroms_SFS$Test, levels = c("Tajima's D: A",
                                                                  "Tajima's D: B",
                                                                  "Fay & Wu's H: A",
                                                                  "Fay & Wu's H: B",
                                                                  "Zeng's E: A",
                                                                  "Zeng's E: B"))
      
      
      
      
        mypalette <- rep(c("grey80", "darkorange1", "black", "firebrick1", "grey60", "tan4"),length.out=40) # chr color palette
      
      
      
        SFS_Stats <- All_Chroms_SFS %>% 
          filter(Chr != "Z") %>% 
          ggplot(.,aes(x=PlotOrder,y = TestStat, color=Chr)) +
          #geom_line(linewidth=0.3)+
          geom_point(size=0.3)+
          theme_bw()+
          geom_hline(yintercept=0,linetype='solid',colour="grey50")+
          scale_colour_manual(values=mypalette) +
          scale_x_continuous(label = MedianChrom$Chr[c(1:12)], breaks= MedianChrom$Middle[c(1:12)]) +
          theme(legend.position ="none") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          xlab("") +
          ylab("") +
          facet_wrap(~Test, ncol=1, scales="free_y")
      
      
        SFS_StatsZ <- ggplot(Z_Chroms_SFS,aes(x=Position/1000000,y = TestStat, color=Chr)) +
          geom_point(size=0.3)+
          theme_bw()+
          geom_hline(yintercept=0,linetype='solid',colour="grey50")+
          scale_colour_manual(values="black") +
          #scale_x_continuous(label = MedianChrom$Chr[c(33)], breaks= MedianChrom$Middle[c(33)]) +
          theme(legend.position ="none") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          geom_vline(xintercept = 6.5, colour="firebrick1", linetype = "longdash") +
          geom_vline(xintercept = 70.1, colour="firebrick1", linetype = "longdash") +
          xlab("Position (Mbp)") +
          ylab("") +
          facet_wrap(~Test, ncol=2, dir="h")
      
      
        SFS_All_Z <- SFS_Stats + SFS_StatsZ +
          plot_layout(ncol=1, nrow=2, widths=c(2.5), heights = c(2, 2)) +
          plot_annotation(tag_levels = 'a')
      
      
        ggsave(filename = paste0("plots/min_seg_sites_", mss, "_",  vcf_type, "_",snptype, "_SFS_All_Z.png"),
               SFS_All_Z,
               width = 8, height = 8, units = "in", device='png')
    
      #####################################
      }
    }
  }
  system(paste0("mv plots/*pdf ", "plots/", vcf_type, "/"))
  system(paste0("mv plots/*png ", "plots/", vcf_type, "/"))
  system(paste0("mv plots/*svg ", "plots/", vcf_type, "/"))
  
  Summary_tables %>% bind_rows() %>%
    mutate(mean = ifelse(stat == "Pi", (mean*1000), mean)) %>%
    mutate(mean = ifelse(stat == "dxy", (mean*1000), mean)) %>%
    mutate(values = paste0(round(mean, 4), " (",round(sd,6), ")")) %>%
    select(-c(mean,sd,n,se)) %>%
    spread(Chr, values) %>%
    filter(snptype %in% c("CorZPHZD", "CorZPHFD", "CorZPHNCXG")) %>%
    mutate(stat=ifelse(stat == "Pi", "pi (*10^-3)", stat),
           stat=ifelse(stat == "dxy", "dxy (*10^-3): A vs B", stat),
           stat = ifelse(stat == "FST", "FST:A vs B", stat)) %>%
    rename(`Variable Site` = snptype) %>%
    mutate(
           `Variable Site` = ifelse(`Variable Site` == "CorZPHZD", "0-fold", `Variable Site`),
           `Variable Site` = ifelse(`Variable Site` == "CorZPHFD", "4-fold", `Variable Site`),
           `Variable Site` = ifelse(`Variable Site` == "CorZPHNCXG", "Non-Coding", `Variable Site`)) %>%
    rename(`Z (outside inversion)` = Z_outside,
           `Z (within inversion)` = Z_inside,
           `A birds (within inversion)` = Z_insideA,
           `B birds (within inversion)` = Z_insideB,
           Measure = stat)  %>%
    .[c(9,7,8,3,1,2,6,4,5),c(1:3,7,4,5,6)] %>%
    write.csv(paste0("indata/PopGenome/", vcf_type, "_Summary_tables.csv"), row.names=FALSE)
  
  
}



