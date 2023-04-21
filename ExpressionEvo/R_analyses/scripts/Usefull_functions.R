vcf_filter <- function(vcf, samples, DP = 1, HETT = 0.9, GQ = 2){
  vcf_out <- vcf
  for (g in 1:nrow(vcf)){ 
    if (str_split(vcf[g,samples], ":", simplify = T)[,1] == "./."){
      vcf_out[g,samples] <- "FAIL"
    } else {
      AD <- which(str_split(vcf[g,]$FORMAT, ":", simplify = T) == "AD")
      if (length(AD) != 0){
        score <- str_split(vcf[g,samples], ":", simplify = T)[AD]
        ADscore <- as.numeric(str_split(score, ",", simplify = T))[1:2]
        if (max(ADscore)/sum(ADscore) > HETT | sum(ADscore) == 0){
          het = "hom"
        } else {
          het = "0/1"
        }
      }
      
      filt <- which(str_split(vcf$FORMAT[g], ":", simplify = T) == "DP")
      score <- as.numeric(str_split(vcf[g,samples], ":", simplify = T)[filt])
      if (score >= DP | is.na(score) == TRUE){
        DP_out <- TRUE
      } else { 
        DP_out <- FALSE
      }
      
      filt <- which(str_split(vcf$FORMAT[g], ":", simplify = T) == "GQ")
      score <- as.numeric(str_split(vcf[g,samples], ":", simplify = T)[filt])
      if (score >= GQ | is.na(score) == TRUE){
        GQ_out <- TRUE
      } else { 
        GQ_out <- FALSE
      }
      
      if (GQ_out == FALSE | DP_out == FALSE){
        vcf_out[g,samples] <- "FAIL"
      } else if (het == "0/1" ){
        temp <- as.vector(str_split(vcf_out[g,samples], ":", simplify = T))
        temp[1] <- het
        temp2 <- paste(temp, collapse = ":")
      }
    }   
    
  }
  return(vcf_out)
}


half_mf <- function(txi, info, half_or_keep = "half", g_or_gs = "G"){
  c=0
  for (t in g_or_gs){
    c=c+1
    ss_info <- info[which(info$tissue == t),]
    if (length(which(info$tissue == t)) == 0) {
      next
    }
    inds <- info$sample[info$tissue == t]
    y <- DGEList(counts = txi$counts[,inds], genes = txi$length[,inds], group = ss_info$Species)
    y <- calcNormFactors(y)
    
    if (half_or_keep == "half"){
      Males = info$sample[info$tissue == t & info$sex == "M"]
      Females = info$sample[info$tissue == t & info$sex == "F"]
      rpkm <- rpkm(y, gene.length = y$genes)
      M_keep = rowSums(rpkm[,Males] >2) >= ceiling(length(Males)/2)
      F_keep = rowSums(rpkm[,Females] >2) >= ceiling(length(Females)/2)
      MF_keep = names(which(M_keep == TRUE | F_keep == TRUE))
      
    } else {
      design <- model.matrix(~sex, data = ss_info)
      k = filterByExpr(y, design = design)
      y <- y[k,,k = FALSE]
      MF_keep = rownames(y$counts)
    }
    
    if (c == 1){
      keep <- MF_keep
    } else {
      keep <- unique(keep, MF_keep)
    }
    
  }
  return((keep))
}


rank_distance <- function(x){
  n <- length(x)
  scores <-c()
  for (i in 1:(n-1)){
    a <- x[i]
    b <- abs(x[(i+1):n] - as.numeric(a))
    bsum <- sum(b)
    scores <- c(scores,bsum)
  }
  f_score <- sum(scores)
  return(f_score)
}


IQR_ss <- function(rpkm){
  genes <- rownames(rpkm)
  for (i in 1:ncol(rpkm)){
    x <- rpkm[,i]
    names(x) <- genes
    s <- summary(x)
    f <- s[2]
    t <- s[5]
    ss <- names(x[which(x > f & x < t)])
    if (i ==1){
      sav_genes = ss
    } else {
      sav_genes <- intersect(ss, sav_genes)
    }
  }
  return(rpkm[sav_genes,])
}

median_scaling_factors <- function(y){

  y$samples$norm.factors <- 1
  rpkm <- rpkm(y, gene.length = y$genes)
  IQ_rpkm <- IQR_ss(rpkm)
  IQ_rpkm_rank <- apply(IQ_rpkm, 2, rank)
  rank_scores <- apply(IQ_rpkm_rank[,1:3], 1, rank_distance)
  top_thousand <- names(rank_scores[order(rank_scores)])[1:1000]
  top_thousand <- top_thousand[-which(is.na(top_thousand))]
  tt_rpkm <- IQ_rpkm[top_thousand,]
  
  medians <- apply(tt_rpkm, 2, median)
  centre <- median(medians)
  scaling_factors <- medians/centre
  
  return(scaling_factors)
}

rpkm_hist <- function(y, info){
  rpkm <- rpkm(y, gene.length = y$genes, log = T, prior.count = 1)
  df <- data.frame("rpkm" = c(), "species" = c() )
  for (i in unique(info$species)){
    s = info$sample[info$species == i]
    r = c(as.matrix(rpkm[,s]))
    l <- rep(i, length(r))
    rl <- cbind(r,l)
    df <- rbind(df, rl)
  }
  names(df) <- c("rpkm", "species")
  df$rpkm <- as.numeric(df$rpkm)
  return(df)
}

