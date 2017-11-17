getExCounts = function(pheno){
      rl = unique(pheno$rls_group)
      if(length(rl)>1)
            stop("Cannot process two different read group lengths at the same time. Please split pheno")
      data(list=paste0("bins_", rl), package = "recountNNLSdata")
      list_totCov= lapply(pheno$bigwig_path, .processSample)
      totCov = do.call(cbind, list_totCov)#/pheno$avg_read_length
      rownames(totCov) = paste0("e", 1:length(bins))
      return(totCov)
}

.processSample = function(sampleFile){
      message("Processing sample ", sampleFile)
      grl = GRangesList()
      for(chr in unique(seqnames(bins))){
            grl[[chr]] = bins[seqnames(bins)==chr]
      }
      regCov <- rtracklayer::import(sampleFile, selection = bins, as = 'RleList')
      regCovKeep = regCov[match(names(grl), names(regCov))]
      list_totCov = sapply(names(grl), .processChr, regCovKeep, grl)
      totCov = do.call(c, list_totCov)
      return(totCov)
}
.processChr = function(chr, rle_cov, grl_bins){
      sum(Views(rle_cov[[chr]], ranges(grl_bins[[chr]])))
}
