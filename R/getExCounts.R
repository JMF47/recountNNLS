#' Obtain Exonic Feature Counts
#'
#' This function returns the coverage of the exonic portion of sufficient
#' features identified based on sequencing read length. It takes as input
#' the phenotype matrix that is created by processPheno().
#' @param pheno The phenotype matrix created by processPheno()
#' @keywords getExCounts
#' @export
getExCounts = function(pheno){
      ## Evaluate the consistency of read lengths supplied in phenotype
      rl = unique(pheno$rls_group)
      if(length(rl)>1)
            stop("Cannot process two different read group lengths at the same time. Please split pheno")

      ## Load appropriate features based on read count
      data(list=paste0("bins_", rl), package = "recountNNLSdata")

      ## Create a GRangesList by chr
      grl = GRangesList()
      for(chr in unique(seqnames(bins)))
            grl[[chr]] = bins[seqnames(bins)==chr]

      ## Counting the coverage of exonic features by file
      list_totCov= lapply(pheno$bigwig_path, .processSample, grl)
      ## Assembling the count information
      totCov = do.call(cbind, list_totCov)
      ## Giving correct exonic feature annotation
      rownames(totCov) = paste0("e", 1:length(bins))

      return(totCov)
}

.processSample = function(sampleFile, grl){
      message("Processing sample ", sampleFile)

      cov_rle = rtracklayer::import(sampleFile, as = 'RleList')
      ## To ensure consistency when some samples have chrs dropped from no mapped reads
      cov_rle_matched = cov_rle[match(names(grl), names(cov_rle), nomatch=0)]
      grl_keep = grl[match(names(cov_rle_matched), names(grl), nomatch=0)]
      cov_binned = sapply(names(grl_keep), .processChr, cov_rle_matched, grl_keep)
            cov_binned = do.call(c, cov_binned)
      id = match(bins, do.call(c, grl_keep))

      cov_out = rep(0, length(bins))
      cov_out[id[!is.na(id)]] = cov_binned[!is.na(id)]
      return(totCov)
}
.processChr = function(chr, rle_cov, grl_bins){
      sum(Views(rle_cov[[chr]], ranges(grl_bins[[chr]])))
}
