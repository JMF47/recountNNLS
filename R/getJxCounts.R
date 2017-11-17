getJxCounts = function(pheno){
      project = unique(pheno$project)
      rse_jx_file = recount::download_study(project, type = 'rse-jx', download = FALSE)
      load(url(rse_jx_file))
      data(gff_jx, package="recountNNLSdata")
      ol = findOverlaps(gff_jx, rowRanges(rse_jx), type="equal")
      rse_jx = rse_jx[subjectHits(ol),]
      counts_jx = assays(rse_jx[, match(pheno$run, colnames(rse_jx))])$counts
      if(dim(counts_jx)[1]>0){
            rownames(counts_jx) = paste0("i", queryHits(ol))
      }
      return(counts_jx)
}
