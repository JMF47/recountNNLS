#' @export
getCounts <- function(pheno, jx_file, cores=1){
      ## Create the appropriate count matrix
      counts_ex = getExCounts(pheno, cores = cores)
      counts_jx = getJxCounts(jx_file, pheno)
      if(length(counts_jx)>0){
            ## Intersect prioritizes the order of the first argument
            ## so that order of the columns matches the order of pheno
            samps = intersect(colnames(counts_ex), colnames(counts_jx))
            counts = rbind(counts_ex[,samps,drop=FALSE],
                     counts_jx[,samps,drop=FALSE])
            return(counts)
      }else{
            return(counts_ex)
      }
}
