#' NNLS for Tx Abundance Calculation
#'
#' This function applies NNLS to calculate the transcript abundance
#' from coverage of exonic and junction coverage statistics for samples
#' annotated in a phenotype matrix created using processPheno().
#' @param pheno The table of phenotype information from processPheno().
#' @param cores The number of processing cores to use.
#' @return Returns an rse object of the estimated number of reads and the associated standard errors.
#' Each row represents a protein-coding gene, and each column represents a sample in the phenotype matrix.
#' @keywords recountNNLS
#' @export
recountNNLS = function(pheno, jx_file=NULL, local=F, counts_ex=NULL, counts_jx=NULL, cores=1){
      rls = unique(pheno$rls_group)
      message(Sys.time(), " ##### There are ", length(rls), " read length groups and ", dim(pheno)[1], " samples")

      if(is.null(counts_jx)){
            if(is.null(jx_file)){
                  counts_jx = getJxCounts(unique(pheno$project), local=local)
            }else{
                  counts_jx = getJxCounts(jx_file)
            }
      }
      rse_list = lapply(rls, processReadLength, pheno, counts_ex, counts_jx, cores)

      message(Sys.time(), " ## Processing all RSEs")
      output = do.call(SummarizedExperiment::cbind, rse_list)

      if(local==T){
            pheno = processPheno(unique(pheno$project))
            colData(output)$bigwig_path = pheno$bigwig_path
      }
      return(output)
}

