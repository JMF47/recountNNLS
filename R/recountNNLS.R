#' NNLS for Tx Abundance Calculation
#'
#' This function applies NNLS to calculate the transcript abundance
#' from coverage of exonic and junction coverage statistics for samples
#' annotated in a phenotype matrix created using processPheno().
#' @param pheno The table of phenotype information from processPheno().
#' @param jx_file The path to the Rail-RNA junction coverage file if processing sample not already in recount2.
#' @param cores The number of processing cores to use.
#' @return Returns an rse object of the estimated number of reads and the associated standard errors.
#' Each row represents a protein-coding gene, and each column represents a sample in the phenotype matrix.
#' @examples
#' project = 'SRP063581'
#' pheno = processPheno(project)
#' @keywords recountNNLS
#' @export
recountNNLS = function(pheno, jx_file=NULL, cores=1){
      rls = unique(pheno$rls_group)
      message(Sys.time(), " ##### There are ", length(rls), " read length groups and ", dim(pheno)[1], " samples")
      if(length(jx_file)==0) jx_file = unique(pheno$project)
      rse_list = lapply(rls, processReadLength, pheno, jx_file, cores)

      message(Sys.time(), " ## Processing all RSEs")
      output = do.call(SummarizedExperiment::cbind, rse_list)

      return(output)
}

