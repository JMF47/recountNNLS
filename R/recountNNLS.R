#' NNLS for Tx Abundance Calculation
#'
#' This function applies NNLS to calculate the transcript abundance
#' from coverage of exonic and junction coverage statistics for samples
#' annotated in a phenotype matrix created using processPheno().
#' @param pheno The table of phenotype information from processPheno().
#' @param jx_file The path to the Rail-RNA junction coverage file if processing sample not already in recount2.
#' @param counts_ex The matrix of exonic feature counts where each row represents a exonic feature and each column a sample.
#' If set to NULL the function will automatically compile the matrix using the information in pheno. Defaults to NULL
#' @param counts_jx The matrix of junction feature counts where each row represents a junction feature and each column a sample.
#' If set to NULL the function will automatically compile the matrix using the information in pheno or the jx_file parameter. Defaults to NULL
#' @param cores The number of processing cores to use.
#' @return Returns an rse object of the estimated number of reads and the associated standard errors.
#' Each row represents a protein-coding gene, and each column represents a sample in the phenotype matrix.
#' @examples
#' '\dontrun{
#' project = 'SRP063581'
#' pheno = processPheno(project)
#' recountNNLS(pheno)
#' }
#' @keywords recountNNLS
#' @export
recountNNLS = function(pheno, jx_file=NULL, counts_ex=NULL, counts_jx=NULL, cores=1){
      rls = unique(pheno$rls_group)
      message(Sys.time(), " ##### There are ", length(rls), " read length groups and ", dim(pheno)[1], " samples")

      if(is.null(counts_jx)){
            if(is.null(jx_file)){
                  counts_jx = getJxCounts(unique(pheno$project), pheno)
            }else{
                  counts_jx = getJxCounts(jx_file, pheno)
            }
      }
      rse_list = lapply(rls, processReadLength, pheno, counts_ex, counts_jx, cores)

      message(Sys.time(), " ## Processing all RSEs")
      output = do.call(SummarizedExperiment::cbind, rse_list)

      # if(local==T){
      #       pheno = processPheno(unique(pheno$project))
      #       colData(output)$bigwig_path = pheno$bigwig_path
      # }
      return(output)
}

