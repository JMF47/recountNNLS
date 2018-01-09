#' Download and/or access quantified SRA projects
#'
#' This function downloads and/or loads into R a
#' \code{rangedSummarizedExperiment} object for transcript level abundance
#' produced by the latest version of recountNNLS.
#' @param project A SRA project name for a project already curated in recount2.
#' @param download_path Where the rse should be downloaded to if not NULL.
#' Defaults to NULL, i.e. no downloading.
#' @return Returns the path to downloaded object.
#' @examples
#' project = 'DRP000366'
#' getRseTx(project)
#' @keywords recountNNLS
#' @export
getRseTx = function(project, download_path=NULL){
      target = paste0('http://duffel.rail.bio/recount/', project, '/rse_tx.RData')
      if(is.null(download_path)){
            error("Please specify a download path.")
      }else{
            download.file(url=paste0('http://duffel.rail.bio/recount/', project, '/rse_tx.RData'), destfile=download_path)
      }
      return(download_path)
}
