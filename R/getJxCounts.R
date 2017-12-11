#' Obtain Exonic Feature Counts
#'
#' This function returns the coverage of the junction portion of the sufficient
#' features identified based on sequencing read length.
#' @param input Either a SRA project name (part of recount2) OR
#' the path to the rail output file containing junction information of the form 'cross_sample_results/junctions.tsv.gz'.
#' @keywords getJxCounts
#' @export
getJxCounts = function(input, local=F){
      message(Sys.time(), " # Getting junction counts")
      data(gff_jx, package="recountNNLSdata")
      url_table <- recount::recount_url

      ## If input is a SRA project name compiled in recount2
      if(input %in% unique(url_table$project)){
            ## Download the junction file if it exists
            if(sum(url_table$project == input & url_table$file_name=='rse_jx.Rdata')>0){
                  if(local==F){
                        rse_jx_file = recount::download_study(input, type = 'rse-jx', download = FALSE)
                        load(url(rse_jx_file))
                  }else{
                        load(url_table$path[url_table$project==input & url_table$file_name=='rse_jx.Rdata'])
                  }

                  ## Look for junctions matching reference junctions
                  ol = findOverlaps(gff_jx, rowRanges(rse_jx), type="equal")
                  rse_jx = rse_jx[subjectHits(ol),]

                  ## Extract counts of reference junctions and annotate with feature name
                  counts_jx = assays(rse_jx[, match(pheno$run, colnames(rse_jx))])$counts
                  if(dim(counts_jx)[1]>0)
                        rownames(counts_jx) = paste0("i", queryHits(ol))
                  rm(rse_jx)
            ## Input is the path to the junction coverage file from rail
            }else{
                  ## Read junction table and parse information
                  tab = read.table(input, header=T, sep="\t")
                  info = stringr::str_split_fixed(tab[,1], ";", n=4)
                  gr = GRanges(info[,1], IRanges(as.numeric(info[,3]), as.numeric(info[,4])))

                  ## Look for junctions matching reference junctions
                  ol = findOverlaps(gff_jx, gr, type="equal")

                  ## Extract counts of reference junctions and annotate with feature name
                  counts_jx = tab[subjectHits(ol),-1]
                  if(dim(counts_jx)[1]>0)
                        rownames(counts_jx) = paste0("i", queryHits(ol))
            }
            return(counts_jx)
      }
      return(NULL)
}
