#' Phenotype Information Creation
#'
#' This function creates a phenotype/meta-information matrix with the information needed
#' for the rest of the analysis. The supplied input is either an SRA
#' project id (from a project compiled on recount2), or a manifest file
#' that contains information where the Rail-RNA outputs are.
#' @param input The project name from SRA or a manifest file for other projects aligned with Rail-RNA.
#' @return A phenotype matrix containing at least the necessary information to run recountNNLS.
#' @keywords processPheno
#'
#' @export
processPheno = function(input, local=F){
      ## If input is of length 1, will interpret as SRA project name
      if(class(input)=="character"){
            project = input
            url_table <- recount::recount_url
            ## Subset URL Data
            url_table <- url_table[url_table$project == project, ]
            if(nrow(url_table) == 0)
                  stop("Invalid 'project' argument. There's no such 'project' in the recount_url data.frame.")
            sampleFiles <- recount::download_study(project = project, type = 'samples', download = FALSE)
            phenoFile <- recount::download_study(project = project, type = 'phenotype',download = FALSE)
            # Read phenotype and process into different rl groups where applicable
            pheno <- .read_pheno(phenoFile, project)
            if(local==F)
                  pheno$bigwig_path = url_table$url[match(pheno$bigwig_file, url_table$file_name)]
            else
                  pheno$bigwig_path = url_table$path[match(pheno$bigwig_file, url_table$file_name)]
            pheno = pheno[!is.na(pheno$bigwig_path),]
            paired = pheno$paired_end*1+1
            pheno$rls = (pheno$avg_read_length/paired)
      }else{
            pheno = input

            required_columns = c("project", "run", "bigwig_path", "rls", "paired_end")
            has_column = sum(required_columns %in% colnames(pheno))
            if(has_column < length(required_columns))
                  stop(paste0("Check the columns of your input meet requirements for function: ", paste0(required_columns, collapse=", ")))

            if(class(pheno$paired_end)!="logical")
                  stop("paired_end must be T/F")

            if(class(pheno$rls)!="numeric")
                  stop("rls should be numeric")

      }
      rls_avail = c(37, 50, 75, 100, 150)
      pheno$rls_group = sapply(pheno$rls, function(x) rls_avail[which.min(abs(rls_avail-x))])
      return(pheno)
}

## Helper function for reading the phenotype files in recount2
.read_pheno <- function(phenoFile, project) {
      if(project %in% c('SRP012682', 'TCGA')) {
            subsets <- c('SRP012682' = 'gtex', 'TCGA' = 'tcga')
            res <- all_metadata(subsets[project], verbose = FALSE)
      } else {
            info <- readLines(phenoFile)
            res <- read.table(text = info[grepl(paste0('^project|^', project),
                                                info)], header = TRUE, stringsAsFactors = FALSE, sep = '\t',
                              comment.char = '', quote = '')
      }
      return(res)
}
