processPheno = function(input){
      ## If input is of length 1, will interpret as SRA project name
      if(length(input)==1){
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
            pheno$bigwig_path = url_table$url[match(pheno$bigwig_file, url_table$file_name)]
      }else{
            pheno = input
      }
      rls_avail = c(37, 50, 75, 100, 150)
      pheno$rls_group = sapply(pheno$avg_read_length, function(x) rls_avail[which.min(abs(rls_avail-x))])
      return(pheno)
}

## Helper function for reading the phenotype files
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
