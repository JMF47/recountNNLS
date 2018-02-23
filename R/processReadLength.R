#' NNLS for Tx Abundance Calculation
#'
#' This function applies runs our NNLS model for a given subset
#' of the pheno matrix that has the same read length specified
#' and reports an rse.
#' @param rl A read length that exists in the pheno table
#' @param pheno The phenotype matrix created by processPheno().
#' @param jx_file The path to the Rail-RNA junction coverage file if processing sample not already in recount2.
#' @param cores The number of processing cores to use.
#' @return A rse of the quantified tx abundances in the samples in pheno.
#' @export
#' @keywords recountNNLS
processReadLength = function(rl, pheno, jx_file, cores){
      message(Sys.time(), paste0(" ### Processing read length group: ", rl))
      pheno = pheno[pheno$rls_group==rl,,drop=FALSE]

      ## Load emission probability matrices
      message(Sys.time(), " # Setting up model covariates")
      data(list=paste0("matrix_", rl), package = "recountNNLSdata")
      matrix_list <- eval(parse(text=paste0("matrix_", rl)))
      genes = names(matrix_list)

      ## Load feature counts
      message(Sys.time(), " # Compiling feature counts")
      counts = getCounts(pheno[pheno$rls_group==rl,], jx_file, cores = cores)

      ## Run the NNLS
      message(Sys.time(), " # Executing model")
      if(rl<=75){
            data(list=paste0("g2l_75"), package = "recountNNLSdata")
            g2l <- eval(parse(text=paste0("g2l_75")))
      }else{
            data(list=paste0("g2l_", rl), package = "recountNNLSdata")
            g2l <- eval(parse(text=paste0("g2l_", rl)))
      }
      loci <- unique(g2l$locus)
      info <- mclapply(loci, .inferReads, g2l, matrix_list, counts, power=1, mc.cores = cores)

      message(Sys.time(), " # Compiling regression information")
      reads = do.call(rbind, sapply(info, function(x)x[[1]]))
      vars = do.call(rbind, sapply(info, function(x)x[[2]]))
      scores = do.call(rbind, sapply(info, function(x)x[[3]]))

      norm_matrix = matrix(rep(as.numeric(pheno$rls), times = dim(reads)[1]), byrow=TRUE, ncol = dim(reads)[2])
      pe_matrix = (matrix(rep(as.numeric(pheno$paired_end), times = dim(reads)[1]), byrow=TRUE, ncol = dim(reads)[2]))*1+1
      norm_matrix = rl/norm_matrix/pe_matrix
      reads = reads*norm_matrix
      vars = vars*norm_matrix^2
      se = sqrt(apply(vars, 2, as.numeric))

      ## Wrap up results in a RSE
      message(Sys.time(), " # Wrap up results in RSE")
      data(tx_grl, package = "recountNNLSdata")
      unquantified = names(tx_grl)[names(tx_grl) %in% rownames(reads)==FALSE]
      uq_info = matrix(NA, ncol = ncol(reads), nrow=length(unquantified))
      reads = rbind(reads, uq_info); se = rbind(se, uq_info)
      ind = match(names(tx_grl), rownames(reads))
      reads = reads[ind,,drop=FALSE]
      se = se[ind,,drop=FALSE]
      scores = scores[ind,,drop=FALSE]

      rownames(se) = NULL; colnames(se) = pheno$run
      rownames(reads) = NULL; colnames(reads) = pheno$run
      rownames(scores) = NULL; colnames(scores) = pheno$run
      assays = list(counts=reads, ses=se, scores = scores)
      data(tx_info, package='recountNNLSdata')
      rse_tx = SummarizedExperiment::SummarizedExperiment(assays=assays, rowRanges=tx_grl, colData=pheno)
      SummarizedExperiment::rowData(rse_tx) = tx_info[match(rownames(rse_tx), tx_info$tx_name),]
      colnames(rse_tx) = pheno$run
      return(rse_tx)
}

.inferReads = function(locus, g2l, ems, counts, power){
      coprocess = as.character(g2l$gene_id[g2l$locus==locus])
      if(length(coprocess)>1){
            ems_sub = ems[coprocess]
            P = .mergeP(ems_sub)
      }else{
            P = ems[[coprocess]]
      }
      ## Proceed to estimate information
      b = NULL; Vb = NULL; scores=NULL
      ## If the emission matrix is not empty
      if(length(P)>0){
            mat = match(rownames(P), rownames(counts))
            ## If there are observed features in the sample
            if(sum(!is.na(mat))>0){
                  P = P[which(!is.na(mat)),,drop=FALSE]
                  counts_sub = counts[mat[!is.na(mat)],, drop=FALSE]
                  # Weighting
                  P_binary = P>0
                  P_bin_sum = apply(P_binary, 1, sum)
                  P_weight = 1/P_bin_sum^power
                  P_weighted = P*P_weight
                  counts_weighted = counts_sub*P_weight
                  ## Use lm to assist in calculating colinearity
                  lm_mods = matrix(apply(counts_weighted, 2, .llm, P_weighted))
                  txs = sapply(lm_mods, function(x) rownames(summary(x)$coefficients))
                  colinear = which(colnames(P_weighted) %in% txs==FALSE)
                  ## If there are colinear transcripts
                  if(length(colinear)>0){
                        colinear_tx = colnames(P_weighted)[colinear]
                        P_weighted = P_weighted[,-colinear,drop=FALSE]
                  }
                  nnls_mods = apply(counts_weighted, 2, .lnnls, P_weighted)
                  for(i in 1:length(nnls_mods)){
                        lm_mods[[i]]$residuals <- nnls_mods[[i]]$residuals
                  }
                  b = matrix(sapply(nnls_mods, function(x) x$x), ncol=dim(counts_sub)[2])
                        rownames(b) = colnames(P_weighted)
                  Vb = matrix(sapply(nnls_mods, .robustVar, P_weighted), ncol=dim(counts_sub)[2])
                  if(dim(P_weighted)[2]>1){
                        scores = .scoreNNLS(P_weighted)
                  }else{
                        scores = 1
                  }
                  ## Adding info for the colinear transcripts
                  if(length(colinear)>0){
                        b_nas = matrix(rep(NA, length(colinear)*dim(b)[2]), ncol = dim(b)[2])
                              rownames(b_nas) = colinear_tx
                        Vb_nas = matrix(rep(NA, length(colinear)*dim(b)[2]), ncol = dim(Vb)[2])
                        scores_nas = rep(0, length(colinear))
                        b = rbind(b, b_nas)
                        Vb = rbind(Vb, Vb_nas)
                        scores = c(scores, scores_nas)
                  }
                  scores = matrix(rep(scores, dim(b)[2]), ncol=dim(b)[2])
            }
      }
      return(list(b, Vb, scores))
}
.llm = function(data, matrix){
      mod=lm(data~.-1, data=matrix)
      return(mod)
}
.robustVar = function(mod, P){
      res = mod$residuals
      X = apply(P, 2, as.numeric)
      XTXinv = solve(t(X) %*% X, tol=0)
      Hat = X %*% XTXinv %*% t(X)
      h = diag(Hat)
            h = sapply(h, function(x) min(x, 0.99))
      n = dim(P)[1]
      p = dim(P)[2]
      d = sapply(n*h/p, function(x) min(x, 4))
      inflate = res^2/(1-h)^d
      if(length(inflate)==1){
            Sigma = inflate
      }else{
            Sigma = diag(as.numeric(inflate))
      }
      bread = XTXinv %*% t(X)
      vcov = bread %*% Sigma %*% t(bread)
      return(diag(vcov))
}
.lnnls = function(data, matrix, boot=100){
      mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
      mod=nnls::nnls(mat, data)
      return(mod)
}
.mergeP = function(P_list){
      empty = which(sapply(P_list, length)==0)
      if(length(empty)>0)
            P_list = P_list[-empty]
      if(length(P_list)>1){
            for(i in 1:length(P_list)){
                  P_list[[i]]$rn = rownames(P_list[[i]])
            }
            P_out = P_list[[1]]
            for(i in 2:length(P_list)){
                  P_out = merge(P_out, P_list[[i]], by="rn", all=TRUE)
            }
            P_out[is.na(P_out)] = 0
            rownames(P_out) = P_out$rn
            P_out$rn = NULL
            return(P_out)
      }
      return(P_list)
}
.scoreNNLS = function(P){
      mod_list = list()
      Pmat = apply(P, 2, as.numeric)
      for(i in 1:dim(P)[2]){
            mod_list[[i]] = nnls::nnls(Pmat[,-i,drop=F], Pmat[,i])
      }
      res = sapply(mod_list, function(x) sum(x$residuals^2))
      base = apply(Pmat^2, 2, sum)
      return(res/base)
}
