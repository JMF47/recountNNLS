#' NNLS for Tx Abundance Calculation
#'
#' This function applies runs our NNLS model for a given subset
#' of the pheno matrix that has the same read length specified
#' and reports an rse.
#' @param rl A read length that exists in the pheno table
#' @param pheno The phenotype matrix created by processPheno().
#' @param counts_ex The matrix of exonic feature counts where each row represents a exonic feature and each column a sample.
#' @param counts_jx The matrix of junction feature counts where each row represents a junction feature and each column a sample.
#' @param cores The number of processing cores to use.
#' @return A rse of the quantified tx abundances in the samples in pheno.
#' @keywords recountNNLS
processReadLength = function(rl, pheno, counts_ex, counts_jx, cores){
      message(Sys.time(), paste0(" ### Processing read length group: ", rl))
      pheno = pheno[pheno$rls_group==rl,,drop=FALSE]

      ## Create the appropriate count matrix
      message(Sys.time(), " # Compiling feature counts")
      if(is.null(counts_ex)) # If we have to calculate counts_ex from bw files
            counts_ex = getExCounts(pheno, cores = cores)
      if(!is.null(counts_jx)){ # If counts_jx is not empty
            counts = rbind(counts_ex, counts_jx[, match(colnames(counts_ex), colnames(counts_jx)), drop=FALSE])
      }else{ # If counts_jx is empty
            counts = counts_ex
      }

      ## Load emission probability matrices
      message(Sys.time(), " # Setting up model covariates")
      data(list=paste0("matrix_", rl), package = "recountNNLSdata")
      matrix_list <- eval(parse(text=paste0("matrix_", rl)))

      ## Run the NNLS
      message(Sys.time(), " # Executing model")
      info = mclapply(unique(g2l$locus), .calculateReads, g2l, matrix_list, counts, power=1, mc.cores = cores)

            message(Sys.time(), " # Compiling regression information")
            reads = do.call(rbind, sapply(info, function(x)x[[1]]))
            norm_matrix = matrix(rep(as.numeric(pheno$rls), times = dim(reads)[1]), byrow=TRUE, ncol = dim(reads)[2])
            norm_matrix = rl/norm_matrix
            reads = reads*norm_matrix
            vars = do.call(rbind, sapply(info, function(x) x[[2]]))
            vars = vars*norm_matrix^2
            colin = do.call(rbind, sapply(info, function(x) x[[3]]))
            se = sqrt(apply(vars, 2, as.numeric))

            ## Process transcripts that are colinear
            colin_mat = matrix(rep(NA, dim(pheno)[1]*dim(colin)[1]), ncol=dim(pheno)[1], nrow=dim(colin)[1])
            rownames(colin_mat) = colin$transcript_id
            reads = rbind(reads, colin_mat)
            se = rbind(se, colin_mat)
            rownames(se) = rownames(reads)

      ## Wrap up results in a RSE
      message(Sys.time(), " # Wrap up results in RSE")
      data(tx_grl, package = "recountNNLSdata")
      unquantified = names(tx_grl)[names(tx_grl) %in% rownames(reads)==FALSE]
      uq_info = matrix(0, ncol = ncol(reads), nrow=length(unquantified))
      reads = rbind(reads, uq_info); se = rbind(se, uq_info)
      ind = match(names(tx_grl), rownames(reads))
      reads = reads[ind,,drop=FALSE]; se = se[ind,,drop=FALSE]

      rownames(se) = NULL; colnames(se) = pheno$run
      rownames(reads) = NULL; colnames(reads) = pheno$run
      rse_tx = SummarizedExperiment::SummarizedExperiment(assays=list(counts=reads, se=se), rowRanges=tx_grl, colData=pheno)
      return(rse_tx)
}

## Gene-wise execution of NNLS
.calculateReads = function(locus, g2l, ems, counts, power, verbose=FALSE){
      if(verbose==TRUE)
            message(locus)
      b = NULL; Vb = NULL; colinear_info = NULL
      genes = as.character(g2l$gene[g2l$locus==locus])
      ems_sub = ems[match(genes, names(ems))]
      ems_sub = ems_sub[sapply(ems_sub, length)>0]
      if(length(ems_sub)>1)
            P = .mergeP(ems_sub)
      if(length(ems_sub)==1)
            P = ems_sub[[1]]
      if(length(ems_sub)==0)
            P = NULL

      ## If the emission matrix is not empty
      if(length(P)>0){
            mat = match(rownames(P), rownames(counts))
            P_size_orig = dim(P)[2]
            if(sum(!is.na(mat))>0){
                  ## Determine weighting of the uniqueness of features
                  P = P[which(!is.na(mat)),,drop=FALSE]
                  P_binary = P>0
                  P_bin_sum = apply(P_binary, 1, sum)
                  P_weight = 1/P_bin_sum^power
                  counts_sub = counts[mat[!is.na(mat)],, drop=FALSE]
                  counts_sub = counts_sub*P_weight
                  P = P*P_weight

                  ## Quantification
                  colinear_info = NULL
                  if(dim(P)[2]>1){ ## More than 1 isoform at this location
                        ## Use lm to assist in calculating colinearity and SE
                        lm_info = matrix(apply(counts_sub, 2, .llm, P))
                        beta = sapply(lm_info, function(x) x[[1]])
                        Vbeta = sapply(lm_info, function(x) x[[2]])^2
                        colinear = which(colnames(P) %in% rownames(beta)==FALSE)
                        if(length(colinear)>0){ ## If there are colinear transcripts
                              colinear_tx = colnames(P)[colinear]
                              P = P[,-colinear]
                              colinear_info = data.frame(locus_id=locus, transcript_id = colinear_tx,
                                                         P_dim=P_size_orig, P_colin=length(colinear))
                        }
                        if(dim(P)[2]>1){ ## If there are non-colinear transcripts
                              b = matrix(apply(counts_sub, 2, .lnnls, P), nrow=dim(P)[2])
                              rownames(b) = colnames(P)
                              d = (1+(b-beta))
                              Vb = d^2*Vbeta
                        }else{ ## If all isoforms are colinear
                              b = NULL
                              Vb = NULL
                        }
                  }else{ ## Just 1 isoform at this location
                        b = matrix(apply(counts_sub, 2, .lnnls, P), nrow=dim(P)[2])
                        rownames(b) = colnames(P)
                        beta = b
                        P = apply(P, 2, as.numeric)
                        fitted = P%*% beta
                        res = counts_sub-fitted
                        res2sum = apply(res, 2, function(x) sum(x^2))
                        Vb = res2sum/(dim(res)[1]-1)
                        Vb = matrix(Vb, nrow=1)
                        rownames(Vb) = colnames(P)
                  }

            }
      }
      return(list(b, Vb, colinear_info))
}

## Wrapper function for lm
.llm = function(data, matrix){
      # mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
      mod=lm(data~.-1, data=matrix)
      mod_out = summary(mod)$coefficients
      return(list(mod_out[,1], mod_out[,2]))
}

## Wrapper function for nnls
.lnnls = function(data, matrix){
      mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
      mod=nnls::nnls(mat, data)
      out = mod$x
      return(out)
}

## Wrapper function to combine emission matrices by locus
.mergeP = function(P_list){
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
