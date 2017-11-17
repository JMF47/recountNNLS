calculateReads = function(gene, ems, counts, junction_weight, power, verbose=TRUE){
      if(verbose==TRUE){
            ind = which(names(ems)==gene)
            if((ind%%1000)==0){
                  print(ind); flush.console()
            }
      }
      b = NULL; Vb = NULL; colinear_info = NULL
      P = ems[[gene]]
      ## If there emission matrix is not empty
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

                  ## Rescale the weight of the junctions
                  junction_ind = stringr::str_detect(rownames(counts_sub), "i")
                  P_weight2 = junction_ind*(junction_weight-1)+1
                  counts_sub = counts_sub*P_weight2
                  P = P*P_weight2

                  ## Quantification
                  colinear_info = NULL
                  if(dim(P)[2]>1){ ## More than 1 isoform at this location
                        ## Use lm to assist in calculating colinearity and SE
                        lm_info = matrix(apply(counts_sub, 2, .llm, P))
                        beta = sapply(lm_info, function(x) x[[1]])
                        Vbeta = sapply(lm_info, function(x) x[[2]])^2
                        colinear = which(colnames(P) %in% rownames(beta)==F)
                        if(length(colinear)>0){ ## If there are colinear transcripts
                              colinear_tx = colnames(P)[colinear]
                              P = P[,-colinear]
                              colinear_info = data.frame(gene_id=gene, transcript_id = colinear_tx,
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
      mod=nnls(mat, data)
      out = mod$x
      return(out)
}
