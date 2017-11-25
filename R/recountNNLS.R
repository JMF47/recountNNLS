#' NNLS for Tx Abundance Calculation
#'
#' This function apply NNLS to calculate the Tx abundance
#' from reduced-representation of coverage information.
#' @param pheno The a table of phenotype information from processPheno().
#' @param cores The number of processing cores to use.
#' @keywords recountNNLS
#' @export
recountNNLS = function(pheno, counts_ex, counts_jx, cores=1){
      ## Stack the counts matrix for input
      counts = rbind(counts_ex, counts_jx)

      ## Load emission probability matrices
      rl = unique(pheno$rls_group)
      if(length(rl)>1)
            stop("Cannot process two different read group lengths at the same time. Please split analysis.")
      data(list=paste0("matrix_", rl), package = "recountNNLSdata")
      genes = names(matrix_list)

      ## NNLS
      info = mclapply(genes, calculateReads, matrix_list, counts, junction_weight=rl, power=1, mc.cores = cores)
      reads = do.call(rbind, sapply(info, function(x)x[[1]]))
            norm_matrix = matrix(rep(as.numeric(pheno$avg_read_length), times = dim(reads)[1]), byrow=T, ncol = dim(reads)[2])
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

      ## Wrap up results in a RSE
      data(tx_grl, package = "recountNNLSdata")
      rowRanges = rng[match(rownames(reads), names(rng))]
      rownames(se) = NULL; colnames(se) = NULL
      rownames(reads) = NULL; colnames(reads) = NULL
      rse_tx = SummarizedExperiment(assays=list(counts=reads, se=se), rowRanges=rowRanges, colData=pheno)
      return(rse_tx)
}
