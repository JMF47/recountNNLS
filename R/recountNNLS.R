recountNNLS = function(pheno, cores=1){
      rl = unique(pheno$rls_group)
      counts_ex = getExCounts(pheno)
      counts_jx = getJxCounts(pheno)
      counts = rbind(counts_ex, counts_jx)
      data(list=paste0("matrix_", rl), package = "recountNNLSdata")
      genes = names(matrix_list)
      # calculateReads(genes[1], matrix_list, counts, power=1)
      info = mclapply(genes, calculateReads, matrix_list, counts, junction_weight=rl, power=1, mc.cores = cores)

      reads = do.call(rbind, sapply(info, function(x)x[[1]]))
            norm_matrix = matrix(rep(as.numeric(pheno$avg_read_length), times = dim(reads)[1]), byrow=T, ncol = dim(reads)[2])
            norm_matrix = rl/norm_matrix
            reads = reads*norm_matrix
      vars = do.call(rbind, sapply(info, function(x) x[[2]]))
            vars = vars*norm_matrix^2
      colin = do.call(rbind, sapply(info, function(x) x[[3]]))
      se = sqrt(apply(vars, 2, as.numeric))

      colin_mat = matrix(rep(NA, dim(pheno)[1]*dim(colin)[1]), ncol=dim(pheno)[1], nrow=dim(colin)[1])
      rownames(colin_mat) = colin$transcript_id

      reads = rbind(reads, colin_mat)
      se = rbind(se, colin_mat)

      data(tx_grl, package = "recountNNLSdata")
      rowRanges = rng[match(rownames(reads), names(rng))]

      rownames(se) = NULL; colnames(se) = NULL
      rownames(reads) = NULL; colnames(reads) = NULL
      rse_tx = SummarizedExperiment(assays=list(counts=reads, se=se), rowRanges=rowRanges, colData=pheno)

      return(rse_tx)
}

# counts = apply(counts, 2, as.numeric)
# # colnames(counts) = NULL; rownames(counts) = NULL
# colnames(se) = NULL; rownames(se) = NULL
# colData = colData(rse_gene)
# colData = colData[match(colnames, rownames(colData)),]
