### Apply hierarchical multiple testing correction
###
### @param cis_assocs data.table of cis-associations from MatrixEQTL
### @param local multiple testing correction method to use at each gene.
### @param global multiple testing correction method to us across all genes.
### @param tests_per_gene a data.table providing the number of cis SNPs / 
###        independent tests for each gene.
###
### @return a data.table containing only significant eGenes or significant
###   eSNPs
hierarchical_correction <- function(cis_assocs, local, global, tests_per_gene=NULL) {
  # Suppress CRAN notes about data.table columns
  local_pval <- NULL
  pvalue <- NULL
  n_tests <- NULL
  gene <- NULL
  statistic <- NULL
  global_pval <- NULL

  # Apply local correction across SNPs at each gene
  if (!is.null(tests_per_gene)) { # i.e. Bonferroni when test only performed on some SNPs
    cis_assocs <- merge(cis_assocs, tests_per_gene, by="gene")
    cis_assocs[, local_pval := adjust_p(pvalue, method=local, N=unique(n_tests)), by=gene]
  } else {
    cis_assocs[, local_pval := adjust_p(pvalue, method=local), by=gene]
  }

  # Apply global correction across genes using the best p-value per gene
  cis_assocs <- cis_assocs[order(abs(statistic), decreasing=TRUE)]
  top_assocs <- cis_assocs[, .SD[which.min(local_pval)], by="gene"]
  top_assocs[,global_pval := adjust_p(local_pval, method=global)]
  cis_assocs <- merge(cis_assocs, top_assocs[,list(gene, global_pval)], by="gene")

  return(cis_assocs)
}

### Adjust p-values for multiple tests
adjust_p <- function(pvals, method, N=NULL) {
  if (method == "qvalue") {
    qvalue::qvalue(pvals)$qvalues
  } else if (method == "eigenMT") {
    pmin(pvals * N, 1)
  } else {
    p.adjust(pvals, method, n=ifelse(is.null(N), length(pvals), N))
  }
}

### Determine the eSNP significance threshold after hiearchical correction
###
### @param cis_assocs a table of hierarchically corrected associations
get_eSNP_threshold <- function(cis_assocs) {
  # Suppress CRAN NOTES
  statistic <- NULL
  local_pval <- NULL
  global_pval <- NULL

  # What is the locally corrected p-value corresponding to the global
  # correction threshold of 0.05?
  cis_assocs <- cis_assocs[order(abs(statistic), decreasing=TRUE)]
  top_assocs <- cis_assocs[, .SD[which.min(local_pval)], by="gene"]
  top_assocs <- top_assocs[order(global_pval)]
  n_sig <- top_assocs[,sum(global_pval < 0.05)]
  eSNP_threshold <- top_assocs[n_sig:(n_sig+1), mean(local_pval)]
  return(eSNP_threshold)
}
