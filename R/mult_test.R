### Summarise MatrixEQTL output to the eGene level
###
### Take the top SNP at each gene and perform multiple testing
### correction across eGenes
###
### @param cis_assocs data.table of cis-associations from MatrixEQTL
### @param local multiple testing correction method to use at each gene.
### @param global multiple testing correction method to us across all genes.
### @param eSNPs optional table specifying eSNPs to use for each gene rather
###        than the top-SNP for global multiple testing correction
###
### @return a data.table containing only significant eGenes and their top SNP
get_eGenes <- function(cis_assocs, local, global,  eSNPs) {
  cis_assocs[, corrected_pval := adjust_p(pvalue, method=local), by=gene]

  if (!missing(eSNPs)) {
    eGenes <- merge(cis_assocs, eSNPs, by=c("gene", "snps"))
  } else {
    # Get SNP with smallest p-value for each gene/ If multiple SNPs have
    # the smallest pvalue, make sure the one with the biggest T-statistic
    # is selected
    cis_assocs <- cis_assocs[order(abs(statistic), decreasing=TRUE)]
    eGenes <- cis_assocs[, .SD[which.min(corrected_pval)], by="gene"] # take row with best p-value per gene
  }
  eGenes[, corrected_pval := adjust_p(corrected_pval, method=global)] # global gene correction


  return(eGenes)
}

### Adjust p-values for multiple tests
adjust_p <- function(pvals, method) {
  if (method == "qvalue") {
    qvalue(pvals)$qvalues
  } else {
    p.adjust(pvals, method)
  }
}
