### Calculate frequency of top SNPs across all significant boostraps
###
### @param boot_eGenes effect size estimates from the bootstrap analysis
### @param eGenes effect size estimates from the eGene discovery analysis
### @param cis_assocs associations for all SNP-gene pairs from the discovery analysis
get_bootstrap_lead_SNPs <- function(boot_eGenes, eGenes, cis_assocs) {
  # CRAN NOTE suppression
  prop_top_eSNP <- NULL
  boots <- NULL
  snp_block <- NULL
  ld_prev <- NULL
  best_boot_eSNP <- NULL

  tryCatch({ # Failure condition: No significant eGenes in any bootstrap
    if (nrow(boot_eGenes) == 0) {
      stop("no significant eGenes")
    }

    # In cases where multiple SNPs occur equally freuqently we need to determine if
    # they are in perfect LD or not.
    top_SNP_count <- boot_eGenes[, .N, by=list(gene, snps)] # In how many bootstraps does each gene-SNP appear?
    sig_boot_count <- boot_eGenes[, list(boots=length(unique(bootstrap))), by=gene] # how many significant bootstraps for each gene?
    top_SNP <- merge(top_SNP_count, sig_boot_count, by="gene")
    top_SNP[, prop_top_eSNP := N/boots] # Divide the occurance of each gene-SNP pair by the number of significant bootstraps for that gene.
    top_SNP <- top_SNP[, .SD[which(prop_top_eSNP == max(prop_top_eSNP))], by=gene] # filter to the most frequent per gene

    # We can use the original cis_assocs table to determine whether top
    # eSNPs are in perfect LD or not. If they are, they will have
    # identical statistics, which we can detect using 'duplicated()'.
    # By sorting the table, we can also pick up independent groups of SNPs
    # in perfect LD that would be more difficult to detect using a
    # correlation matrix on minor allele dosage
    top_SNP <- merge(top_SNP[,list(gene, snps, prop_top_eSNP)],
                     cis_assocs[,list(gene, snps, statistic, beta, pvalue)],
                     by=c("gene", "snps"))
    top_SNP <- top_SNP[order(beta)][order(statistic)][order(pvalue)][order(gene)]
    top_SNP <- top_SNP[,list(snps, prop_top_eSNP, ld_prev=duplicated(.SD, by=c("statistic", "beta", "pvalue"))), by=gene]
    # Assign each distinct block a number so we can collapse information
    top_SNP[, snp_block := NA_real_]
    block <- 0
    for (ii in top_SNP[,.I]) {
      block <- ifelse(top_SNP[ii, ld_prev], block, block + 1)
      top_SNP[ii, snp_block := block]
    }

    # Now summarise at gene level: where there are multiple SNPs that
    # occurred with equal frequency across all bootstraps then SNPs in
    # perfect LD will be separated by a "/" while independent SNPs
    # or SNP blocks will be separated with a ";".
    top_SNP <- top_SNP[, list(gene=unique(gene), prop_top_eSNP=unique(prop_top_eSNP),
                              best_boot_eSNP=paste(snps, collapse="/")), by=snp_block]
    top_SNP <- top_SNP[, list(best_boot_eSNP=paste(best_boot_eSNP, collapse=";"),
                              prop_top_eSNP=unique(prop_top_eSNP)), by=gene]
    eGenes <- merge(eGenes, top_SNP, by="gene", all.x=TRUE)
  }, error=function(e) {
    warning("No significant eGenes were significant in any bootstrap", immediate.=TRUE)
    eGenes[, best_boot_eSNP := NA_real_]
    eGenes[, prop_top_eSNP := NA_real_]
  })
  return(eGenes)
}
