### Get the top eSNP(s) for each significant eGene
###
### In some cases there may be multiple top eSNPs with identical statistics
### e.g. due to perfect LD.
###
### @param cis_assocs data.table of cis-associations from MatrixEQTL
### @param eGenes data.table of significant eGenes.
### @param collapse logical; should multiple SNPs be collapsed into one row?
###
### @return a table containing one or more rows for each eGene with the top eSNPs
get_eSNPs <- function(cis_assocs, eGenes, collapse=FALSE) {
  top_blocks <- merge(cis_assocs[,list(gene, snps, statistic, beta, pvalue)],
                      eGenes[,list(gene, statistic, beta, pvalue)],
                       by=c("gene", "statistic", "beta", "pvalue"))
  if (collapse) {
    top_blocks <- top_blocks[, list(snps=paste(snps, collapse="/")), by="gene"]
  }
  eGenes[, snps := NULL]
  eGenes <- merge(eGenes, top_blocks[,list(gene, snps)], by="gene")
  return(eGenes)
}
