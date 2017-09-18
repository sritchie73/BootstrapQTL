### Calculate the effect of the winner's curse across
###
### @param boot_eGenes effect size estimates from the bootstrap analysis
### @param eGenes effect size estimates from the eGene discovery analysis
### @param estimator estimator type
### @param force_sign logical; should the bootstrap effect sizes be forced
###        to always be in the same direction as the nominal effect size?
###        Must be TRUE when bootstrap_eSNP is "top".
correct_winners_curse <- function(boot_eGenes, eGenes, estimator="shrinkage", force_sign=FALSE) {
  # CRAN NOTE suppression:
  gene <- NULL
  detection_beta <- NULL
  estimation_beta <- NULL
  bootstrap <- NULL
  nominal_beta <- NULL
  winners_curse <- NULL
  corrected_beta <- NULL
  correction_boots <- NULL
  corrected_pval <- NULL

  effect_sizes <- merge(boot_eGenes[,list(gene, detection_beta, estimation_beta, bootstrap)],
                        eGenes[,list(gene, nominal_beta)], by="gene")

  if (force_sign) {
    effect_sizes[, detection_beta := abs(detection_beta) * sign(nominal_beta)]
    effect_sizes[, estimation_beta := abs(estimation_beta) * sign(nominal_beta)]
  }

  if (estimator == "shrinkage") {
    effect_sizes <- effect_sizes[, list(
      corrected_beta=unique(nominal_beta) - mean(detection_beta - estimation_beta),
      correction_boots=.N, nominal_beta=unique(nominal_beta)
    ), by=gene]
  } else if (estimator == "out_of_sample") {
    effect_sizes <- effect_sizes[, list(
      corrected_beta=mean(estimation_beta),
      correction_boots=.N, nominal_beta=unique(nominal_beta)
    ), by=gene]
  } else if (estimator == "weighted") {
    effect_sizes <- effect_sizes[, list(
      corrected_beta=0.368*unique(nominal_beta) + 0.632 * mean(estimation_beta),
      correction_boots=.N, nominal_beta=unique(nominal_beta)
    ), by=gene]
  }
  effect_sizes[, winners_curse := nominal_beta - corrected_beta]

  eGenes <- merge(eGenes, effect_sizes[,list(gene, winners_curse, corrected_beta, correction_boots)], by="gene", all.x=TRUE)
  eGenes[corrected_pval < 0.05 & is.na(correction_boots), correction_boots := 0] # In case no bootstraps are significant

  return(eGenes)
}
