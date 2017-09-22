### Calculate the effect of the winner's curse across
###
### @param boot_eGenes effect size estimates from the bootstrap analysis
### @param sig_assocs effect size estimates from the mapping analysis
### @param estimator estimator type
correct_winners_curse <- function(boot_eGenes, sig_assocs, estimator="shrinkage") {
  # CRAN NOTE suppression:
  gene <- NULL
  snps <- NULL
  detection_beta <- NULL
  estimation_beta <- NULL
  bootstrap <- NULL
  nominal_beta <- NULL
  winners_curse <- NULL
  corrected_beta <- NULL
  correction_boots <- NULL
  winners_curse <- NULL
  corrected_beta <- NULL
  correction_boots <- NULL

  tryCatch({ # Failure condition: no significant bootstraps
    effect_sizes <- merge(boot_eGenes[,list(gene, snps, detection_beta, estimation_beta, bootstrap)],
                          sig_assocs[,list(gene, snps, nominal_beta=beta)], by=c("gene", "snps"))

    # Force the effect size sign to always be in the same direction
    effect_sizes[, detection_beta := abs(detection_beta) * sign(nominal_beta)]
    effect_sizes[, estimation_beta := abs(estimation_beta) * sign(nominal_beta)]

    if (estimator == "shrinkage") {
      effect_sizes <- effect_sizes[, list(
        corrected_beta=unique(nominal_beta) - mean(detection_beta - estimation_beta),
        correction_boots=.N, nominal_beta=unique(nominal_beta)
      ), by=list(gene, snps)]
    } else if (estimator == "out_of_sample") {
      effect_sizes <- effect_sizes[, list(
        corrected_beta=mean(estimation_beta),
        correction_boots=.N, nominal_beta=unique(nominal_beta)
      ), by=list(gene, snps)]
    } else if (estimator == "weighted") {
      effect_sizes <- effect_sizes[, list(
        corrected_beta=0.368*unique(nominal_beta) + 0.632 * mean(estimation_beta),
        correction_boots=.N, nominal_beta=unique(nominal_beta)
      ), by=list(gene, snps)]
    }
    effect_sizes[, winners_curse := nominal_beta - corrected_beta]

    sig_assocs <- merge(sig_assocs,
                        effect_sizes[,list(gene, snps, winners_curse, corrected_beta, correction_boots)],
                        by=c("gene", "snps"))
  }, error=function(e) {
    warning("No significant eSNPs were significant in any bootstrap", immediate.=TRUE)
    sig_assocs[, winners_curse := NA_real_]
    sig_assocs[, corrected_beta := NA_real_]
    sig_assocs[, correction_boots := NA_real_]
  })
  sig_assocs[is.na(correction_boots), correction_boots := 0] # For genes where no bootstraps are significant
  return(sig_assocs)
}
