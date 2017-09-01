#'
#'
#' @param snps see MatrixEQTL
#' @param gene see MatrixEQTL
#' @param cvrt see MatrixEQTL
#' @param n_bootstraps number of bootstraps to run
#' @param ... arguments for MatrixEQTL

bootstrapQTL <- function(snps, gene, cvrt=SlicedData$new(), n_bootstraps=200, ...) {

  # Run Matrix eQTL to determine significant eGenes and get nominal
  # estimates for their effect sizes
  eQTLs <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    pvOutputThreshold = 0, # we don't need the Trans eQTL pvalues
    useModel = useModel,
    errorCovariance = errorCovariance,
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    output_file_name=NULL,
    output_file_name.cis=NULL,
    noFDRsaveMemory=FALSE
  );

  # Apply hierarchical multiple testing correction
  cis_assocs <- as.data.table(eQTLs$cis$eqtls) # Note need to read in file if output_file_name.cis != NULL
  cis_assocs[, corrected_pval := p.adjust(pvalue, method="bonferroni"), by=gene] # local SNP correction
  eGenes <- cis_assocs[, .SD[which.min(corrected_pval)], by="gene"] # take row with best p-value per gene
  eGenes[, corrected_pval := p.adjust(pvalue, method="BH")] # global gene correction
  eGenes <- eGenes[,list(gene, top_snp=snps, statistic, nominal_beta=beta,
                         nominal_pval=pvalue, corrected_pval)]

  # Run MatrixEQTL in each bootstrap detection group and estimation
  # group
  boot_eGenes <- foreach(id_boot = seq_len(n_bootstraps),
                             .inorder = FALSE,
                             .packages = c("MatrixEQTL", "data.table"),
                             .combine=rbind
  ) %dopar% {

    tryCatch({
      sample_size <- snps$nCols()
      id_detection <- sample(seq_len(sample_size), sample_size, replace = TRUE)
      id_estimation <- setdiff(seq_len(sample_size), id_detection)

      # Copy and subset the gene and snp data for the bootstrap
      # detection group
      gene_detection <- gene$Clone()
      gene_detection$ColumnSubsample(id_detection)
      snps_detection <- snps$Clone()
      snps_detection$ColumnSubsample(id_detection)
      cvrt_detection <- cvrt$Clone()
      cvrt_detection$ColumnSubsample(id_detection)

      # Run MatrixEQTL on the detection group
      # Revisit to remove unneccesary parameters
      eQTL_detection <- Matrix_eQTL_main(
        snps = snps_detection,
        gene = gene_detection,
        cvrt = cvrt_detection,
        pvOutputThreshold = 0, # we don't need the Trans eQTL pvalues
        useModel = useModel,
        errorCovariance = errorCovariance,
        pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        output_file_name=NULL,
        output_file_name.cis=NULL,
        noFDRsaveMemory=FALSE,
        verbose=FALSE
      )

      # Do hierarchical multiple testing correction to determine
      # significant eGenes
      detection_cis_assocs <- as.data.table(eQTL_detection$cis$eqtls)
      detection_cis_assocs[,pvalue := p.adjust(pvalue, method="bonferroni"), by=gene] # local SNP correction
      detection_eGene_assocs <- detection_cis_assocs[, .SD[which.min(pvalue)], by="gene"] # take row with best p-value per gene
      detection_eGene_assocs[, pvalue := p.adjust(pvalue, method="BH")] # global gene correction

      # Filter to only the significant genes and the columns we care about
      sig_boot_assocs <- detection_eGene_assocs[pvalue < 0.05]
      sig_boot_assocs <- sig_boot_assocs[, list(gene, snps, detection_beta=beta)]

      # If there are no significant eGenes in this bootstrap we can move
      # to the next one
      if(nrow(sig_boot_assocs) == 0) {
        return(NULL)
      }

      # Run MatrixEQTL in the estimation group. We only care about the
      # significant SNP-Gene pairs from the detection group, so we can
      # save computation time by filtering their respective matrices.
      # Note: where multiple significant genes or snps are located within
      # the same cis window the filtering will not be perfect.
      sig_boot_genes <- sig_boot_assocs$gene
      sig_boot_snps <- sig_boot_assocs$snps

      gene_estimation <- gene$Clone()
      gene_estimation$RowReorder(which(gene_estimation$GetAllRowNames() %in% sig_boot_genes))
      gene_estimation$ColumnSubsample(id_estimation)
      snps_estimation <- snps$Clone()
      snps_estimation$RowReorder(which(snps_estimation$GetAllRowNames() %in% sig_boot_snps))
      snps_estimation$ColumnSubsample(id_estimation)
      cvrt_estimation <- cvrt$Clone()
      cvrt_estimation$ColumnSubsample(id_estimation)

      eQTL_estimation <- Matrix_eQTL_main(
        snps = snps_estimation,
        gene = gene_estimation,
        cvrt = cvrt_estimation,
        pvOutputThreshold = 0, # we don't need the Trans eQTL pvalues
        useModel = useModel,
        errorCovariance = errorCovariance,
        pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        output_file_name=NULL,
        output_file_name.cis=NULL,
        noFDRsaveMemory=FALSE,
        verbose=FALSE
      )

      # Add the estimation beta to the significant bootstrap association
      # table
      estimation_cis_assocs <- as.data.table(eQTL_estimation$cis$eqtls)
      estimation_assocs <- estimation_cis_assocs[, list(gene, snps, estimation_beta=beta)]

      # Note the merge step filters out cases where the SNP/Gene filtering
      # lead to multiple SNPs for one gene
      sig_boot_assocs <- merge(sig_boot_assocs, estimation_assocs, by=c("gene", "snps"))
      sig_boot_assocs[,bootstrap := id_boot]

      return(sig_boot_assocs)
    }, error=function(e) {
      warning(e$message)
      return(data.table(gene=NA, snps=NA, detection_beta=NA, estimation_beta=NA, bootstrap=id_boot))
    })
  }
  failed_bootstraps <- boot_eGenes[is.na(gene), length(unique(bootstrap))]
  warning(failed_bootstraps, "/", n_bootstraps, " bootstraps failed ",
          " due to insufficient study sample size (see warnings).")

  # Adjust for winner's curse

  return(eGenes)
}
