#' Bootstrap eQTL analysis for accurate effect size estimation
#'
#' Performs bootstrap analysis on significant eQTL eGenes to obtain
#' accurate eQTL effect size estimates by correcting for the "Winner's
#' Curse".
#'
#' @param snps see MatrixEQTL
#' @param gene see MatrixEQTL
#' @param cvrt see MatrixEQTL
#' @param n_bootstraps number of bootstraps to run
#' @param n_cores number of cores to parallise the bootstrap procedure
#'  over
#' @param ... arguments for MatrixEQTL
#'
#' @import foreach
#' @import data.table
#' @import MatrixEQTL
#'
#' @export
#'
#' @examples
#' # Locations for example data from the MatrixEQTL package
#' base.dir = find.package('MatrixEQTL');
#' SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
#' snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="");
#' expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
#' gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");
#' covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");
#'
#' # Load the SNP data
#' snps = SlicedData$new();
#' snps$fileDelimiter = "\t";      # the TAB character
#' snps$fileOmitCharacters = "NA"; # denote missing values;
#' snps$fileSkipRows = 1;          # one row of column labels
#' snps$fileSkipColumns = 1;       # one column of row labels
#' snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
#' snps$LoadFile(SNP_file_name);
#'
#' # Load the data detailing the position of each SNP
#' snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
#'
#' # Load the gene expression data
#' gene = SlicedData$new();
#' gene$fileDelimiter = "\t";      # the TAB character
#' gene$fileOmitCharacters = "NA"; # denote missing values;
#' gene$fileSkipRows = 1;          # one row of column labels
#' gene$fileSkipColumns = 1;       # one column of row labels
#' gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
#' gene$LoadFile(expression_file_name);
#'
#' # Load the data detailing the position of each gene
#' genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
#'
#' # Load the covariates data
#' cvrt = SlicedData$new();
#' cvrt$fileDelimiter = "\t";      # the TAB character
#' cvrt$fileOmitCharacters = "NA"; # denote missing values;
#' cvrt$fileSkipRows = 1;          # one row of column labels
#' cvrt$fileSkipColumns = 1;       # one column of row labels
#' if(length(covariates_file_name)>0) {
#'   cvrt$LoadFile(covariates_file_name);
#' }
#'
#' # Run the BootstrapEQTL analysis
#' eGenes <- BootstrapEQTL(snps, gene, cvrt, snpspos, genespos)
#'
BootstrapEQTL <- function(snps, gene, cvrt=SlicedData$new(),
                          snpspos, genespos, n_bootstraps=200, n_cores=1) {

  # Set up parallel computing environment
  par_setup <- setupParallel(n_cores, verbose=TRUE, reporterCore=FALSE)
  on.exit({
    cleanupCluster(par_setup$cluster, par_setup$predef)
  })

  # Check if the user has already loaded data.table: if not, load it and
  # make sure we return the table as a data.frame
  has.data.table <- "data.table" %in% names(sessionInfo()$otherPkgs)
  if (!has.data.table) {
    pkgReqCheck("data.table") # silently load without tutorial message
  }

  # Run Matrix eQTL to determine significant eGenes and get nominal
  # estimates for their effect sizes
  eQTLs <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    pvOutputThreshold = 0, # we don't need the Trans eQTL pvalues
    errorCovariance = numeric(),
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = 1e6,
    output_file_name=NULL,
    output_file_name.cis=NULL,
    noFDRsaveMemory=FALSE,
    useModel=modelLINEAR,
    verbose=FALSE
  );

  # Apply hierarchical multiple testing correction
  cis_assocs <- as.data.table(eQTLs$cis$eqtls) # Note need to read in file if output_file_name.cis != NULL
  cis_assocs[, corrected_pval := p.adjust(pvalue, method="bonferroni"), by=gene] # local SNP correction
  cis_assocs <- cis_assocs[order(abs(statistic), decreasing=TRUE)] # If multiple SNPs have the smallest pvalue, make sure the one with the biggest T-statistic is selected
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
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = 1e6,
        output_file_name=NULL,
        output_file_name.cis=NULL,
        noFDRsaveMemory=FALSE,
        verbose=FALSE
      )

      # Do hierarchical multiple testing correction to determine
      # significant eGenes
      detection_cis_assocs <- as.data.table(eQTL_detection$cis$eqtls)
      detection_cis_assocs[,pvalue := p.adjust(pvalue, method="bonferroni"), by=gene] # local SNP correction

      # We want to keep two SNPs: the top SNP in this bootstrap, and the
      # SNP that was detected by the initial eQTL analysis.
      # The top SNP is used to determine the most likely causal SNP, and
      # the eQTL SNP is used to correct the effect size estimate.
      detection_cis_assocs <- detection_cis_assocs[order(abs(statistic), decreasing=TRUE)] # If multiple SNPs have the smallest pvalue, make sure the one with the biggest T-statistic is selected
      detection_top_SNPs <- detection_cis_assocs[, .SD[which.min(pvalue)], by="gene"] # take row with best p-value per gene
      detection_top_SNPs[, pvalue := p.adjust(pvalue, method="BH")] # global gene correction
      detection_top_SNPs[, snp_type := "top"]

      detection_eQTL_SNPs <- merge(detection_cis_assocs, eGenes[,list(gene, snps=top_snp)], by=c("gene", "snps"))
      detection_eQTL_SNPs[, pvalue := p.adjust(pvalue, method="BH")] # global gene correction
      detection_eQTL_SNPs[, snp_type := "eQTL"]

      detection_eGene_assocs <- rbind(detection_top_SNPs, detection_eQTL_SNPs)

      # Filter to only the genes that were significant in the original
      # eQTL analysis and are also significant in the bootstrap
      # detection group
      sig_boot_assocs <- detection_eGene_assocs[gene %in% eGenes[corrected_pval < 0.05, gene] & pvalue < 0.05]
      sig_boot_assocs <- sig_boot_assocs[, list(gene, snps, snp_type, detection_beta=beta)]

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
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = 1e6,
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
      sig_boot_assocs[,error := NA_character_]

      return(sig_boot_assocs)
    }, error=function(e) {
      return(data.table(gene=NA, snps=NA, snp_type=NA, detection_beta=NA,
                        estimation_beta=NA, bootstrap=id_boot,
                        error=gsub("\n", " ", e$message)))
    })
  }

  # Report failed bootstraps if any
  if (boot_eGenes[,sum(!is.na(error))] > 0) {
    errors <- boot_eGenes[!is.na(error), .N, by=error]
    for (ii in errors[,.I]) {
      warning(errors[ii, N], "/", n_bootstraps, " bootstraps failed with",
              " error '", errors[ii, error], "' in the bootstrap",
              " estimation group.")
    }
    n_bootstraps <- n_bootstraps - boot_eGenes[,sum(!is.na(error))]
  }
  boot_eGenes <- boot_eGenes[is.na(error)]

  # Adjust for winner's curse - shrinkage estimator
  correction <- boot_eGenes[snp_type == "eQTL",
                            list(winners_curse=mean(detection_beta-estimation_beta)),
                            by=gene]
  eGenes <- merge(eGenes, correction, by="gene", all.x=TRUE)
  eGenes[, winners_curse := winners_curse]
  eGenes[, corrected_beta := nominal_beta - winners_curse]

  # Report the total number of significant bootstraps for each eGene
  correction_boots <- boot_eGenes[snp_type == "eQTL", list(correction_boots=.N), by=gene]
  eGenes <- merge(eGenes, correction_boots, by="gene", all.x=TRUE)

  # Take the SNP which was most frequently the top SNP across all
  # significant bootstraps as the probable causal SNP (best_eSNP).
  top_SNP <- boot_eGenes[snp_type == "top", .N, by=list(gene, snps)]
  top_SNP <- top_SNP[, .SD[which(N == max(N))], by=gene]
  top_SNP <- top_SNP[, list(best_boot_eSNP=paste(snps, collapse="; ")), by=gene]
  eGenes <- merge(eGenes, top_SNP, by="gene", all.x=TRUE)

  # Report proportion of significant bootstraps where the best eSNP was
  # the top SNP
  top_snp_props <- boot_eGenes[snp_type == "top", .N, by=list(gene, snps)]
  top_snp_props <- top_snp_props[, list(snps, prop_top_eSNP=N/sum(N)), by=gene]
  top_snp_props <- top_snp_props[, .SD[which.max(prop_top_eSNP)], by=gene]
  eGenes <- merge(eGenes, top_snp_props[, list(gene, prop_top_eSNP)], by="gene", all.x=TRUE)

  # Sort eGenes by significance
  eGenes <- eGenes[order(abs(statistic), decreasing=TRUE)][order(corrected_pval)]

  # If the user has not loaded data.table themselves cast back to a
  # data frame to avoid confusion
  if (!has.data.table) {
    eGenes <- as.data.frame(eGenes)
  }
  return(eGenes)
}
