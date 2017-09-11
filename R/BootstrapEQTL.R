#' Bootstrap eQTL analysis for accurate effect size estimation
#'
#' Performs a cis-eQTL analysis using MatrixEQTL with an additional
#' bootstrap analysis on significant eQTL eGenes to obtain accurate
#' eQTL effect size estimates by correcting for the "Winner's
#' Curse".
#'
#' @details
#'  \code{BootstrapEQTL} extends the \code{\link{MatrixEQTL}} package by
#'  performing hierarchical multiple testing correction that controls
#'  the eGene false discovery rate at 5\% and provides accurate eQTL
#'  effect sizes for significant eGenes by performing a bootstrap
#'  procedure that corrects for the overestimation of effect sizes
#'  ("The Winner's Curse effect").
#'
#'  \subsection{EGene detection:}{
#'  Detection of significant eGenes is performed using a hieararchical
#'  multiple testing correction procedure. At each gene, local multiple
#'  testing adjustment is performed using the method specified by the
#'  \code{'local_correction'} argument (Bonferroni correction by default)
#'  to correct for the number of \emph{cis}-SNPs tested at each gene.
#'  Global multiple testing adjustment is subsequently performed across
#'  all genes by taking the minimum locally corrected p-value at each
#'  gene, then performing the \code{'global_correction'} adjustment across
#'  all genes (FDR correction by default). The default settings best
#'  control ethe eGene false discovery rate without sacrificing
#'  sensitivity (see citation).
#'  }
#'  \subsection{Winner's curse adjustment:}{
#'  EQTL effect sizes of significant eGenes are typically overestimated
#'  when compared to replication datasets due to the selection of the
#'  best \emph{cis}-SNP at each gene. \code{BootstrapEQTL} removes this
#'  overestimation by performing a bootstrap procedure after significant
#'  eGene discovery by \code{\link{MatrixEQTL}}.
#'
#'  The bootstrap procedure provides an estimate of the effect size
#'  overestimation ("Winner's Curse") by partitioning the discovery
#'  dataset into two groups: an eQTL detection group, and an effect size
#'  estimation group. The eQTL detection group is determined by
#'  bootstrapping (sampling with replacement) the dataset samples, with
#'  the remaining samples (those not selected via the sampling with
#'  replacement) composing the eQTL estimation group.
#'  \code{\link{MatrixEQTL}} is then performed in the detection group in
#'  the same way as described in the "EGene discovery:" section above to
#'  determine whether significant eGenes remain significant in the
#'  bootstrapped detection group when considering the same SNP-gene
#'  pairs detected in the initial eGene discovery analysis. If the
#'  significant eGenes remain significant in the bootstrapped detection
#'  group, then the effect size of the SNP-Gene pair is calculated in
#'  the left-out estimation group. The difference between the effect size
#'  in the bootstrap detection group and the left-out estimation group
#'  is then used as an estimate of the "Winner's Curse" effect for that
#'  SNP-gene pair. This procedure is repeated several hundred times (see
#'  the \code{'n_bootstraps'} argument) to obtain a robust "Winner's
#'  Curse" effect for each SNP-Gene pair by averaging over significant
#'  bootstraps. This is reported in the \code{winners_curse} column of
#'  the returned table. The eQTL effect size corrected for the "Winner's
#'  Curse" is reported in the \code{corrected_beta} column.
#'
#'  Note that SNP-gene pairs may not remain significant in all
#'  bootstraps, so the effective number of bootstraps used obtain the
#'  "Winner's Curse" estimate will typically be lower than the number of
#'  bootstraps specified in \code{'n_bootstraps'}. The number of
#'  bootstraps that were significant for each eGene is reported in the
#'  \code{'correction_boots'} column of the returned table.
#'
#'  In addition to correcting for the Winner's Curse, the bootstrap
#'  procedure also collects statistics about which SNP(s) are the top
#'  SNP for each tested eGene in the bootstrapped detection group. This
#'  may differ from the top SNP from the initial eGene detection
#'  analysis in any given bootstrap. The \code{'best_boot_eSNP'} column
#'  reports the SNP that was the top SNP most frequently across the
#'  bootstrap procedure. There is some evidence from simulation studies
#'  that this is more likely to be the causal SNP than the
#'  \code{'top_snp'} from the eGene detection analysis (see citation),
#'  although the \code{'top_snp'} and \code{'best_boot_eSNP'} are
#'  typically the same. Multiple SNPs may be reported in the
#'  \code{'best_boot_eSNP'} column. SNPs separated by a "/" indicate
#'  multiple SNPs in perfect linkage disequilbrium (LD), \emph{i.e.}
#'  with identical test statistics/effect sizes/p-values in all
#'  bootstraps, while SNPs separated by a ";" indicate multiple SNPs (or
#'  SNP groups in perfect LD) that were the top SNP(s) with equal
#'  frequency across the bootstrap procedure (\emph{e.g.}
#'  "SNP1;SNP2/SNP3" could indicate SNP1 was the top SNP in 50\% of
#'  significant bootstraps and SNP2 and SNP3 are in perfect LD and were
#'  was the top SNP in the other 50\% of significant bootstraps). The
#'  \code{'prop_top_eSNP'} column reports the proportion of significant
#'  bootstraps for which the \code{'best_boot_eSNP'} was the top SNP.
#'  A low \code{'prop_top_eSNP'} (\emph{e.g.} < 0.5) can indicate
#'  complex LD structure around the top eSNP or multiple causal eSNPs.
#'  Note that here the number of significant bootstraps may be higher
#'  than that reported in \code{'correction_boots'} as an eGene may be
#'  significant when performing global FDR correction across the top
#'  bootstrap detection group SNP but not when performing global FDR
#'  correction on the top SNP from the origina eGene detection analysis.
#'  }
#'  \subsection{Warning messages:}{
#'  It is possible for bootstrap analyses to fail due to the reduced
#'  sample sizes of the bootstrap detection and bootstrap estimation
#'  groups. For example, the bootstrap resampling may lead to an
#'  estimation group in which all individuals are homozygous for the
#'  eQTL SNP or are all the same age or sex. This causes the eQTL
#'  analysis to crash since the linear models cannot be fit. These
#'  errors are collated at the end of the bootstrap procedure, grouped
#'  together, and reported as a series of warnings detailing the number
#'  of bootstraps in which the error occurred as well as whether they
#'  occurred in the detection or estimation groups.
#'  }
#'
#' @param snps \code{\link[MatrixEQTL]{SlicedData}} object containing genotype
#'   information used as input into \code{\link[MatrixEQTL]{Matrix_eQTL_main}}.
#' @param gene \code{\link[MatrixEQTL]{SlicedData}} object containing gene expression
#'   information used as input into \code{\link[MatrixEQTL]{Matrix_eQTL_main}}.
#' @param snpspos \code{data.frame} object with information about SNP locations.
#'   Used in conjunction with \code{'genespos'} and \code{'cisDist'} to
#'   determine SNPs in \emph{cis} of each gene. Must have three columns: \enumerate{
#'     \item 'snpid' describing the name of the SNP and corresponding to rows in
#'     the 'snps' matrix.
#'     \item 'chr' describing the chromosome for each SNP.
#'     \item 'pos' describing the position of the SNP on the chromosome.
#'   }
#' @param genepos \code{data.frame} object with information about transcript locations.
#'   Used in conjunction with \code{'snpspos'} and \code{'cisDist'} to
#'   determine SNPs in \emph{cis} of each gene. Must have four columns: \enumerate{
#'     \item 'geneid' describing the name of the gene and corresponding to rows in
#'     the 'gene' matrix.
#'     \item 'chr' describing the chromosome for each SNP.
#'     \item 'left' describing the start position of the transcript.
#'     \item 'right' describing the end position of the transcript.
#'   }
#' @param cvrt \code{\link[MatrixEQTL]{SlicedData}} object containing covariate
#'   information used as input into \code{\link[MatrixEQTL]{Matrix_eQTL_main}}.
#'   Argument can be ignored in the case of no covariates.
#' @param n_bootstraps number of bootstraps to run.
#' @param n_cores number of cores to parallise the bootstrap procedure
#'  over.
#' @param eGene_detection_file_name \code{character}, \code{connection} or \code{NULL}.
#'  File to save local \emph{cis} associations to in the eGene detection analysis. Corresponds
#'  to the \code{output_file_name.cis} argument in \code{\link[MatrixEQTL]{Matrix_eQTL_main}}.
#'  If a file with this name exists it is overwritten, if \code{NULL} output is not saved
#'  to file.
#' @param bootstrap_file_directory \code{character} or \code{NULL}. If not \code{NULL},
#'  files will be saved in this directory storing local \emph{cis} associations for
#'  the bootstrap eGene detection group (detection_bootstrapnumber.txt) and local
#'  \emph{cis} associations the bootstrap left-out eGene effect size estimation
#'  group (estimation_bootstrapnumber.txt). Estimation group files will only be saved
#'  where signficant eGenes are also significant in the bootstrap detection group
#'  (see Details). Corresponds to the \code{output_file_name.cis} argument in the
#'  respective calls to \code{\link[MatrixEQTL]{Matrix_eQTL_main}}. Files in this
#'  directory will be overwritten if they already exist.
#' @param cisDist \code{numeric}. Argument to \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
#'  controlling maximum distance SNP-gene pairs are considered local. The distance is
#'  measured to the nearest end of the gene.
#' @param local_correction multiple testing correction method to use when
#'  correcting p-values across all SNPs at each eGene. Must be a method
#'  specified in \code{\link[stats]{p.adjust.methods}} or "qvalue" for
#'  the \code{\link[qvalue]{qvalue}} package.
#' @param global_correction multiple testing correction method to use when
#'  correcting p-values across all eGenes. Must be a method specified in
#'  \code{\link[stats]{p.adjust.methods}} or "qvalue" for the
#'  \code{\link[qvalue]{qvalue}} package.
#'
#' @return
#'  A \code{data.frame} (or \code{\link[data.table]{data.table}} if the
#'  user has the library loaded) containing the results for each significant eGene:
#'  \describe{
#'    \item{\code{'gene':}}{The eGene.}
#'    \item{\code{'top_snp':}}{The SNP with the smallest p-value for the eGene in the MatrixEQTL analysis.
#'      Multiple SNPs separated by a "/" indicate multiple SNPs in perfect LD with identical
#'      test statistics, effect sizes, and p-values.}
#'    \item{\code{'statistic':}}{The test statistic for the \code{top_snp}-\code{gene} pair.}
#'    \item{\code{'nominal_beta':}}{The eQTL effect size for the \code{top_snp}-\code{gene} pair in the
#'      MatrixEQTL analysis.}
#'    \item{\code{'nominal_pval':}}{The p-value for the \code{top_snp}-\code{gene} pair from the MatrixEQTL
#'      analysis prior to multiple testing correction.}
#'    \item{\code{'corrected_pval':}}{The hierarchical multiple-testing adjusted p-value for the
#'    \code{top_snp}-\code{gene} pair}
#'    \item{\code{'winners_curse':}}{The amount of effect size overestimation determined by the
#'      bootstrap analysis (See Details).}
#'    \item{\code{'corrected_beta':}}{The eGene effect size after adjustment for the \code{winners_curse}.}
#'    \item{\code{'correction_boots':}}{The number of bootstraps that contributed to the estimation of
#'      the \code{winners_curse}, \emph{i.e.} the number of bootstraps in which the
#'      \code{top_snp}-\code{gene} pair was significant (see details)}
#'    \item{\code{'best_boot_eSNP':}}{The SNP or SNPs that were most frequently the top eSNP in
#'      the bootstrap procedure (see details).}
#'    \item{\code{'prop_top_eSNP':}}{The proportion of bootstraps in which the \code{'best_boot_eSNP'}
#'      was the top eSNP (see details).}
#'  }
#'
#' @import foreach
#' @import data.table
#' @import MatrixEQTL
#' @importFrom stats p.adjust p.adjust.methods
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
#' eGenes <- BootstrapEQTL(snps, gene, snpspos, genepos, cvrt, n_bootstraps=10, n_cores=2)
#'
BootstrapEQTL <- function(
  snps, gene, snpspos, genepos, cvrt=SlicedData$new(),
  n_bootstraps=500, n_cores=1, eGene_detection_file_name=NULL,
  bootstrap_file_directory=NULL, cisDist=1e6, local_correction="bonferroni",
  global_correction="fdr"
) {

  # R CMD check complains about data.table columns and foreach iterators
  # being undefined global variables or functions. The following
  # statements are functionally useless other than to suppress the R CMD
  # check warnings
  N <- NULL
  bootstrap <- NULL
  corrected_beta <- NULL
  corrected_pval <- NULL
  detection_beta <- NULL
  error <- NULL
  error_type <- NULL
  estimation_beta <- NULL
  id_boot <- NULL
  nominal_beta <- NULL
  # nominal_pval <- NULL
  prop_top_eSNP <- NULL
  pvalue <- NULL
  snp_type <- NULL
  statistic <- NULL
  top_snp <- NULL
  winners_curse <- NULL
  snp_block <- NULL
  ld_prev <- NULL
  best_boot_eSNP <- NULL

  # stringsAsFactors=TRUE causes crashes here
  saf <- options("stringsAsFactors")[[1]]
  options(stringsAsFactors=FALSE)
  on.exit(options(stringsAsFactors=saf))

  # Check column names are in same order
  if (!(all(snps$columnNames == gene$columnNames))) {
    stop("'snps' and 'genes' column names must be in same order")
  }
  if (cvrt$nCols() != 0 && !(all(cvrt$columnNames == gene$columnNames))) {
    stop("'cvrt' column names must be in same order as 'snps' and 'genes'")
  }

  # Check for special characters in SNP ids
  if (any(grepl("/", rownames(snps)) || any(grepl(";", rownames(snps))))) {
    stop('special characters ";" and "/" not allowed in SNP identifiers')
  }
  if (any(grepl("/", snpspos[,1]) || any(grepl(";", snpspos[,1])))) {
    # I.e. the user will rename the problem variants in 'snps' but
    # forget to also fix 'snpspos'.
    stop('special characters ";" and "/" found in \'snpspos\' SNP identifiers')
  }

  # Check multiple testing adjustment methods are ok
  mult.test.methods <- c(p.adjust.methods, "qvalue")
  if (!(local_correction %in% mult.test.methods ) ||
      !(global_correction %in% mult.test.methods )) {
    stop("'local_correction' and 'global_correction' must be one of ",
         paste(paste0('"', p.adjust.methods, '"'), collapse=", "), ", or ",
         '"qvalue".')
  }
  if ((local_correction == "qvalue" || global_correction == "qvalue") &&
      !pkgReqCheck("qvalue")) {
    stop("'qvalue' package not installed")
  }

  # Force cis-eQTL analysis
  if (missing(snpspos) || is.na(snpspos) || is.null(snpspos)) {
    stop("'snpspos' must be provided")
  }
  if (missing(genepos) || is.na(genepos) || is.null(genepos)) {
    stop("'genepos' must be provided")
  }

  # Create bootstrap file directory
  if (!is.null(bootstrap_file_directory)) {
    dir.create(bootstrap_file_directory, showWarnings=FALSE)
  }

  # Check if the user has already loaded data.table: if not, load it and
  # make sure we return the table as a data.frame
  has.data.table <- isNamespaceLoaded("data.table")
  if (!has.data.table) {
    suppressMessages(requireNamespace("data.table")) # silently load without tutorial message
  }

  cat("Testing cis-eQTL associations...\n")
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
    cisDist = cisDist,
    output_file_name=NULL,
    output_file_name.cis=eGene_detection_file_name,
    noFDRsaveMemory=ifelse(is.null(eGene_detection_file_name), FALSE, TRUE),
    useModel=modelLINEAR,
    verbose=FALSE
  );

  cat("Identifying significant eGenes..\n")
  cis_assocs <- get_cis_assocs(eQTLs, eGene_detection_file_name)
  eGenes <- get_eGenes(cis_assocs, local_correction, global_correction)

  cat("Identifying SNPs in perfect LD with lead SNPs...\n")
  eGenes <- get_eSNPs(cis_assocs, eGenes, collapse=TRUE)

  # Relabel columns
  eGenes <- eGenes[,list(gene, top_snp=snps, statistic, nominal_beta=beta,
                         nominal_pval=pvalue, corrected_pval)]

  if (eGenes[corrected_pval < 0.05, .N] == 0) {
    warning("No significant eGenes detected")
    return(eGenes)
  }

  # Make sure only necessary objects are ported to each worker
  boot_objs <-  c("cvrt", "snps", "gene", "cis_assocs", "eGenes", "snpspos",
                  "genepos", "bootstrap_file_directory")
  other_objs <- ls()[!(ls() %in% boot_objs)]

  # Set up parallel computing environment
  par_setup <- setupParallel(n_cores, verbose=TRUE, reporterCore=FALSE)
  on.exit({
    cleanupCluster(par_setup$cluster, par_setup$predef)
  }, add = TRUE)

  cat("Running bootstrap procedure for", n_bootstraps, "bootstraps.\n")
  # Run MatrixEQTL in each bootstrap detection group and estimation
  # group
  boot_eGenes <- foreach(id_boot = seq_len(n_bootstraps),
                         .inorder = FALSE,
                         .export = boot_objs, .noexport=other_objs,
                         .combine=rbind
  ) %dopar% {
    # Silently load packages on parallel workers
    suppressMessages(requireNamespace("data.table"))
    suppressMessages(requireNamespace("MatrixEQTL"))

    tryCatch({
      n_samples <- ncol(snps)
      id_detection <- sample(seq_len(n_samples), n_samples, replace = TRUE)
      id_estimation <- setdiff(seq_len(n_samples), id_detection)

      # Copy and subset the gene and snp data for the bootstrap
      # detection group
      gene_detection <- gene$Clone()
      gene_detection$ColumnSubsample(id_detection)
      snps_detection <- snps$Clone()
      snps_detection$ColumnSubsample(id_detection)
      cvrt_detection <- cvrt$Clone()
      cvrt_detection$ColumnSubsample(id_detection)

      # File to save detection group cis associations
      if (is.null(bootstrap_file_directory)) {
        detection_file <- NULL
      } else {
        detection_file <- file.path(bootstrap_file_directory, paste0("detection_", id_boot, ".txt"))
      }

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
        cisDist = cisDist,
        output_file_name=NULL,
        output_file_name.cis=detection_file,
        noFDRsaveMemory=ifelse(is.null(detection_file), FALSE, TRUE),
        verbose=FALSE
      )
      detection_cis_assocs <- get_cis_assocs(eQTL_detection, detection_file)

      # For the Winner's Curse correction we will examine the SNP-gene
      # pairs for significant eGenes from the eGene detection analysis.
      # If there were multiple SNPs in perfect LD, take the first --
      # we must have only 1 pvalue for the subsequent FDR correction.
      detection_eQTL_SNPs <- get_eGenes(detection_cis_assocs, local_correction, global_correction,
                                        eSNPs=eGenes[,list(gene, snps=gsub("/.*", "", top_snp))])

      # Filter to significant eGenes
      detection_eQTL_SNPs <- detection_eQTL_SNPs[
        corrected_pval < 0.05 & # Is the gene significant in this bootstrap?
        gene %in% eGenes[corrected_pval < 0.05, gene] # And in the original eGene detection analysis?
      ]

      # We also want to collect statistics about the top eSNP in the
      # detection group
      detection_top_SNPs <- get_eGenes(detection_cis_assocs, local_correction, global_correction)
      # Filter to significant eGenes
      detection_top_SNPs <- detection_top_SNPs[
        corrected_pval < 0.05 & # Is the gene significant in this bootstrap?
        gene %in% eGenes[corrected_pval < 0.05, gene] # And in the original eGene detection analysis?
      ]
      # Get eSNPs with identical statisticis due to perfect LD in
      # bootstrap detection group
      detection_top_SNPs <- get_eSNPs(detection_cis_assocs, detection_top_SNPs, collapse=FALSE)

      # These are the columns we want to keep out the other end
      detection_eQTL_SNPs <- detection_eQTL_SNPs[, list(gene, snps, detection_beta=beta)]
      detection_top_SNPs <- detection_top_SNPs[,list(gene, snps, detection_beta=NA_real_,
                                                      estimation_beta=NA_real_, snp_type="top")]

      # If there are no significant eGenes in this bootstrap we can move
      # to the next one
      if(nrow(detection_eQTL_SNPs) == 0) {
        if (nrow(detection_top_SNPs) == 0) {
          return(NULL)
        } else {
          sig_boot_assocs <- detection_top_SNPs
          sig_boot_assocs[, bootstrap := id_boot]
          sig_boot_assocs[, error := NA_character_]
          sig_boot_assocs[, error_type := NA_character_]
          return(sig_boot_assocs)
        }
      }

      tryCatch({
        # Run MatrixEQTL in the estimation group. We only care about the
        # significant SNP-Gene pairs from the detection group, so we can
        # save computation time by filtering their respective matrices.
        # Note: where multiple significant genes or snps are located within
        # the same cis window the filtering will not be perfect.
        sig_boot_genes <- detection_eQTL_SNPs[, gene]
        sig_boot_snps <- detection_eQTL_SNPs[, snps]

        gene_estimation <- gene$Clone()
        gene_estimation$RowReorder(which(gene_estimation$GetAllRowNames() %in% sig_boot_genes))
        gene_estimation$ColumnSubsample(id_estimation)
        snps_estimation <- snps$Clone()
        snps_estimation$RowReorder(which(snps_estimation$GetAllRowNames() %in% sig_boot_snps))
        snps_estimation$ColumnSubsample(id_estimation)
        cvrt_estimation <- cvrt$Clone()
        cvrt_estimation$ColumnSubsample(id_estimation)

        # File to save estimation group associations
        if (is.null(bootstrap_file_directory)) {
          estimation_file <- NULL
        } else {
          estimation_file <- file.path(bootstrap_file_directory, paste0("detection_", id_boot, ".txt"))
        }

        # Run MatrixEQTl in the estimation group
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
          cisDist = cisDist,
          output_file_name=NULL,
          output_file_name.cis=estimation_file,
          noFDRsaveMemory=ifelse(is.null(estimation_file), FALSE, TRUE),
          verbose=FALSE
        )
        estimation_cis_assocs <- get_cis_assocs(eQTL_estimation, estimation_file)

        # Create final table to return
        sig_boot_assocs <- merge(detection_eQTL_SNPs,
                                 estimation_cis_assocs[, list(gene, snps, estimation_beta=beta)],
                                 by=c("gene", "snps"))

        # Add the top detection SNPs table
        sig_boot_assocs[, snp_type := "eQTL"]
        sig_boot_assocs <- rbind(sig_boot_assocs, detection_top_SNPs)

        # Other housekeeping info
        sig_boot_assocs[,bootstrap := id_boot]
        sig_boot_assocs[,error := NA_character_]
        sig_boot_assocs[, error_type := NA_character_]

        return(sig_boot_assocs)
      }, warning=function(w) {
        warning(w$message) # pass through warnings -- shouldnt happen unless there's bugs
      }, error=function(e) {
        # If the estimation group eQTL analysis fails, return just the
        # top SNPs
        sig_boot_assocs <- detection_top_SNPs
        sig_boot_assocs[, bootstrap := id_boot]
        sig_boot_assocs[, error := NA_character_]
        sig_boot_assocs[, error_type := NA_character_]

        sig_boot_assocs <- rbind(sig_boot_assocs,
           data.table(gene=NA, snps=NA, snp_type=NA, detection_beta=NA,
                      estimation_beta=NA, bootstrap=id_boot, error=gsub("\n", " ", e$message),
                      error_type="estimation"))
        return(sig_boot_assocs)
      })
    }, warning=function(w) {
      warning(w$message) # pass through warnings -- shouldnt happen unless there's bugs
    }, error=function(e) {
      # If the detection group eQTL analysis fails return just the errors
      return(data.table(gene=NA, snps=NA, snp_type=NA, detection_beta=NA,
                        estimation_beta=NA, bootstrap=id_boot, error=gsub("\n", " ", e$message),
                        error_type="detection"))
    })
  }

  # Report failed bootstraps if any
  if (boot_eGenes[,sum(!is.na(error))] > 0) {
    errors <- boot_eGenes[!is.na(error), .N, by=list(error, error_type)]
    for (ii in errors[,.I]) {
      warning(errors[ii, N], "/", n_bootstraps, " bootstraps failed with",
              " error '", errors[ii, error], "' in the bootstrap ",
              errors[ii, error_type], " group.")
    }
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
  eGenes[corrected_pval < 0.05 & is.na(correction_boots), correction_boots := 0] # In case no bootstraps are significant

  # Calculate frequency of top SNPs across all significant boostraps to
  # report the best bootstrap eSNP (probable causal eSNP). In cases
  # where multiple SNPs occur equally freuqently we need to determine if
  # they are in perfect LD or not.
  top_SNP <- boot_eGenes[snp_type == "top", list(prop_top_eSNP=.N/length(unique(bootstrap))), by=list(gene, snps)]
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

  # Sort eGenes by significance
  eGenes <- eGenes[order(abs(statistic), decreasing=TRUE)][order(corrected_pval)]

  # If the user has not loaded data.table themselves cast back to a
  # data frame to avoid confusion
  if (!has.data.table) {
    eGenes <- as.data.frame(eGenes)
  }
  return(eGenes)
}
