#' Bootstrap QTL analysis for accurate effect size estimation
#'
#' Performs cis-QTL mapping using MatrixEQTL then performs a bootstrap
#' analysis to obtain unbiased effect size estimates for traits with
#' significant evidence of genetic regulation correcting for the
#' "Winner's Curse" arising from lead-SNP selection.
#'
#' @details
#'  Although the package interface and documentation describe the use of
#'  \code{BootstrapQTL} for \emph{cis}-eQTL mapping, the package can be
#'  applied to any QTL study of quantitative traits with chromosomal
#'  positions, for example \emph{cis}-QTL mapping of epigenetic
#'  modifications. Any matrix of molecular trait data can be provided
#'  to the \code{'gene'} argument provided a corresponding \code{'genepos'}
#'  'data.frame' detailing the chromosomal positions of each trait is
#'  provided.
#'
#'  \subsection{Cis-eQTL mapping:}{
#'  EQTL mapping is performed using the
#'  \code{\link[MatrixEQTL:MatrixEQTL-package]{MatrixEQTL}} package. A three step
#'  hieararchical multiple testing correction procedure is used to
#'  determine significant eGenes and eSNPs. At the first step, nominal
#'  p-values from \code{\link[MatrixEQTL:MatrixEQTL-package]{MatrixEQTL}} for all
#'  \emph{cis}-SNPs are adjusted for each gene separately using the
#'  method specified in the \code{'local_correction'} argument
#'  (Bonferroni correction by default). In the second step, the best
#'  adjusted p-value is taken for each gene, and this set of locally
#'  adjusted p-values is corrected for multiple testing across all genes
#'  using the methods pecified in the \code{'global_correction'} argument
#'  (FDR correction by default). In the third step, an eSNP significance
#'  threshold on the locally corrected p-values is determined as the
#'  locally corrected p-value corresponding to the globally corrected
#'  p-value threshold of 0.05.
#'
#'  A gene is considered a significant eGene if its globally corrected
#'  p-value is < 0.05, and a SNP is considered a significant eSNP for
#'  that eGene if its locally corrected p-value < the eSNP significance
#'  threshold.
#'
#'  The default settings for \code{'local_correction'} and
#'  \code{'global_correction'} were found to best control eGene false
#'  discovery rate without sacrificing sensitivity (see citation).
#'  }
#'  \subsection{Winner's Curse correction:}{
#'  EQTL effect sizes of significant eSNPs on significant eGenes are
#'  typically overestimated when compared to replication datasets
#'  (see citation). \code{BootstrapEQTL} removes this overestimation by
#'  performing a bootstrap procedure after eQTL mapping.
#'
#'  Three Winner's Curse correction methods are available: the Shrinkage
#'  method, the Out of Sample method, and the Weighted Estimator method.
#'  All three methods work on the same basic principle of performing
#'  repeated sample bootstrapping to partition the dataset into two
#'  groups: an eQTL detection group comprising study samples select via
#'  random sampling with replacement, and an eQTL effect size estimation
#'  group comprising the remaining samples not selected via the random
#'  sampling. The default estimator, \code{'correction_type = "shrinkage"'},
#'  provided the most accurate corrected effect sizes in our simulation
#'  study (see citation).
#'
#'  The \strong{shrinkage method} ("shrinkage" in
#'  \code{'correction_type'}) corrects for the Winner's Curse by
#'  measuring the average difference between the eQTL effect size
#'  in the bootstrap detection group and the bootstrap estimation group,
#'  then subtracting this difference from the naive eQTL effect size
#'  estimate obtained from the eGene detection analysis prior to the
#'  bootstrap procedure.
#'
#'  The \strong{out of sample method} ("out_of_sample" in
#'  \code{'correction_type'}) corrects for the Winner's Curse by taking
#'  the average eQTL effect size across bootstrap estimation groups as
#'  an unbiased effect size estimate.
#'
#'  The \strong{weighted estimator method} ("weighted" in
#'  \code{'correction_type'}) corrects for the Winner's Curse by taking
#'  a weighted average of the nominal estimate of the eQTL effect size
#'  and the average of eQTL effect sizes across the bootstrap estimation
#'  groups: \eqn{0.368 * naive_estimate + 0.632 *
#'  mean(bootstrap estimation group effect sizes)}.
#'
#'  In all three methods bootstrap effect sizes only contribute to
#'  the Winner's Curse correction if the corresponding eSNP is
#'  significantly associated with the eGene in the bootstrap detection
#'  group (locally corrected bootstrap P-value < eSNP significance
#'  threshold determing in the eQTL mapping step).
#'
#'  Note that eQTLs may not remain significant in all bootstraps, so the
#'  effective number of bootstraps used to obtain the Winner's Curse
#'  estimate will typically be lower than the number of bootstraps
#'  specified in \code{'n_bootstraps'}. The number of bootstraps that
#'  were significant for each eQTL are reported in the
#'  \code{'correction_boots'} column of the returned table.
#'  }
#'  \subsection{Winner's Curse corrected effect sizes}{
#'  The user should be aware that ability to correct for Winner's Curse
#'  can vary across significant eQTLs depending on their statistical
#'  power (\emph{i.e. minor allele frequency, true effect size, and
#'  study sample size}). Users should be skeptical of corrected effect
#'  sizes that are larger than the nominal effect sizes estimated by
#'  \code{\link[MatrixEQTL:MatrixEQTL-package]{MatrixEQTL}}, which likely reflects low
#'  power for eQTL detection rather than an underestimated effect size.
#'  }
#'  \subsection{Bootstrap warning messages:}{
#'  It is possible for bootstrap analyses to fail due to the reduced
#'  sample sizes of the bootstrap detection and bootstrap estimation
#'  groups. For example, the bootstrap resampling may lead to an
#'  detection or estimation groups in which all individuals are
#'  homozygous for an eSNP or have no variance in their supplied
#'  covariates (\emph{e.g.} the estimation group may comprise
#'  individuals all of the same sex). In this case the bootstrap will
#'  fail for all eQTLs since \code{\link[MatrixEQTL:MatrixEQTL-package]{MatrixEQTL}} will
#'  be unable to perform the model fitting.
#'
#'  Failed bootstraps are reported after the bootstrap procedure in
#'  a series of warning messages indicating the number of bootstrap
#'  failures grouped by the reason for the bootstrap failure.
#'  }
#'
#' @param snps \code{\link[MatrixEQTL:SlicedData-class]{SlicedData}} object containing genotype
#'   information used as input into \code{\link[MatrixEQTL]{Matrix_eQTL_main}}.
#' @param gene \code{\link[MatrixEQTL:SlicedData-class]{SlicedData}} object containing gene expression
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
#'   Note that \code{\link[MatrixEQTL]{Matrix_eQTL_main}} tests all
#'   variants within \code{cisDist} of the start or end of the gene.
#'   If you wish instead to test all variants within \code{cisDist} of
#'   the transcription start site, you should specify this location in
#'   both the 'left' and 'right' columns of the \code{genepos}
#'   data.frame. Similarly, when analysing a molecular phenotype that
#'   have a single chromosomal position then the 'left' and 'right'
#'   columns should both contain the same position.
#' @param cvrt \code{\link[MatrixEQTL:SlicedData-class]{SlicedData}} object containing covariate
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
#'  controlling the maximum distance from a gene to consider tests for
#'  eQTL mapping.
#' @param local_correction multiple testing correction method to use when
#'  correcting p-values across all SNPs at each gene (see EQTL mapping
#'  section in Details). Can be a method specified in \code{\link[stats:p.adjust]{p.adjust.methods}},
#'  "qvalue" for the \code{\link[qvalue]{qvalue}} package, or "eigenMT"
#'  if EigenMT has been used to estimate the number effective independent 
#'  tests (see \code{eigenMT_tests_per_gene}).
#' @param global_correction multiple testing correction method to use when
#'  correcting p-values across all genes after performing local correction
#'  (see EQTL mapping section in Details). Must be a method specified in
#'  \code{\link[stats:p.adjust]{p.adjust.methods}} or "qvalue" for the
#'  \code{\link[qvalue]{qvalue}} package.
#' @param correction_type \code{character}. One of "shrinkage", "out_of_sample"
#'  or "weighted". Determines which Winner's Curse correction method is
#'  used (see Details).
#' @param errorCovariance \code{numeric matrix} argument to \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
#'  specifying the error covariance.
#' @param useModel \code{integer} argument to \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
#'  specifying the type of model to fit between each SNP and gene. Should be one of
#'  \code{\link[MatrixEQTL]{modelLINEAR}}, \code{\link[MatrixEQTL]{modelANOVA}}, or
#'  \code{\link[MatrixEQTL]{modelLINEAR_CROSS}}.
#' @param eigenMT_tests_per_gene \code{data.frame} containing the number of effective
#'  independent tests for each gene estimated by the EigenMT (\url{https://github.com/joed3/eigenMT}).
#'  Ignore unless \code{'local_correction="eigenMT"'}.
#'
#' @return
#'  A \code{data.frame} (or \code{\link[data.table]{data.table}} if the
#'  user has the library loaded) containing the results for each significant eQTL:
#'  \describe{
#'    \item{\code{'eGene':}}{The eQTL eGene.}
#'    \item{\code{'eSNP':}}{The eQTL eSNP.}
#'    \item{\code{'statistic':}}{The test statistic for the association between the eGene and eSNP.}
#'    \item{\code{'nominal_beta':}}{The eQTL effect size for the \code{eGene}-\code{eSNP}
#'      pair estimated by \code{\link[MatrixEQTL:MatrixEQTL-package]{MatrixEQTL}}.}
#'    \item{\code{'corrected_beta':}}{The eQTL effect size after adjustment for the \code{winners_curse}.}
#'    \item{\code{'winners_curse':}}{The amount of effect size overestimation determined by the
#'      bootstrap analysis (See Details).}
#'    \item{\code{'correction_boots':}}{The number of bootstraps that contributed to the estimation of
#'      the \code{winners_curse}, \emph{i.e.} the number of bootstraps in which the \code{eSNP}
#'      remained significantly associated with the \code{eGene} (see Details).}
#'    \item{\code{'nominal_pval':}}{The p-value for the \code{eGene}-\code{eSNP} pair
#'      from the \code{\link[MatrixEQTL:MatrixEQTL-package]{MatrixEQTL}} analysis.}
#'    \item{\code{'eSNP_pval':}}{The locally corrected p-value for the \code{eGene}-\code{eSNP} pair (see Details).}
#'    \item{\code{'eGene_pval':}}{The globally corrected p-value for the \code{eGene} based on its top eSNP (see Details).}
#'  }
#'
#' @import foreach
#' @import data.table
#' @import MatrixEQTL
#' @importFrom stats p.adjust p.adjust.methods
#' @importFrom utils sessionInfo
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
#' # Run the BootstrapQTL analysis
#' eQTLs <- BootstrapQTL(snps, gene, snpspos, genepos, cvrt, n_bootstraps=10, n_cores=2)
#'
BootstrapQTL <- function(
  snps, gene, snpspos, genepos, cvrt=SlicedData$new(),
  n_bootstraps=200, n_cores=1, eGene_detection_file_name=NULL,
  bootstrap_file_directory=NULL, cisDist=1e6,
  local_correction="bonferroni", global_correction="fdr",
  correction_type="shrinkage", errorCovariance=numeric(),
  useModel=modelLINEAR, eigenMT_tests_per_gene=NULL
) {

  # R CMD check complains about data.table columns and foreach iterators
  # being undefined global variables or functions. The following
  # statements are functionally useless other than to suppress the R CMD
  # check warnings
  global_pval <- NULL
  local_pval <- NULL
  id_boot <- NULL
  bootstrap <- NULL
  id_boot <- NULL
  error <- NULL
  error_type <- NULL
  N <- NULL
  statistic <- NULL
  corrected_beta <- NULL
  winners_curse <- NULL
  correction_boots <- NULL
  pvalue <- NULL
  eSNP_pval <- NULL
  eGene_pval <- NULL

  ##--------------------------------------------------------------------
  ## Check inputs will be OK first - want function to crash at the start
  ## instead of after several hours running MatrixEQTL if inputs
  ## are bad.
  ##--------------------------------------------------------------------

  # Check class of 'snps' and 'gene'
  if (!("SlicedData" %in% class(snps)) || !("SlicedData" %in% class(gene)) ||
      !("SlicedData" %in% class(cvrt))) {
    stop("'snps', 'gene', and 'cvrt' must be objects of type \"SlicedData\"")
  }
  # Check column names are in same order
  if (!(all(snps$columnNames == gene$columnNames))) {
    stop("'snps' and 'gene' column names must be in same order")
  }
  if (cvrt$nCols() != 0 && !(all(cvrt$columnNames == gene$columnNames))) {
    stop("'cvrt' column names must be in same order as 'snps' and 'gene'")
  }

  # Check multiple testing adjustment methods are ok
  mult.test.methods <- c(p.adjust.methods, "qvalue")
  if (length(local_correction) > 1 || !(local_correction %in% c(mult.test.methods, "eigenMT"))) {
    stop("'local_correction' must be one of ", paste(paste0('"',  mult.test.methods, '"'), ", or \"eigenMT\""))
  }
  if (length(global_correction) > 1 || !(global_correction %in% mult.test.methods)) {
    stop("'local_correction' must be one of ", paste(paste0('"', p.adjust.methods, '"'), ", or \"qvalue\""))
  }
  if ((local_correction == "qvalue" || global_correction == "qvalue") &&
      !pkgReqCheck("qvalue")) {
    stop("'qvalue' package not installed")
  }
  if (local_correction == "eigenMT" && (
    is.null(eigenMT_tests_per_gene) ||  
    !is.data.frame(eigenMT_tests_per_gene) ||
    nrow(eigenMT_tests_per_gene) != gene$nRows() ||
    ncol(eigenMT_tests_per_gene) != 2 || 
    !is.numeric(eigenMT_tests_per_gene[,2]) || 
    any(is.na(eigenMT_tests_per_gene[,2])) ||
    !all(eigenMT_tests_per_gene[,1] %in% rownames(gene)) ||
    !all(rownames(gene) %in% eigenMT_tests_per_gene[,1])
  )) {
    stop("'eigenMT_tests_per_gene' must be a data.frame containing the number of effective independent tests per gene")
  }
  if (!is.null(eigenMT_tests_per_gene)) {
    eigenMT_tests_per_gene <- as.data.table(eigenMT_tests_per_gene)
    setnames(eigenMT_tests_per_gene, c("gene", "n_tests"))
  }

  # Check correction_type input is ok
  if (length(correction_type) > 1 || !(correction_type %in% c("shrinkage", "out_of_sample", "weighted"))) {
    stop('\'correction_type\' must be one of "shrinkage", "out_of_sample", or "weighted"')
  }

  # Force cis-eQTL analysis
  if (missing(snpspos) || is.na(snpspos) || is.null(snpspos)) {
    stop("'snpspos' must be provided")
  }
  if (missing(genepos) || is.na(genepos) || is.null(genepos)) {
    stop("'genepos' must be provided")
  }

  # Check snpspos and genepos inputs
  if (!is.data.frame(snpspos) && ncol(snpspos) != 3) {
    stop("'snpspos' must be a data.frame with 3 columns")
  }
  if (!is.data.frame(genepos) && ncol(genepos) != 4) {
    stop("'genepos' must be a data.frame with 4 columns")
  }

  # Check files and directories are ok
  if (!is.null(eGene_detection_file_name) && length(eGene_detection_file_name) > 1) {
    stop("Only 1 file may be specified in 'eGene_detection_file_name'")
  }
  if (!is.null(bootstrap_file_directory) && length(bootstrap_file_directory) > 1) {
    stop("Only 1 directory may be specified in 'bootstrap_file_directory'")
  }
  if (!is.null(bootstrap_file_directory)) {
    tryCatch({
      dir.create(bootstrap_file_directory, showWarnings=FALSE)
    }, error=function(e) {
      stop("Could not create or access directory specified in 'bootstrap_file_directory'")
    })
  }

  # Check number of bootstraps
  if (n_bootstraps < 0) {
    stop("'n_bootstraps' must be >= 0")
  }

  # Check cisDist
  if (length(cisDist) != 1 || !is.numeric(cisDist) || !is.finite(cisDist) || cisDist < 1) {
    stop("'cisDist' must be a number > 0")
  }

  ##--------------------------------------------------------------------
  ## Set up and check options that need to be reset at end of the
  ## function
  ##--------------------------------------------------------------------
  on.exit({ cat("Restoring global environment...\n") })
  # stringsAsFactors=TRUE causes crashes here
  saf <- options("stringsAsFactors")[[1]]
  options(stringsAsFactors=FALSE)
  on.exit(options(stringsAsFactors=saf), add=TRUE)

  # Check if the user has already loaded data.table: if not, load it and
  # make sure we return the table as a data.frame
  has.data.table <- !("data.table" %in% names(sessionInfo()$loadedOnly))
  if (!has.data.table) {
    suppressMessages(requireNamespace("data.table")) # silently load without tutorial message
  }

  # Set up parallel computing environment
  par_setup <- setupParallel(n_cores, verbose=TRUE, reporterCore=FALSE)
  on.exit({
    cleanupCluster(par_setup$cluster, par_setup$predef)
  }, add = TRUE)

  ##--------------------------------------------------------------------
  ## Perform QTL mapping and nominal effect size estimation
  ##--------------------------------------------------------------------

  cat("Mapping QTLs and obtaining nominal effect size estimates...\n")
  # Run Matrix eQTL to determine significant eGenes and get nominal
  # estimates for their effect sizes
  eQTLs <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    pvOutputThreshold = 0, # we don't need the Trans eQTL pvalues
    errorCovariance = errorCovariance,
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    output_file_name=NULL,
    output_file_name.cis=eGene_detection_file_name,
#   noFDRsaveMemory=ifelse(is.null(eGene_detection_file_name), FALSE, TRUE),
    noFDRsaveMemory=FALSE, # See: https://github.com/andreyshabalin/MatrixEQTL/issues/8
    useModel=modelLINEAR,
    verbose=FALSE
  );
  # Load the table of all cis-Assocations
  cis_assocs <- get_cis_assocs(eQTLs, eGene_detection_file_name)
  # Perform hierarchcial correction
  cis_assocs <- hierarchical_correction(cis_assocs, local_correction, global_correction, eigenMT_tests_per_gene)
  # Determine significance threshold for eSNPs
  eSNP_threshold <- get_eSNP_threshold(cis_assocs)
  # Filter to significant associations
  sig_assocs <- cis_assocs[global_pval < 0.05 & local_pval < eSNP_threshold]

  if (nrow(sig_assocs) == 0) {
    stop("No significant associations detected")
  }

  if (n_bootstraps > 0) {
    ##--------------------------------------------------------------------
    ## Apply filters prior to bootstrapping to speed up calculations and
    ## minimise memory usage
    ##--------------------------------------------------------------------
  
    cat("Optimising data for bootstrap procedure...\n")
  
    # Filter 'genes' and 'snps' to just the significant assocations
    sig_pairs <- sig_assocs[,list(gene, snps)]
  
    gene_boot <- gene$Clone()
    gene_boot$RowReorder(which(rownames(gene_boot) %in% sig_assocs[,gene]))
    genepos <- genepos[genepos[,1] %in% sig_assocs[,gene], ]
  
    # If we don't need to run the bootstrap analysis for all cis-SNPs, we
    # can filter to just the significant eSNPs
    if (local_correction %in% c("bonferroni", "none")) {
      tests_per_gene <- cis_assocs[gene %in% sig_assocs[,gene], list(n_tests=.N), by=gene]
    } else if (local_correction == "eigenMT") {
      tests_per_gene <- eigenMT_tests_per_gene
    } else {
      tests_per_gene <- NULL
    }
    
    if (!is.null(tests_per_gene)) {
      snps_boot <- snps$Clone()
      snps_boot$RowReorder(which(rownames(snps_boot) %in% sig_assocs[,snps]))
    } else {
      snps_boot <- snps
    }
  
    snpspos <- snpspos[snpspos[,1] %in% sig_assocs[,snps],]
  
    # Make sure only necessary objects are ported to each worker
    boot_objs <-  c("cvrt", "snps_boot", "gene_boot", "snpspos",
                    "genepos", "bootstrap_file_directory", "sig_pairs",
                    "tests_per_gene", "local_correction", "eSNP_threshold",
                    "errorCovariance", "cisDist", "useModel")
    other_objs <- ls()[!(ls() %in% boot_objs)]
  
    ##--------------------------------------------------------------------
    ## Run the bootstrap procedure
    ##--------------------------------------------------------------------
  
    cat("Running bootstrap procedure for", n_bootstraps, "bootstraps on", getDoParWorkers(), "cores.\n")
    # Run MatrixEQTL in each bootstrap detection group and estimation
    # group
    boot_eGenes <- foreach(id_boot = seq_len(n_bootstraps),
                           .inorder = FALSE,
                           .export = boot_objs, .noexport=other_objs,
                           .init = data.table(error=character(0), error_type=character(0)),
                           .combine = rbind_dt  # rbind with FILL = TRUE as default
    ) %dopar% {
      ##--------------------------------------------------------------------
      ## Detection group gene detection
      ##--------------------------------------------------------------------
      tryCatch({ # WHY IS EXCEPTION HANDLING SO HARD IN THIS TERRIBLE LANGUAGE?
        # Silently load packages on parallel workers
        suppressMessages(requireNamespace("data.table"))
        suppressMessages(requireNamespace("MatrixEQTL"))
  
        n_samples <- ncol(snps_boot)
        id_detection <- sample(seq_len(n_samples), n_samples, replace = TRUE)
        id_estimation <- setdiff(seq_len(n_samples), id_detection)
  
  
        # Copy and subset the gene and snp data for the bootstrap
        # detection group
        gene_detection <- gene_boot$Clone()
        gene_detection$ColumnSubsample(id_detection)
        snps_detection <- snps_boot$Clone()
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
          useModel = useModel,
          errorCovariance = errorCovariance,
          pvOutputThreshold.cis = 1,
          snpspos = snpspos,
          genepos = genepos,
          cisDist = cisDist,
          output_file_name=NULL,
          output_file_name.cis=detection_file,
          #noFDRsaveMemory=ifelse(is.null(detection_file), FALSE, TRUE),
          noFDRsaveMemory=FALSE, # See: https://github.com/andreyshabalin/MatrixEQTL/issues/8
          verbose=FALSE
        )
        # Load the table of all cis-Assocations
        detection_cis_assocs <- get_cis_assocs(eQTL_detection, detection_file)
        # Perform local SNP correction
        detection_cis_assocs <- hierarchical_correction(detection_cis_assocs, local_correction, "none", tests_per_gene)
        # Filter to just significant eSNP-eGene pairs (i.e. some SNPs may
        # be in cis with multiple genes, but only associated with 1)
        detection_cis_assocs <- detection_cis_assocs[sig_pairs, on=c("gene", "snps")]
        # Filter to associations that are significant in this bootstrap
        detection_sig_assocs <- detection_cis_assocs[local_pval < eSNP_threshold]
        # Filter columns
        detection_sig_assocs <- detection_sig_assocs[, list(gene, snps, detection_beta=beta)]
        # Add identifier columns
        detection_sig_assocs[, bootstrap := id_boot]
  
        # If there are no significant associations in this bootstrap we
        # can move to the next one
        if(nrow(detection_sig_assocs) == 0) {
          return(NULL)
        }
  
        ##----------------------------------------------------------------
        ## Estimation group effect size re-estimation
        ##----------------------------------------------------------------
        tryCatch({
          # Filter to significant genes and SNPs in the detection group
          gene_estimation <- gene_boot$Clone()
          gene_estimation$RowReorder(which(gene_estimation$GetAllRowNames() %in% detection_sig_assocs[,gene]))
          gene_estimation$ColumnSubsample(id_estimation)
          snps_estimation <- snps_boot$Clone()
          snps_estimation$RowReorder(which(snps_estimation$GetAllRowNames() %in% detection_sig_assocs[,snps]))
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
            useModel = useModel,
            errorCovariance = errorCovariance,
            pvOutputThreshold.cis = 1,
            snpspos = snpspos,
            genepos = genepos,
            cisDist = cisDist,
            output_file_name=NULL,
            output_file_name.cis=estimation_file,
            #noFDRsaveMemory=ifelse(is.null(detection_file), FALSE, TRUE),
            noFDRsaveMemory=FALSE, # See: https://github.com/andreyshabalin/MatrixEQTL/issues/8
            verbose=FALSE
          )
          # Load the table of all cis-assocations
          estimation_cis_assocs <- get_cis_assocs(eQTL_estimation, estimation_file)
  
          # Create final table to return
          sig_boot_assocs <- merge(detection_sig_assocs,
                                   estimation_cis_assocs[, list(gene, snps, estimation_beta=beta)],
                                   by=c("gene", "snps"))
          return(sig_boot_assocs)
        }, error=function(e) {
          return(data.table(bootstrap=id_boot, error=e$message, error_type="estimation"))
        })
      }, error=function(e) {
        return(data.table(bootstrap=id_boot, error=e$message, error_type="detection"))
      })
    }
  
    ##--------------------------------------------------------------------
    ## Collate bootstrap errors
    ##--------------------------------------------------------------------
    # Report failed bootstraps if any
    if (boot_eGenes[,sum(!is.na(error))] > 0) {
      errors <- boot_eGenes[!is.na(error), .N, by=list(error, error_type)]
      for (ii in errors[,.I]) {
        warning(errors[ii, N], "/", n_bootstraps, " bootstraps failed with",
                " error '", errors[ii, error], "' in the bootstrap ",
                errors[ii, error_type], " group.", immediate.=TRUE)
      }
    }
    boot_eGenes <- boot_eGenes[is.na(error)]
  
    ##--------------------------------------------------------------------
    ## Remove effect of the winners curse based on bootstrap estimates
    ##--------------------------------------------------------------------
  
    cat("Removing winner's curse...\n")
    # Calculate the effect of the winner's curse. If effect sizes have been
    # estimated using the top bootstrap eSNP then we need to account for the
    # fact that different eSNPs with anti-correlated minor allele frequencies
    # may have been the top SNP at each bootstrap
    sig_assocs  <- correct_winners_curse(boot_eGenes, sig_assocs, correction_type)
    
    ##--------------------------------------------------------------------
    ##  Prepare table to be returned
    ##--------------------------------------------------------------------
    
    # Relabel columns
    sig_assocs <- sig_assocs[,list(eGene=gene, eSNPs=snps, statistic, nominal_beta=beta,
                                   corrected_beta, winners_curse, correction_boots,
                                   nominal_pval=pvalue, eSNP_pval=local_pval,
                                   eGene_pval=global_pval)]
  } else {
    sig_assocs <- sig_assocs[,list(eGene=gene, eSNPs=snps, statistic, nominal_beta=beta,
                                   nominal_pval=pvalue, eSNP_pval=local_pval,
                                   eGene_pval=global_pval)]
  }

  # Sort table by significance
  sig_assocs <- sig_assocs[order(abs(statistic), decreasing=TRUE)][order(eSNP_pval)][order(eGene_pval)]

  # If the user has not loaded data.table themselves cast back to a
  # data frame to avoid confusion
  if (!has.data.table) {
    sig_assocs  <- as.data.frame(sig_assocs)
  }

  on.exit({ cat("Done!\n") }, add=TRUE)
  return(sig_assocs)
}
