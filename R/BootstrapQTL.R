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
#'  \subsection{EGene detection:}{
#'  Detection of significant eGenes is performed using a hieararchical
#'  multiple testing correction procedure. At each gene, local multiple
#'  testing adjustment is performed using the method specified by the
#'  \code{'local_correction'} argument (Bonferroni correction by default)
#'  to correct for the number of \emph{cis}-SNPs tested at each gene.
#'  Global multiple testing adjustment is subsequently performed across
#'  all genes by taking the minimum locally corrected p-value at each
#'  gene, then performing the \code{'global_correction'} adjustment across
#'  all genes (FDR correction by default). These default settings best
#'  control the eGene false discovery rate without sacrificing
#'  sensitivity (see citation).
#'  }
#'  \subsection{Winner's Curse correction:}{
#'  EQTL effect sizes of significant eGenes are typically overestimated
#'  when compared to replication datasets due to the selection of the
#'  best \emph{cis}-SNP at each gene. \code{BootstrapEQTL} removes this
#'  overestimation by performing a bootstrap procedure after significant
#'  eGene discovery by \code{\link{MatrixEQTL}}.
#'
#'  Three Winner's Curse correction methods are available: the Shrinkage
#'  method, the Out of Sample method, and the Weighted Estimator method.
#'  All three methods work on the same basic principle of performing
#'  repeated sample bootstrapping to partition the dataset into two
#'  groups: an eQTL detection group comprising study samples select via
#'  random sampling with replacement, and an eQTL effect size estimation
#'  group comprising the remaining samples not selected via the random
#'  sampling.
#'
#'  The \strong{shrinkage method} ("shrinkage" in
#'  \code{'correction_type'}) corrects for the winners curse by
#'  measuring the average difference between the eGene effect size in
#'  the bootstrap detection group and the bootstrap estimation group,
#'  then subtracting this difference from the naive eQTL effect size
#'  estimate obtained from the eGene detection analysis prior to the
#'  bootstrap procedure.
#'
#'  The \strong{out of sample method} ("out_of_sample" in
#'  \code{'correction_type'}) corrects for the winners curse by taking
#'  the average eGene effect size across bootstrap estimation groups as
#'  an unbiased effect size estimate.
#'
#'  The \strong{weighted estimator method} ("weighted" in
#'  \code{'correction_type'}) corrects for the winners curse by taking a
#'  weighted average of the naive estimate of the effect size and the
#'  average of eGene effect sizes across the bootstrap estimation
#'  groups: \eqn{0.368 * naive_estimate + 0.632 *
#'  mean(bootstrap_estimation_group_effect_sizes)}.
#'
#'  In all three methods bootstrap group effect sizes only contribute to
#'  the winner's curse correction if the corresponding eGene is
#'  significant in the bootstrap detection group. Two approaches are
#'  provided for measuring eGene significance and effect sizes in the
#'  bootstrap groups: (1) using the SNP in each bootstrap detection
#'  group with the smallest p-value ("top" in \code{'bootstrap_eSNPs'}),
#'  and (2) using the SNP with the smallest p-value in the eGene
#'  detection analysis performed prior to the bootstrap procedure
#'  ("discovery" in \code{'bootstrap_eSNPs'}). In both cases bootstraps
#'  are only considered in the final winner's curse calculation where
#'  the gene-SNP assocation is significant in the bootstrap detection
#'  group. Each gene-SNP pair is considered significant if its locally
#'  corrected bootstrap detection group P-value is smaller than the
#'  locally corrected p-value in the eGene disocvery analysis
#'  corresponding to the global multiple testing correction threshold of
#'  0.05 in the eGene discovery analysis.
#'
#'  The default settings, \code{'correction_type = "shrinkage"'} and
#'  \code{'bootstrap_eSNPs = "discovery"'} provided the most accurate
#'  corrected effect sizes in our simulation study (see citation).
#'
#'  Note that use of \code{'bootstrap_eSNPs = "top"'} will
#'  dramatically slow down the bootstrap procedure as associaitons will
#'  need to be calculated between all \emph{cis}-SNPs and significant
#'  eGenes rather than just between the significant eGenes and their
#'  lead SNPs.
#'
#'  Note that SNP-gene pairs may not remain significant in all
#'  bootstraps, so the effective number of bootstraps used obtain the
#'  "Winner's Curse" estimate will typically be lower than the number of
#'  bootstraps specified in \code{'n_bootstraps'}. The number of
#'  bootstraps that were significant for each eGene is reported in the
#'  \code{'correction_boots'} column of the returned table.
#'  }
#'  \subsection{Collating the most frequent top bootstrap SNP}{
#'  In addition to correcting for the Winner's Curse, the bootstrap
#'  procedure can also collects statistics about which SNP(s) are the top
#'  SNP for each tested eGene in the bootstrapped detection group. This
#'  can differ from the top SNP from the initial eGene detection
#'  analysis in any given bootstrap. The \code{'best_boot_eSNP'} column
#'  reports the SNP that was the top SNP most frequently across the
#'  bootstrap procedure. There is some evidence from simulation studies
#'  that this is more likely to be the causal SNP than the
#'  \code{'top_snp'} from the eGene detection analysis (see citation),
#'  although the \code{'top_snp'} and \code{'best_boot_eSNP'} are
#'  typically the same.
#'
#'  Multiple SNPs may be reported in the \code{'best_boot_eSNP'} column.
#'  SNPs separated by a "/" indicate multiple SNPs in perfect linkage
#'  disequilbrium (LD), \emph{i.e.} with identical test
#'  statistics/effect sizes/p-values in all bootstraps, while SNPs
#'  separated by a ";" indicate multiple SNPs (or SNP groups in perfect
#'  LD) that were the top SNP(s) with equal frequency across the
#'  bootstrap procedure (\emph{e.g.} "SNP1;SNP2/SNP3" could indicate
#'  SNP1 was the top SNP in 50\% of significant bootstraps and SNP2 and
#'  SNP3 are in perfect LD and were was the top SNP in the other 50\% of
#'  significant bootstraps). The \code{'prop_top_eSNP'} column reports
#'  the proportion of significant bootstraps for which the
#'  \code{'best_boot_eSNP'} was the top SNP. A low
#'  \code{'prop_top_eSNP'} (\emph{e.g.} < 0.5) can indicate complex LD
#'  structure around the top eSNP or multiple causal eSNPs. Note that
#'  here the number of significant bootstraps may be higher than that
#'  reported in \code{'correction_boots'} if
#'  \code{'bootstrap_eSNPs="discovery"'}, as an eGene may be significant
#'  when performing global multiple testing correction across the top
#'  bootstrap detection group SNP but not when performing global
#'  multiple testing correction on the top SNP from the original eGene
#'  detection analysis.
#'  }
#'  \subsection{Bootstrap warning messages:}{
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
#'   If analysing a molecular phenotype that have a single chromosomal
#'   position then the 'left' and 'right' columns should both contain
#'   the same position.
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
#' @param bootstrap_eSNPs \code{character}. One of "discovery" or "top".
#'  Controls which SNPs are used for effect size estimation in the bootstrap
#'  procedure. If "discovery" effect sizes are estimated based on the top
#'  eSNPs from the eGene detection analysis performed prior to the bootstrap
#'  procedure. If "top", effect sizes are estimated based on the top eSNP
#'  in each bootstrap (note: this makes the bootstrap procedure very slow;
#'  see Details).
#' @param correction_type \code{character}. One of "shrinkage", "out_of_sample"
#'  or "weighted". Determines which Winner's Curse correction method is
#'  used (see Details).
#' @param collate_top_snps \code{logical}. Should the collation of the
#'  most frequent top SNP per eGene be skipped? Allows fast bootstrap
#'  by allowing each bootstrap to test only the lead SNP-gene pairs from
#'  the eGene detection analysis. Only works given the default package
#'  settings: \code{'bootstrap_eSNPs="discovery"'} and \code{'local_correction="bonferroni"'}
#'  or (\code{'local_correction="none"'}).
#' @param errorCovariance \code{numeric matrix} argument to \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
#'  specifying the error covariance.
#' @param useModel \code{integer} argument to \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
#'  specifying the type of model to fit between each SNP and gene. Should be one of
#'  \code{\link[MatrixEQTL]{modelLINEAR}}, \code{\link[MatrixEQTL]{modelANOVA}}, or
#'  \code{\link[MatrixEQTL]{modelLINEAR_CROSS}}.
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
#' # Run the BootstrapQTL analysis
#' eGenes <- BootstrapQTL(snps, gene, snpspos, genepos, cvrt, n_bootstraps=10, n_cores=2)
#'
BootstrapQTL <- function(
  snps, gene, snpspos, genepos, cvrt=SlicedData$new(),
  n_bootstraps=200, n_cores=1, eGene_detection_file_name=NULL,
  bootstrap_file_directory=NULL, cisDist=1e6,
  local_correction="bonferroni", global_correction="fdr",
  bootstrap_eSNPs="discovery", correction_type="shrinkage",
  collate_top_snps=ifelse(bootstrap_eSNPs=="discovery" && local_correction == "bonferroni", FALSE, TRUE),
  errorCovariance=numeric(), useModel=modelLINEAR
) {

  # R CMD check complains about data.table columns and foreach iterators
  # being undefined global variables or functions. The following
  # statements are functionally useless other than to suppress the R CMD
  # check warnings
  global_pval <- NULL
  local_bonf <- NULL
  pvalue <- NULL
  n_snps <- NULL
  id_boot <- NULL
  local_pval <- NULL
  corrected_pval <- NULL
  error <- NULL
  error_type <- NULL
  N <- NULL
  statistic <- NULL
  snp_type <- NULL
  bootstrap <- NULL

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
  if (length(local_correction) > 1 || length(global_correction) > 1 ||
      !(local_correction %in% mult.test.methods ) ||
      !(global_correction %in% mult.test.methods )) {
    stop("'local_correction' and 'global_correction' must be one of ",
         paste(paste0('"', p.adjust.methods, '"'), collapse=", "), ", or ",
         '"qvalue"')
  }
  if ((local_correction == "qvalue" || global_correction == "qvalue") &&
      !pkgReqCheck("qvalue")) {
    stop("'qvalue' package not installed")
  }

  # Check bootstrap_eSNPs input is ok
  if (length(bootstrap_eSNPs) > 1 || !(bootstrap_eSNPs %in% c("top", "discovery"))) {
    stop("'bootstrap_eSNPs' must be either \"top\" or \"discovery\"")
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
  if (n_bootstraps < 1) {
    stop("'n_bootstraps' must be larger than 0")
  }

  # Check cisDist
  if (length(cisDist) != 1 || !is.numeric(cisDist) || !is.finite(cisDist) || cisDist < 1) {
    stop("'cisDist' must be a number > 0")
  }

  # Check collate_top_snps
  if (length(collate_top_snps) > 1 || !is.finite(collate_top_snps) || !is.logical(collate_top_snps)) {
    stop("'collate_top_snps' must either be TRUE or FALSE")
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
  has.data.table <- isNamespaceLoaded("data.table")
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
    noFDRsaveMemory=ifelse(is.null(eGene_detection_file_name), FALSE, TRUE),
    useModel=modelLINEAR,
    verbose=FALSE
  );

  cis_assocs <- get_cis_assocs(eQTLs, eGene_detection_file_name)
  eGenes <- get_eGenes(cis_assocs, local_correction, global_correction)

  ##--------------------------------------------------------------------
  ## Identify cases where multiple SNPs are the top eSNP
  ##--------------------------------------------------------------------
  cat("Identifying SNPs in perfect LD with lead eSNPs...\n")
  eGenes <- get_eSNPs(cis_assocs, eGenes, collapse=TRUE)

  if (eGenes[global_pval < 0.05, .N] == 0) {
    warning("No significant eGenes detected")
    return(eGenes)
  }

  ##--------------------------------------------------------------------
  ## Apply filters prior to bootstrapping to speed up calculations and
  ## minimise memory usage
  ##--------------------------------------------------------------------

  cat("Optimising data for bootstrap procedure...\n")

  # The threshold for eGene significance is always based on the bonferroni
  # corrected local pvalues - regardless of the multiple correction
  # testing method used for eGene detection.
  #
  # eGenes <- merge(eGenes, snps_per_gene, by="gene")
  # eGenes[, local_bonf := pmin(pvalue * n_snps, 1)]
  n_sig <- eGenes[, sum(global_pval < 0.05)]
  bootstrap_threshold <- eGenes[order(global_pval)][n_sig:(n_sig+1), mean(local_pval)]

  # Get list of significant eGenes and their eSNPs for later
  significant_eGenes <- eGenes[global_pval < 0.05, gene]
  lead_eSNPs <- eGenes[global_pval < 0.05, gsub("/.*", "", snps)] # just take first if in perfect ld

  # Filter to significant eGenes
  gene_boot <- gene$Clone()
  gene_boot$RowReorder(which(rownames(gene_boot) %in% significant_eGenes))
  genepos_boot <- genepos[genepos[,1] %in% significant_eGenes, ]

  # Do we need to run the analysis on all cis-SNPs?
  collate_top_snps <- collate_top_snps || bootstrap_eSNPs == "top" || !(local_correction %in% c("bonferroni", "none"))

  # If we don't need to determine the top SNP at each bootstrap, we can
  # massively speed up the analysis by filtering the snps dataset to
  # just the lead eSNPs.
  if (!collate_top_snps) {
    snps_per_gene <- cis_assocs[gene %in% significant_eGenes, list(n_snps=.N), by=gene]
    snps_boot <- snps$Clone()
    snps_boot$RowReorder(which(rownames(snps_boot) %in% lead_eSNPs))
  } else {
    snps_per_gene <- NULL
    snps_boot <- snps
  }

  # Make sure only necessary objects are ported to each worker
  boot_objs <-  c("cvrt", "snps_boot", "gene_boot", "snpspos",
                  "genepos_boot", "bootstrap_file_directory",
                  "snps_per_gene", "significant_eGenes", "lead_eSNPs",
                  "bootstrap_threshold", "collate_top_snps",
                  "errorCovariance", "cisDist", "useModel",
                  "bootstrap_eSNPs", "local_correction")
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
        genepos = genepos_boot,
        cisDist = cisDist,
        output_file_name=NULL,
        output_file_name.cis=detection_file,
        noFDRsaveMemory=ifelse(is.null(detection_file), FALSE, TRUE),
        verbose=FALSE
      )
      detection_cis_assocs <- get_cis_assocs(eQTL_detection, detection_file)

      # Do eGene detection and effect size estimation in this bootstrap
      if (bootstrap_eSNPs == "discovery") {
        eSNPs <- data.table(gene=significant_eGenes, snps=lead_eSNPs)
        detection_eQTL_SNPs <- get_eGenes(detection_cis_assocs, local_correction, "none", eSNPs, snps_per_gene)
      } else {
        detection_eQTL_SNPs <- get_eGenes(detection_cis_assocs, local_correction, "none")
      }
      # Filter to genes that pass the nominal local correction threshold
      detection_eQTL_SNPs <- detection_eQTL_SNPs[local_pval < bootstrap_threshold]
      # Filter columns
      detection_eQTL_SNPs <- detection_eQTL_SNPs[, list(gene, snps, detection_beta=beta)]
      # Add identifier columns
      detection_eQTL_SNPs[, snp_type := "mapping"]
      detection_eQTL_SNPs[, bootstrap := id_boot]

      if (collate_top_snps) {
        # We also want to collect statistics about the top eSNP in the detection group
        detection_top_SNPs <- get_eGenes(detection_cis_assocs, local_correction, "none", snps_per_gene = snps_per_gene)
        # Filter to significant eGenes
        detection_top_SNPs <- detection_top_SNPs[local_pval < bootstrap_threshold]
        # Get eSNPs with identical statisticis due to perfect LD in bootstrap detection group
        detection_top_SNPs <- get_eSNPs(detection_cis_assocs, detection_top_SNPs, collapse=FALSE)
        # Filter columns
        detection_top_SNPs <- detection_top_SNPs[, list(gene, snps)]
        # Add identifier columns
        detection_top_SNPs[, snp_type := "top"]
        detection_top_SNPs[, bootstrap := id_boot]
      }

      # If there are no significant eGenes in this bootstrap we can move
      # to the next one
      if(nrow(detection_eQTL_SNPs) == 0) {
        if (!(collate_top_snps) || nrow(detection_top_SNPs) == 0) {
          return(NULL)
        } else {
          return(detection_top_SNPs)
        }
      }

      ##----------------------------------------------------------------
      ## Estimation group effect size re-estimation
      ##----------------------------------------------------------------
      tryCatch({
        # Run MatrixEQTL in the estimation group. We only care about the
        # significant SNP-Gene pairs from the detection group, so we can
        # save computation time by filtering their respective matrices.
        # Note: where multiple significant genes or snps are located within
        # the same cis window the filtering will not be perfect.
        sig_boot_genes <- detection_eQTL_SNPs[, gene]
        sig_boot_snps <- detection_eQTL_SNPs[, snps]

        gene_estimation <- gene_boot$Clone()
        gene_estimation$RowReorder(which(gene_estimation$GetAllRowNames() %in% sig_boot_genes))
        gene_estimation$ColumnSubsample(id_estimation)
        snps_estimation <- snps_boot$Clone()
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
          useModel = useModel,
          errorCovariance = errorCovariance,
          pvOutputThreshold.cis = 1,
          snpspos = snpspos,
          genepos = genepos_boot,
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

        if (collate_top_snps) {
          # Add the top detection SNPs table
          sig_boot_assocs <- rbind(sig_boot_assocs, detection_top_SNPs, fill=TRUE)
        }

        return(sig_boot_assocs)
      }, error=function(e) {
        if (exists("detection_top_SNPs") && nrow(detection_top_SNPs) > 0 && detection_top_SNPs[,bootstrap] == id_boot) {
          return(rbind_dt(detection_top_SNPs, data.table(bootstrap=id_boot, error=e$message, error_type="estimation")))
        } else {
          return(data.table(bootstrap=id_boot, error=e$message, error_type="estimation"))
        }
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
  # Relabel columns
  eGenes <- eGenes[,list(gene, top_snp=snps, statistic, nominal_beta=beta,
                         nominal_pval=pvalue, corrected_pval=global_pval)]

  # Calculate the effect of the winner's curse. If effect sizes have been
  # estimated using the top bootstrap eSNP then we need to account for the
  # fact that different eSNPs with anti-correlated minor allele frequencies
  # may have been the top SNP at each bootstrap
  eGenes <- correct_winners_curse(
    boot_eGenes[snp_type == "mapping"], eGenes, correction_type,
    force_sign=ifelse(bootstrap_eSNPs == "top", TRUE, FALSE))

  ##--------------------------------------------------------------------
  ## Collate information about the most frequent eSNPs in bootstrap
  ##--------------------------------------------------------------------

  if (collate_top_snps) {
    cat("Identifying lead SNPs in bootstrap analysis...\n")
    eGenes <- get_bootstrap_lead_SNPs(boot_eGenes[snp_type == "top"], eGenes, cis_assocs)
  }

  # Sort eGenes by significance
  eGenes <- eGenes[order(abs(statistic), decreasing=TRUE)][order(corrected_pval)]

  # If the user has not loaded data.table themselves cast back to a
  # data frame to avoid confusion
  if (!has.data.table) {
    eGenes <- as.data.frame(eGenes)
  }

  on.exit({ cat("Done!\n") }, add=TRUE)
  return(eGenes)
}
