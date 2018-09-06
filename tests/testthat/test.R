context("BootstraQTL")

test_that("Run BootstrapQTL", {
  # Toy data locations
  base.dir = find.package('MatrixEQTL');
  SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
  snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="");
  expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
  gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");
  covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");

  # Load toy data as per MatrixEQTL vignette
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);

  # Add a SNP in perfect LD with the significant eGene-eSNP
  snps$RowReorder(c(1:15,5))
  rownames(snps)[16] <- "Snp_16"  # duplated of Snp_05
  snpspos <- rbind(snpspos, data.frame(snpid="Snp_16", chr="chr1", pos=792481))

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }

  outfile <- tempfile()
  outdir <- tempdir()
  eQTLs <- BootstrapQTL(snps, gene, snpspos, genepos, cvrt, n_bootstraps=10, n_cores=2,
                        eGene_detection_file_name = outfile,
                        bootstrap_file_directory = outdir)
  unlink(outfile)
  unlink(outdir, recursive=TRUE)

  eQTLs <- BootstrapQTL(snps, gene, snpspos, genepos, cvrt, n_bootstraps=10, n_cores=2)
  
  eQTLs <- BootstrapQTL(snps, gene, snpspos, genepos, cvrt, n_bootstraps=0)
})




