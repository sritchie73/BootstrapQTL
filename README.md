# BootstrapQTL
## A *cis*-QTL mapping method that corrects for the Winner's Curse

BootstrapQTL is a *cis*-QTL mapping tool that uses a fast bootstrap procedure 
to correct for the overestimation of effect sizes present in  *cis*-QTL 
mapping studies ("The Winner's Curse effect").

For more information see the associated publication in *Nucleic Acids Research*; [Power, false discovery rate and Winner’s Curse in eQTL studies](https://doi.org/10.1093/nar/gky780).

## Installation

The latest stable version of BootstrapQTL can be installed directly from
from either this GitHub repository or from CRAN:

```{r}
# GitHub Install:
library(devtools)
install_github("InouyeLab/BootstrapQTL")

# CRAN Install:
install.packages("BootstrapQTL")
```

## Package tutorial

BootstrapQTL makes use of the [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/)
package and therefore requires data to be loaded into R as per the 
[MatrixEQTL tutorial](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html).

The following script, package interface, and package documentation all
describe the use of BootstrapQTL for cis-eQTL mapping. However, this 
approach and package can be applied more generally to any QTL study of
quantitative traits with chromosomal positions, for example cis-QTL 
studies of epigenetic modifications (described below). 

The following code shows an example of loading genotype data, gene 
expression data, covariates data, snp position data, and gene position
data using the MatrixEQTL package example dataset:

```{r}
library(BootstrapQTL)

# Locations for example data from the MatrixEQTL package
base.dir = find.package('MatrixEQTL');
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");

# Load the SNP data - columns must be samples and rows genotypes
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# Load the gene expression data - columns must be samples and rows genes
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# Load the covariates data - columns must be samples and rows covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariates_file_name);

# Note, the ordering of columns in 'snps', 'genes', and 'cvrt' must be
# identical.

# Load the data detailing the position of each SNP - this should contain
# three columns: 
#  (1) 'snpid' describing the name of the SNP and corresponding to rows 
#       in the 'snps' matrix.
#  (2) 'chr' describing the chromosome for each SNP.
#  (3) 'pos' describing the position of the SNP on the chromosome.
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);


# Load the data detailing the position of each gene - this should contain
# four columns:
#  (1) 'geneid' describing the name of the gene and corresponding to rows 
#       in the 'gene' matrix.
#  (2) 'chr' describing the chromosome for each SNP.
#  (3) 'left' describing the start position of the transcript.
#  (4) 'right' describing the end position of the transcript.
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
```

Once the data has been loaded into R the BootstrapEQTL analysis can be
run with a single command:

```{r}
# Run the BootstrapQTL analysis
eGenes <- BootstrapQTL(snps, gene, snpspos, genepos,
                        n_bootstraps=200, n_cores=2,
                        eGene_detection_file_name = "cis_eQTL_associations.txt",
                        bootstrap_file_directory = "bootstrap_analyses/")
```

Analysis of different types of  molecular data procedes in a similar 
fashion as shown above. Simply load in the quantitative trait data into 
R as shown above for the `gene` variable, then create and load a 
`data.frame` of positional information as shown above for the `genepos`
variable. If the molecular traits have a single genomic position, rather
than spanning a range of base pairs, then you will need to duplicate 
the trait position in the `genepos` `data.frame` (*i.e.* the `"left"`,
and `"right"` columns will contain the same location).
