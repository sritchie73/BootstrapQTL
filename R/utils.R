### Set up a parallel backend
###
### Set up a backend with the requested number of cores, or use existing backend
### if the user has set one up already.
###
### @param nCores number of cores requested.
### @param verbose logical; Is verbose printing on?
### @param reporterCore logical; is there a reporter core?
###
### @details
###  If \code{verbose = TRUE} and \code{reporterCore = TRUE} then an extra core
###  will be registered for reporting: This core simply sits and prints out
###  progress bars for the permutation procedure.
###
### @return
###  A list containing the total number of cores registered, the registered
###  cluster if on a Windows machine, whether an existing parallel backend
###  is being used, and the number of OMP and BLAS threads set in the current
###  R session.
setupParallel <- function(nCores, verbose, reporterCore) {
  if (is.null(nCores)) {
    if (pkgReqCheck("parallel")) {
      if (parallel::detectCores() > 1) {
        nCores <- parallel::detectCores() - 1
      } else {
        nCores <- 1
      }
    } else {
      nCores <- 1
      if (verbose) cat("Unable to find 'parallel' package, running on 1 core.\n")
    }
  }
  if (!is.numeric(nCores) || length(nCores) > 1 || nCores < 1)
    stop("'n_cores' must be a single number greater than 0")

  # Defaults to return
  cl <- NULL
  predef <- FALSE

  # First, check whether the user has already set up a parallel backend. In this
  # case, we can ignore the `nCores` argument.
  if (foreach::getDoParWorkers() > 1) {
    if (verbose) cat("Ignoring 'n_cores': parallel backend detected.\n")
    if (reporterCore) {
      if (verbose) {
        cat("Reserving 1 core for progress reporting.", foreach::getDoParWorkers() - 1,
            "cores will be used for computation.\n")
      }
    }
    nCores <- foreach::getDoParWorkers()
    predef <- TRUE
  }
  # If the user is on a Windows machine, we have to use the `doParallel` package
  else if (.Platform$OS.type == "windows") {
    # Quietly load parallel backend packages. Throw our own warning and
    # continue
    if(pkgReqCheck("doParallel")) {
      # we need an additional thread to monitor and report progress
      workerCores <- nCores
      if (verbose && reporterCore) {
        nCores <- nCores + 1
      }

      cl <- parallel::makeCluster(nCores, type="PSOCK")
      doParallel::registerDoParallel(cl)

      if (verbose) cat("Registering", workerCores, "cores for bootstrap procedure.\n")
      if (workerCores > parallel::detectCores()) {
        stop(
          "Requested number of cores (", workerCores, ") is higher than the ",
          "number of available cores (", parallel::detectCores(),
          "). Using too many cores may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end
      # once the analysis has finished.
      cat("Warning: unable to utilise multiple cores. Please install the 'doParallel' package",
          "to enable parallel computation.\n", file=stderr())
      warning("Package required for parallel computation not installed")
    }
  } else if (.Platform$OS.type == "unix" && nCores > 1) {
    # Quietly load parallel backend packages. Throw our own warning and
    # continue
    if (pkgReqCheck("doMC")) {
      # we need an additional thread to monitor and report progress
      workerCores <- nCores
      if (verbose && reporterCore) {
        nCores <- nCores + 1
      }
      doMC::registerDoMC(nCores)
      if (verbose) cat("Registering", workerCores, "cores for bootstrap procedure.\n")
      if ((nCores - 1) > parallel::detectCores()) {
        stop(
          "Requested number of cores (", workerCores, ") is higher than the ",
          "number of available cores (", parallel::detectCores(),
          "). Using too many cores may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end
      # once the analysis has finished.
      cat("Warning: unable to utilise multiple cores. Please install the 'doMC' package",
          "to enable parallel computation.\n", file=stderr())
      warning("Package required for parallel computation not installed")
    }
  } else {
    if (verbose) cat("Registering 1 core for bootstrap procedure.\n")
  }

  # Suppress annoying foreach warning generated when using %dopar% and running
  # in serial
  if (nCores == 1) {
    suppressWarnings({
      ii <- 0 # suppress R CMD check note
      foreach(ii = 1:2) %dopar% { ii }
    })
  }

  return(list(nCores=nCores, cluster=cl, predef=predef))
}

### De-register a parallel backend
###
### @param cluster registered cluster on a Windows machine
### @param predef logical; was a pre-existing parallel backend used?
cleanupCluster <- function(cluster, predef) {
  if (!is.null(cluster)) {
    if (pkgReqCheck("parallel")) {
      # Clobber the backend
      parallel::stopCluster(cluster)
      cl <- parallel::makeCluster(1, type="PSOCK")
      doParallel::registerDoParallel(cl)
      closeAllConnections()
    }
  } else if (!predef) {
    if (pkgReqCheck("doMC")) {
      doMC::registerDoMC(1)
    }
  }
}

### Silently check and load a package into the namespace
###
### @param pkg name of the package to check
###
### @return logical; \code{TRUE} if the package is installed and can be loaded.
pkgReqCheck <- function(pkg) {
  suppressMessages(suppressWarnings(requireNamespace(pkg, quietly=TRUE)))
}

### Get table of cis associations from MatrixEQTL output
###
### @param meqtl_obj object returend by Matrix_eQTL_main
### @param output_file file associations were saved to
###
### @return data.table of all cis-associations
get_cis_assocs <- function(meqtl_obj, output_file) {
  if (is.null(output_file)) {
    cis_assocs <- as.data.table(meqtl_obj$cis$eqtls)
  } else {
    cis_assocs <- fread(output_file)
    setnames(cis_assocs, c("SNP", "t-stat", "p-value"), c("snps", "statistic", "pvalue"))
  }
  return(cis_assocs)
}

### rbind two data tables while filling in missing columns with NA
###
### used as an argument to foreach's .combine
rbind_dt <- function(...) {
  rbind(..., fill=TRUE)
}
