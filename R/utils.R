
### Verbose Concatenate and Print with Indentation
###
### Concatenate and output the objects only if the \code{verbose} flag is set
### to \code{TRUE}. Allows for indentation, adding a series of spaces to the
### beginning of each line, 2 for every increment in \code{ind}.
###
### @details
###  \code{vCat} is slightly more intelligent than regulat \code{cat} in the way
###  it formats the output and breaks it into lines. As a result, \code{fill} is
###  set to \code{TRUE} by default. Another notable difference from \code{cat} is
###  the way newline objects are handled. For example, the call:
###  \code{cat("hello", "world", "\\n", "foo", "bar")} won't wrap the newline
###  character with spaces. This avoids the need to set \code{sep} to \code{""}
###  and embed multiple \code{paste} calls. Finally, a newline character is
###  appended to the end of the whole message, avoiding the need to manually
###  specify this when calling \code{vCat}.
###
### @seealso \code{\link[base]{cat}}
### @param verbose logical. If \code{TRUE}, passes the rest of the arguments to
###   \code{\link{cat}}
### @param ind an integer corresponding to the level of indentation. Each
###   indentation level corresponds to two spaces.
### @param ... Arguments to pass to \code{\link[base]{cat}}
### @param sep a character vector of strings to append after each element.
### @param fill a logical or (positive) numeric controlling how the output is
###   broken into successive lines. If \code{FALSE}, only newlines created
###   explicitly by "\\n" are printed. Otherwise, the output is broken into lines
###   with print width equal to the option width if fill is \code{TRUE}
###   (default), or the value of fill if this is numeric. Non-positive fill
###   values are ignored, with a warning.
### @param labels character vector of labels for the lines printed. Ignored if
###   fill is \code{FALSE}.
###
### @keywords internal
vCat <- function(verbose, ind=0,  ..., sep=" ", fill=TRUE, labels=NULL) {
  if (!(is.vector(verbose) && !is.list(verbose) && is.logical(verbose) &&
        length(verbose) == 1 && !is.na(verbose))) {
    stop("'verbose' must be one of 'TRUE' or 'FALSE'")
  }

  if(verbose) {
    # Put a timestamp at the start of the message:
    if (is.null(labels))
      labels <- paste0("[", format(Sys.time(), usetz=TRUE), "] ")

    # We need to format each line with the indendation level
    if (ind > 0) {
      indent <- paste(rep("  ", ind), collapse="")
    } else {
      indent = ""
    }

    args <- list(...)
    if (is.null(names(args))) {
      str <- paste(args, collapse=sep)
      named <- NULL
    } else {
      str <- paste(args[names(args) == ""], collapse=sep)
      named <- args[names(args) != ""]
    }
    str <- gsub(" \n ", "\n", str) # make it easier to insert newlines
    lines <- strsplit(str, "\n")[[1]]
    # Handle automatic line wrapping
    if (fill) {
      if (is.logical(fill)) {
        fillWidth <- options("width")
      } else if (!is.numeric(fill)) {
        stop("invalid 'fill' argument")
      } else if (fill < 1) {
        warning("non-positive 'fill' argument will be ignored")
        fillWidth <- options("width")
      } else {
        fillWidth <- fill
      }
      words <- strsplit(lines, " ")
      # Create new lines by accumulating words in each line until the addition
      # of the next word would exceed the fillWidth.
      formatted <- lapply(words, function(lw) {
        newlines <- c("")
        curnl <- 1
        for (w in lw) {
          if (newlines[curnl] == "") {
            newlines[curnl] <- paste0(labels, indent, w)
          } else if(nchar(newlines[curnl]) + 1 + nchar(w) < fillWidth) {
            newlines[curnl] <- paste(newlines[curnl], w)
          } else {
            curnl <- curnl + 1
            labels <- paste(rep(" ", nchar(labels)), collapse="")
            newlines[curnl] <- paste0(labels, indent, w)
          }
        }
        paste(newlines, collapse="\n")
      })
      lines <- strsplit(paste(formatted, collapse="\n"), "\n")[[1]]
    }
    str <- paste0(lines, "\n", collapse="")
    if (is.null(named)) {
      cat(str)
    } else {
      # build expression from remaining named arguments
      args = sapply(seq_along(named), function(n) {
        if (is.character(named[[n]])) {
          paste0(names(named)[n], "=", "'", named[[n]], "'")
        } else {
          paste0(names(named)[n], "=", named[[n]])
        }
      })
      eval(parse(text=paste0(paste0(c("cat(str", args), collapse=", "), ")")))
    }
  }
}

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
      vCat(verbose, 1, "Unable to find 'parallel' package, running on 1 core")
    }
  }
  if (!is.numeric(nCores) || length(nCores) > 1 || nCores < 1)
    stop("'nCores' must be a single number greater than 0")

  # Defaults to return
  cl <- NULL
  predef <- FALSE

  # First, check whether the user has already set up a parallel backend. In this
  # case, we can ignore the `nCores` argument.
  if (foreach::getDoParWorkers() > 1) {
    vCat(verbose, 1, "Ignoring 'nCores': parallel backend detected.")
    if (reporterCore) {
      vCat(
        verbose, 1, "Reserving 1 core for progress reporting.",
        foreach::getDoParWorkers() - 1, "cores will be used for computation"
      )
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

      vCat(verbose, 1, "Running on", workerCores, "cores.")
      if (workerCores > parallel::detectCores()) {
        stop(
          "Requested number of threads (", workerCores, ") is higher than the ",
          "number of available cores (", parallel::detectCores(),
          "). Using too many threads may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end
      # once the analysis has finished.
      vCat(
        TRUE, 1, file=stderr(),
        "Warning: running on 1 core. Please install the 'doParallel' package",
        "to enable parallel computation"
      )
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
      vCat(verbose, 1, "Running on", workerCores, "cores.")
      if ((nCores - 1) > parallel::detectCores()) {
        stop(
          "Requested number of threads (", workerCores, ") is higher than the ",
          "number of available cores (", parallel::detectCores(),
          "). Using too many threads may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end
      # once the analysis has finished.
      vCat(
        TRUE, 1, file=stderr(),
        "Warning: running on 1 core. Please install the 'doMC' package",
        "to enable parallel computation"
      )
      warning("Package required for parallel computation not installed")
    }
  } else {
    vCat(verbose, 1, "Running on 1 cores.")
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
