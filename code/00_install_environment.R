# Install the locked project environment from the repository root.
# Run in a fresh R 4.5.1 session before the analysis scripts.

local({
  required_r <- "4.5.1"
  required_renv <- "1.1.5"
  required_bioc <- "3.21"
  project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  lockfile <- file.path(project, "renv.lock")

  if (!file.exists(lockfile) || !file.exists("GeoMx_revised.Rproj")) {
    stop("Open GeoMx_revised.Rproj and run this script from the repository root.",
         call. = FALSE)
  }
  if (getRversion() != package_version(required_r)) {
    stop("This environment requires R ", required_r, "; current R is ",
         getRversion(), ". Install and select R 4.5.1, restart R, and rerun.",
         call. = FALSE)
  }
  if (!isTRUE(capabilities("cairo"))) {
    stop("This R installation needs Cairo support for the analysis PDFs.",
         call. = FALSE)
  }
  if ("renv" %in% loadedNamespaces() &&
      getNamespaceVersion("renv") != package_version(required_renv)) {
    stop("A different renv version is loaded. Restart R without loading renv, ",
         "then rerun this script.", call. = FALSE)
  }

  # Bootstrap the recorded renv version in a project-local library.
  bootstrap <- file.path(project, "renv", "bootstrap", required_renv)
  dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(bootstrap) || file.access(bootstrap, 2L) != 0L) {
    stop("The project directory is not writable: ", project, call. = FALSE)
  }
  options(timeout = max(600, getOption("timeout", 60)), renv.consent = TRUE)
  if (!("renv" %in% loadedNamespaces())) {
    installed <- tryCatch(utils::packageVersion("renv", lib.loc = bootstrap),
                          error = function(e) NULL)
    if (is.null(installed) || installed != package_version(required_renv)) {
      archive <- tempfile("renv-", fileext = ".tar.gz")
      on.exit(unlink(archive), add = TRUE)
      url <- paste0("https://cran.r-project.org/src/contrib/Archive/renv/renv_",
                    required_renv, ".tar.gz")
      message("Installing renv ", required_renv, "...")
      utils::download.file(url, archive, mode = "wb")
      utils::install.packages(archive, repos = NULL, type = "source", lib = bootstrap)
    }
    .libPaths(c(bootstrap, .libPaths()))
    loadNamespace("renv", lib.loc = bootstrap)
  }
  stopifnot(getNamespaceVersion("renv") == package_version(required_renv))

  lock <- renv::lockfile_read(lockfile)
  if (!identical(lock$R$Version, required_r) ||
      !identical(lock$Bioconductor$Version, required_bioc) ||
      !identical(lock$Packages$renv$Version, required_renv)) {
    stop("The setup script and renv.lock specify different environment versions.",
         call. = FALSE)
  }
  expected <- vapply(lock$Packages, function(x) x$Version, character(1))
  loaded <- intersect(names(expected), loadedNamespaces())
  mismatched <- loaded[vapply(loaded, function(pkg) {
    getNamespaceVersion(pkg) != package_version(expected[[pkg]])
  }, logical(1))]
  if (length(mismatched)) {
    stop("Restart R before installation; different package versions are loaded: ",
         paste(mismatched, collapse = ", "), call. = FALSE)
  }

  # Bootstrap BiocManager before activating the project library.
  required_biocmanager <- lock$Packages$BiocManager$Version
  if (!("BiocManager" %in% loadedNamespaces())) {
    installed <- tryCatch(utils::packageVersion("BiocManager", lib.loc = bootstrap),
                          error = function(e) NULL)
    if (is.null(installed) || installed != package_version(required_biocmanager)) {
      bioc_archive <- tempfile("BiocManager-", fileext = ".tar.gz")
      on.exit(unlink(bioc_archive), add = TRUE)
      url <- paste0("https://cran.r-project.org/src/contrib/Archive/BiocManager/",
                    "BiocManager_", required_biocmanager, ".tar.gz")
      utils::download.file(url, bioc_archive, mode = "wb")
      utils::install.packages(bioc_archive, repos = NULL, type = "source",
                              lib = bootstrap)
    }
    loadNamespace("BiocManager", lib.loc = bootstrap)
  }
  stopifnot(getNamespaceVersion("BiocManager") ==
              package_version(required_biocmanager))

  # Install packages sequentially in dependency order.
  previous_install_options <- options(Ncpus = 1L, renv.config.cache.symlinks = FALSE)
  on.exit(options(previous_install_options), add = TRUE)

  # Restore every recorded package without updating the lockfile.
  renv::load(project = project)
  project_library <- renv::paths$library(project = project)
  dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
  biocmanager <- paste0("BiocManager@", lock$Packages$BiocManager$Version)
  renv::install(biocmanager, project = project, library = project_library,
                prompt = FALSE)
  renv::restore(project = project, library = project_library,
                lockfile = lockfile, prompt = FALSE, transactional = TRUE)
  renv::load(project = project)
  project_library <- renv::paths$library(project = project)

  # Include matching recommended packages supplied with R.
  verification_libraries <- unique(c(project_library, .Library))
  actual <- vapply(names(expected), function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg, lib.loc = verification_libraries)),
             error = function(e) NA_character_)
  }, character(1))
  matched <- vapply(names(expected), function(pkg) {
    !is.na(actual[[pkg]]) &&
      package_version(actual[[pkg]]) == package_version(expected[[pkg]])
  }, logical(1))
  if (!all(matched)) {
    stop("Package verification failed: ",
         paste(names(expected)[!matched], collapse = ", "), call. = FALSE)
  }

  analysis_packages <- c(
    "BiocParallel", "clusterProfiler", "dplyr", "edgeR", "fgsea", "ggplot2",
    "limma", "matrixStats", "pheatmap", "purrr", "RColorBrewer", "readr",
    "readxl", "scales", "SpatialDecon", "statmod", "stringr", "tibble", "tidyr"
  )
  for (pkg in analysis_packages) {
    tryCatch(loadNamespace(pkg, lib.loc = project_library), error = function(e) {
      stop("Cannot load ", pkg, ": ", conditionMessage(e), call. = FALSE)
    })
  }

  if (as.character(BiocManager::version()) != required_bioc) {
    stop("The installed Bioconductor release does not match renv.lock.", call. = FALSE)
  }
  loaded <- intersect(names(expected), loadedNamespaces())
  mismatched <- loaded[vapply(loaded, function(pkg) {
    getNamespaceVersion(pkg) != package_version(expected[[pkg]])
  }, logical(1))]
  if (length(mismatched)) {
    stop("Restart R and rerun setup; loaded package versions differ: ",
         paste(mismatched, collapse = ", "), call. = FALSE)
  }

  # Keep a local installation record outside the analysis outputs.
  report_dir <- file.path(project, "renv", "setup-report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(Package = names(expected), Expected = expected,
                             Installed = actual, row.names = NULL),
                   file.path(report_dir, "package_versions.csv"), row.names = FALSE)
  writeLines(c(paste("Checked:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
               paste("Lockfile MD5:", unname(tools::md5sum(lockfile))),
               paste("Bioconductor:", required_bioc),
               paste("Project library:", project_library),
               capture.output(sessionInfo())),
             file.path(report_dir, "sessionInfo.txt"))
  message("Environment ready: R ", required_r, "; Bioconductor ", required_bioc,
          "; ", length(expected), " package versions verified.")
  message("Continue with the analysis scripts. Installation record: renv/setup-report/")
})
