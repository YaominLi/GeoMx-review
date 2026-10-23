# Exercise the released scripts with synthetic inputs on Windows.
# These fixtures are software tests and are not study observations.

stopifnot(.Platform$OS.type == "windows")
source("code/00_install_environment.R")

local({
  project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  fixture <- file.path(tempdir(), "\u6587\u6863", "GeoMx test")
  dir.create(fixture, recursive = TRUE, showWarnings = FALSE)
  setwd(fixture)
  on.exit(setwd(project), add = TRUE)
  set.seed(6193)

  genes <- sprintf("GENE%04d", seq_len(300))
  meta <- expand.grid(Patient_number = seq_len(4),
                      infiltration_category = c("desert", "excluded", "inflamed"),
                      AOI_code = c("tumor", "fibroblast", "immune"),
                      stringsAsFactors = FALSE)
  meta$Slide_number <- meta$Patient_number
  meta$aoi_id <- sprintf("AOI%03d", seq_len(nrow(meta)))
  counts <- matrix(0, 301, nrow(meta),
                   dimnames = list(c(genes, "NegProbe-WTX"), meta$aoi_id))
  baseline <- runif(300, 100, 300)
  for (j in seq_len(nrow(meta))) {
    g <- match(meta$infiltration_category[j], c("desert", "excluded", "inflamed"))
    mu <- baseline * exp(0.08 * meta$Patient_number[j])
    marker <- seq.int((g - 1) * 35 + 1, g * 35)
    mu[marker] <- mu[marker] * 2.5
    counts[seq_len(300), j] <- rnbinom(300, mu = mu, size = 8)
    counts[301, j] <- 5
  }
  q3 <- apply(counts, 2, quantile, probs = 0.75)
  norm <- sweep(counts, 2, q3 / median(q3), "/")
  write.csv(counts, "counts.csv")
  write.csv(norm, "q3norm_nonlog.csv")
  openxlsx::write.xlsx(meta, "meta_cleaned_v5.xlsx")
  make_reference <- function(types) {
    x <- matrix(runif(300 * length(types), 1, 5), 300,
                dimnames = list(genes, types))
    for (i in seq_along(types)) {
      idx <- seq.int((i - 1) * 20 + 1, i * 20)
      x[idx, i] <- x[idx, i] + 20
    }
    x
  }
  saveRDS(make_reference(c("ITGAE+ Tex", "GZMK+ Tex", "S100A11+ Tex",
                           "ZNF683+ Teff", "Tn", "Cycling T", "Tpex")),
          "ref_scAtlasVAE.rds")
  saveRDS(make_reference(c("iCAF", "mCAF", "vCAF", "Pericyte", "tpCAF",
                           "hsp_tpCAF", "IDO_CAF", "apCAF", "rCAF", "dCAF")),
          "CAF_reference.rds")

  scripts <- c("01_tumor_TMM_voom_heatmap_PCA.R",
               "02_fibroblast_TMM_voom_heatmap_PCA.R",
               "03_immune_TMM_voom_heatmap_PCA.R", "05_tumor_gene_overlap.R",
               "06_prepare_GSEA_rankings.R", "11_Tcell_CAF_Q3_deconvolution.R")
  for (script in scripts) {
    message("WINDOWS TEST: ", script)
    source(file.path(project, "code", script), local = .GlobalEnv)
  }
  stopifnot(length(list.files(pattern = "\\.pdf$")) == 9,
            length(list.files("enrichment_inputs", pattern = "\\.csv$")) == 13,
            all(file.info(list.files(pattern = "\\.pdf$"))$size > 1000))
  for (cell in c("tumor", "fibroblast", "immune")) {
    d <- read.csv(paste0("all_genes_DE_with_expr_and_fdr_", cell, ".csv"))
    stopifnot(nrow(d) == 301, all(is.finite(d$AveExpr)))
  }

  # Use the Windows socket backend selected by script 07.
  d <- read.csv("enrichment_inputs/subtype_tumor_excluded.csv")
  terms <- data.frame(term = rep(c("TEST_SET_A", "TEST_SET_B"), each = 50),
                      gene = c(genes[1:50], genes[51:100]))
  bp <- BiocParallel::SnowParam(4, type = "SOCK", progressbar = FALSE)
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  set.seed(2001)
  result <- clusterProfiler::GSEA(setNames(d$score, d$Gene), TERM2GENE = terms,
                                  minGSSize = 10, maxGSSize = 500,
                                  pvalueCutoff = 1, pAdjustMethod = "BH",
                                  verbose = FALSE, seed = TRUE, BPPARAM = bp)
  result <- as.data.frame(result)
  stopifnot(nrow(result) == 2, all(is.finite(result$NES)))
  record <- file.path(project, "renv", "setup-report", "windows_smoke.txt")
  writeLines(c("WINDOWS SMOKE PASSED",
               "Synthetic inputs; working path includes Chinese characters and a space.",
               "Six released scripts; nine PDFs; thirteen rankings; socket GSEA.",
               capture.output(sessionInfo())), record)
  message("WINDOWS SMOKE PASSED")
})
