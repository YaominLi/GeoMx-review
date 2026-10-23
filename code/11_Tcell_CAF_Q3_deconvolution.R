# GeoMx T-cell and CAF deconvolution
# Inputs in the working directory:
#   counts.csv, q3norm_nonlog.csv, meta_cleaned_v5.xlsx,
#   ref_scAtlasVAE.rds, CAF_reference.rds
# q3norm_nonlog.csv contains linear Q3-normalized expression.
# Outputs: T-cell IE/II and CAF Desert/Excluded/Inflamed PDF panels only.
# Separate reference matrices define separate compositional denominators.
# Estimates are reference-relative proportions, not absolute tissue cell counts.

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2001)

# 1. Packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(SpatialDecon)

# Validated with R 4.5.1, SpatialDecon 1.18.0 and ggplot2 3.5.2.
# Install SpatialDecon with BiocManager::install("SpatialDecon").
if (as.character(packageVersion("SpatialDecon")) != "1.18.0") {
  warning("Validated SpatialDecon version: 1.18.0; installed version: ",
          as.character(packageVersion("SpatialDecon")))
}

# 2. Input data
getwd()

input_files <- c("counts.csv", "q3norm_nonlog.csv", "meta_cleaned_v5.xlsx",
                 "ref_scAtlasVAE.rds", "CAF_reference.rds")
if (any(!file.exists(input_files))) {
  stop("Missing input files in the working directory: ",
       paste(input_files[!file.exists(input_files)], collapse = ", "))
}

meta <- read_excel("meta_cleaned_v5.xlsx") |> as.data.frame()
cts  <- read.csv("counts.csv", row.names = 1, check.names = FALSE) |> as.matrix()
norm <- read.csv("q3norm_nonlog.csv", row.names = 1, check.names = FALSE) |> as.matrix()

ref_tcell <- readRDS("ref_scAtlasVAE.rds")
ref_caf   <- readRDS("CAF_reference.rds")

# 3. AOI alignment
# Match AOIs by canonical identifiers, then restore metadata identifiers.

stopifnot(all(c("aoi_id", "AOI_code", "Patient_number",
                "infiltration_category") %in% names(meta)))
stopifnot(!anyNA(meta$aoi_id), !anyNA(meta$Patient_number))

clean_aoi <- function(x) gsub("[^a-z0-9]", "", tolower(x))
meta_ids <- clean_aoi(meta$aoi_id)
cts_ids  <- clean_aoi(colnames(cts))
norm_ids <- clean_aoi(colnames(norm))

stopifnot(!anyDuplicated(meta_ids), !anyDuplicated(cts_ids),
          !anyDuplicated(norm_ids), ncol(cts) == nrow(meta),
          ncol(norm) == nrow(meta))
cts_index  <- match(meta_ids, cts_ids)
norm_index <- match(meta_ids, norm_ids)
stopifnot(!anyNA(cts_index), !anyNA(norm_index))

cts  <- cts[, cts_index, drop = FALSE]
norm <- norm[, norm_index, drop = FALSE]
colnames(cts) <- colnames(norm) <- meta$aoi_id
rownames(meta) <- meta$aoi_id

stopifnot(identical(rownames(cts), rownames(norm)),
          all(is.finite(cts)), all(cts >= 0),
          all(is.finite(norm)), all(norm >= 0),
          !anyDuplicated(toupper(rownames(norm))))

# Assign infiltration groups from metadata.
# T-cell iib and iic are both classified as Inflamed (II).
group_map <- c(desert = "ID", excluded = "IE", inflamed = "II")
meta$group <- unname(group_map[tolower(trimws(meta$infiltration_category))])
stopifnot(!anyNA(meta$group))

# 4. Shared SpatialDecon procedure
# Estimate background from the complete linear-Q3 matrix before gene intersection.
# Fit all reference cell types without a fallback estimator.
run_deconvolution <- function(norm_matrix, count_matrix, reference) {

  stopifnot(is.matrix(reference), !is.null(rownames(reference)),
            !is.null(colnames(reference)),
            !anyDuplicated(toupper(rownames(reference))),
            !anyDuplicated(colnames(reference)),
            all(is.finite(reference)), all(reference >= 0))

  bg <- derive_GeoMx_background(
    norm = norm_matrix,
    probepool = rep(1, nrow(norm_matrix)),
    negnames = "NegProbe-WTX"
  )
  stopifnot(identical(dim(bg), dim(norm_matrix)),
            all(is.finite(bg)), all(bg > 0))

  rownames(bg) <- rownames(norm_matrix) <- rownames(count_matrix) <-
    toupper(rownames(norm_matrix))
  rownames(reference) <- toupper(rownames(reference))

  common_genes <- Reduce(intersect, list(rownames(reference),
                                        rownames(norm_matrix),
                                        rownames(count_matrix)))
  stopifnot(length(common_genes) >= 50)

  reference    <- reference[common_genes, , drop = FALSE]
  norm_matrix  <- norm_matrix[common_genes, , drop = FALSE]
  count_matrix <- count_matrix[common_genes, , drop = FALSE]
  bg           <- bg[common_genes, , drop = FALSE]

  # Use syntactically valid cell-type names during fitting.
  original_names <- colnames(reference)
  safe_names <- make.names(original_names, unique = TRUE)
  colnames(reference) <- safe_names

  set.seed(2001)
  fit <- spatialdecon(
    norm = norm_matrix,
    bg = bg,
    X = reference,
    raw = count_matrix,
    resid_thresh = 3,
    lower_thresh = 0.5,
    align_genes = TRUE,
    n_tumor_clusters = 10,
    maxit = 1000
  )

  # Restore reference names in the proportion matrix; retain the original fit structure.
  prop <- fit$prop_of_all
  rownames(prop) <- original_names[match(rownames(prop), safe_names)]
  stopifnot(!anyNA(rownames(prop)),
            setequal(rownames(prop), original_names),
            setequal(colnames(prop), colnames(norm_matrix)),
            all(is.finite(prop)), all(prop >= 0),
            max(abs(colSums(prop) - 1)) < 1e-8)

  list(fit = fit, prop = prop, genes = common_genes,
       name_mapping = data.frame(internal = safe_names, reference = original_names))
}

# 5. T-cell deconvolution
immmeta <- meta[meta$AOI_code == "immune", , drop = FALSE]
immnorm <- norm[, immmeta$aoi_id, drop = FALSE]
immcts  <- cts[, immmeta$aoi_id, drop = FALSE]

message("Running T-cell deconvolution...")
immune_result <- run_deconvolution(immnorm, immcts, ref_tcell)
immdcvn <- immune_result$fit
immprop <- immune_result$prop

# Aggregate the three Tex subtypes. Retain Tpex in the fitted denominator.
tcell_mapping <- list(
  Tex = c("ITGAE+ Tex", "GZMK+ Tex", "S100A11+ Tex"),
  Teff = "ZNF683+ Teff",
  Tn = "Tn",
  `Cycling T` = "Cycling T",
  Tpex = "Tpex"
)
stopifnot(setequal(unlist(tcell_mapping), rownames(immprop)))
tcell_prop <- do.call(rbind, lapply(tcell_mapping, function(x) {
  colSums(immprop[x, , drop = FALSE])
}))
rownames(tcell_prop) <- names(tcell_mapping)
stopifnot(max(abs(colSums(tcell_prop) - 1)) < 1e-8)

# 6. CAF deconvolution
fibmeta <- meta[meta$AOI_code == "fibroblast", , drop = FALSE]
fibnorm <- norm[, fibmeta$aoi_id, drop = FALSE]
fibcts  <- cts[, fibmeta$aoi_id, drop = FALSE]

message("Running CAF deconvolution...")
caf_result <- run_deconvolution(fibnorm, fibcts, ref_caf)
fibdcvn <- caf_result$fit
fibprop <- caf_result$prop

# 7. Group summaries
# Size: AOI-weighted group mean proportion.
# Color: standardize patient-group mean proportions within each cell type
# using the sample SD, then average z-scores within each group.
# Paired patients contribute one mean per observed group.
# Constant types have undefined z-scores. No significance tests are performed.
summarise_bubbles <- function(prop_matrix, aoi_meta, groups) {

  aoi_meta <- aoi_meta[aoi_meta$group %in% groups, , drop = FALSE]
  stopifnot(nrow(aoi_meta) > 0, all(groups %in% aoi_meta$group))
  aoi_meta$group <- factor(aoi_meta$group, levels = groups)

  df_long <- as.data.frame(t(prop_matrix[, aoi_meta$aoi_id, drop = FALSE]))
  df_long$aoi_id <- rownames(df_long)
  df_long <- df_long %>%
    pivot_longer(-aoi_id, names_to = "celltype", values_to = "prop") %>%
    left_join(aoi_meta[, c("aoi_id", "group", "Patient_number")], by = "aoi_id")

  # Average AOIs within each patient and group before standardization.
  patient <- df_long %>%
    group_by(Patient_number, group, celltype) %>%
    summarise(prop = mean(prop), AOIs = n(), .groups = "drop") %>%
    group_by(celltype) %>%
    mutate(z_center = mean(prop), z_sd = sd(prop),
           patient_z = if (sd(prop) > 0) (prop - mean(prop)) / sd(prop) else NA_real_) %>%
    ungroup()

  color_data <- patient %>%
    group_by(group, celltype) %>%
    summarise(mean_z = if (all(is.na(patient_z))) NA_real_ else mean(patient_z),
              patient_mean = mean(prop), .groups = "drop")

  plot_data <- df_long %>%
    group_by(group, celltype) %>%
    summarise(mean_prop = mean(prop), AOIs = n(),
              patients = n_distinct(Patient_number), .groups = "drop") %>%
    left_join(color_data, by = c("group", "celltype"))

  sample_sizes <- aoi_meta %>% group_by(group) %>%
    summarise(AOIs = n(), patients = n_distinct(Patient_number), .groups = "drop")

  list(AOI_data = df_long, patient_data = patient,
       plot_data = plot_data, sample_sizes = sample_sizes)
}

tcell_summary <- summarise_bubbles(tcell_prop, immmeta, c("IE", "II"))
caf_summary_3 <- summarise_bubbles(fibprop, fibmeta, c("ID", "IE", "II"))

# 8. Bubble plots
tcell_order <- c("Tex", "Teff", "Tn", "Cycling T")
caf_order <- c("iCAF", "mCAF", "vCAF", "Pericyte", "tpCAF", "hsp_tpCAF",
               "IDO_CAF", "apCAF", "rCAF", "dCAF")
stopifnot(all(caf_order %in% rownames(fibprop)))

plot_bubbles <- function(summary_data, cell_order, groups, group_labels, title) {

  df_plot <- summary_data$plot_data %>%
    filter(celltype %in% cell_order) %>%
    mutate(group = factor(group, levels = groups),
           celltype = factor(celltype, levels = rev(cell_order)))

  ggplot(df_plot, aes(group, celltype)) +
    geom_point(aes(size = mean_prop, fill = mean_z),
               shape = 21, color = "grey40", stroke = 0.3) +
    scale_size(range = c(5, 9), limits = c(0, 1), breaks = c(0, 0.2, 0.75),
               labels = label_percent(accuracy = 1),
               name = "Mean proportion\n(compressed size)") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
      limits = c(-0.5, 0.5), oob = scales::squish,
      breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5"),
      na.value = "grey80", name = "Scaled abundance\n(mean z-score)"
    ) +
    guides(
      size = guide_legend(order = 1,
                          override.aes = list(fill = "grey55", colour = "grey55")),
      fill = guide_colorbar(order = 2, barheight = grid::unit(22, "mm"),
                            barwidth = grid::unit(4, "mm"))
    ) +
    scale_x_discrete(limits = groups, labels = group_labels,
                     expand = expansion(add = 0.35)) +
    scale_y_discrete(limits = rev(cell_order), expand = expansion(add = 0.45)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = 13),
      axis.text.x = element_text(margin = margin(t = 3)),
      axis.text.y = element_text(margin = margin(r = 3)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 15, color = "black",
                                margin = margin(b = 4)),
      legend.position = "right",
      legend.title = element_text(size = 11, color = "black"),
      legend.text = element_text(size = 11, color = "black"),
      legend.key.height = grid::unit(8, "mm"),
      legend.key.width = grid::unit(9, "mm"),
      legend.spacing.y = grid::unit(2, "mm"),
      legend.box.spacing = grid::unit(1.5, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(3, 4, 3, 3)
    )
}

p_tcell <- plot_bubbles(tcell_summary, tcell_order, c("IE", "II"),
                        c(IE = "Excluded", II = "Inflamed"), "T-cell") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9, color = "black"),
        legend.location = "plot",
        legend.justification = "top",
        legend.key.height = grid::unit(6, "mm"),
        legend.spacing.y = grid::unit(1, "mm"))
p_caf_3 <- plot_bubbles(caf_summary_3, caf_order, c("ID", "IE", "II"),
                        c(ID = "Desert", IE = "Excluded", II = "Inflamed"), "CAF")

# 9. Export the two final panels to the working directory.
# Point sizes use a compressed 5-9 scale; area is not proportional to abundance.
# Zero estimates retain a minimum-size bubble; undefined z-scores are gray.
# Colors saturate at -0.5 and +0.5 without modifying the underlying values.
# Tpex is omitted from the selected T-cell panel without renormalization.
ggsave("Tcell_IE_II_scaled_abundance_bubble.pdf", p_tcell, device = cairo_pdf,
       width = 3.0, height = 2.9, bg = "white")
ggsave("CAF_ID_IE_II_scaled_abundance_bubble.pdf", p_caf_3, device = cairo_pdf,
       width = 4.7, height = 4.2, bg = "white")

message("Saved Tcell_IE_II_scaled_abundance_bubble.pdf and CAF_ID_IE_II_scaled_abundance_bubble.pdf to: ", getwd())
