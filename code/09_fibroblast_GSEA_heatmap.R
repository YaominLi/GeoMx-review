# Selected fibroblast pathways. Color represents NES; missing estimates remain NA.
# Run the GSEA calculation first.

rm(list = ls()); set.seed(2001)
library(tidyr); library(tibble); library(purrr); library(ggplot2)
library(pheatmap)
library(dplyr)
library(stringr)

df <- read.csv("GSEA_MSigDB_all_subtypes_combined_fibroblast.csv", check.names = FALSE)

selected_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
  "GOBP_CARTILAGE_DEVELOPMENT",
  "GOBP_CARBOHYDRATE_PHOSPHORYLATION",
  "KEGG_GALACTOSE_METABOLISM",
  "GOMF_FATTY_ACID_BINDING",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
  "GOBP_DENDRITIC_CELL_CHEMOTAXIS"
)

sel <- df %>%
  dplyr::filter(
    PathwayID %in% selected_pathways |
      Description %in% selected_pathways |
      PathwayKey %in% paste0("H::", selected_pathways)
  ) %>%
  dplyr::distinct(Description, .keep_all = TRUE)

found_names <- unique(c(sel$PathwayID, sel$Description))
not_found <- setdiff(selected_pathways, found_names)
if (length(not_found)) {
  message("Pathways not found: ", paste(not_found, collapse = ", "))
}

sel <- sel %>%
  mutate(.order = match(Description, selected_pathways)) %>%
  arrange(.order, Description) %>%
  select(-.order)

cols <- c("fibroblast_desert_NES", "fibroblast_excluded_NES", "fibroblast_inflamed_NES")

missing <- setdiff(cols, names(sel))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

mat <- sel %>%
  dplyr::select(dplyr::all_of(cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  as.matrix()

rn <- if ("Description" %in% names(sel)) sel$Description else sel$PathwayID
rownames(mat) <- stringr::str_wrap(rn, width = 35)

colnames(mat) <- c("Desert", "Excluded", "Inflamed")

lim <- max(abs(mat), na.rm = TRUE)

center_low  <- -1
center_high <-  1
if (lim <= 1) { center_low <- -lim; center_high <- lim }

target  <- 1.5
qtarget <- 0.70

if (lim > center_high) {
  ratio <- (min(target, lim) - center_high) / (lim - center_high + 1e-9)
  ratio <- max(min(ratio, 0.999), 1e-6)
  gamma_pos <- log(ratio) / log(qtarget)
  gamma_pos <- max(1.2, min(gamma_pos, 8))
} else {
  gamma_pos <- 1.2
}

n_center <- 81; n_pos <- 40

t_pos      <- seq(0, 1, length.out = n_pos)
pos_breaks <- center_high + (lim - center_high) * (t_pos ^ gamma_pos)
neg_breaks <- -rev(pos_breaks)
center_breaks <- seq(center_low, center_high, length.out = n_center)

bk  <- c(neg_breaks, center_breaks[-1], pos_breaks[-1])

neg_cols    <- colorRampPalette(c("#2b6cb0", "#90cdf4", "#eef2f7"))(length(neg_breaks) - 1)

center_cols <- colorRampPalette(
  c("#dbeafe",
    "#e5e7eb",
    "#fecaca")
)(length(center_breaks) - 1)

pos_cols    <- colorRampPalette(c("#fee2e2", "#fb7185", "#b91c1c"))(length(pos_breaks) - 1)
pal <- c(neg_cols, center_cols, pos_cols)

pheatmap(
  mat,
  color           = pal,
  breaks          = bk,
  na_col          = "grey90",
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = FALSE,
  fontsize        = 10,
  cellheight      = 18,
  cellwidth       = 45,
  angle_col       = 0,
  border_color    = "grey85",
  main            = "Fibroblast pathways",
  filename        = "gsea_selected_pathways_fibroblast.pdf"
)

message("Saved: gsea_selected_pathways_fibroblast.pdf")

# Retain the complete statistics for every displayed pathway.
write.csv(sel, "selected_GSEA_fibroblast_statistics.csv", row.names = FALSE)
