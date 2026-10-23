# Selected composite-ranking GSEA bars.
# Bar length and fill represent NES; BH FDR resolves duplicate pathway entries.
rm(list = ls()); set.seed(2001)
suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr); library(ggplot2)})

input <- "gsea_results/GSEA_three_sets_long.csv"
if (!file.exists(input)) stop("Run the GSEA calculation before this script.")
x <- read_csv(input, show_col_types = FALSE)
stopifnot(all(c("SubtypeShort", "Collection", "Description", "NES", "pvalue", "p.adjust") %in% names(x)))
dir.create("quick_plots", showWarnings = FALSE)

selected <- list(
  EXC = c("GOBP_EMBRYONIC_SKELETAL_SYSTEM_MORPHOGENESIS", "GOBP_NEURON_CELL_CELL_ADHESION",
          "GOMF_WNT_PROTEIN_BINDING", "THEODOROU_MAMMARY_TUMORIGENESIS", "GOMF_GAP_JUNCTION_CHANNEL_ACTIVITY"),
  INT = c("GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT", "FARDIN_HYPOXIA_11",
          "GOCC_COMPLEX_OF_COLLAGEN_TRIMERS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_APICAL_JUNCTION"),
  DES = c("FARDIN_HYPOXIA_11", "WP_AEROBIC_GLYCOLYSIS_AUGMENTED", "GOCC_DESMOSOME",
          "REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS", "GOMF_INSULIN_LIKE_GROWTH_FACTOR_BINDING")
)
colours <- list(EXC = c("#D7ECFF", "#1E64B7"), INT = c("#EBD4FF", "#5B2ABF"), DES = c("#DDDDDD", "#000000"))
for (set_code in names(selected)) {
  d <- x %>% filter(SubtypeShort == set_code, Description %in% selected[[set_code]])
  if (!nrow(d)) stop("Selected pathways are missing: ", set_code)
  p_col <- if (any(is.finite(d$p.adjust))) "p.adjust" else "pvalue"
  d <- d %>% mutate(p_use = .data[[p_col]]) %>% group_by(Description) %>%
    slice_min(p_use, n = 1, with_ties = FALSE) %>% ungroup() %>%
    filter(is.finite(NES), NES > 1) %>% arrange(desc(NES)) %>%
    mutate(label = str_trunc(paste0(Collection, "::", Description), width = 70))
  if (!nrow(d)) stop("No selected pathways have NES > 1: ", set_code)
  write_csv(d, file.path("quick_plots", paste0("Selected_", set_code, "_statistics.csv")))
  d$label_f <- factor(d$label, levels = rev(d$label))
  p <- ggplot(d, aes(NES, label_f, fill = NES)) +
    geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", NES)), hjust = -0.15, size = 3) +
    scale_fill_gradient(low = colours[[set_code]][1], high = colours[[set_code]][2], name = "NES") +
    coord_cartesian(xlim = c(0, max(d$NES) * 1.15)) +
    labs(title = paste0("Selected terms (", set_code, ") - NES>1"), x = "NES", y = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold"), legend.position = "right")
  ggsave(file.path("quick_plots", paste0("Selected_", set_code, "_horizontal.pdf")),
         p, width = 8, height = max(3, 0.45 * nrow(d) + 1.4))
}
