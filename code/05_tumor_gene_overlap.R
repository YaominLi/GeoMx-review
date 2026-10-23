# Tumor genes upregulated in Desert or Excluded relative to Inflamed.
# Set membership uses nominal P < 0.05 without a Top-N restriction.
# DI-only and EI-only do not establish a Desert-versus-Excluded difference.
library(grid)

master <- read.csv("all_genes_DE_with_expr_and_fdr_tumor.csv", check.names = FALSE)
di <- master$Gene[master$logFC_desert_vs_inflamed > 0 & master$P.desert_vs_inflamed < 0.05]
ei <- master$Gene[master$logFC_excluded_vs_inflamed > 0 & master$P.excluded_vs_inflamed < 0.05]
gene_sets <- list(DI_only = setdiff(di, ei), shared = intersect(di, ei), EI_only = setdiff(ei, di))
for (nm in names(gene_sets)) {
  write.csv(data.frame(Gene = gene_sets[[nm]]), paste0("tumor_", nm, "_genes.csv"), row.names = FALSE)
}
write.csv(data.frame(set = names(gene_sets), genes = lengths(gene_sets)),
          "tumor_intersection_counts.csv", row.names = FALSE)

pdf("tumor_DI_EI_venn.pdf", width = 6, height = 4)
grid.newpage()
# Equal-area schematic; circle area does not encode gene counts.
pushViewport(viewport(width = unit(3.7, "inches"), height = unit(3.7, "inches")))
grid.circle(x = 0.35, y = 0.48, r = 0.29,
            gp = gpar(fill = adjustcolor("grey60", alpha.f = 0.5), col = "grey30"))
grid.circle(x = 0.65, y = 0.48, r = 0.29,
            gp = gpar(fill = adjustcolor("#7AA6DC", alpha.f = 0.5), col = "grey30"))
grid.text(lengths(gene_sets), x = c(0.23, 0.5, 0.77), y = 0.48,
          gp = gpar(fontsize = 16, col = "black"))
grid.text(c("Desert > Inflamed", "Excluded > Inflamed"), x = c(0.23, 0.77), y = 0.84,
          gp = gpar(fontsize = 11, col = "black"))
popViewport()
dev.off()
