# Export full-gene rankings for MSigDB GSEA.
# Run the three marker scripts first.
# Subtype scores use signed -log10(P); three-group scores average two contrasts.

out <- "enrichment_inputs"
dir.create(out, showWarnings = FALSE)
save_rank <- function(x, name) {
  x <- x[is.finite(x)]
  if (anyDuplicated(names(x))) {
    x <- x[order(abs(x), decreasing = TRUE)]
    x <- x[!duplicated(names(x))]
  }
  x <- sort(x, decreasing = TRUE)
  write.csv(data.frame(Gene = names(x), score = unname(x)), file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

signed_score <- function(lf, p) {
  s <- rep(0, length(lf))
  ok <- !(is.na(lf) | is.na(p) | p <= 0)
  s[ok] <- sign(lf[ok]) * -log10(p[ok])
  s
}

three_group_scores <- function(d) {
  de <- signed_score(d$logFC_desert_vs_excluded, d$P.desert_vs_excluded)
  di <- signed_score(d$logFC_desert_vs_inflamed, d$P.desert_vs_inflamed)
  ei <- signed_score(d$logFC_excluded_vs_inflamed, d$P.excluded_vs_inflamed)
  list(desert = setNames(rowMeans(cbind(de, di), na.rm = TRUE), d$Gene),
       excluded = setNames(rowMeans(cbind(ei, -de), na.rm = TRUE), d$Gene),
       inflamed = setNames(rowMeans(cbind(-di, -ei), na.rm = TRUE), d$Gene))
}

for (cell in c("tumor", "fibroblast")) {
  d <- read.csv(paste0("all_genes_DE_with_expr_and_fdr_", cell, ".csv"), check.names = FALSE)
  ranks <- three_group_scores(d)
  for (group in names(ranks)) save_rank(ranks[[group]], paste0("subtype_", cell, "_", group))
}

d <- read.csv("all_genes_DE_with_expr_and_fdr_immune.csv", check.names = FALSE)
ei <- signed_score(d$logFC_excluded_vs_inflamed, d$P.excluded_vs_inflamed)
save_rank(setNames(ei, d$Gene), "subtype_immune_excluded")
save_rank(setNames(-ei, d$Gene), "subtype_immune_inflamed")

# Composite tumor scores are ranking statistics, not independent Z tests.
d <- read.csv("all_genes_DE_with_expr_and_fdr_tumor.csv", check.names = FALSE)
pDI <- d$P.desert_vs_inflamed; pEI <- d$P.excluded_vs_inflamed
stopifnot(all(is.finite(pDI)), all(is.finite(pEI)), all(pDI > 0), all(pEI > 0))

zDI <- sign(d$logFC_desert_vs_inflamed) * qnorm(1 - pDI/2)
zEI <- sign(d$logFC_excluded_vs_inflamed) * qnorm(1 - pEI/2)
stopifnot(all(is.finite(zDI)), all(is.finite(zEI)))
zDI <- as.numeric(scale(zDI)); zEI <- as.numeric(scale(zEI))
save_rank(setNames((zDI + zEI)/sqrt(2), d$Gene), "composite_INT")
save_rank(setNames(zDI - pmax(zEI, 0), d$Gene), "composite_DES")
save_rank(setNames(zEI - pmax(zDI, 0), d$Gene), "composite_EXC")
save_rank(setNames(sign(d$logFC_desert_vs_inflamed) * -log10(pDI), d$Gene), "contrast_DI_all")
save_rank(setNames(sign(d$logFC_excluded_vs_inflamed) * -log10(pEI), d$Gene), "contrast_EI_all")
