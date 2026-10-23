# Tumor differential expression, marker heatmap and PCA.
# Inputs: counts.csv and meta_cleaned_v5.xlsx in the working directory.
# Nominal P selects exploratory markers; BH FDR is retained in the full table.

rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(edgeR); library(limma)
  library(pheatmap); library(RColorBrewer); library(ggplot2); library(grid)
})

set.seed(2001)

# Analysis parameters.
keep_groups   <- c("tumor_desert","tumor_excluded","tumor_inflamed")
p_thresh      <- 0.05
topN_each     <- 40
zcap          <- 2
out_master_csv   <- "all_genes_DE_with_expr_and_fdr_tumor.csv"
out_selected_csv <- "selected_subtype_specific_genes_tumor.csv"
out_heatmap_pdf  <- "tumor_specific_modules_voom_modeled_heatmap.pdf"

# Read counts and align AOIs.
counts <- read.csv("counts.csv", row.names=1, check.names=FALSE)
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl before running this script.")
meta <- as.data.frame(readxl::read_xlsx("meta_cleaned_v5.xlsx"))
stopifnot(all(c("aoi_id","AOI_code","infiltration_category","Slide_number") %in% names(meta)))
if(!"group" %in% colnames(meta)){
  meta$group <- paste(meta$AOI_code, meta$infiltration_category, sep="_")
}
meta$aoi_id <- as.character(meta$aoi_id)

meta0 <- meta %>% filter(group %in% keep_groups)
mm <- function(x) gsub("[^A-Za-z0-9]", "", x)
stopifnot(!anyDuplicated(mm(meta$aoi_id)), !anyDuplicated(mm(colnames(counts))),
          !anyDuplicated(rownames(counts)), all(is.finite(as.matrix(counts))),
          all(as.matrix(counts) >= 0), !anyNA(meta$Slide_number))
map_idx <- match(mm(meta0$aoi_id), mm(colnames(counts)))
stopifnot(!anyNA(map_idx))
meta0$counts_colname <- colnames(counts)[map_idx]
meta2 <- meta0 %>% filter(!is.na(counts_colname))
stopifnot(nrow(meta2)>0)
meta2$group <- factor(meta2$group, levels=keep_groups)

counts2 <- counts[, meta2$counts_colname, drop=FALSE]
stopifnot(all(colnames(counts2) == meta2$counts_colname))

# TMM factors are estimated within this compartment from counts.csv.
# Voom and duplicateCorrelation model repeated AOIs using Slide_number as the block.
dge <- DGEList(counts=counts2)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~0 + meta2$group); colnames(design) <- levels(meta2$group)
slide_factor <- factor(meta2$Slide_number)

v0 <- voom(dge, design=design, plot=FALSE)
dupcor <- duplicateCorrelation(v0, design=design, block=slide_factor)
v  <- voom(dge, design=design, block=slide_factor, correlation=dupcor$consensus.correlation, plot=FALSE)
fit <- lmFit(v, design, block=slide_factor, correlation=dupcor$consensus.correlation)

# Pairwise subtype contrasts.
contr <- makeContrasts(
  desert_vs_excluded   = tumor_desert - tumor_excluded,
  desert_vs_inflamed   = tumor_desert - tumor_inflamed,
  excluded_vs_inflamed = tumor_excluded - tumor_inflamed,
  levels = design
)
contrast_names <- colnames(contr)
fit2 <- eBayes(contrasts.fit(fit, contr), trend=TRUE, robust=TRUE)

# Slide adjustment is used for visualization only.
expr_batch_removed <- removeBatchEffect(v$E, batch=slide_factor, design=design)
all_genes <- rownames(v$E)
master <- data.frame(Gene=all_genes, AveExpr=rowMeans(expr_batch_removed), stringsAsFactors=FALSE)

for(co in contrast_names){
  tt <- topTable(fit2, coef=co, number=Inf, sort.by="none")
  master[[paste0("logFC_", co)]] <- tt$logFC[match(master$Gene, rownames(tt))]
  master[[paste0("P.",    co)]]  <- tt$P.Value[match(master$Gene, rownames(tt))]
  master[[paste0("FDR.",  co)]]  <- tt$adj.P.Val[match(master$Gene, rownames(tt))]
}

avg_by_grp <- t(apply(expr_batch_removed, 1, function(x) tapply(x, meta2$group, mean, na.rm=TRUE)))
avg_by_grp <- as.data.frame(avg_by_grp); colnames(avg_by_grp) <- paste0("meanExpr_", colnames(avg_by_grp))
master <- cbind(master, avg_by_grp)
write.csv(master, out_master_csv, row.names=FALSE)

# Select genes by nominal P and consistent contrast direction.
specific_genes <- list(
  tumor_desert = master %>% filter(logFC_desert_vs_excluded > 0,  `P.desert_vs_excluded`  < p_thresh,
                                   logFC_desert_vs_inflamed > 0,  `P.desert_vs_inflamed`  < p_thresh) %>% pull(Gene),
  tumor_excluded = master %>% filter(logFC_excluded_vs_inflamed > 0, `P.excluded_vs_inflamed` < p_thresh,
                                     logFC_desert_vs_excluded  < 0, `P.desert_vs_excluded`   < p_thresh) %>% pull(Gene),
  tumor_inflamed = master %>% filter(logFC_desert_vs_inflamed   < 0, `P.desert_vs_inflamed`   < p_thresh,
                                     logFC_excluded_vs_inflamed < 0, `P.excluded_vs_inflamed` < p_thresh) %>% pull(Gene)
)
if(!is.null(topN_each) && is.numeric(topN_each)){
  specific_genes <- lapply(names(specific_genes), function(st){
    g <- unique(specific_genes[[st]]); if(length(g)==0) return(character(0))
    ref <- switch(st,
                  tumor_desert   = c("logFC_desert_vs_excluded","logFC_desert_vs_inflamed"),
                  tumor_excluded = c("logFC_excluded_vs_inflamed","logFC_desert_vs_excluded"),
                  tumor_inflamed = c("logFC_desert_vs_inflamed","logFC_excluded_vs_inflamed"))
    s <- rowMeans(abs(master[match(g, master$Gene), ref, drop=FALSE]), na.rm=TRUE)
    g[order(s, decreasing=TRUE)][seq_len(min(topN_each, length(g)))]
  })
  names(specific_genes) <- keep_groups
}


expr_all <- expr_batch_removed
heat_genes <- unique(unlist(specific_genes))
if(length(heat_genes)==0) stop("No genes meet the heatmap selection criteria.")

# Cluster rows within subtype modules.
ordered_genes <- c(); module_sizes <- c()
for(st in keep_groups){
  gs <- intersect(specific_genes[[st]], heat_genes)
  if(length(gs)==0) next
  mat <- expr_all[gs, , drop=FALSE]
  ord <- if(nrow(mat)>1){
    cm <- cor(t(mat), use="pairwise.complete.obs"); cm[is.na(cm)] <- 0
    hclust(as.dist(1-cm), method="complete")$order
  } else 1
  ordered_genes <- c(ordered_genes, gs[ord])
  module_sizes  <- c(module_sizes, length(gs))
}
gaps_row <- if(length(module_sizes)>1) head(cumsum(module_sizes), -1) else NULL

# Scale expression by gene and cap the display at +/-2.
heat_mat <- expr_all[ordered_genes, , drop=FALSE]
heat_z <- t(scale(t(heat_mat)))
heat_z[!is.finite(heat_z)] <- 0
heat_z[heat_z>zcap] <- zcap; heat_z[heat_z< -zcap] <- -zcap

# Cluster AOIs within subtype blocks.
samples_cols <- colnames(heat_z)
map_idx <- match(samples_cols, meta2$counts_colname)
meta_ord <- meta2[map_idx, , drop=FALSE]; rownames(meta_ord) <- samples_cols
col_order <- unlist(lapply(keep_groups, function(g){
  idx <- which(as.character(meta_ord$group)==g)
  if(length(idx)<=1) return(idx)
  cm <- cor(heat_z[,idx,drop=FALSE], use="pairwise.complete.obs"); cm[is.na(cm)] <- 0
  idx[hclust(as.dist(1-cm), method="complete")$order]
}))
heat_z   <- heat_z[, col_order, drop=FALSE]
meta_ord <- meta_ord[col_order, , drop=FALSE]
group_counts <- sapply(keep_groups, function(g) sum(meta_ord$group==g))
gaps_col <- if(length(group_counts)>1){ tmp <- cumsum(group_counts); tmp[-length(tmp)] } else NULL

# Heatmap annotations and palette.
annotation_col <- data.frame(
  slide_number = factor(meta_ord$Slide_number),
  subtype      = factor(meta_ord$group, levels=keep_groups)
)
rownames(annotation_col) <- colnames(heat_z)
pal <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(201)
brk <- seq(-zcap, zcap, length.out=length(pal)+1)
ann_colors <- list(
  slide_number = setNames(colorRampPalette(brewer.pal(8,"Set3"))(nlevels(annotation_col$slide_number)),
                          levels(annotation_col$slide_number)),
  subtype      = setNames(c("#000000","#0000FF","#FF0000"), keep_groups)
)

pdf(out_heatmap_pdf, width=12, height=max(6, 0.12*nrow(heat_z)+4))
pheatmap(heat_z,
         annotation_col=annotation_col, annotation_colors=ann_colors,
         cluster_rows=FALSE, cluster_cols=FALSE,
         gaps_row=gaps_row, gaps_col=gaps_col,
         show_rownames=TRUE,  show_colnames=FALSE,
         color=pal, breaks=brk, border_color=NA,
         main=sprintf("Tumor subtype-specific modules (voom+dupCor)\n%s (cap=%.1f)", "row_z", zcap),
         fontsize_row=if(nrow(heat_z)>80) 7 else 8
)
dev.off()
message("Heatmap saved: ", normalizePath(out_heatmap_pdf))

# Export selected genes and their contrast statistics.
out_tab <- data.frame(Gene = ordered_genes, stringsAsFactors = FALSE)

out_tab$Subtype <- sapply(out_tab$Gene, function(g){
  found <- sapply(names(specific_genes), function(st) g %in% specific_genes[[st]])
  if(any(found)) return(names(specific_genes)[which(found)[1]]) else return(NA_character_)
})

for(cn in contrast_names){
  pcol <- paste0("P.", cn); lcol <- paste0("logFC_", cn)
  if(pcol %in% colnames(master)) out_tab[[pcol]] <- master[[pcol]][match(out_tab$Gene, master$Gene)] else out_tab[[pcol]] <- NA
  if(lcol %in% colnames(master)) out_tab[[lcol]] <- master[[lcol]][match(out_tab$Gene, master$Gene)] else out_tab[[lcol]] <- NA
}
write.csv(out_tab, out_selected_csv, row.names = FALSE)
message("Selected gene table saved: ", out_selected_csv)

# PCA of the 1,000 most variable genes after slide adjustment.
vars <- matrixStats::rowVars(expr_batch_removed)
hvgs <- order(vars, decreasing = TRUE)[seq_len(min(1000, sum(vars > 0)))]
pca_z <- t(scale(t(expr_batch_removed[hvgs, , drop = FALSE])))
pca_z[!is.finite(pca_z)] <- 0
pc <- prcomp(t(pca_z), center = FALSE, scale. = FALSE)
explained <- pc$sdev^2 / sum(pc$sdev^2)
pca_df <- data.frame(aoi_id = meta2$aoi_id, subtype = meta2$infiltration_category,
                     patient = meta2$Patient_number, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
pca_df$subtype <- factor(pca_df$subtype, levels = c("desert", "excluded", "inflamed"))
p <- ggplot(pca_df, aes(PC1, PC2, colour = subtype)) +
  geom_point(size = 2.2) + stat_ellipse(type = "t", linewidth = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = c(desert = "#000000", excluded = "#0000FF", inflamed = "#FF0000")) +
  labs(title = "PCA: tumor", subtitle = "voom_slide_removed_HVG1000",
       x = sprintf("PC1 (%.1f%%)", 100 * explained[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * explained[2])) + theme_bw(base_size = 11)
ggsave("PCA_tumor_HVG1000.pdf", p, width = 5.2, height = 5.6)
write.csv(pca_df, "PCA_tumor_HVG1000_coordinates.csv", row.names = FALSE)
