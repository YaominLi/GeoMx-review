# Full-gene MSigDB enrichment using the retained ranking definitions.
# Checkpoints include the RNG state so a resumed run follows the same sequence.
suppressPackageStartupMessages({library(clusterProfiler); library(BiocParallel); library(ggplot2)})
stopifnot(as.character(packageVersion("fgsea"))=="1.34.0",
          as.character(packageVersion("clusterProfiler"))=="4.16.0")
out <- "gsea_results"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
writeLines("RUNNING", file.path(out, "RUN_STATUS.txt"))
collections <- c(H="h",C2="c2",C5="c5",C7="c7")
t2g <- lapply(collections,function(x) {
  d <- clusterProfiler::read.gmt(file.path("gene_sets",paste0(x,".all.v2025.1.Hs.symbols.gmt")))
  n <- table(d$term); d[d$term %in% names(n)[n>=10 & n<=500],,drop=FALSE]
})
bp <- if (.Platform$OS.type=="windows") SnowParam(4,type="SOCK",progressbar=FALSE) else MulticoreParam(4,progressbar=FALSE)
register(bp,default=TRUE)
families <- list(
  tumor=c("subtype_tumor_desert","subtype_tumor_excluded","subtype_tumor_inflamed"),
  fibroblast=c("subtype_fibroblast_desert","subtype_fibroblast_excluded","subtype_fibroblast_inflamed"),
  immune=c("subtype_immune_excluded","subtype_immune_inflamed"),
  intersection=c("composite_INT","composite_DES","composite_EXC","contrast_DI_all","contrast_EI_all"))
summary <- list(); all <- list()
for (family in names(families)) {
  set.seed(2001)
  for (rank_name in families[[family]]) {
    d <- read.csv(file.path("enrichment_inputs",paste0(rank_name,".csv")))
    ranked <- setNames(d$score,d$Gene)
    stopifnot(!anyDuplicated(names(ranked)),all(is.finite(ranked)),all(diff(ranked)<=0))
    od <- file.path(out,rank_name); dir.create(od,showWarnings=FALSE)
    for (coll in names(collections)) {
      message(format(Sys.time())," ",rank_name," / ",coll)
      cache <- file.path(od,paste0(coll,"_checkpoint.rds"))
      context <- list(rank_md5=unname(tools::md5sum(file.path("enrichment_inputs",paste0(rank_name,".csv")))),
                      gmt_md5=unname(tools::md5sum(file.path("gene_sets",paste0(collections[[coll]],".all.v2025.1.Hs.symbols.gmt")))),
                      fgsea=as.character(packageVersion("fgsea")),clusterProfiler=as.character(packageVersion("clusterProfiler")))
      if (file.exists(cache)) {
        z <- readRDS(cache)
        if (!identical(z$context,context)) stop("Checkpoint inputs/software changed: ",cache)
        result <- z$table; assign(".Random.seed",z$rng,envir=.GlobalEnv)
        writeLines(unique(z$warnings),file.path(od,paste0(coll,"_warnings.txt")))
      } else {
        warnings <- character()
        result_obj <- withCallingHandlers(clusterProfiler::GSEA(
          geneList=ranked,TERM2GENE=t2g[[coll]],minGSSize=10,maxGSSize=500,
          pvalueCutoff=1,pAdjustMethod="BH",verbose=FALSE,seed=TRUE,BPPARAM=bp),
          warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
        result <- as.data.frame(result_obj)
        saveRDS(list(table=result,rng=.Random.seed,warnings=warnings,context=context),cache)
        writeLines(unique(warnings),file.path(od,paste0(coll,"_warnings.txt")))
      }
      stopifnot(nrow(result)>0, all(is.finite(result$NES)))
      write.csv(result,file.path(od,paste0(coll,"_full.csv")),row.names=FALSE)
      summary[[length(summary)+1]] <- data.frame(rank=rank_name,collection=coll,n_ranked=length(ranked),n_terms=nrow(result),FDR05=sum(result$p.adjust<.05,na.rm=TRUE),nominal05=sum(result$pvalue<.05,na.rm=TRUE))
      result$Rank <- rank_name; result$Collection <- coll
      all[[length(all)+1]] <- result
      write.csv(do.call(rbind,summary),file.path(out,"run_summary.csv"),row.names=FALSE)
    }
  }
}
long <- do.call(rbind,all)
# Match the schema expected by the original selected-pathway heatmap script.
for (cell in c("tumor","fibroblast","immune")) {
  ranks <- families[[cell]]; x <- long[long$Rank %in% ranks,]
  base <- unique(data.frame(PathwayKey=paste(x$Collection,x$ID,sep="::"),Collection=x$Collection,PathwayID=x$ID,Description=x$Description))
  for (rank in ranks) {
    d <- x[x$Rank==rank,]; idx <- match(base$PathwayKey,paste(d$Collection,d$ID,sep="::"))
    nm <- sub("^subtype_","",rank)
    for (v in c(NES="NES",pvalue="pvalue",padj="p.adjust",genes="core_enrichment")) {
      suffix <- names(c(NES="NES",pvalue="pvalue",padj="p.adjust",genes="core_enrichment"))[match(v,c("NES","pvalue","p.adjust","core_enrichment"))]
      base[[paste0(nm,"_",suffix)]] <- d[[v]][idx]
    }
  }
  write.csv(base,paste0("GSEA_MSigDB_all_subtypes_combined_",cell,".csv"),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

# Composite ranks for the selected NES bars.
three <- long[long$Rank %in% c("composite_INT", "composite_DES", "composite_EXC"), ]
three$SubtypeShort <- sub("^composite_", "", three$Rank)
write.csv(three, "gsea_results/GSEA_three_sets_long.csv", row.names = FALSE)
writeLines("COMPLETED: 13 ranking definitions x 4 GMT collections.",file.path(out,"RUN_STATUS.txt"))
