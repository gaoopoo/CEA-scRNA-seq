# Project root: the folder containing SMC/, macro/, Fig3/, Fig5_CellChat/, etc.
# Every path in this script is built from this one setting, so the code no longer
# depends on a single hard-coded drive letter. Either start R/RStudio in the
# project root, or edit the one line below to point at it directly.
# Stored as an R option so it survives the rm(list = ls()) calls between sections.
options(cea.project_root = normalizePath(".", mustWork = FALSE))

#-----Figure 5. TREM2??? mature lipid-associated macrophages (LAMs) occupy a late repair-associated state along the monocyte-to-macrophage continuum
rm(list = ls())
library(Seurat)
library(tidyverse)
library(patchwork)

setwd(paste0(getOption("cea.project_root"), "/macro"))
macro <- readRDS(file = "macro.rds")
-----------------------------------------------------------------------------
  # ============================================================================
# 共用准备代码：原路径定位与 00_common.R 完整原文
# ==============================================================================


# ==============================================================================
source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE))
  source_path <- rstudioapi::getActiveDocumentContext()$path
if (!nzchar(source_path)) stop("Use source('.../01_donor_panels.R', encoding='UTF-8').")
FIG4_HOME <- dirname(source_path)

# 
# Shared settings/helpers. Source from 01/02/03; do not modify original objects.
options(stringsAsFactors = FALSE)
if (dir.exists(file.path(FIG4_HOME,"R-library"))) .libPaths(c(file.path(FIG4_HOME,"R-library"),.libPaths()))
cfg <- list(
  input = paste0(getOption("cea.project_root"), "/macro/macro.rds"),
  existing_pb = paste0(getOption("cea.project_root"), "/macro/Fig4D/Fig4D_donor_pseudobulk_results"),
  old_trajectory = paste0(getOption("cea.project_root"), "/macro/Fig4E/macro_with_pseudotime_3grp.rds"),
  populations = c("Monocytes", "Macrophages1", "Macrophages2"),
  groups = c("ASYM", "SYM"), min_cells = 20L, seed = 20260901L,
  bootstrap = 5000L, bins = 15L, bin_cells = 10L, bin_donors = 4L,
  trajectory_seeds = c(20260901L, 20260902L),
  # Not a known technical-batch variable. Do not regress out clinical group.
  alignment_group = NULL
)
programs <- list(
  Inflammatory = c("IL1B","TNF","IL6","NLRP3","CCL3","CXCL2","CCL4",
                   "S100A8","S100A9","FCN1","VCAN"),
  `Lipid/LAM-associated` = c("TREM2","APOE","CD9","SPP1","GPNMB","LPL",
                             "FABP5","LIPA","CTSD","CTSB"),
  `Efferocytosis-associated` = c("MERTK","CD163","C1QA"),
  `Antigen presentation` = c("CIITA","HLA-DRA","CD74"),
  `Legacy LAM (sensitivity)` = c("TREM2","APOE","CD9","SPP1","GPNMB","LPL",
                                 "FABP5","LIPA","CTSD","CTSB","MERTK","CD163")
)
primary_programs <- names(programs)[1:4]
focused_genes <- c("TREM2","APOE","CD9","LIPA","CTSD","MERTK","CD163","C1QA","C1QB")
group_colors <- c(ASYM = "#F9C80E", SYM = "#C71000")
pop_colors <- c(Monocytes = "#F9AB3C", Macrophages1 = "#CA480D", Macrophages2 = "#C86C18")

need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages, then restart: ", paste(missing, collapse=", "))
}
need(c("Matrix","SeuratObject","edgeR","dplyr","tidyr","ggplot2","patchwork","ragg","svglite"))
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
theme_fig <- theme_classic(base_size=10, base_family="sans") + theme(
  text=element_text(family="sans"), axis.text=element_text(color="black"),
  plot.title=element_text(face="bold",size=11), plot.subtitle=element_text(size=8,color="grey30"),
  strip.background=element_rect(fill="grey95",color="grey75"), strip.text=element_text(face="bold"))
program_labeller <- labeller(program=c("Inflammatory"="Inflammatory",
                                       "Lipid/LAM-associated"="Lipid /\nLAM", "Efferocytosis-associated"="Efferocytosis\nassociated",
                                       "Antigen presentation"="Antigen\npresentation"))

new_run <- function(stage) {
  parent <- file.path(FIG4_HOME,"results")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)
  out <- file.path(parent,paste0(stage,"_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if (dir.exists(out)) stop("Run folder already exists; wait one second and rerun.")
  dir.create(out)
  for (x in c("figures","tables","objects")) dir.create(file.path(out,x))
  saveRDS(cfg,file.path(out,"settings.rds"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo_start.txt"))
  out
}
csv <- function(x,name,out) write.csv(x,file.path(out,"tables",paste0(name,".csv")),row.names=FALSE,na="",fileEncoding="UTF-8")
save_fig <- function(p,name,out,w=7.2,h=4.5) {
  # Never use an unregistered font in base pdf(). Cairo/SVG and raster agree.
  if (requireNamespace("showtext",quietly=TRUE)) showtext::showtext_auto(FALSE)
  if (!capabilities("cairo")) stop("Cairo unavailable. Use SVG/PNG or install a Cairo-enabled R build.")
  p <- p & theme(text=element_text(family="sans"))
  # Use short relative device paths: some Windows graphics devices truncate long
  # UTF-8 absolute paths even when the path is shorter than 260 Unicode characters.
  previous_wd <- getwd()
  setwd(file.path(out,"figures"))
  on.exit(setwd(previous_wd),add=TRUE)
  base <- name
  ggsave(paste0(base,".pdf"),p,width=w,height=h,device=grDevices::cairo_pdf,bg="white")
  ggsave(paste0(base,".svg"),p,width=w,height=h,device=svglite::svglite,bg="white")
  ggsave(paste0(base,".png"),p,width=w,height=h,dpi=300,device=ragg::agg_png,bg="white")
  paths <- paste0(base,c(".pdf",".svg",".png"))
  if(any(!file.exists(paths)) || any(file.info(paths)$size<100)) stop("Incomplete figure export: ",name)
  saveRDS(p,file.path(out,"objects",paste0(name,"_plot.rds")))
}
read_source <- function() {
  x <- readRDS(cfg$input)
  md <- x@meta.data
  stopifnot(all(c("orig.ident","sample","group","celltype") %in% names(md)))
  a <- x@assays[["RNA"]]
  if (inherits(a,"Assay5")) {
    layers <- SeuratObject::Layers(a,search="^counts")
    if (length(layers)!=1L) stop("RNA has multiple count layers. JoinLayers on a COPY, then specify that input.")
    counts <- SeuratObject::LayerData(a,layer=layers)
  } else counts <- methods::slot(a,"counts")
  stopifnot(nrow(counts)>0L,ncol(counts)>0L,!anyDuplicated(rownames(counts)),
            !anyDuplicated(colnames(counts)),all(colnames(counts)%in%rownames(md)))
  values <- if (inherits(counts,"sparseMatrix")) counts@x else as.vector(counts)
  if (any(!is.finite(values)|values<0|abs(values-round(values))>1e-6)) stop("RNA counts are not raw nonnegative integers.")
  md <- md[colnames(counts),,drop=FALSE]
  md$cell_id <- rownames(md); md$donor <- as.character(md$orig.ident)
  md$sample_label <- as.character(md$sample); md$group <- as.character(md$group)
  md$population <- as.character(md$celltype)
  lookup <- distinct(md,donor,sample_label,group)
  if(anyNA(lookup) || anyDuplicated(lookup$donor) || anyDuplicated(lookup$sample_label))
    stop("Donor/sample/group is missing or not one-to-one. Resolve metadata before continuing.")
  list(counts=counts,md=md)
}
scope_table <- function(md,out) {
  tab <- md %>% count(population,group,sample_label,name="n_cells") %>%
    mutate(in_focused_analysis=population%in%cfg$populations & group%in%cfg$groups,
           reason=case_when(!group%in%cfg$groups ~ "Outside ASYM/SYM donor comparison",
                            !population%in%cfg$populations ~ "Outside focused three-population scope; not proof of contamination",
                            TRUE ~ "Included before per-donor cell threshold"))
  csv(tab,"scope_all_cells_in_input",out)
  writeLines(c("Focused analysis: Monocytes, Macrophages1, Macrophages2 only.",
               "Macrophages3 is not deleted or reclassified; annotation/QC remains unresolved.",
               "The full atlas (including legacy Neutrophils/Macrophages4) is not reannotated here.",
               "This figure cannot establish completeness of the monocyte/macrophage compartment.",
               "Absence of a population from macro.rds is not evidence of its absence in the full atlas."),file.path(out,"SCOPE.txt"))
}
make_pb <- function(counts,md) {
  use <- md$population%in%cfg$populations & md$group%in%cfg$groups
  m <- md[use,,drop=FALSE]
  m$pb_id <- paste(m$donor,m$population,sep="___")
  ids <- unique(m$pb_id)
  aggregation <- Matrix::sparseMatrix(i=seq_len(nrow(m)),j=match(m$pb_id,ids),x=1,
                                      dims=c(nrow(m),length(ids)),dimnames=list(m$cell_id,ids))
  cnt <- counts[,m$cell_id,drop=FALSE] %*% aggregation
  meta <- m %>% count(pb_id,donor,sample_label,group,population,name="n_cells")
  meta <- meta[match(ids,meta$pb_id),]
  meta$keep <- meta$n_cells>=cfg$min_cells
  list(counts=cnt[,meta$keep,drop=FALSE],meta=as.data.frame(meta[meta$keep,]),qc=meta)
}
row_z <- function(m) {
  s <- apply(m,1,sd,na.rm=TRUE)
  ans <- sweep(sweep(m,1,rowMeans(m,na.rm=TRUE),"-"),1,s,"/")
  ans[!is.finite(ans)] <- NA_real_
  ans
}
long_matrix <- function(m,value="value") {
  d <- as.data.frame(as.table(m),stringsAsFactors=FALSE)
  names(d) <- c("gene","pb_id",value); d
}
paired_result <- function(d,a,b,value="score") {
  a0 <- d[d$population==a,c("donor",value)]; b0 <- d[d$population==b,c("donor",value)]
  names(a0)[2] <- "a"; names(b0)[2] <- "b"
  z <- merge(a0,b0,by="donor"); z <- z[is.finite(z$a)&is.finite(z$b),]
  n <- nrow(z)
  if(n<3L) return(data.frame(contrast=paste(b,"-",a),n_pairs=n,p=NA,mean_difference=NA,
                             ci_low=NA,ci_high=NA,method="Insufficient donors"))
  delta <- z$b-z$a
  # Exact conditional signed-rank enumeration, including ties; zero differences omitted.
  nz <- delta[delta!=0]; ranks <- rank(abs(nz),ties.method="average")
  if(length(nz)==0L) p <- 1 else {
    if(length(nz)>20L) stop("Exact enumeration limited to 20 donor pairs.")
    signs <- as.matrix(expand.grid(rep(list(c(-1,1)),length(nz))))
    observed <- abs(sum(sign(nz)*ranks))
    p <- mean(abs(drop(signs%*%ranks)) >= observed-1e-10)
  }
  draws <- replicate(cfg$bootstrap,mean(sample(delta,n,replace=TRUE)))
  ci <- quantile(draws,c(.025,.975),names=FALSE)
  data.frame(contrast=paste(b,"-",a),n_pairs=n,p=p,mean_difference=mean(delta),
             ci_low=ci[1],ci_high=ci[2],method="Two-sided exact signed-rank sign enumeration; symmetry assumption")
}
latest_run <- function(stage) {
  dirs <- list.dirs(file.path(FIG4_HOME,"results"),recursive=FALSE,full.names=TRUE)
  dirs <- dirs[startsWith(basename(dirs),paste0(stage,"_")) & file.exists(file.path(dirs,"COMPLETE.txt"))]
  if (!length(dirs)) stop("Run and complete stage ",stage," first.")
  sort(dirs,decreasing=TRUE)[1]
}

----------------------------------------------------------
# ==========================================================================================
# ---------FigA | Donor-balanced functional-marker heatmap

local({
  source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
  if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE))
    source_path <- rstudioapi::getActiveDocumentContext()$path
  if (!nzchar(source_path)) stop("Use source('.../01_donor_panels.R', encoding='UTF-8').")
  FIG4_HOME <- dirname(source_path)
  source(file.path(FIG4_HOME,"00_common.R"),local=TRUE,encoding="UTF-8")
  out <- new_run("01_donor"); set.seed(cfg$seed)
  message("Read-only source: ",cfg$input)
  dat <- read_source(); scope_table(dat$md,out)
  pb <- make_pb(dat$counts,dat$md); csv(pb$qc,"pseudobulk_cell_count_QC",out)
  rm(dat); gc()
  dge <- edgeR::DGEList(as.matrix(pb$counts[Matrix::rowSums(pb$counts)>0,,drop=FALSE]))
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  lcpm <- edgeR::cpm(dge,log=TRUE,prior.count=2)
  meta <- pb$meta; meta$population <- factor(meta$population,levels=cfg$populations)
  meta$group <- factor(meta$group,levels=cfg$groups)
  csv(cbind(meta,library_size=dge$samples$lib.size,TMM_factor=dge$samples$norm.factors),"pseudobulk_metadata",out)
  membership <- bind_rows(lapply(names(programs),function(nm) data.frame(program=nm,gene=programs[[nm]])))
  membership$origin <- "Investigator-defined; partitioned from existing Fig4 marker panels; not a functional assay"
  membership$in_RNA <- membership$gene%in%rownames(lcpm)
  sds <- apply(lcpm,1,sd)
  membership$nonconstant <- unname(sds[membership$gene]>0)
  csv(membership,"program_definitions_and_QC",out)
  z <- row_z(lcpm)
  scores <- bind_rows(lapply(names(programs),function(nm) {
    genes <- intersect(programs[[nm]],rownames(z)); genes <- genes[is.finite(z[genes,1])]
    if (length(genes)<max(2,ceiling(.7*length(programs[[nm]]))))
      stop("Low gene-set coverage: ",nm,". Inspect program_definitions_and_QC.csv.")
    cbind(meta,program=nm,score=colMeans(z[genes,,drop=FALSE]),n_genes=length(genes))
  }))
  csv(scores,"all_donor_program_scores",out)
  stats <- bind_rows(lapply(names(programs),function(nm) {
    d <- scores[scores$program==nm,]
    bind_rows(lapply(combn(cfg$populations,2,simplify=FALSE),function(pair)
      cbind(program=nm,paired_result(d,pair[1],pair[2]))))
  }))
  stats$p_Holm_15 <- p.adjust(stats$p,"holm")
  stats$multiplicity_family <- "All five fixed programs x all three population pairs (15 tests)"
  csv(stats,"paired_program_statistics",out)
  # A: equal-donor expression summary. No cell-level test; row-z uses all retained PBs.
  genes_main <- membership %>% filter(program%in%primary_programs,in_RNA,nonconstant) %>% select(program,gene)
  hm <- long_matrix(z[unique(genes_main$gene),,drop=FALSE],"gene_z") %>%
    left_join(meta,by="pb_id") %>% left_join(genes_main,by="gene")
  csv(hm,"A_marker_donor_source",out)
  hmean <- hm %>% group_by(program,gene,population) %>%
    summarise(mean_z=mean(gene_z),n_donors=n_distinct(donor),.groups="drop")
  hmean$gene <- factor(hmean$gene,levels=unique(genes_main$gene))
  hmean$population <- factor(hmean$population,levels=rev(cfg$populations))
  hmean$program <- factor(hmean$program,levels=primary_programs)
  limit <- max(abs(hmean$mean_z)); limit <- max(1,ceiling(limit*2)/2)
  pA <- ggplot(hmean,aes(gene,population,fill=mean_z)) + geom_tile(color="white",linewidth=.3) +
    facet_grid(~program,scales="free_x",space="free_x",labeller=labeller(program=c(
      "Inflammatory"="Inflammatory", "Lipid/LAM-associated"="Lipid/LAM-associated",
      "Efferocytosis-associated"="Efferocytosis\nassociated", "Antigen presentation"="Antigen\npresentation"))) +
    scale_fill_gradient2(low="#326A9C",mid="white",high="#B62932",limits=c(-limit,limit),name="Mean gene z") +
    labs(x=NULL,y=NULL,title="Marker programs in the focused three-population analysis",
         subtitle="RNA-count pseudobulks; equal weight per donor; not a lineage or functional assay") +
    theme_fig + theme(axis.text.x=element_text(angle=60,hjust=1,face="italic",size=8),axis.ticks=element_blank())
  save_fig(pA,"Fig4A_marker_program_heatmap",out,11,3.5)
})

---------------------------------------------------------------------------------
# ==========================================================================================
# Figure B : GSEA: pooled Macrophages1/Macrophages2 versus Monocytes
rm(list = ls())
gc()

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(fgsea)
  library(ggplot2)
  library(patchwork)
})

macro <- readRDS(paste0(getOption("cea.project_root"), "/macro/macro.rds"))
out_dir <- paste0(getOption("cea.project_root"), "/macro/Fig4C_v3/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. 准备 subset ──────────────────────────────────────────────────
macro$lam_vs_mono <- case_when(
  macro$celltype %in% c("Macrophages1", "Macrophages2") ~ "LAM_lineage",
  macro$celltype == "Monocytes" ~ "Monocytes",
  TRUE ~ NA_character_
)

macro_sub <- subset(macro, subset = !is.na(lam_vs_mono))
Idents(macro_sub) <- macro_sub$lam_vs_mono
cat("Subset cells:\n"); print(table(Idents(macro_sub)))

# ── 2. 全基因组 DE（用于构建 ranked list）───────────────────────────
cat("\nRunning FindMarkers on all genes (this takes ~2-5 min)...\n")
de_all <- FindMarkers(macro_sub,
                      ident.1 = "LAM_lineage",
                      ident.2 = "Monocytes",
                      logfc.threshold = 0,
                      min.pct = 0.01,
                      test.use = "wilcox") %>%
  rownames_to_column("gene") %>%
  filter(!is.na(avg_log2FC), is.finite(avg_log2FC))

cat("Total genes in ranking:", nrow(de_all), "\n")

# Ranked gene list
ranks <- de_all$avg_log2FC
names(ranks) <- de_all$gene
ranks <- sort(ranks, decreasing = TRUE)

# 简单诊断
cat("Top 10 ranked genes (LAM lineage up):\n")
print(head(ranks, 10))
cat("\nBottom 10 ranked genes (Monocytes up):\n")
print(tail(ranks, 10))

# ── 3. 3 个发表的 TREM2??? LAM signature ──────────────────────────────
LAM_signatures <- list(
  
  Cochain_2018 = c("TREM2","CD9","SPP1","GPNMB","LPL","FABP5",
                   "FABP4","PLIN2","CTSD","CTSB","CTSL","CTSZ",
                   "ACP5","MMP9","APOE","APOC1","LIPA","CD63",
                   "PSAP","SLC11A1","ANXA1"),
  
  Patterson_2023 = c("TREM2","APOE","APOC1","SPP1","GPNMB","LPL",
                     "FABP5","FABP4","PLIN2","LIPA","CTSD","CTSB",
                     "CTSL","CD9","CD63","CD68","LGALS3","MMP9",
                     "ITGAX","CHIT1","CHI3L1"),
  
  Willemsen_2024 = c("TREM2","CD9","APOE","SPP1","GPNMB","LPL",
                     "FABP5","LIPA","CTSD","CTSB","MERTK","CD163",
                     "PLIN2","FABP4","CD36","MRC1","STAB1","C1QA",
                     "C1QB","C1QC")
)

# 过滤到 ranked list 里存在的基因
LAM_signatures <- lapply(LAM_signatures, function(g) intersect(g, names(ranks)))

cat("\nSignature coverage:\n")
for (sig_name in names(LAM_signatures)) {
  cat("  ", sig_name, ":", length(LAM_signatures[[sig_name]]), "genes\n")
}

# ── 4. 跑 fgsea ─────────────────────────────────────────────────────
set.seed(42)
gsea_result <- fgseaMultilevel(
  pathways = LAM_signatures,
  stats    = ranks,
  minSize  = 5,
  maxSize  = 200,
  eps      = 0
) %>%
  as_tibble() %>%
  arrange(desc(NES))

cat("\n=== GSEA Results ===\n")
print(gsea_result %>% select(pathway, NES, pval, padj, size))

# 保存结果
write.csv(gsea_result %>%
            select(pathway, NES, pval, padj, size, leadingEdge) %>%
            mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")),
          paste0(out_dir, "GSEA_results.csv"),
          row.names = FALSE)

# ── 5. 美化版 GSEA enrichment plot ──────────────────────────────────
plot_gsea_curve <- function(pathway_name, ranks, pathway_genes, gsea_row) {
  
  ranks_sorted <- sort(ranks, decreasing = TRUE)
  hits <- names(ranks_sorted) %in% pathway_genes
  
  N <- length(ranks_sorted)
  N_hit <- sum(hits)
  N_miss <- N - N_hit
  
  weight <- abs(ranks_sorted) * hits
  hit_score <- weight / sum(weight)
  miss_score <- (!hits) / N_miss
  running_ES <- cumsum(hit_score - miss_score)
  
  curve_color <- ifelse(gsea_row$NES > 0, "#C04040", "#3A6FA0")
  
  df_curve <- data.frame(rank = seq_along(running_ES), ES = running_ES)
  
  # ES 曲线
  p1 <- ggplot(df_curve, aes(x = rank, y = ES)) +
    geom_line(color = curve_color, linewidth = 1.2) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    annotate("text",
             x = N * 0.55,
             y = max(running_ES) * 0.85,
             label = sprintf("NES = %.2f\nP = %s\nP.adj = %s",
                             gsea_row$NES,
                             format.pval(gsea_row$pval, digits = 2, eps = 1e-50),
                             format.pval(gsea_row$padj, digits = 2, eps = 1e-50)),
             hjust = 0, size = 3.6, fontface = "bold",
             color = curve_color) +
    labs(title = gsub("_", " ", pathway_name),
         y = "Enrichment score", x = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.text = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.4),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # 命中位置 ticks
  hit_ranks <- which(hits)
  df_hits <- data.frame(rank = hit_ranks)
  
  p2 <- ggplot(df_hits, aes(x = rank)) +
    geom_segment(aes(xend = rank, y = 0, yend = 1),
                 color = curve_color, linewidth = 0.4) +
    xlim(0, N) +
    labs(x = "Gene rank (← LAM lineage up | Monocytes up →)") +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth = 0.4),
          axis.text.x = element_text(color = "black", size = 8),
          plot.margin = margin(0, 0, 5, 0))
  
  p1 / p2 + plot_layout(heights = c(4, 0.6))
}

# 生成 3 张
gsea_plots <- list()
for (i in seq_len(nrow(gsea_result))) {
  sig_name <- gsea_result$pathway[i]
  gsea_plots[[sig_name]] <- plot_gsea_curve(
    pathway_name  = sig_name,
    ranks         = ranks,
    pathway_genes = LAM_signatures[[sig_name]],
    gsea_row      = gsea_result[i, ]
  )
}

# 拼接
p_C <- gsea_plots[[1]] | gsea_plots[[2]] | gsea_plots[[3]]
p_C <- p_C + plot_annotation(
  title = "Figure 4C. TREM2??? LAM lineage validated by three independent published signatures",
  subtitle = "GSEA on genes ranked by LAM lineage (foamy + mature) vs Monocytes",
  theme = theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )
)

ggsave(paste0(out_dir, "Fig4C_v2.pdf"), p_C,
       width = 15, height = 5.5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig4C_v2.png"), p_C,
       width = 15, height = 5.5, dpi = 600)

cat("\n??? Fig 4C v2 saved to:", out_dir, "\n")

--------------------------------------------------------------------------------
  
# ==========================================================================================
# Figure C : Inflammatory and Lipid/LAM-associated programs

local({
  source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
  if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE))
    source_path <- rstudioapi::getActiveDocumentContext()$path
  if (!nzchar(source_path)) stop("Use source('.../01_donor_panels.R', encoding='UTF-8').")
  FIG4_HOME <- dirname(source_path)
  source(file.path(FIG4_HOME,"00_common.R"),local=TRUE,encoding="UTF-8")
  out <- new_run("01_donor"); set.seed(cfg$seed)
  message("Read-only source: ",cfg$input)
  dat <- read_source(); scope_table(dat$md,out)
  pb <- make_pb(dat$counts,dat$md); csv(pb$qc,"pseudobulk_cell_count_QC",out)
  rm(dat); gc()
  dge <- edgeR::DGEList(as.matrix(pb$counts[Matrix::rowSums(pb$counts)>0,,drop=FALSE]))
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  lcpm <- edgeR::cpm(dge,log=TRUE,prior.count=2)
  meta <- pb$meta; meta$population <- factor(meta$population,levels=cfg$populations)
  meta$group <- factor(meta$group,levels=cfg$groups)
  csv(cbind(meta,library_size=dge$samples$lib.size,TMM_factor=dge$samples$norm.factors),"pseudobulk_metadata",out)
  membership <- bind_rows(lapply(names(programs),function(nm) data.frame(program=nm,gene=programs[[nm]])))
  membership$origin <- "Investigator-defined; partitioned from existing Fig4 marker panels; not a functional assay"
  membership$in_RNA <- membership$gene%in%rownames(lcpm)
  sds <- apply(lcpm,1,sd)
  membership$nonconstant <- unname(sds[membership$gene]>0)
  csv(membership,"program_definitions_and_QC",out)
  z <- row_z(lcpm)
  scores <- bind_rows(lapply(names(programs),function(nm) {
    genes <- intersect(programs[[nm]],rownames(z)); genes <- genes[is.finite(z[genes,1])]
    if (length(genes)<max(2,ceiling(.7*length(programs[[nm]]))))
      stop("Low gene-set coverage: ",nm,". Inspect program_definitions_and_QC.csv.")
    cbind(meta,program=nm,score=colMeans(z[genes,,drop=FALSE]),n_genes=length(genes))
  }))
  csv(scores,"all_donor_program_scores",out)
  stats <- bind_rows(lapply(names(programs),function(nm) {
    d <- scores[scores$program==nm,]
    bind_rows(lapply(combn(cfg$populations,2,simplify=FALSE),function(pair)
      cbind(program=nm,paired_result(d,pair[1],pair[2]))))
  }))
  stats$p_Holm_15 <- p.adjust(stats$p,"holm")
  stats$multiplicity_family <- "All five fixed programs x all three population pairs (15 tests)"
  csv(stats,"paired_program_statistics",out)
  
  # B: donor pairing is shown as grey lines, not a group-average developmental arrow.
  ds <- scores %>% filter(program%in%primary_programs)
  ds$program <- factor(ds$program,levels=primary_programs)
  pB <- ggplot(ds,aes(population,score)) +
    geom_line(aes(group=donor),color="grey75",linewidth=.35) +
    geom_boxplot(aes(fill=population),outlier.shape=NA,width=.5,alpha=.25) +
    geom_point(aes(fill=group),shape=21,size=1.9,stroke=.3) +
    facet_wrap(~program,ncol=2,scales="free_y") +
    scale_fill_manual(values=c(pop_colors,group_colors),breaks=cfg$groups,name="Clinical group") +
    labs(x=NULL,y="Mean gene z (donor pseudobulk)",title="Donor-level program comparison",
         subtitle=paste0("One point per donor/population; n = ",n_distinct(meta$donor)," donors; paired lines; all 15 adjusted tests in source table")) +
    theme_fig + theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="bottom")
  save_fig(pB,"Fig4B_donor_programs",out,8,6)
})

--------------------------------------------------------------------------------
# ==========================================================================================
# ------Figure D :Efferocytosis and antigen-presentation programs
local({
  source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
  if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE))
    source_path <- rstudioapi::getActiveDocumentContext()$path
  if (!nzchar(source_path)) stop("Use source('.../01_donor_panels.R', encoding='UTF-8').")
  FIG4_HOME <- dirname(source_path)
  source(file.path(FIG4_HOME,"00_common.R"),local=TRUE,encoding="UTF-8")
  out <- new_run("01_donor"); set.seed(cfg$seed)
  message("Read-only source: ",cfg$input)
  dat <- read_source(); scope_table(dat$md,out)
  pb <- make_pb(dat$counts,dat$md); csv(pb$qc,"pseudobulk_cell_count_QC",out)
  rm(dat); gc()
  dge <- edgeR::DGEList(as.matrix(pb$counts[Matrix::rowSums(pb$counts)>0,,drop=FALSE]))
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  lcpm <- edgeR::cpm(dge,log=TRUE,prior.count=2)
  meta <- pb$meta; meta$population <- factor(meta$population,levels=cfg$populations)
  meta$group <- factor(meta$group,levels=cfg$groups)
  csv(cbind(meta,library_size=dge$samples$lib.size,TMM_factor=dge$samples$norm.factors),"pseudobulk_metadata",out)
  membership <- bind_rows(lapply(names(programs),function(nm) data.frame(program=nm,gene=programs[[nm]])))
  membership$origin <- "Investigator-defined; partitioned from existing Fig4 marker panels; not a functional assay"
  membership$in_RNA <- membership$gene%in%rownames(lcpm)
  sds <- apply(lcpm,1,sd)
  membership$nonconstant <- unname(sds[membership$gene]>0)
  csv(membership,"program_definitions_and_QC",out)
  z <- row_z(lcpm)
  scores <- bind_rows(lapply(names(programs),function(nm) {
    genes <- intersect(programs[[nm]],rownames(z)); genes <- genes[is.finite(z[genes,1])]
    if (length(genes)<max(2,ceiling(.7*length(programs[[nm]]))))
      stop("Low gene-set coverage: ",nm,". Inspect program_definitions_and_QC.csv.")
    cbind(meta,program=nm,score=colMeans(z[genes,,drop=FALSE]),n_genes=length(genes))
  }))
  csv(scores,"all_donor_program_scores",out)
  stats <- bind_rows(lapply(names(programs),function(nm) {
    d <- scores[scores$program==nm,]
    bind_rows(lapply(combn(cfg$populations,2,simplify=FALSE),function(pair)
      cbind(program=nm,paired_result(d,pair[1],pair[2]))))
  }))
  stats$p_Holm_15 <- p.adjust(stats$p,"holm")
  stats$multiplicity_family <- "All five fixed programs x all three population pairs (15 tests)"
  csv(stats,"paired_program_statistics",out)
  # B: donor pairing is shown as grey lines, not a group-average developmental arrow.
  ds <- scores %>% filter(program%in%primary_programs)
  ds$program <- factor(ds$program,levels=primary_programs)
  pB <- ggplot(ds,aes(population,score)) +
    geom_line(aes(group=donor),color="grey75",linewidth=.35) +
    geom_boxplot(aes(fill=population),outlier.shape=NA,width=.5,alpha=.25) +
    geom_point(aes(fill=group),shape=21,size=1.9,stroke=.3) +
    facet_wrap(~program,ncol=2,scales="free_y") +
    scale_fill_manual(values=c(pop_colors,group_colors),breaks=cfg$groups,name="Clinical group") +
    labs(x=NULL,y="Mean gene z (donor pseudobulk)",title="Donor-level program comparison",
         subtitle=paste0("One point per donor/population; n = ",n_distinct(meta$donor)," donors; paired lines; all 15 adjusted tests in source table")) +
    theme_fig + theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="bottom")
  save_fig(pB,"Fig4B_donor_programs",out,8,6)
})


-------------------------------------------------------------------------------
# ==============================================================================
#------------ Figure E : Focused 9-gene LAM score in Macrophages2

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# 1. Paths and prespecified settings
# ==============================================================================

input_rds <- paste0(getOption("cea.project_root"), "/macro/macro.rds")

out_dir <- file.path(
  paste0(getOption("cea.project_root"), "/A 颈动脉单细胞文章/TREM2??? macrophage–reparative fibroblast healing axis"),
  "Communications_Biology_submission/Communications_Biology_submission",
  "Fig4D_donor_pseudobulk_results"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260831)
bootstrap_iterations <- 10000L
minimum_cells_per_pseudobulk <- 20L

mature_lam_signature <- c(
  "TREM2", "APOE", "CD9",
  "LIPA", "CTSD",
  "MERTK", "CD163", "C1QA", "C1QB"
)

group_colors <- c(
  "ASYM" = "#F9C80EFF",
  "SYM"  = "#C71000FF"
)

theme_fig4 <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey35"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4)
  )

# ==============================================================================
# 2. Load object and validate donor metadata
# ==============================================================================

macro <- readRDS(input_rds)

required_metadata <- c("orig.ident", "sample", "group", "celltype")
missing_metadata <- setdiff(required_metadata, colnames(macro@meta.data))
if (length(missing_metadata) > 0L) {
  stop("Missing metadata column(s): ", paste(missing_metadata, collapse = ", "))
}

if (!"RNA" %in% Assays(macro)) {
  stop("RNA assay is absent; raw-count pseudobulk cannot be generated.")
}

# The saved macro object has integrated as its active assay and an empty RNA
# normalized-data slot. This analysis intentionally uses RNA raw counts, not the
# integrated assay, because integrated values are not appropriate pseudobulk
# expression measurements.
rna_counts <- GetAssayData(macro, assay = "RNA", layer = "counts")
if (nrow(rna_counts) == 0L || ncol(rna_counts) == 0L) {
  stop("RNA counts are empty; raw-count pseudobulk cannot be generated.")
}

md <- macro@meta.data %>%
  mutate(
    cell_id = rownames(macro@meta.data),
    donor_id = as.character(orig.ident),
    sample_label = as.character(sample),
    group = as.character(group),
    celltype = as.character(celltype)
  ) %>%
  filter(
    celltype == "Macrophages2",
    group %in% c("ASYM", "SYM"),
    !is.na(donor_id)
  )

donor_check <- md %>% distinct(donor_id, sample_label, group)
if (any(count(donor_check, donor_id)$n != 1L)) {
  stop("At least one donor is assigned to more than one ASYM/SYM group.")
}

group_n <- donor_check %>% count(group)
if (!all(c("ASYM", "SYM") %in% group_n$group)) {
  stop("ASYM or SYM donors are missing from Macrophages2.")
}

# ==============================================================================
# 3. Construct one Macrophages2 raw-count pseudobulk per donor
# ==============================================================================

md <- md %>%
  mutate(
    pseudobulk_id = paste(donor_id, group, sep = "___"),
    pseudobulk_id = factor(pseudobulk_id, levels = unique(pseudobulk_id))
  )

if (!all(md$cell_id %in% colnames(rna_counts))) {
  stop("Some selected metadata cell IDs are absent from the RNA count matrix.")
}
rna_counts <- rna_counts[, md$cell_id, drop = FALSE]

design_pb <- Matrix::sparse.model.matrix(~ 0 + pseudobulk_id, data = md)
colnames(design_pb) <- levels(md$pseudobulk_id)
pseudobulk_counts <- rna_counts %*% design_pb

pseudobulk_meta <- md %>%
  count(donor_id, sample_label, group, pseudobulk_id, name = "n_cells") %>%
  mutate(pseudobulk_id = as.character(pseudobulk_id))

pseudobulk_meta <- pseudobulk_meta[
  match(colnames(pseudobulk_counts), pseudobulk_meta$pseudobulk_id),
  , drop = FALSE
]

if (any(is.na(pseudobulk_meta$pseudobulk_id))) {
  stop("Pseudobulk count columns could not be matched to donor metadata.")
}

excluded_low_cell <- pseudobulk_meta %>%
  filter(n_cells < minimum_cells_per_pseudobulk)

keep_pb <- pseudobulk_meta$n_cells >= minimum_cells_per_pseudobulk
pseudobulk_counts <- pseudobulk_counts[, keep_pb, drop = FALSE]
pseudobulk_meta <- pseudobulk_meta[keep_pb, , drop = FALSE]

if (nrow(excluded_low_cell) > 0L) {
  warning(
    "Excluded donor pseudobulks with fewer than ",
    minimum_cells_per_pseudobulk, " Macrophages2 cells: ",
    paste(excluded_low_cell$donor_id, collapse = ", ")
  )
}

if (n_distinct(pseudobulk_meta$donor_id) != nrow(pseudobulk_meta)) {
  stop("Each retained donor must contribute exactly one Macrophages2 pseudobulk.")
}

# ==============================================================================
# 4. edgeR TMM normalization and mature-LAM signature score
# ==============================================================================

nonzero_gene <- Matrix::rowSums(pseudobulk_counts) > 0
pseudobulk_counts <- pseudobulk_counts[nonzero_gene, , drop = FALSE]

dge <- edgeR::DGEList(counts = pseudobulk_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

signature_present <- intersect(mature_lam_signature, rownames(log_cpm))
signature_missing <- setdiff(mature_lam_signature, signature_present)

if (length(signature_present) < 5L) {
  stop(
    "Fewer than five mature-LAM signature genes are present: ",
    paste(signature_present, collapse = ", ")
  )
}

signature_log_cpm <- log_cpm[signature_present, , drop = FALSE]
gene_sd <- apply(signature_log_cpm, 1, sd)
constant_genes <- names(gene_sd)[!is.finite(gene_sd) | gene_sd == 0]

if (length(constant_genes) > 0L) {
  warning(
    "Constant signature genes were excluded from z scoring: ",
    paste(constant_genes, collapse = ", ")
  )
  signature_log_cpm <- signature_log_cpm[
    setdiff(rownames(signature_log_cpm), constant_genes),
    , drop = FALSE
  ]
}

signature_gene_z <- t(scale(t(signature_log_cpm)))
mature_lam_score <- colMeans(signature_gene_z, na.rm = TRUE)

pseudobulk_meta <- pseudobulk_meta %>%
  mutate(
    mature_lam_pb_score = mature_lam_score[pseudobulk_id],
    # Sensitivity metric without per-gene standardization.
    mature_lam_mean_logCPM = colMeans(signature_log_cpm)[pseudobulk_id],
    group = factor(group, levels = c("ASYM", "SYM"))
  )

if (any(!is.finite(pseudobulk_meta$mature_lam_pb_score))) {
  stop("At least one mature-LAM pseudobulk score is non-finite.")
}

# ==============================================================================
# 5. Donor-level inference, effect size and confidence interval
# ==============================================================================

bootstrap_independent_difference <- function(x1, x0,
                                             iterations = bootstrap_iterations) {
  boot <- replicate(
    iterations,
    mean(sample(x1, length(x1), replace = TRUE)) -
      mean(sample(x0, length(x0), replace = TRUE))
  )
  ci <- unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
  c(
    estimate = mean(x1) - mean(x0),
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

run_donor_test <- function(data, score_col, analysis_label) {
  x0 <- data[[score_col]][data$group == "ASYM"]
  x1 <- data[[score_col]][data$group == "SYM"]
  exact_ok <- !anyDuplicated(c(x0, x1))
  
  fit <- wilcox.test(
    x1, x0,
    paired = FALSE,
    alternative = "two.sided",
    exact = exact_ok,
    correct = !exact_ok
  )
  
  u <- unname(fit$statistic)
  rank_biserial <- 2 * u / (length(x1) * length(x0)) - 1
  boot <- bootstrap_independent_difference(x1, x0)
  
  tibble(
    analysis = analysis_label,
    score_column = score_col,
    comparison = "SYM - ASYM",
    n_ASYM = length(x0),
    n_SYM = length(x1),
    test = ifelse(exact_ok, "Exact Wilcoxon rank-sum (Mann-Whitney)",
                  "Wilcoxon rank-sum with continuity correction"),
    statistic = unname(fit$statistic),
    p_value = fit$p.value,
    rank_biserial = rank_biserial,
    mean_difference = unname(boot["estimate"]),
    mean_difference_ci_low = unname(boot["ci_low"]),
    mean_difference_ci_high = unname(boot["ci_high"])
  )
}

primary_test <- run_donor_test(
  pseudobulk_meta,
  "mature_lam_pb_score",
  "Primary: mean gene-wise z score"
)

sensitivity_test <- run_donor_test(
  pseudobulk_meta,
  "mature_lam_mean_logCPM",
  "Sensitivity: mean signature-gene logCPM"
)

descriptive_stats <- pseudobulk_meta %>%
  group_by(group) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    total_macrophages2_cells = sum(n_cells),
    mean_score = mean(mature_lam_pb_score),
    sd_score = sd(mature_lam_pb_score),
    median_score = median(mature_lam_pb_score),
    q1_score = unname(quantile(mature_lam_pb_score, 0.25)),
    q3_score = unname(quantile(mature_lam_pb_score, 0.75)),
    .groups = "drop"
  )
write.csv(
  pseudobulk_meta,
  file.path(out_dir, "Fig4D_donor_pseudobulk_scores.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  as.data.frame(signature_log_cpm),
  file.path(out_dir, "Fig4D_signature_gene_logCPM.csv"),
  row.names = TRUE, fileEncoding = "UTF-8"
)

write.csv(
  as.data.frame(signature_gene_z),
  file.path(out_dir, "Fig4D_signature_gene_z_scores.csv"),
  row.names = TRUE, fileEncoding = "UTF-8"
)

write.csv(
  primary_test,
  file.path(out_dir, "Fig4D_primary_donor_Wilcoxon_test.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  sensitivity_test,
  file.path(out_dir, "Fig4D_sensitivity_mean_logCPM_test.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  descriptive_stats,
  file.path(out_dir, "Fig4D_descriptive_statistics.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  excluded_low_cell,
  file.path(out_dir, "Fig4D_excluded_low_cell_pseudobulks.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# ==========================================================================================
# E：当前对应评分图（原 pC）
# ==========================================================================================

local({
  source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
  if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE))
    source_path <- rstudioapi::getActiveDocumentContext()$path
  if (!nzchar(source_path)) stop("Use source('.../01_donor_panels.R', encoding='UTF-8').")
  FIG4_HOME <- dirname(source_path)
  source(file.path(FIG4_HOME,"00_common.R"),local=TRUE,encoding="UTF-8")
  out <- new_run("01_donor"); set.seed(cfg$seed)
  message("Read-only source: ",cfg$input)
  dat <- read_source(); scope_table(dat$md,out)
  pb <- make_pb(dat$counts,dat$md); csv(pb$qc,"pseudobulk_cell_count_QC",out)
  rm(dat); gc()
  dge <- edgeR::DGEList(as.matrix(pb$counts[Matrix::rowSums(pb$counts)>0,,drop=FALSE]))
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  lcpm <- edgeR::cpm(dge,log=TRUE,prior.count=2)
  meta <- pb$meta; meta$population <- factor(meta$population,levels=cfg$populations)
  meta$group <- factor(meta$group,levels=cfg$groups)
  csv(cbind(meta,library_size=dge$samples$lib.size,TMM_factor=dge$samples$norm.factors),"pseudobulk_metadata",out)
  membership <- bind_rows(lapply(names(programs),function(nm) data.frame(program=nm,gene=programs[[nm]])))
  membership$origin <- "Investigator-defined; partitioned from existing Fig4 marker panels; not a functional assay"
  membership$in_RNA <- membership$gene%in%rownames(lcpm)
  sds <- apply(lcpm,1,sd)
  membership$nonconstant <- unname(sds[membership$gene]>0)
  csv(membership,"program_definitions_and_QC",out)
  z <- row_z(lcpm)
  existing <- function(name) read.csv(file.path(cfg$existing_pb,name),check.names=FALSE)
  m2 <- existing("Fig4D_donor_pseudobulk_scores.csv")
  test <- existing("Fig4D_primary_donor_Wilcoxon_test.csv")
  oldz <- as.matrix(read.csv(file.path(cfg$existing_pb,"Fig4D_signature_gene_z_scores.csv"),row.names=1,check.names=FALSE))
  stopifnot(!anyDuplicated(m2$donor_id),all(m2$group%in%cfg$groups),setequal(rownames(oldz),focused_genes),
            setequal(m2$pseudobulk_id,colnames(oldz)))
  oldz <- oldz[focused_genes,m2$pseudobulk_id,drop=FALSE]
  if(max(abs(colMeans(oldz)-m2$mature_lam_pb_score))>1e-8) stop("Existing score and gene-z source mismatch.")
  current <- meta[as.character(meta$population)=="Macrophages2",]
  if(!setequal(m2$donor_id,current$donor)) stop("Existing Macro2 donors differ from current input.")
  current <- current[match(m2$donor_id,current$donor),]
  stopifnot(all(m2$n_cells==current$n_cells),all(m2$sample_label==current$sample_label),
            all(m2$group==as.character(current$group)))
  # Recompute the Macro2-only normalization to verify provenance, not to choose P.
  cnt2 <- pb$counts[,current$pb_id,drop=FALSE]
  dy2 <- edgeR::calcNormFactors(edgeR::DGEList(as.matrix(cnt2[Matrix::rowSums(cnt2)>0,,drop=FALSE])))
  newz <- row_z(edgeR::cpm(dy2,log=TRUE,prior.count=2)[focused_genes,,drop=FALSE])
  difference <- max(abs(unname(newz)-unname(oldz)))
  csv(data.frame(max_absolute_gene_z_difference=difference,tolerance=1e-7),"existing_Macro2_reproduction_check",out)
  if(difference>1e-7) stop("Existing Macro2 scores do not reproduce from current RNA counts. Do not reuse P.")
  pcheck <- wilcox.test(m2$mature_lam_pb_score[m2$group=="SYM"],m2$mature_lam_pb_score[m2$group=="ASYM"],exact=TRUE)$p.value
  stopifnot(abs(pcheck-test$p_value)<1e-10)
  m2$group <- factor(m2$group,levels=cfg$groups)
  csv(m2,"C_existing_primary_Macro2_scores",out); csv(test,"C_existing_primary_Macro2_statistics",out)
  csv(existing("Fig4D_sensitivity_mean_logCPM_test.csv"),"C_existing_logCPM_sensitivity_statistics",out)
  pC <- ggplot(m2,aes(group,mature_lam_pb_score,fill=group)) +
    geom_boxplot(width=.5,outlier.shape=NA,alpha=.5) +
    geom_point(shape=21,size=2.8,position=position_jitter(width=.06,height=0,seed=cfg$seed)) +
    scale_fill_manual(values=group_colors) +
    scale_x_discrete(labels=function(x)paste0(x,"\n(n=",vapply(x,function(g)sum(m2$group==g),integer(1))," donors)")) +
    labs(x=NULL,y="Focused LAM score (mean gene z)",title="Focused LAM program in Macrophages2",
         subtitle=sprintf("Existing donor-pseudobulk analysis: two-sided P = %.3f; r_rb = %.2f",test$p_value,test$rank_biserial)) +
    theme_fig + theme(legend.position="none")
  save_fig(pC,"Fig4C_Macro2_ASYM_SYM",out,5.8,4.5)
})

-----------------------------------------------------------------------------------
# ==========================================================================================
# Figure F : Donor-level 9-gene expression heatmap in Macrophages2
local({
  source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
  if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE))
    source_path <- rstudioapi::getActiveDocumentContext()$path
  if (!nzchar(source_path)) stop("Use source('.../01_donor_panels.R', encoding='UTF-8').")
  FIG4_HOME <- dirname(source_path)
  source(file.path(FIG4_HOME,"00_common.R"),local=TRUE,encoding="UTF-8")
  out <- new_run("01_donor"); set.seed(cfg$seed)
  message("Read-only source: ",cfg$input)
  dat <- read_source(); scope_table(dat$md,out)
  pb <- make_pb(dat$counts,dat$md); csv(pb$qc,"pseudobulk_cell_count_QC",out)
  rm(dat); gc()
  dge <- edgeR::DGEList(as.matrix(pb$counts[Matrix::rowSums(pb$counts)>0,,drop=FALSE]))
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  lcpm <- edgeR::cpm(dge,log=TRUE,prior.count=2)
  meta <- pb$meta; meta$population <- factor(meta$population,levels=cfg$populations)
  meta$group <- factor(meta$group,levels=cfg$groups)
  csv(cbind(meta,library_size=dge$samples$lib.size,TMM_factor=dge$samples$norm.factors),"pseudobulk_metadata",out)
  membership <- bind_rows(lapply(names(programs),function(nm) data.frame(program=nm,gene=programs[[nm]])))
  membership$origin <- "Investigator-defined; partitioned from existing Fig4 marker panels; not a functional assay"
  membership$in_RNA <- membership$gene%in%rownames(lcpm)
  sds <- apply(lcpm,1,sd)
  membership$nonconstant <- unname(sds[membership$gene]>0)
  csv(membership,"program_definitions_and_QC",out)
  z <- row_z(lcpm)
  existing <- function(name) read.csv(file.path(cfg$existing_pb,name),check.names=FALSE)
  m2 <- existing("Fig4D_donor_pseudobulk_scores.csv")
  test <- existing("Fig4D_primary_donor_Wilcoxon_test.csv")
  oldz <- as.matrix(read.csv(file.path(cfg$existing_pb,"Fig4D_signature_gene_z_scores.csv"),row.names=1,check.names=FALSE))
  stopifnot(!anyDuplicated(m2$donor_id),all(m2$group%in%cfg$groups),setequal(rownames(oldz),focused_genes),
            setequal(m2$pseudobulk_id,colnames(oldz)))
  oldz <- oldz[focused_genes,m2$pseudobulk_id,drop=FALSE]
  if(max(abs(colMeans(oldz)-m2$mature_lam_pb_score))>1e-8) stop("Existing score and gene-z source mismatch.")
  current <- meta[as.character(meta$population)=="Macrophages2",]
  if(!setequal(m2$donor_id,current$donor)) stop("Existing Macro2 donors differ from current input.")
  current <- current[match(m2$donor_id,current$donor),]
  stopifnot(all(m2$n_cells==current$n_cells),all(m2$sample_label==current$sample_label),
            all(m2$group==as.character(current$group)))
  # Recompute the Macro2-only normalization to verify provenance, not to choose P.
  cnt2 <- pb$counts[,current$pb_id,drop=FALSE]
  dy2 <- edgeR::calcNormFactors(edgeR::DGEList(as.matrix(cnt2[Matrix::rowSums(cnt2)>0,,drop=FALSE])))
  newz <- row_z(edgeR::cpm(dy2,log=TRUE,prior.count=2)[focused_genes,,drop=FALSE])
  difference <- max(abs(unname(newz)-unname(oldz)))
  csv(data.frame(max_absolute_gene_z_difference=difference,tolerance=1e-7),"existing_Macro2_reproduction_check",out)
  if(difference>1e-7) stop("Existing Macro2 scores do not reproduce from current RNA counts. Do not reuse P.")
  pcheck <- wilcox.test(m2$mature_lam_pb_score[m2$group=="SYM"],m2$mature_lam_pb_score[m2$group=="ASYM"],exact=TRUE)$p.value
  stopifnot(abs(pcheck-test$p_value)<1e-10)
  m2$group <- factor(m2$group,levels=cfg$groups)
  csv(m2,"C_existing_primary_Macro2_scores",out); csv(test,"C_existing_primary_Macro2_statistics",out)
  csv(existing("Fig4D_sensitivity_mean_logCPM_test.csv"),"C_existing_logCPM_sensitivity_statistics",out)
  
  order <- m2$sample_label[order(m2$group,m2$sample_label)] # fixed label order, not sorted by effect
  h2 <- long_matrix(oldz,"gene_z") %>% left_join(m2,by=c("pb_id"="pseudobulk_id"))
  h2$gene <- factor(h2$gene,levels=rev(focused_genes)); h2$sample_label <- factor(h2$sample_label,levels=order)
  csv(h2,"D_Macro2_gene_expression_source",out)
  zlimit <- ceiling(max(abs(h2$gene_z))*2)/2
  pD <- ggplot(h2,aes(sample_label,gene,fill=gene_z)) + geom_tile(color="white",linewidth=.3) +
    facet_grid(~group,scales="free_x",space="free_x") +
    scale_fill_gradient2(low="#326A9C",mid="white",high="#B62932",limits=c(-zlimit,zlimit),name="Gene z") +
    labs(x="Donor",y=NULL,title="Signature genes in Macrophages2",subtitle="Same nine genes and normalization as the score panel") +
    theme_fig + theme(axis.text.x=element_text(angle=45,hjust=1),axis.text.y=element_text(face="italic"),axis.ticks=element_blank())
  save_fig(pD,"Fig4D_Macro2_gene_heatmap",out,7.6,4.5)
})
