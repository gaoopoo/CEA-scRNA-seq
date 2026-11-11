#-Supplemental Figure 5
# Figure 5 — 网络中心性分析
suppressMessages({
  library(CellChat); library(tidyverse); library(ggplot2)
  library(patchwork); library(ggrepel)
})

setwd("G:/13sample")
out_dir <- "G:/13sample/Fig5_centrality/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cellchat_sym <- readRDS("cellchat/cellchat_sym.rds")
cellchat_sym <- netAnalysis_computeCentrality(cellchat_sym, slot.name = "netP")

core_celltypes <- c("Fibroblasts1","Fibroblasts2","Proliferating cells2",
                    "SMCs1","SMCs3","Macrophages2","ECs","T cells1")
celltype_short <- c("Fibroblasts1"="Fibro1","Fibroblasts2"="Fibro2",
                    "Proliferating cells2"="PRFs","SMCs1"="SMC1","SMCs3"="SMC3",
                    "Macrophages2"="mLAM","ECs"="ECs","T cells1"="T cells")

# ── 提取 Sender / Receiver (聚合通路) ────────────────────────
centr <- cellchat_sym@netP$centr
all_cells <- levels(cellchat_sym@idents)
n_cells <- length(all_cells)
agg_role <- function(role) {
  mat <- sapply(centr, function(x) {
    v <- x[[role]]
    out <- setNames(rep(0, n_cells), all_cells); out[names(v)] <- v; out
  })
  rowSums(mat)
}
role_core <- data.frame(cell = all_cells,
                        Sender = agg_role("outdeg"),
                        Receiver = agg_role("indeg")) %>%
  filter(cell %in% core_celltypes) %>%
  mutate(short = celltype_short[cell])

norm01 <- function(x) (x-min(x))/(max(x)-min(x)+1e-9)
role_norm <- role_core %>% mutate(across(c(Sender,Receiver), norm01))

cat("\n=== Sender / Receiver 得分 ===\n")
print(role_core %>% mutate(across(where(is.numeric),~round(.,3))) %>% arrange(desc(Receiver)))
write.csv(role_core %>% mutate(across(where(is.numeric),~round(.,4))),
          paste0(out_dir,"centrality_scores_SR.csv"), row.names=FALSE)

# ── 图1: heatmap (细胞 × Sender/Receiver) ─────────────────
hm <- role_norm %>% dplyr::select(short, Sender, Receiver) %>%
  pivot_longer(-short, names_to="role", values_to="score")
hm$role <- factor(hm$role, levels=c("Sender","Receiver"))
cell_order <- role_norm %>% arrange(Receiver) %>% pull(short)
hm$short <- factor(hm$short, levels=cell_order)

p_hm <- ggplot(hm, aes(role, short, fill=score)) +
  geom_tile(color="white", linewidth=0.9) +
  geom_text(aes(label=sprintf("%.2f",score)), size=3.3, color="grey15") +
  scale_fill_gradientn(colors=c("#F7FBFF","#9ECAE1","#3A6FA0","#08306B"),
                       name="Normalized\nstrength") +
  labs(title="Network centrality (SYM)",
       subtitle="Sender and receiver strength, aggregated across pathways",
       x=NULL, y=NULL) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=9.5,color="grey40"),
        axis.text.x=element_text(face="bold",size=11),
        axis.text.y=element_text(face="bold",size=11),
        panel.grid=element_blank(),
        legend.key.width=unit(0.4,"cm"))

ggsave(paste0(out_dir,"centrality_heatmap.pdf"), p_hm, width=4.8, height=5, device=cairo_pdf)
ggsave(paste0(out_dir,"centrality_heatmap.png"), p_hm, width=4.8, height=5, dpi=600)


# =============================================================================
# Supplementary Figure 5 - Step 1
# Fixed gene sets and deterministic module scoring

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
})

set.seed(20260901)

# -----------------------------------------------------------------------------
# 1. Paths and fixed settings
# -----------------------------------------------------------------------------

mesenchymal_rds <- "G:/13sample/Fig3/Fig3_v3/fib_prf_v3_pseudotime.rds"
macrophage_rds  <- "G:/13sample/macro/Fig4E/macro_with_pseudotime_3grp.rds"

output_dir <- "G:/13sample/Fig5_CellChat/Fig5_temporal_donor_validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_groups <- c("ASYM", "SYM")

# HGNC symbols are deduplicated deliberately:
# CCN2 is used instead of the duplicate synonym CTGF;
# CCN1 is used instead of the duplicate synonym CYR61.
gene_sets <- list(
  Layer1_Matricellular = c(
    "POSTN", "TNC", "COMP", "THBS1", "THBS2", "SPARC", "CCN2", "CCN1"
  ),
  Layer2_ECM = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1",
    "COL6A1", "COL6A2", "COL6A3", "FN1", "LAMA4", "LAMB1",
    "LAMC1", "SPP1"
  ),
  Layer3_ImmuneIntegration = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1", "CD74", "APP", "CXCL9", "CXCL10",
    "CXCL12", "CXCL16"
  ),
  MatureLAM = c(
    "TREM2", "APOE", "CD9", "LIPA", "CTSD", "MERTK", "CD163",
    "C1QA", "C1QB"
  )
)

lineage_gene_sets <- list(
  Mesenchymal = gene_sets[c("Layer1_Matricellular", "Layer2_ECM")],
  Macrophage  = gene_sets[c("Layer3_ImmuneIntegration", "MatureLAM")]
)

# -----------------------------------------------------------------------------
# 2. Compatibility helpers
# -----------------------------------------------------------------------------

get_rna_data <- function(object) {
  # Seurat v5 uses layer=; older Seurat versions use slot=.
  tryCatch(
    GetAssayData(object, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(object, assay = "RNA", slot = "data")
  )
}

safe_minmax <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) {
    stop("Pseudotime cannot be normalized because its range is zero or invalid.")
  }
  (x - r[1]) / diff(r)
}

choose_celltype_column <- function(metadata, lineage) {
  if (lineage == "Mesenchymal") {
    candidates <- c("celltype_v3", "celltype")
  } else {
    candidates <- c("celltype_v4", "celltype")
  }
  out <- candidates[candidates %in% colnames(metadata)][1]
  if (is.na(out)) stop("No compatible cell-type column was found for ", lineage)
  out
}

# -----------------------------------------------------------------------------
# 3. Deterministic primary module score
# -----------------------------------------------------------------------------

calculate_gene_z_scores <- function(object, lineage, modules) {
  md <- object@meta.data
  required <- c("orig.ident", "group", "pseudotime")
  missing_required <- setdiff(required, colnames(md))
  if (length(missing_required) > 0L) {
    stop(
      lineage, " object is missing metadata columns: ",
      paste(missing_required, collapse = ", ")
    )
  }
  
  celltype_col <- choose_celltype_column(md, lineage)
  eligible <- rownames(md)[
    md$group %in% analysis_groups & is.finite(as.numeric(md$pseudotime))
  ]
  if (length(eligible) == 0L) stop("No eligible plaque cells for ", lineage)
  
  expr <- get_rna_data(object)
  eligible <- intersect(eligible, colnames(expr))
  md <- md[eligible, , drop = FALSE]
  expr <- expr[, eligible, drop = FALSE]
  
  score_list <- vector("list", length(modules))
  qc_list <- vector("list", length(modules))
  names(score_list) <- names(modules)
  
  for (i in seq_along(modules)) {
    module_name <- names(modules)[i]
    requested <- unique(modules[[i]])
    available <- intersect(requested, rownames(expr))
    missing_genes <- setdiff(requested, available)
    
    if (length(available) < 3L) {
      stop(module_name, " has fewer than three available genes in ", lineage)
    }
    
    # Only a small gene-set x cell matrix is made dense.
    small <- as.matrix(expr[available, , drop = FALSE])
    gene_mean <- rowMeans(small, na.rm = TRUE)
    gene_sd <- apply(small, 1, sd, na.rm = TRUE)
    informative <- is.finite(gene_sd) & gene_sd > 0
    
    if (sum(informative) < 3L) {
      stop(module_name, " has fewer than three variable genes in ", lineage)
    }
    
    z <- sweep(
      small[informative, , drop = FALSE],
      1,
      gene_mean[informative],
      FUN = "-"
    )
    z <- sweep(z, 1, gene_sd[informative], FUN = "/")
    
    score_list[[i]] <- colMeans(z, na.rm = TRUE)
    qc_list[[i]] <- data.frame(
      lineage = lineage,
      module = module_name,
      requested_n = length(requested),
      available_n = length(available),
      variable_n = sum(informative),
      available_genes = paste(available, collapse = ";"),
      missing_genes = paste(missing_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  
  score_wide <- as.data.frame(score_list, check.names = FALSE)
  score_long <- score_wide %>%
    mutate(cell_id = eligible) %>%
    pivot_longer(
      cols = all_of(names(modules)),
      names_to = "module",
      values_to = "score_gene_z"
    )
  
  cell_metadata <- data.frame(
    cell_id = eligible,
    donor_id = as.character(md$orig.ident),
    group = as.character(md$group),
    lineage = lineage,
    celltype = as.character(md[[celltype_col]]),
    pseudotime = as.numeric(md$pseudotime),
    stringsAsFactors = FALSE
  ) %>%
    mutate(pt_norm = safe_minmax(pseudotime))
  
  score_long <- score_long %>%
    left_join(cell_metadata, by = "cell_id") %>%
    relocate(
      cell_id, donor_id, group, lineage, celltype,
      pseudotime, pt_norm, module, score_gene_z
    )
  
  list(
    scores = score_long,
    qc = bind_rows(qc_list),
    eligible_cells = eligible,
    metadata = cell_metadata
  )
}

# -----------------------------------------------------------------------------
# 4. Optional UCell sensitivity score (only if already installed)
# -----------------------------------------------------------------------------

calculate_ucell_sensitivity <- function(object, eligible_cells, modules, lineage) {
  if (!requireNamespace("UCell", quietly = TRUE)) return(NULL)
  
  message("UCell detected; calculating optional sensitivity scores for ", lineage)
  object_small <- subset(object, cells = eligible_cells)
  DefaultAssay(object_small) <- "RNA"
  
  object_small <- UCell::AddModuleScore_UCell(
    object_small,
    features = modules,
    assay = "RNA"
  )
  
  expected <- paste0(names(modules), "_UCell")
  found <- intersect(expected, colnames(object_small@meta.data))
  if (length(found) != length(expected)) {
    warning("UCell columns were not found as expected; sensitivity scores skipped.")
    return(NULL)
  }
  
  object_small@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    select(cell_id, all_of(found)) %>%
    pivot_longer(
      cols = all_of(found),
      names_to = "module_ucell",
      values_to = "score_ucell"
    ) %>%
    mutate(
      module = sub("_UCell$", "", module_ucell),
      lineage = lineage
    ) %>%
    select(cell_id, lineage, module, score_ucell)
}

# -----------------------------------------------------------------------------
# 5. Run one lineage at a time to limit memory usage
# -----------------------------------------------------------------------------

stopifnot(file.exists(mesenchymal_rds), file.exists(macrophage_rds))

message("Loading mesenchymal trajectory object...")
mes_obj <- readRDS(mesenchymal_rds)
mes_result <- calculate_gene_z_scores(
  mes_obj,
  lineage = "Mesenchymal",
  modules = lineage_gene_sets$Mesenchymal
)
mes_ucell <- calculate_ucell_sensitivity(
  mes_obj,
  eligible_cells = mes_result$eligible_cells,
  modules = lineage_gene_sets$Mesenchymal,
  lineage = "Mesenchymal"
)
rm(mes_obj)
gc()

message("Loading macrophage trajectory object...")
mac_obj <- readRDS(macrophage_rds)
mac_result <- calculate_gene_z_scores(
  mac_obj,
  lineage = "Macrophage",
  modules = lineage_gene_sets$Macrophage
)
mac_ucell <- calculate_ucell_sensitivity(
  mac_obj,
  eligible_cells = mac_result$eligible_cells,
  modules = lineage_gene_sets$Macrophage,
  lineage = "Macrophage"
)
rm(mac_obj)
gc()

cell_scores <- bind_rows(mes_result$scores, mac_result$scores)
gene_set_qc <- bind_rows(mes_result$qc, mac_result$qc)

ucell_scores <- bind_rows(mes_ucell, mac_ucell)
if (nrow(ucell_scores) > 0L) {
  cell_scores <- cell_scores %>%
    left_join(
      ucell_scores,
      by = c("cell_id", "lineage", "module")
    )
} else {
  cell_scores$score_ucell <- NA_real_
}

# Score-method agreement is descriptive sensitivity information only.
score_method_agreement <- cell_scores %>%
  filter(is.finite(score_ucell)) %>%
  group_by(lineage, module) %>%
  summarise(
    n_cells = n(),
    spearman_rho = cor(score_gene_z, score_ucell, method = "spearman"),
    .groups = "drop"
  )

cell_count_qc <- cell_scores %>%
  distinct(cell_id, donor_id, group, lineage, celltype) %>%
  count(lineage, group, donor_id, celltype, name = "n_cells")

donor_module_means <- cell_scores %>%
  group_by(lineage, group, donor_id, module) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(score_gene_z, na.rm = TRUE),
    median_score = median(score_gene_z, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 6. Export auditable intermediate files
# -----------------------------------------------------------------------------

saveRDS(
  cell_scores,
  file.path(output_dir, "Fig5S_cell_level_module_scores.rds")
)
write_csv(
  gene_set_qc,
  file.path(output_dir, "Fig5S_gene_set_QC.csv")
)
write_csv(
  cell_count_qc,
  file.path(output_dir, "Fig5S_cell_count_QC.csv")
)
write_csv(
  donor_module_means,
  file.path(output_dir, "Fig5S_donor_module_means.csv")
)
write_csv(
  score_method_agreement,
  file.path(output_dir, "Fig5S_score_method_agreement.csv")
)

methods_record <- c(
  "Primary module score: mean gene-wise z score from log-normalized RNA expression.",
  "Gene symbols were fixed before analysis and duplicate synonyms were removed.",
  "Only plaque donors (ASYM and SYM) were included; Normal n=1 was excluded from inference.",
  "Pseudotime was normalized separately within the mesenchymal and macrophage trajectories.",
  "Absolute module-score or pseudotime values must not be compared between lineages.",
  if (nrow(ucell_scores) > 0L) {
    "UCell was available and calculated as an optional sensitivity score."
  } else {
    "UCell was not installed; no UCell sensitivity score was calculated."
  }
)
writeLines(
  methods_record,
  file.path(output_dir, "Fig5S_step1_methods_record.txt"),
  useBytes = TRUE
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "Fig5S_step1_sessionInfo.txt")
)

message("Step 1 completed. Output directory: ", output_dir)
message("Next run: Fig5S_step2_donor_GAM_and_statistics.R")

-------------------------------------------------------------------------
# =============================================================================
# Supplementary Figure 5 - Step 2

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(mgcv)
})

set.seed(20260901)

# -----------------------------------------------------------------------------
# 1. Paths and prespecified settings
# -----------------------------------------------------------------------------

output_dir <- "G:/13sample/Fig5_CellChat/Fig5_temporal_donor_validation"
score_rds <- file.path(output_dir, "Fig5S_cell_level_module_scores.rds")
donor_mean_csv <- file.path(output_dir, "Fig5S_donor_module_means.csv")

stopifnot(file.exists(score_rds), file.exists(donor_mean_csv))

n_pseudotime_bins <- 10L
min_cells_per_donor_bin <- 20L
smooth_spar <- 0.60
bootstrap_repetitions <- 10000L

module_levels <- c(
  "Layer1_Matricellular",
  "Layer2_ECM",
  "Layer3_ImmuneIntegration",
  "MatureLAM"
)

# -----------------------------------------------------------------------------
# 2. One equal-weight observation per donor x bin x module
# -----------------------------------------------------------------------------

cell_scores <- readRDS(score_rds) %>%
  filter(
    group %in% c("ASYM", "SYM"),
    is.finite(pt_norm),
    is.finite(score_gene_z)
  ) %>%
  mutate(
    pt_bin = pmin(
      n_pseudotime_bins,
      floor(pt_norm * n_pseudotime_bins) + 1L
    ),
    pt_mid = (pt_bin - 0.5) / n_pseudotime_bins,
    module = factor(module, levels = module_levels)
  )

donor_bin_all <- cell_scores %>%
  group_by(lineage, group, donor_id, module, pt_bin, pt_mid) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(score_gene_z, na.rm = TRUE),
    median_score = median(score_gene_z, na.rm = TRUE),
    .groups = "drop"
  )

excluded_donor_bins <- donor_bin_all %>%
  filter(n_cells < min_cells_per_donor_bin)

donor_bin <- donor_bin_all %>%
  filter(n_cells >= min_cells_per_donor_bin) %>%
  group_by(lineage, module) %>%
  mutate(
    # Standardization changes only the within-module y scale. It does not align
    # the two independently inferred lineage pseudotimes.
    score_std = as.numeric(scale(mean_score))
  ) %>%
  ungroup()

if (any(!is.finite(donor_bin$score_std))) {
  stop("Non-finite standardized donor-bin scores were generated.")
}

write_csv(
  donor_bin,
  file.path(output_dir, "Fig5S_donor_pseudotime_bins.csv")
)
write_csv(
  excluded_donor_bins,
  file.path(output_dir, "Fig5S_excluded_donor_bins.csv")
)

bin_coverage <- donor_bin %>%
  group_by(lineage, module, group, donor_id) %>%
  summarise(
    retained_bins = n_distinct(pt_bin),
    total_cells = sum(n_cells),
    min_cells_in_bin = min(n_cells),
    max_cells_in_bin = max(n_cells),
    .groups = "drop"
  )

write_csv(
  bin_coverage,
  file.path(output_dir, "Fig5S_donor_bin_coverage.csv")
)

# -----------------------------------------------------------------------------
# 3. Donor-balanced lineage-specific GAM curves
# -----------------------------------------------------------------------------

fit_lineage_gam <- function(data, lineage_name, requested_modules) {
  dat <- data %>%
    filter(lineage == lineage_name, module %in% requested_modules) %>%
    mutate(
      module = factor(as.character(module), levels = requested_modules),
      donor_id = factor(donor_id)
    )
  
  if (n_distinct(dat$donor_id) < 6L) {
    stop("Too few donors for ", lineage_name, " GAM.")
  }
  
  # Donor-bin rows receive equal model weight. n_cells is an inclusion
  # threshold only, so large donors cannot dominate the fitted curve.
  model <- mgcv::gam(
    score_std ~ module +
      s(pt_mid, by = module, k = 5, bs = "cr") +
      s(donor_id, bs = "re"),
    data = dat,
    method = "REML"
  )
  
  prediction_grid <- expand_grid(
    module = factor(requested_modules, levels = requested_modules),
    pt_mid = seq(0.001, 0.999, length.out = 500),
    donor_id = factor(levels(dat$donor_id)[1], levels = levels(dat$donor_id))
  )
  
  prediction <- predict(
    model,
    newdata = prediction_grid,
    type = "link",
    se.fit = TRUE,
    exclude = "s(donor_id)"
  )
  
  prediction_df <- prediction_grid %>%
    mutate(
      lineage = lineage_name,
      fit = as.numeric(prediction$fit),
      se = as.numeric(prediction$se.fit),
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se
    )
  
  list(model = model, prediction = prediction_df)
}

mes_gam <- fit_lineage_gam(
  donor_bin,
  lineage_name = "Mesenchymal",
  requested_modules = c("Layer1_Matricellular", "Layer2_ECM")
)

mac_gam <- fit_lineage_gam(
  donor_bin,
  lineage_name = "Macrophage",
  requested_modules = c("Layer3_ImmuneIntegration", "MatureLAM")
)

gam_predictions <- bind_rows(mes_gam$prediction, mac_gam$prediction)
write_csv(
  gam_predictions,
  file.path(output_dir, "Fig5S_GAM_marginal_predictions.csv")
)

capture.output(
  summary(mes_gam$model),
  file = file.path(output_dir, "Fig5S_mesenchymal_GAM_summary.txt")
)
capture.output(
  summary(mac_gam$model),
  file = file.path(output_dir, "Fig5S_macrophage_GAM_summary.txt")
)

# -----------------------------------------------------------------------------
# 4. Donor-specific smoothed peak positions
# -----------------------------------------------------------------------------

estimate_donor_peak <- function(one_donor_module) {
  z <- one_donor_module %>%
    arrange(pt_mid) %>%
    distinct(pt_mid, .keep_all = TRUE)
  
  if (nrow(z) < 5L) {
    return(data.frame(
      peak_position = NA_real_,
      observed_min = min(z$pt_mid),
      observed_max = max(z$pt_mid),
      retained_bins = nrow(z),
      boundary_peak = NA
    ))
  }
  
  fit <- smooth.spline(
    x = z$pt_mid,
    y = z$score_std,
    spar = smooth_spar
  )
  grid <- seq(min(z$pt_mid), max(z$pt_mid), length.out = 901)
  fitted <- predict(fit, x = grid)$y
  peak <- grid[which.max(fitted)]
  
  data.frame(
    peak_position = peak,
    observed_min = min(z$pt_mid),
    observed_max = max(z$pt_mid),
    retained_bins = nrow(z),
    boundary_peak = peak %in% range(grid)
  )
}

donor_peaks <- donor_bin %>%
  group_by(lineage, group, donor_id, module) %>%
  group_modify(~ estimate_donor_peak(.x)) %>%
  ungroup()

write_csv(
  donor_peaks,
  file.path(output_dir, "Fig5S_donor_module_peak_positions.csv")
)

# -----------------------------------------------------------------------------
# 5. Prespecified donor-level metrics
# -----------------------------------------------------------------------------

mes_peak_difference <- donor_peaks %>%
  filter(
    lineage == "Mesenchymal",
    module %in% c("Layer1_Matricellular", "Layer2_ECM")
  ) %>%
  select(group, donor_id, module, peak_position) %>%
  pivot_wider(names_from = module, values_from = peak_position) %>%
  drop_na(Layer1_Matricellular, Layer2_ECM) %>%
  mutate(peak_difference = Layer2_ECM - Layer1_Matricellular)

macro_bin_wide <- donor_bin %>%
  filter(
    lineage == "Macrophage",
    module %in% c("Layer3_ImmuneIntegration", "MatureLAM")
  ) %>%
  select(group, donor_id, pt_mid, module, score_std) %>%
  pivot_wider(names_from = module, values_from = score_std) %>%
  drop_na(Layer3_ImmuneIntegration, MatureLAM)

macro_donor_metrics <- macro_bin_wide %>%
  group_by(group, donor_id) %>%
  summarise(
    retained_bins = n(),
    rho_Layer3_vs_pseudotime = cor(
      pt_mid,
      Layer3_ImmuneIntegration,
      method = "spearman",
      use = "complete.obs"
    ),
    rho_Layer3_vs_MatureLAM = cor(
      Layer3_ImmuneIntegration,
      MatureLAM,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

write_csv(
  mes_peak_difference,
  file.path(output_dir, "Fig5S_mesenchymal_paired_peak_difference.csv")
)
write_csv(
  macro_donor_metrics,
  file.path(output_dir, "Fig5S_macrophage_donor_consistency_metrics.csv")
)

# -----------------------------------------------------------------------------
# 6. Effect sizes, donor bootstrap, and primary tests
# -----------------------------------------------------------------------------

bootstrap_median_ci <- function(x, repetitions = bootstrap_repetitions) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  boot <- replicate(
    repetitions,
    median(sample(x, size = length(x), replace = TRUE))
  )
  unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
}

bootstrap_independent_median_difference_ci <- function(
    x,
    y,
    repetitions = bootstrap_repetitions) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) {
    return(c(NA_real_, NA_real_))
  }
  boot <- replicate(
    repetitions,
    median(sample(x, length(x), replace = TRUE)) -
      median(sample(y, length(y), replace = TRUE))
  )
  unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
}

signed_rank_biserial <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (length(x) == 0L) return(NA_real_)
  ranks <- rank(abs(x), ties.method = "average")
  (sum(ranks[x > 0]) - sum(ranks[x < 0])) / sum(ranks)
}

one_sample_signed_summary <- function(x, analysis_name) {
  x <- x[is.finite(x)]
  nonzero <- x[x != 0]
  exact_ok <- length(nonzero) == length(x) && !anyDuplicated(abs(nonzero))
  
  if (length(nonzero) == 0L) {
    p_value <- 1
    statistic <- 0
  } else {
    test <- wilcox.test(
      x,
      mu = 0,
      alternative = "two.sided",
      exact = exact_ok,
      correct = !exact_ok
    )
    p_value <- test$p.value
    statistic <- unname(test$statistic)
  }
  
  ci <- bootstrap_median_ci(x)
  data.frame(
    analysis = analysis_name,
    test = "Two-sided one-sample Wilcoxon signed-rank test",
    n_donors = length(x),
    estimate_median = median(x),
    ci95_lower = ci[1],
    ci95_upper = ci[2],
    statistic = statistic,
    p_value = p_value,
    rank_biserial = signed_rank_biserial(x),
    positive_n = sum(x > 0),
    zero_n = sum(x == 0),
    negative_n = sum(x < 0),
    exact_test = exact_ok,
    stringsAsFactors = FALSE
  )
}

primary_statistics <- bind_rows(
  one_sample_signed_summary(
    mes_peak_difference$peak_difference,
    "Mesenchymal peak: Layer2_ECM minus Layer1_Matricellular"
  ),
  one_sample_signed_summary(
    macro_donor_metrics$rho_Layer3_vs_pseudotime,
    "Macrophage donor rho: Layer3 versus pseudotime"
  ),
  one_sample_signed_summary(
    macro_donor_metrics$rho_Layer3_vs_MatureLAM,
    "Macrophage donor rho: Layer3 versus MatureLAM"
  )
) %>%
  mutate(holm_q_value = p.adjust(p_value, method = "holm"))

write_csv(
  primary_statistics,
  file.path(output_dir, "Fig5S_primary_donor_statistics.csv")
)

# -----------------------------------------------------------------------------
# 7. Donor-level ASYM versus SYM sensitivity comparisons
# -----------------------------------------------------------------------------

donor_means <- read_csv(donor_mean_csv, show_col_types = FALSE) %>%
  mutate(module = factor(module, levels = module_levels))

independent_group_summary <- function(one_module) {
  x <- one_module$mean_score[one_module$group == "SYM"]
  y <- one_module$mean_score[one_module$group == "ASYM"]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  exact_ok <- !anyDuplicated(c(x, y))
  test <- wilcox.test(
    x,
    y,
    paired = FALSE,
    alternative = "two.sided",
    exact = exact_ok,
    correct = !exact_ok
  )
  
  # R's two-sample W equals the Mann-Whitney U for the first group (SYM).
  u_statistic <- unname(test$statistic)
  rank_biserial <- 2 * u_statistic / (length(x) * length(y)) - 1
  ci <- bootstrap_independent_median_difference_ci(x, y)
  
  data.frame(
    lineage = unique(one_module$lineage),
    module = as.character(unique(one_module$module)),
    comparison = "SYM minus ASYM",
    n_SYM = length(x),
    n_ASYM = length(y),
    median_SYM = median(x),
    median_ASYM = median(y),
    median_difference = median(x) - median(y),
    ci95_lower = ci[1],
    ci95_upper = ci[2],
    wilcoxon_U = u_statistic,
    p_value = test$p.value,
    rank_biserial = rank_biserial,
    exact_test = exact_ok,
    stringsAsFactors = FALSE
  )
}

group_statistics <- donor_means %>%
  group_split(lineage, module) %>%
  map_dfr(independent_group_summary) %>%
  mutate(holm_q_value = p.adjust(p_value, method = "holm"))

write_csv(
  group_statistics,
  file.path(output_dir, "Fig5S_ASYM_SYM_donor_module_statistics.csv")
)

# -----------------------------------------------------------------------------
# 8. Save all inputs for plotting in Step 3
# -----------------------------------------------------------------------------

analysis_results <- list(
  settings = list(
    n_pseudotime_bins = n_pseudotime_bins,
    min_cells_per_donor_bin = min_cells_per_donor_bin,
    smooth_spar = smooth_spar,
    bootstrap_repetitions = bootstrap_repetitions
  ),
  donor_bin = donor_bin,
  excluded_donor_bins = excluded_donor_bins,
  bin_coverage = bin_coverage,
  gam_predictions = gam_predictions,
  mesenchymal_gam = mes_gam$model,
  macrophage_gam = mac_gam$model,
  donor_peaks = donor_peaks,
  mes_peak_difference = mes_peak_difference,
  macro_donor_metrics = macro_donor_metrics,
  primary_statistics = primary_statistics,
  donor_means = donor_means,
  group_statistics = group_statistics
)

saveRDS(
  analysis_results,
  file.path(output_dir, "Fig5S_analysis_results.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "Fig5S_step2_sessionInfo.txt")
)

message("Step 2 completed. Primary statistics:")
print(primary_statistics)
message("Next run: Fig5S_step3_make_figures.R")

--------------------------------------------------------------------------------
  
  # =============================================================================
# Supplementary Figure 5 - Step 3
# Result-adapted, donor-balanced, font-safe version
# =============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260901)


# =============================================================================
# 1. Paths and input
# =============================================================================

output_dir <- paste0(
  "G:/13sample/Fig5_CellChat/",
  "Fig5_temporal_donor_validation"
)

result_rds <- file.path(
  output_dir,
  "Fig5S_analysis_results.rds"
)

if (!file.exists(result_rds)) {
  stop(
    "Cannot find Step 2 analysis object: ",
    result_rds
  )
}

res <- readRDS(result_rds)

donor_bin <- res$donor_bin
gam_predictions <- res$gam_predictions
macro_donor_metrics <- res$macro_donor_metrics
primary_statistics <- res$primary_statistics
donor_peaks <- res$donor_peaks
donor_means <- res$donor_means
group_statistics <- res$group_statistics


# =============================================================================
# 2. Colours and labels
# =============================================================================

module_colours <- c(
  Layer1_Matricellular = "#D55E00",
  Layer2_ECM = "#E69F00",
  Layer3_ImmuneIntegration = "#0072B2",
  MatureLAM = "#56B4E9"
)

group_colours <- c(
  ASYM = "#F0C419",
  SYM = "#D7301F"
)

module_labels <- c(
  Layer1_Matricellular = "Layer 1: Matricellular",
  Layer2_ECM = "Layer 2: ECM",
  Layer3_ImmuneIntegration = "Layer 3: Immune integration",
  MatureLAM = "Mature-LAM signature"
)


# =============================================================================
# 3. Font-safe publication theme
# =============================================================================

theme_publication <- theme_classic(
  base_size = 9,
  base_family = "sans"
) +
  theme(
    text = element_text(
      family = "sans",
      colour = "black"
    ),
    
    axis.title = element_text(
      family = "sans",
      colour = "black",
      size = 9
    ),
    
    axis.text = element_text(
      family = "sans",
      colour = "black",
      size = 8
    ),
    
    axis.line = element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    plot.title = element_text(
      family = "sans",
      face = "bold",
      size = 10,
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      family = "sans",
      size = 8,
      colour = "grey30",
      hjust = 0
    ),
    
    plot.tag = element_text(
      family = "sans",
      face = "bold",
      size = 12
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey75"
    ),
    
    strip.text = element_text(
      family = "sans",
      face = "bold",
      size = 8.5
    ),
    
    legend.title = element_text(
      family = "sans",
      size = 8
    ),
    
    legend.text = element_text(
      family = "sans",
      size = 8
    ),
    
    legend.key.width = grid::unit(
      10,
      "mm"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 6,
      b = 5,
      l = 5
    )
  )


# =============================================================================
# 4. Font-safe figure export
# =============================================================================

save_figure <- function(
    plot,
    filename_stem,
    width,
    height) {
  
  pdf_file <- file.path(
    output_dir,
    paste0(filename_stem, ".pdf")
  )
  
  png_file <- file.path(
    output_dir,
    paste0(filename_stem, ".png")
  )
  
  # Cairo PDF supports Windows TrueType/sans-serif fonts reliably.
  if (isTRUE(capabilities("cairo"))) {
    
    ggsave(
      filename = pdf_file,
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )
    
  } else {
    
    # Fallback to the built-in Helvetica/sans PDF font.
    ggsave(
      filename = pdf_file,
      plot = plot,
      device = grDevices::pdf,
      width = width,
      height = height,
      units = "in",
      family = "sans",
      useDingbats = FALSE,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  # Save 600-dpi raster output.
  if (isTRUE(capabilities("cairo"))) {
    
    ggsave(
      filename = png_file,
      plot = plot,
      device = "png",
      type = "cairo",
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    
  } else {
    
    ggsave(
      filename = png_file,
      plot = plot,
      device = "png",
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  message("Saved: ", pdf_file)
  message("Saved: ", png_file)
}


# =============================================================================
# 5. Donor coverage and peak-position QC
# =============================================================================

donor_support <- donor_bin %>%
  group_by(
    lineage,
    module,
    pt_mid
  ) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    
    n_ASYM = n_distinct(
      donor_id[group == "ASYM"]
    ),
    
    n_SYM = n_distinct(
      donor_id[group == "SYM"]
    ),
    
    .groups = "drop"
  )


peak_quality <- donor_peaks %>%
  group_by(
    lineage,
    module
  ) %>%
  summarise(
    n_donors = n(),
    
    boundary_peak_n = sum(
      boundary_peak %in% TRUE,
      na.rm = TRUE
    ),
    
    boundary_peak_fraction = mean(
      boundary_peak %in% TRUE,
      na.rm = TRUE
    ),
    
    median_peak_position = median(
      peak_position,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


write_csv(
  donor_support,
  file.path(
    output_dir,
    "Fig5S_donor_support_per_bin.csv"
  )
)

write_csv(
  peak_quality,
  file.path(
    output_dir,
    "Fig5S_peak_position_quality_summary.csv"
  )
)


# =============================================================================
# 6. Donor-level Layer 1-Layer 2 coupling
# =============================================================================
# Peak differences remain available in the Step 2 outputs.
# Because most peaks occurred at the observed boundary, peak ordering is not
# used as the principal display. Instead, one within-donor correlation is
# calculated for Layer 1 versus Layer 2.

mes_coupling <- donor_bin %>%
  filter(
    lineage == "Mesenchymal",
    module %in% c(
      "Layer1_Matricellular",
      "Layer2_ECM"
    )
  ) %>%
  select(
    group,
    donor_id,
    pt_mid,
    module,
    score_std
  ) %>%
  pivot_wider(
    names_from = module,
    values_from = score_std
  ) %>%
  drop_na(
    Layer1_Matricellular,
    Layer2_ECM
  ) %>%
  group_by(
    group,
    donor_id
  ) %>%
  summarise(
    retained_bins = n(),
    
    rho_Layer1_vs_Layer2 = cor(
      Layer1_Matricellular,
      Layer2_ECM,
      method = "spearman",
      use = "complete.obs"
    ),
    
    .groups = "drop"
  )


# =============================================================================
# 7. Statistical helper functions
# =============================================================================

bootstrap_median_ci <- function(
    x,
    repetitions = 10000L) {
  
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    return(c(
      NA_real_,
      NA_real_
    ))
  }
  
  boot_values <- replicate(
    repetitions,
    median(
      sample(
        x,
        size = length(x),
        replace = TRUE
      )
    )
  )
  
  unname(
    quantile(
      boot_values,
      probs = c(
        0.025,
        0.975
      ),
      na.rm = TRUE
    )
  )
}


signed_rank_biserial <- function(x) {
  
  x <- x[
    is.finite(x) &
      x != 0
  ]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  ranks <- rank(
    abs(x),
    ties.method = "average"
  )
  
  (
    sum(ranks[x > 0]) -
      sum(ranks[x < 0])
  ) / sum(ranks)
}


one_sample_rho_summary <- function(
    x,
    analysis_name) {
  
  x <- x[is.finite(x)]
  
  exact_ok <- (
    all(x != 0) &&
      !anyDuplicated(abs(x))
  )
  
  test_result <- wilcox.test(
    x,
    mu = 0,
    alternative = "two.sided",
    exact = exact_ok,
    correct = !exact_ok
  )
  
  ci <- bootstrap_median_ci(x)
  
  data.frame(
    analysis = analysis_name,
    n_donors = length(x),
    estimate_median = median(x),
    ci95_lower = ci[1],
    ci95_upper = ci[2],
    statistic = unname(
      test_result$statistic
    ),
    p_value = test_result$p.value,
    rank_biserial = signed_rank_biserial(x),
    positive_n = sum(x > 0),
    zero_n = sum(x == 0),
    negative_n = sum(x < 0),
    exact_test = exact_ok,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 8. Donor-level correlation statistics
# =============================================================================

correlation_statistics <- bind_rows(
  
  one_sample_rho_summary(
    mes_coupling$rho_Layer1_vs_Layer2,
    paste0(
      "Mesenchymal donor rho: ",
      "Layer1 versus Layer2"
    )
  ),
  
  one_sample_rho_summary(
    macro_donor_metrics$
      rho_Layer3_vs_pseudotime,
    paste0(
      "Macrophage donor rho: ",
      "Layer3 versus pseudotime"
    )
  ),
  
  one_sample_rho_summary(
    macro_donor_metrics$
      rho_Layer3_vs_MatureLAM,
    paste0(
      "Macrophage donor rho: ",
      "Layer3 versus MatureLAM"
    )
  )
) %>%
  mutate(
    holm_q_value = p.adjust(
      p_value,
      method = "holm"
    )
  )


write_csv(
  mes_coupling,
  file.path(
    output_dir,
    "Fig5S_mesenchymal_donor_coupling_metrics.csv"
  )
)

write_csv(
  correlation_statistics,
  file.path(
    output_dir,
    "Fig5S_donor_correlation_statistics.csv"
  )
)


# =============================================================================
# 9. Curve plotting function
# =============================================================================

make_curve_panel <- function(
    lineage_name,
    requested_modules,
    panel_title,
    x_limits,
    panel_subtitle) {
  
  observed <- donor_bin %>%
    filter(
      lineage == lineage_name,
      module %in% requested_modules,
      pt_mid >= x_limits[1],
      pt_mid <= x_limits[2]
    ) %>%
    mutate(
      module = factor(
        as.character(module),
        levels = requested_modules
      )
    )
  
  
  predicted <- gam_predictions %>%
    filter(
      lineage == lineage_name,
      module %in% requested_modules,
      pt_mid >= x_limits[1],
      pt_mid <= x_limits[2]
    ) %>%
    mutate(
      module = factor(
        as.character(module),
        levels = requested_modules
      )
    )
  
  
  ggplot() +
    
    # Thin lines represent individual donors.
    geom_line(
      data = observed,
      aes(
        x = pt_mid,
        y = score_std,
        colour = module,
        group = interaction(
          donor_id,
          module
        )
      ),
      linewidth = 0.32,
      alpha = 0.18,
      na.rm = TRUE
    ) +
    
    # Marginal 95% confidence intervals.
    geom_ribbon(
      data = predicted,
      aes(
        x = pt_mid,
        ymin = lower,
        ymax = upper,
        fill = module
      ),
      alpha = 0.16,
      colour = NA,
      na.rm = TRUE
    ) +
    
    # Donor-balanced marginal GAM curves.
    geom_line(
      data = predicted,
      aes(
        x = pt_mid,
        y = fit,
        colour = module
      ),
      linewidth = 1.05,
      na.rm = TRUE
    ) +
    
    scale_colour_manual(
      name = "Module",
      values = module_colours,
      breaks = requested_modules,
      labels = module_labels[
        requested_modules
      ]
    ) +
    
    scale_fill_manual(
      name = "Module",
      values = module_colours,
      breaks = requested_modules,
      labels = module_labels[
        requested_modules
      ]
    ) +
    
    scale_x_continuous(
      limits = x_limits,
      breaks = seq(
        0,
        1,
        by = 0.25
      ),
      expand = expansion(
        mult = c(
          0.01,
          0.01
        )
      )
    ) +
    
    labs(
      title = panel_title,
      subtitle = panel_subtitle,
      x = paste0(
        "Within-lineage pseudotime ",
        "(normalized)"
      ),
      y = "Standardized module activity"
    ) +
    
    theme_publication +
    
    theme(
      legend.position = "bottom"
    )
}


# =============================================================================
# 10. Mesenchymal trajectory panel
# =============================================================================
# The 0.05 and 0.15 mesenchymal bins had only two and one retained donors.
# Therefore, the donor-supported display starts at pseudotime 0.25.

p_mes_curves <- make_curve_panel(
  lineage_name = "Mesenchymal",
  
  requested_modules = c(
    "Layer1_Matricellular",
    "Layer2_ECM"
  ),
  
  panel_title = "Mesenchymal trajectory",
  
  x_limits = c(
    0.25,
    0.95
  ),
  
  panel_subtitle = paste0(
    "Donor-balanced curves; ",
    "sparse leading bins omitted"
  )
)


# =============================================================================
# 11. Macrophage trajectory panel
# =============================================================================

p_mac_curves <- make_curve_panel(
  lineage_name = "Macrophage",
  
  requested_modules = c(
    "Layer3_ImmuneIntegration",
    "MatureLAM"
  ),
  
  panel_title = "Macrophage trajectory",
  
  x_limits = c(
    0.05,
    0.95
  ),
  
  panel_subtitle = "Donor-balanced curves"
)


# =============================================================================
# 12. Layer 1-Layer 2 coupling panel
# =============================================================================

p_mes_coupling <- ggplot(
  mes_coupling,
  aes(
    x = "Layer 1 vs Layer 2",
    y = rho_Layer1_vs_Layer2
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.45
  ) +
  
  geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.55
  ) +
  
  geom_point(
    aes(
      fill = group
    ),
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    position = position_jitter(
      width = 0.07,
      height = 0
    )
  ) +
  
  scale_fill_manual(
    name = "Clinical group",
    values = group_colours
  ) +
  
  scale_y_continuous(
    limits = c(
      -1,
      1
    ),
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  
  labs(
    title = "Layer 1-Layer 2 coupling",
    subtitle = "One correlation per donor",
    x = NULL,
    y = "Donor-level Spearman rho"
  ) +
  
  theme_publication +
  
  theme(
    legend.position = "bottom"
  )


# =============================================================================
# 13. Macrophage donor-correlation panel
# =============================================================================

macro_rho_long <- macro_donor_metrics %>%
  select(
    group,
    donor_id,
    rho_Layer3_vs_pseudotime,
    rho_Layer3_vs_MatureLAM
  ) %>%
  pivot_longer(
    cols = starts_with("rho_"),
    names_to = "metric",
    values_to = "rho"
  ) %>%
  mutate(
    metric = factor(
      metric,
      
      levels = c(
        "rho_Layer3_vs_pseudotime",
        "rho_Layer3_vs_MatureLAM"
      ),
      
      labels = c(
        "Layer 3 vs\npseudotime",
        "Layer 3 vs\nMature-LAM"
      )
    )
  )


p_macro_rho <- ggplot(
  macro_rho_long,
  aes(
    x = metric,
    y = rho
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.45
  ) +
  
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.55
  ) +
  
  geom_point(
    aes(
      fill = group
    ),
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  
  scale_fill_manual(
    name = "Clinical group",
    values = group_colours
  ) +
  
  scale_y_continuous(
    limits = c(
      -1,
      1
    ),
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  
  labs(
    title = "Within-donor immune-state alignment",
    subtitle = "One correlation per donor",
    x = NULL,
    y = "Donor-level Spearman rho"
  ) +
  
  theme_publication +
  
  theme(
    legend.position = "bottom"
  )


# =============================================================================
# 14. Complete Supplementary Figure
# =============================================================================

fig5s_validation <- (
  p_mes_curves +
    p_mes_coupling +
    p_mac_curves +
    p_macro_rho
) +
  
  plot_layout(
    ncol = 2,
    guides = "collect",
    widths = c(
      1.55,
      1
    )
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    legend.position = "bottom"
  )


save_figure(
  plot = fig5s_validation,
  filename_stem =
    "Fig5S_donor_balanced_temporal_validation",
  width = 7.2,
  height = 6.7
)


# =============================================================================
# 15. Candidate replacement for current Fig. 5G
# =============================================================================

fig5g_candidate <- (
  p_mes_curves +
    p_mac_curves
) +
  
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    legend.position = "bottom"
  )


save_figure(
  plot = fig5g_candidate,
  filename_stem =
    "Fig5G_candidate_donor_balanced_lineage_curves",
  width = 7.2,
  height = 3.4
)


# =============================================================================
# 16. Donor-level ASYM/SYM sensitivity panel
# =============================================================================

donor_means_plot <- donor_means %>%
  mutate(
    module = factor(
      as.character(module),
      levels = names(module_labels),
      labels = unname(module_labels)
    ),
    
    group = factor(
      group,
      levels = c(
        "ASYM",
        "SYM"
      )
    )
  )


p_group_sensitivity <- ggplot(
  donor_means_plot,
  aes(
    x = group,
    y = mean_score,
    fill = group
  )
) +
  
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.62,
    linewidth = 0.55
  ) +
  
  geom_point(
    shape = 21,
    size = 2.25,
    stroke = 0.45,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  
  facet_wrap(
    ~ module,
    scales = "free_y",
    ncol = 2
  ) +
  
  scale_fill_manual(
    values = group_colours
  ) +
  
  labs(
    x = NULL,
    
    y = paste0(
      "Donor-aggregated module score\n",
      "(mean gene-wise z score)"
    )
  ) +
  
  theme_publication +
  
  theme(
    legend.position = "none"
  )


save_figure(
  plot = p_group_sensitivity,
  filename_stem =
    "Fig5S_ASYM_SYM_donor_module_sensitivity",
  width = 6.3,
  height = 5.2
)


# =============================================================================
# 17. Formatting helpers for caption
# =============================================================================

fmt_num <- function(
    x,
    digits = 3) {
  
  ifelse(
    is.na(x),
    "NA",
    formatC(
      x,
      digits = digits,
      format = "f"
    )
  )
}


fmt_p <- function(x) {
  
  ifelse(
    is.na(x),
    "NA",
    format.pval(
      x,
      digits = 3,
      eps = 0.001
    )
  )
}


get_corr <- function(pattern) {
  
  correlation_statistics %>%
    filter(
      grepl(
        pattern,
        analysis,
        fixed = TRUE
      )
    ) %>%
    slice(1)
}


mes_corr <- get_corr(
  "Layer1 versus Layer2"
)

pt_corr <- get_corr(
  "Layer3 versus pseudotime"
)

lam_corr <- get_corr(
  "Layer3 versus MatureLAM"
)


peak_stat <- primary_statistics %>%
  filter(
    grepl(
      "Mesenchymal peak",
      analysis,
      fixed = TRUE
    )
  ) %>%
  slice(1)


# =============================================================================
# 18. Generate caption draft
# =============================================================================

caption_lines <- c(
  
  paste0(
    "Supplementary Figure 5 | ",
    "Donor-balanced validation of lineage-specific ",
    "temporal module activity."
  ),
  
  paste0(
    "(A) Layer 1 matricellular and Layer 2 ECM activity ",
    "along the mesenchymal trajectory. Thin lines represent ",
    "donor-bin summaries; thick lines and ribbons show marginal ",
    "GAM fits and 95% confidence intervals. The sparsely supported ",
    "0.05 and 0.15 bins were omitted from display."
  ),
  
  paste0(
    "(B) Within-donor correlation between Layer 1 and Layer 2 ",
    "module activity. Median Spearman rho = ",
    fmt_num(
      mes_corr$estimate_median
    ),
    "; donor-bootstrap 95% CI ",
    fmt_num(
      mes_corr$ci95_lower
    ),
    " to ",
    fmt_num(
      mes_corr$ci95_upper
    ),
    "; two-sided Wilcoxon P = ",
    fmt_p(
      mes_corr$p_value
    ),
    ", Holm q = ",
    fmt_p(
      mes_corr$holm_q_value
    ),
    "; positive in ",
    mes_corr$positive_n,
    "/",
    mes_corr$n_donors,
    " donors. This secondary analysis evaluates temporal ",
    "coupling, not causal or sequential ordering."
  ),
  
  paste0(
    "(C) Layer 3 immune-integration and Mature-LAM module ",
    "activity along macrophage pseudotime, displayed using ",
    "the same donor-balanced procedure."
  ),
  
  paste0(
    "(D) Within-donor correlations. For Layer 3 versus ",
    "pseudotime, median rho = ",
    fmt_num(
      pt_corr$estimate_median
    ),
    ", 95% CI ",
    fmt_num(
      pt_corr$ci95_lower
    ),
    " to ",
    fmt_num(
      pt_corr$ci95_upper
    ),
    ", P = ",
    fmt_p(
      pt_corr$p_value
    ),
    ", Holm q = ",
    fmt_p(
      pt_corr$holm_q_value
    ),
    "; positive in ",
    pt_corr$positive_n,
    "/",
    pt_corr$n_donors,
    " donors. For Layer 3 versus Mature-LAM, median rho = ",
    fmt_num(
      lam_corr$estimate_median
    ),
    ", 95% CI ",
    fmt_num(
      lam_corr$ci95_lower
    ),
    " to ",
    fmt_num(
      lam_corr$ci95_upper
    ),
    ", P = ",
    fmt_p(
      lam_corr$p_value
    ),
    ", Holm q = ",
    fmt_p(
      lam_corr$holm_q_value
    ),
    "; positive in ",
    lam_corr$positive_n,
    "/",
    lam_corr$n_donors,
    " donors."
  ),
  
  paste0(
    "Layer 1 and Layer 2 peak positions were not temporally ",
    "separable (median difference = ",
    fmt_num(
      peak_stat$estimate_median
    ),
    "; P = ",
    fmt_p(
      peak_stat$p_value
    ),
    "), and most peaks occurred at the observed boundary. ",
    "Accordingly, the modules are interpreted as overlapping ",
    "components of an early structural phase rather than ",
    "strictly ordered steps."
  ),
  
  paste0(
    "Pseudotime was normalized independently within each lineage; ",
    "mesenchymal and macrophage x-axis positions are therefore ",
    "not interpreted as a shared absolute clock. Donor was the ",
    "unit of inference."
  )
)


writeLines(
  caption_lines,
  con = file.path(
    output_dir,
    "Fig5S_caption_draft.txt"
  ),
  useBytes = TRUE
)


# =============================================================================
# 19. Generate interpretation guardrails
# =============================================================================

guardrail_lines <- c(
  
  paste0(
    "FIGURE 5 INTERPRETATION RULES BASED ON ",
    "THE OBSERVED DONOR-LEVEL RESULTS"
  ),
  
  "",
  
  paste0(
    "1. Retain the three-layer functional architecture: ",
    "matricellular, ECM, and immune-integration modules remain ",
    "biologically and topologically distinct."
  ),
  
  paste0(
    "2. Do not claim a resolved Layer 1-to-Layer 2 temporal ",
    "sequence. Their peak difference was centered at zero and ",
    "21/24 mesenchymal module peaks were boundary peaks."
  ),
  
  paste0(
    "3. Describe Layer 1 and Layer 2 as overlapping or coupled ",
    "components of an early structural phase. Their activities ",
    "were positively correlated in 12/12 donors."
  ),
  
  paste0(
    "4. Layer 3 increased along macrophage pseudotime in 12/12 ",
    "donors and covaried with the Mature-LAM program in 12/12 ",
    "donors. It is defensible to describe immune integration as ",
    "late within macrophage maturation."
  ),
  
  paste0(
    "5. Do not claim a cross-lineage absolute temporal axis: ",
    "the two pseudotimes were inferred and normalized separately."
  ),
  
  paste0(
    "6. The ASYM/SYM donor-module plot is a sensitivity analysis. ",
    "Layer 2 showed a positive trend, whereas Layer 3 donor means ",
    "did not increase in SYM."
  ),
  
  paste0(
    "7. Group-pooled CellChat strengths remain descriptive unless ",
    "CellChat is recomputed as one network per donor."
  )
)


writeLines(
  guardrail_lines,
  con = file.path(
    output_dir,
    "Fig5S_interpretation_guardrails.txt"
  ),
  useBytes = TRUE
)


# =============================================================================
# 20. Export plotting tables and session information
# =============================================================================

write_csv(
  macro_rho_long,
  file.path(
    output_dir,
    "Fig5S_plotdata_macro_donor_rho.csv"
  )
)

write_csv(
  group_statistics,
  file.path(
    output_dir,
    "Fig5S_group_statistics_for_legend.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "Fig5S_step3_sessionInfo.txt"
  )
)


# =============================================================================
# 21. Completion message
# =============================================================================

message("")
message("Step 3 completed successfully.")
message("")
message("Open these files first:")
message(
  file.path(
    output_dir,
    "Fig5S_donor_balanced_temporal_validation.pdf"
  )
)
message(
  file.path(
    output_dir,
    "Fig5G_candidate_donor_balanced_lineage_curves.pdf"
  )
)
message(
  file.path(
    output_dir,
    "Fig5S_ASYM_SYM_donor_module_sensitivity.pdf"
  )
)
message(
  file.path(
    output_dir,
    "Fig5S_interpretation_guardrails.txt"
  )
)

--------------------------------------------------------------------------------
  # =============================================================================
# Supplementary Figure 5 - Step 3 (adapted to the observed Step 2 results)
# =============================================================================
# Main decisions encoded here:
# 1. Layer 1/Layer 2 ordering is not inferred from boundary peak positions.
# 2. Their temporal coupling is summarized by one Spearman rho per donor.
# 3. Layer 3 timing is evaluated within the macrophage trajectory only.
# 4. Mesenchymal and macrophage pseudotimes are not treated as one clock.
# =============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260901)

output_dir <- "G:/13sample/Fig5_CellChat/Fig5_temporal_donor_validation"
result_rds <- file.path(output_dir, "Fig5S_analysis_results.rds")
stopifnot(file.exists(result_rds))

res <- readRDS(result_rds)

donor_bin <- res$donor_bin
gam_predictions <- res$gam_predictions
macro_donor_metrics <- res$macro_donor_metrics
primary_statistics <- res$primary_statistics
donor_peaks <- res$donor_peaks
donor_means <- res$donor_means
group_statistics <- res$group_statistics

module_colours <- c(
  Layer1_Matricellular = "#D55E00",
  Layer2_ECM = "#E69F00",
  Layer3_ImmuneIntegration = "#0072B2",
  MatureLAM = "#56B4E9"
)

group_colours <- c(ASYM = "#F0C419", SYM = "#D7301F")

module_labels <- c(
  Layer1_Matricellular = "Layer 1: Matricellular",
  Layer2_ECM = "Layer 2: ECM",
  Layer3_ImmuneIntegration = "Layer 3: Immune integration",
  MatureLAM = "Mature-LAM signature"
)

theme_publication <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.title = element_text(colour = "black", size = 9),
    axis.text = element_text(colour = "black", size = 8),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.subtitle = element_text(size = 8, colour = "grey30", hjust = 0),
    strip.background = element_rect(fill = "grey95", colour = "grey75"),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(10, "mm"),
    plot.margin = margin(5, 6, 5, 5)
  )

save_figure <- function(plot, filename_stem, width, height) {
  ggsave(
    filename = file.path(output_dir, paste0(filename_stem, ".pdf")),
    plot = plot,
    device = grDevices::pdf,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(output_dir, paste0(filename_stem, ".png")),
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1. Donor coverage and peak-position QC
# -----------------------------------------------------------------------------

donor_support <- donor_bin %>%
  group_by(lineage, module, pt_mid) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    n_ASYM = n_distinct(donor_id[group == "ASYM"]),
    n_SYM = n_distinct(donor_id[group == "SYM"]),
    .groups = "drop"
  )

peak_quality <- donor_peaks %>%
  group_by(lineage, module) %>%
  summarise(
    n_donors = n(),
    boundary_peak_n = sum(boundary_peak %in% TRUE, na.rm = TRUE),
    boundary_peak_fraction = mean(boundary_peak %in% TRUE, na.rm = TRUE),
    median_peak_position = median(peak_position, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  donor_support,
  file.path(output_dir, "Fig5S_donor_support_per_bin.csv")
)
write_csv(
  peak_quality,
  file.path(output_dir, "Fig5S_peak_position_quality_summary.csv")
)

# -----------------------------------------------------------------------------
# 2. Secondary donor-level coupling metric for Layer 1 versus Layer 2
# -----------------------------------------------------------------------------
# Peak differences remain reported in the Step 2 output. Because 21/24
# mesenchymal module peaks occurred at the observed boundary, peak ordering is
# not used as the principal display. The following correlation quantifies the
# observed temporal coupling without claiming equivalence or causality.

mes_coupling <- donor_bin %>%
  filter(
    lineage == "Mesenchymal",
    module %in% c("Layer1_Matricellular", "Layer2_ECM")
  ) %>%
  select(group, donor_id, pt_mid, module, score_std) %>%
  pivot_wider(names_from = module, values_from = score_std) %>%
  drop_na(Layer1_Matricellular, Layer2_ECM) %>%
  group_by(group, donor_id) %>%
  summarise(
    retained_bins = n(),
    rho_Layer1_vs_Layer2 = cor(
      Layer1_Matricellular,
      Layer2_ECM,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

bootstrap_median_ci <- function(x, repetitions = 10000L) {
  x <- x[is.finite(x)]
  boot <- replicate(
    repetitions,
    median(sample(x, size = length(x), replace = TRUE))
  )
  unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
}

signed_rank_biserial <- function(x) {
  x <- x[is.finite(x) & x != 0]
  ranks <- rank(abs(x), ties.method = "average")
  (sum(ranks[x > 0]) - sum(ranks[x < 0])) / sum(ranks)
}

one_sample_rho_summary <- function(x, analysis_name) {
  x <- x[is.finite(x)]
  exact_ok <- all(x != 0) && !anyDuplicated(abs(x))
  test <- wilcox.test(
    x,
    mu = 0,
    alternative = "two.sided",
    exact = exact_ok,
    correct = !exact_ok
  )
  ci <- bootstrap_median_ci(x)
  
  data.frame(
    analysis = analysis_name,
    n_donors = length(x),
    estimate_median = median(x),
    ci95_lower = ci[1],
    ci95_upper = ci[2],
    statistic = unname(test$statistic),
    p_value = test$p.value,
    rank_biserial = signed_rank_biserial(x),
    positive_n = sum(x > 0),
    zero_n = sum(x == 0),
    negative_n = sum(x < 0),
    exact_test = exact_ok,
    stringsAsFactors = FALSE
  )
}

correlation_statistics <- bind_rows(
  one_sample_rho_summary(
    mes_coupling$rho_Layer1_vs_Layer2,
    "Mesenchymal donor rho: Layer1 versus Layer2"
  ),
  one_sample_rho_summary(
    macro_donor_metrics$rho_Layer3_vs_pseudotime,
    "Macrophage donor rho: Layer3 versus pseudotime"
  ),
  one_sample_rho_summary(
    macro_donor_metrics$rho_Layer3_vs_MatureLAM,
    "Macrophage donor rho: Layer3 versus MatureLAM"
  )
) %>%
  mutate(holm_q_value = p.adjust(p_value, method = "holm"))

write_csv(
  mes_coupling,
  file.path(output_dir, "Fig5S_mesenchymal_donor_coupling_metrics.csv")
)
write_csv(
  correlation_statistics,
  file.path(output_dir, "Fig5S_donor_correlation_statistics.csv")
)

# -----------------------------------------------------------------------------
# 3. Curve plotting helper
# -----------------------------------------------------------------------------

make_curve_panel <- function(
    lineage_name,
    requested_modules,
    panel_title,
    x_limits,
    panel_subtitle) {
  observed <- donor_bin %>%
    filter(
      lineage == lineage_name,
      module %in% requested_modules,
      pt_mid >= x_limits[1],
      pt_mid <= x_limits[2]
    ) %>%
    mutate(module = factor(as.character(module), levels = requested_modules))
  
  predicted <- gam_predictions %>%
    filter(
      lineage == lineage_name,
      module %in% requested_modules,
      pt_mid >= x_limits[1],
      pt_mid <= x_limits[2]
    ) %>%
    mutate(module = factor(as.character(module), levels = requested_modules))
  
  ggplot() +
    geom_line(
      data = observed,
      aes(
        x = pt_mid,
        y = score_std,
        colour = module,
        group = interaction(donor_id, module)
      ),
      linewidth = 0.32,
      alpha = 0.18,
      na.rm = TRUE
    ) +
    geom_ribbon(
      data = predicted,
      aes(x = pt_mid, ymin = lower, ymax = upper, fill = module),
      alpha = 0.16,
      colour = NA
    ) +
    geom_line(
      data = predicted,
      aes(x = pt_mid, y = fit, colour = module),
      linewidth = 1.05
    ) +
    scale_colour_manual(
      name = "Module",
      values = module_colours,
      breaks = requested_modules,
      labels = module_labels[requested_modules]
    ) +
    scale_fill_manual(
      name = "Module",
      values = module_colours,
      breaks = requested_modules,
      labels = module_labels[requested_modules]
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = panel_title,
      subtitle = panel_subtitle,
      x = "Within-lineage pseudotime (normalized)",
      y = "Standardized module activity"
    ) +
    theme_publication +
    theme(legend.position = "bottom")
}

# Mesenchymal 0.05 and 0.15 bins had only two and one retained donors,
# respectively. The displayed donor-supported range therefore starts at 0.25.
p_mes_curves <- make_curve_panel(
  lineage_name = "Mesenchymal",
  requested_modules = c("Layer1_Matricellular", "Layer2_ECM"),
  panel_title = "Mesenchymal trajectory",
  x_limits = c(0.25, 0.95),
  panel_subtitle = "Donor-balanced curves; sparse leading bins omitted"
)

p_mac_curves <- make_curve_panel(
  lineage_name = "Macrophage",
  requested_modules = c("Layer3_ImmuneIntegration", "MatureLAM"),
  panel_title = "Macrophage trajectory",
  x_limits = c(0.05, 0.95),
  panel_subtitle = "Donor-balanced curves"
)

# -----------------------------------------------------------------------------
# 4. Donor-level raw-data panels
# -----------------------------------------------------------------------------

p_mes_coupling <- ggplot(
  mes_coupling,
  aes(x = "Layer 1 vs Layer 2", y = rho_Layer1_vs_Layer2)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.55
  ) +
  geom_point(
    aes(fill = group),
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    position = position_jitter(width = 0.07, height = 0)
  ) +
  scale_fill_manual(name = "Clinical group", values = group_colours) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(
    title = "Layer 1-Layer 2 coupling",
    subtitle = "One correlation per donor",
    x = NULL,
    y = "Donor-level Spearman rho"
  ) +
  theme_publication +
  theme(legend.position = "bottom")

macro_rho_long <- macro_donor_metrics %>%
  select(
    group,
    donor_id,
    rho_Layer3_vs_pseudotime,
    rho_Layer3_vs_MatureLAM
  ) %>%
  pivot_longer(
    cols = starts_with("rho_"),
    names_to = "metric",
    values_to = "rho"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "rho_Layer3_vs_pseudotime",
        "rho_Layer3_vs_MatureLAM"
      ),
      labels = c(
        "Layer 3 vs\npseudotime",
        "Layer 3 vs\nMature-LAM"
      )
    )
  )

p_macro_rho <- ggplot(macro_rho_long, aes(x = metric, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.55
  ) +
  geom_point(
    aes(fill = group),
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  scale_fill_manual(name = "Clinical group", values = group_colours) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(
    title = "Within-donor immune-state alignment",
    subtitle = "One correlation per donor",
    x = NULL,
    y = "Donor-level Spearman rho"
  ) +
  theme_publication +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# 5. Supplementary panel and candidate replacement for Fig. 5G
# -----------------------------------------------------------------------------

fig5s_validation <- (
  p_mes_curves + p_mes_coupling +
    p_mac_curves + p_macro_rho
) +
  plot_layout(ncol = 2, guides = "collect", widths = c(1.55, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_figure(
  fig5s_validation,
  "Fig5S_donor_balanced_temporal_validation",
  width = 7.2,
  height = 6.7
)

fig5g_candidate <- (
  p_mes_curves + p_mac_curves
) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_figure(
  fig5g_candidate,
  "Fig5G_candidate_donor_balanced_lineage_curves",
  width = 7.2,
  height = 3.4
)

# -----------------------------------------------------------------------------
# 6. Donor-level ASYM/SYM sensitivity panel
# -----------------------------------------------------------------------------

donor_means_plot <- donor_means %>%
  mutate(
    module = factor(
      as.character(module),
      levels = names(module_labels),
      labels = unname(module_labels)
    ),
    group = factor(group, levels = c("ASYM", "SYM"))
  )

p_group_sensitivity <- ggplot(
  donor_means_plot,
  aes(x = group, y = mean_score, fill = group)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.62,
    linewidth = 0.55
  ) +
  geom_point(
    shape = 21,
    size = 2.25,
    stroke = 0.45,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  facet_wrap(~ module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colours) +
  labs(
    x = NULL,
    y = "Donor-aggregated module score\n(mean gene-wise z score)"
  ) +
  theme_publication +
  theme(legend.position = "none")

save_figure(
  p_group_sensitivity,
  "Fig5S_ASYM_SYM_donor_module_sensitivity",
  width = 6.3,
  height = 5.2
)

# -----------------------------------------------------------------------------
# 7. Result-aware caption and interpretation guardrails
# -----------------------------------------------------------------------------

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", format.pval(x, digits = 3, eps = 0.001))
}

get_corr <- function(pattern) {
  correlation_statistics %>%
    filter(grepl(pattern, analysis, fixed = TRUE)) %>%
    slice(1)
}

mes_corr <- get_corr("Layer1 versus Layer2")
pt_corr <- get_corr("Layer3 versus pseudotime")
lam_corr <- get_corr("Layer3 versus MatureLAM")
peak_stat <- primary_statistics %>%
  filter(grepl("Mesenchymal peak", analysis, fixed = TRUE)) %>%
  slice(1)

caption_lines <- c(
  "Supplementary Figure 5 | Donor-balanced validation of lineage-specific temporal module activity.",
  paste0(
    "(A) Layer 1 matricellular and Layer 2 ECM activity along the mesenchymal trajectory. ",
    "Thin lines represent donor-bin summaries; thick lines and ribbons show marginal GAM fits and 95% confidence intervals. ",
    "The sparsely supported 0.05 and 0.15 bins were omitted from display."
  ),
  paste0(
    "(B) Within-donor correlation between Layer 1 and Layer 2 module activity. Median Spearman rho = ",
    fmt_num(mes_corr$estimate_median), "; donor-bootstrap 95% CI ",
    fmt_num(mes_corr$ci95_lower), " to ", fmt_num(mes_corr$ci95_upper),
    "; two-sided Wilcoxon P = ", fmt_p(mes_corr$p_value),
    ", Holm q = ", fmt_p(mes_corr$holm_q_value),
    "; positive in ", mes_corr$positive_n, "/", mes_corr$n_donors, " donors. ",
    "This secondary analysis evaluates temporal coupling, not causal or sequential ordering."
  ),
  paste0(
    "(C) Layer 3 immune-integration and Mature-LAM module activity along macrophage pseudotime, ",
    "displayed using the same donor-balanced procedure."
  ),
  paste0(
    "(D) Within-donor correlations. For Layer 3 versus pseudotime, median rho = ",
    fmt_num(pt_corr$estimate_median), ", 95% CI ", fmt_num(pt_corr$ci95_lower),
    " to ", fmt_num(pt_corr$ci95_upper), ", P = ", fmt_p(pt_corr$p_value),
    ", Holm q = ", fmt_p(pt_corr$holm_q_value), "; positive in ",
    pt_corr$positive_n, "/", pt_corr$n_donors, " donors. For Layer 3 versus Mature-LAM, ",
    "median rho = ", fmt_num(lam_corr$estimate_median), ", 95% CI ",
    fmt_num(lam_corr$ci95_lower), " to ", fmt_num(lam_corr$ci95_upper),
    ", P = ", fmt_p(lam_corr$p_value), ", Holm q = ",
    fmt_p(lam_corr$holm_q_value), "; positive in ", lam_corr$positive_n,
    "/", lam_corr$n_donors, " donors."
  ),
  paste0(
    "Layer 1 and Layer 2 peak positions were not temporally separable (median difference = ",
    fmt_num(peak_stat$estimate_median), "; P = ", fmt_p(peak_stat$p_value),
    "), and most peaks occurred at the observed boundary. Accordingly, the modules are interpreted as ",
    "overlapping components of an early structural phase rather than strictly ordered steps."
  ),
  paste0(
    "Pseudotime was normalized independently within each lineage; mesenchymal and macrophage x-axis ",
    "positions are therefore not interpreted as a shared absolute clock. Donor was the unit of inference."
  )
)

writeLines(
  caption_lines,
  con = file.path(output_dir, "Fig5S_caption_draft.txt"),
  useBytes = TRUE
)

guardrail_lines <- c(
  "FIGURE 5 INTERPRETATION RULES BASED ON THE OBSERVED DONOR-LEVEL RESULTS",
  "",
  paste0(
    "1. Retain the three-layer functional architecture: matricellular, ECM, and immune-integration modules remain biologically and topologically distinct."
  ),
  paste0(
    "2. Do not claim a resolved Layer 1-to-Layer 2 temporal sequence. Their peak difference was centered at zero and 21/24 mesenchymal module peaks were boundary peaks."
  ),
  paste0(
    "3. Describe Layer 1 and Layer 2 as overlapping or coupled components of an early structural phase. Their activities were positively correlated in 12/12 donors."
  ),
  paste0(
    "4. Layer 3 increased along macrophage pseudotime in 12/12 donors and covaried with the Mature-LAM program in 12/12 donors. It is defensible to describe immune integration as late within macrophage maturation."
  ),
  paste0(
    "5. Do not claim a cross-lineage absolute temporal axis: the two pseudotimes were inferred and normalized separately."
  ),
  paste0(
    "6. The ASYM/SYM donor-module plot is a sensitivity analysis. Layer 2 showed a positive trend, whereas Layer 3 donor means did not increase in SYM."
  ),
  paste0(
    "7. Group-pooled CellChat strengths remain descriptive unless CellChat is recomputed as one network per donor."
  )
)

writeLines(
  guardrail_lines,
  con = file.path(output_dir, "Fig5S_interpretation_guardrails.txt"),
  useBytes = TRUE
)

write_csv(
  macro_rho_long,
  file.path(output_dir, "Fig5S_plotdata_macro_donor_rho.csv")
)
write_csv(
  group_statistics,
  file.path(output_dir, "Fig5S_group_statistics_for_legend.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "Fig5S_step3_sessionInfo.txt")
)

message("Step 3 completed. Open these first:")
message(file.path(output_dir, "Fig5S_donor_balanced_temporal_validation.pdf"))
message(file.path(output_dir, "Fig5G_candidate_donor_balanced_lineage_curves.pdf"))
message(file.path(output_dir, "Fig5S_interpretation_guardrails.txt"))



--------------------------------------------------------------------------------
  # FigureE — 时序可视化 v3 
suppressMessages({
  library(Seurat); library(tidyverse); library(ggplot2); library(ggridges); library(patchwork)
})

setwd("G:/13sample")
out_dir <- "G:/13sample/Fig5_temporal/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 三模块基因 ─────────────────
mod_initiate <- c("POSTN","TNC","COMP","THBS1","THBS2","SPARC","CCN2","CCN1","CTGF","CYR61")
mod_execute  <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL6A1","COL6A2","COL6A3",
                  "FN1","LAMA4","LAMB1","LAMC1","SPP1")
mod_integrate<- c("HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
                  "CD74","APP","CXCL9","CXCL10","CXCL12","CXCL16")
modules_all <- list(Module1_Initiate=mod_initiate, Module2_Execute=mod_execute, Module3_Integrate=mod_integrate)
mod_cols <- c("Module1_Initiate"="#C0392B","Module2_Execute"="#E8772E","Module3_Integrate"="#3A6FA0")
mod_labs <- c("Module1_Initiate"="Initiate (matricellular)",
              "Module2_Execute"="Execute (ECM)",
              "Module3_Integrate"="Integrate (immune)")

# ── 提取模块score + pseudotime, 并 min-max 标准化 ────────────────
get_score_df <- function(obj, pt_col, show_modules) {
  mods <- lapply(modules_all, function(g) intersect(g, rownames(obj)))
  obj <- AddModuleScore(obj, features=mods, name="Mod_", assay="RNA")
  sc <- paste0("Mod_", seq_along(mods))
  df <- obj@meta.data[, c(pt_col, sc)]; colnames(df) <- c("pseudotime", names(mods))
  df <- df[is.finite(df$pseudotime), c("pseudotime", show_modules)]
  # 标准化 pseudotime 到 0-1 (便于两谱系对齐展示)
  df$pt_norm <- (df$pseudotime - min(df$pseudotime))/(max(df$pseudotime)-min(df$pseudotime))
  df
}

# min-max 标准化每个模块的score (让幅度可比)
norm01 <- function(x){ r<-range(x,na.rm=TRUE); (x-r[1])/(r[2]-r[1]+1e-9) }

# ── 绘图函数: 标准化曲线 + 达峰竖线 ─────────────────────────────
plot_temporal_clean <- function(df, show_modules, title, subtitle) {
  long <- df %>% dplyr::select(pt_norm, all_of(show_modules)) %>%
    pivot_longer(-pt_norm, names_to="module", values_to="score")
  # 每模块标准化
  long <- long %>% group_by(module) %>% mutate(score_n = norm01(score)) %>% ungroup()
  long$module <- factor(long$module, levels=show_modules)
  
  # loess拟合求达峰
  peaks <- data.frame()
  for (m in show_modules) {
    sub <- long %>% filter(module==m) %>% arrange(pt_norm)
    lo <- loess(score_n ~ pt_norm, data=sub, span=0.5)
    px <- sub$pt_norm[which.max(predict(lo))]
    peaks <- rbind(peaks, data.frame(module=m, peak=px))
  }
  peaks$module <- factor(peaks$module, levels=show_modules)
  
  ggplot(long, aes(pt_norm, score_n, color=module, fill=module)) +
    geom_smooth(method="loess", span=0.5, se=TRUE, linewidth=1.4, alpha=0.13) +
    geom_vline(data=peaks, aes(xintercept=peak, color=module),
               linetype="dashed", linewidth=0.7, show.legend=FALSE) +
    geom_point(data=peaks, aes(x=peak, y=1.02, color=module), size=3, show.legend=FALSE) +
    scale_color_manual(values=mod_cols[show_modules], labels=mod_labs[show_modules], name=NULL) +
    scale_fill_manual(values=mod_cols[show_modules], guide="none") +
    scale_y_continuous(limits=c(-0.05,1.08), expand=c(0,0)) +
    labs(title=title, subtitle=subtitle, x="Pseudotime (normalized)", y="Module activity (scaled)") +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold",size=12.5),
          plot.subtitle=element_text(size=9,color="grey40"),
          axis.text=element_text(color="black"), axis.line=element_line(linewidth=0.4),
          legend.position="top", legend.text=element_text(size=9))
}

# ============================================================================
# 谱系1: 间质 — 展示 Initiate & Execute (早期)
# ============================================================================
mes <- readRDS("G:/13sample/Fig3_v3/fib_prf_v3_pseudotime.rds")  
pt_mes <- "pseudotime"                         
df_mes <- get_score_df(mes, pt_mes, c("Module1_Initiate","Module2_Execute"))
p_mes <- plot_temporal_clean(df_mes, c("Module1_Initiate","Module2_Execute"),
                             "Early structural-repair phase (mesenchymal lineage)",
                             "Matricellular initiation and ECM execution peak early and are tightly coupled")

# ============================================================================
# 谱系2: 巨噬 — 展示 Execute & Integrate (免疫后)
# ============================================================================
macro <- readRDS("macro/Fig4E/macro_with_pseudotime_3grp.rds")
pt_macro <- grep("pseudotime|Pseudotime|monocle", colnames(macro@meta.data), value=TRUE)[1]
cat("巨噬pseudotime列:", pt_macro, "\n")
df_mac <- get_score_df(macro, pt_macro, c("Module2_Execute","Module3_Integrate"))
p_mac <- plot_temporal_clean(df_mac, c("Module2_Execute","Module3_Integrate"),
                             "Late immune-integration phase (macrophage lineage)",
                             "Immune integration peaks markedly later, at TREM2\u207a LAM maturation")

# ── 组合图 (两谱系并排, 讲清两阶段) ──────────────────────────────
p_combined <- p_mes + p_mac + plot_layout(ncol=2) +
  plot_annotation(title="Two-phase temporal organization of the healing program",
                  theme=theme(plot.title=element_text(face="bold",size=14)))
ggsave(paste0(out_dir,"temporal_clean_combined.pdf"), p_combined, width=13, height=5.5, device=cairo_pdf)
ggsave(paste0(out_dir,"temporal_clean_combined.png"), p_combined, width=13, height=5.5, dpi=600)

# ============================================================================
# 山脊图 (各模块"高活性细胞"沿pseudotime的密度分布)
# ============================================================================
make_ridge_df <- function(obj, pt_col, lineage){
  mods <- lapply(modules_all, function(g) intersect(g, rownames(obj)))
  obj <- AddModuleScore(obj, features=mods, name="Mod_", assay="RNA")
  sc <- paste0("Mod_", seq_along(mods))
  df <- obj@meta.data[, c(pt_col, sc)]; colnames(df) <- c("pseudotime", names(mods))
  df <- df[is.finite(df$pseudotime),]
  df$pt_norm <- (df$pseudotime-min(df$pseudotime))/(max(df$pseudotime)-min(df$pseudotime))
  df$lineage <- lineage
  df
}
# 各模块取"活性最高的前40%细胞", 看它们的pseudotime分布
ridge_long <- bind_rows(
  make_ridge_df(mes, pt_mes, "mes") %>% pivot_longer(c(Module1_Initiate,Module2_Execute),
                                                     names_to="module", values_to="score"),
  make_ridge_df(macro, pt_macro, "mac") %>% pivot_longer(c(Module2_Execute,Module3_Integrate),
                                                         names_to="module", values_to="score")
) %>% group_by(lineage,module) %>%
  filter(score > quantile(score, 0.6)) %>% ungroup()
ridge_long$module <- factor(ridge_long$module,
                            levels=c("Module1_Initiate","Module2_Execute","Module3_Integrate"))

p_ridge <- ggplot(ridge_long, aes(x=pt_norm, y=module, fill=module)) +
  geom_density_ridges(alpha=0.75, scale=1.3, color="white", linewidth=0.4) +
  scale_fill_manual(values=mod_cols, labels=mod_labs, name=NULL) +
  scale_y_discrete(labels=mod_labs) +
  labs(title="Pseudotime distribution of high-activity cells per module",
       subtitle="Initiate/Execute cells concentrate early; Integrate cells concentrate late",
       x="Pseudotime (normalized)", y=NULL) +
  theme_ridges(font_size=11) +
  theme(plot.title=element_text(face="bold",size=12.5),
        plot.subtitle=element_text(size=9,color="grey40"), legend.position="none")

ggsave(paste0(out_dir,"temporal_ridge.pdf"), p_ridge, width=8, height=4.5, device=cairo_pdf)
ggsave(paste0(out_dir,"temporal_ridge.png"), p_ridge, width=8, height=4.5, dpi=600)

cat("\n完成. 输出:\n")
cat("  temporal_clean_combined — 两谱系标准化曲线+达峰标记(主推)\n")
cat("  temporal_ridge — 山脊图(备选, 最直观)\n")
print(p_combined); print(p_ridge)




