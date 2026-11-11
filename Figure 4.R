#Figure 4. Rare proliferating reparative fibroblasts (PRFs) represent a potential source of reparative fibroblast states

rm(list = ls())
gc()

suppressMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggpubr)
  library(ggsignif)
  library(scales)
})
library(ggrastr)
setwd("G:/13sample/")
out_dir <- "G:/13sample/Fig3_v3/"

# 加载对象（已含 pseudotime + Healing1）
fib_prf <- readRDS(paste0(out_dir, "fib_prf_v3_pseudotime.rds"))
cds     <- readRDS(paste0(out_dir, "fib_prf_v3_cds.rds"))

# ── 统一配色 ─────────────────────────────────────────────────────────
celltype_colors_v3 <- c(
  "Fibroblasts1" = "#D4A89AFF",   # ★ 改：浅粉棕（让出视觉焦点给 PRFs）
  "Fibroblasts2" = "#E2ED50FF",   # 黄
  "PRFs"         = "#524559FF"    # 深红 ★ 主角色
)

# 三组配色（保留供其他 panel 用）
group_colors_full <- c(
  "Normal" = "#8BEAFFFF",   # 绿 (与 Fig1/2 一致)
  "ASYM"   = "#FFE705FF",   # 黄
  "SYM"    = "#C71000FF"    # 红
)
# 仅 ASYM/SYM（用于 Panel E 趋势图）
group_colors_2 <- c("ASYM" = "#FFE705FF",
                    "SYM"  = "#C71000FF")

# 通用主题
theme_fig3 <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text     = element_text(color = "black"),
    axis.line     = element_line(linewidth = 0.4),
    axis.ticks    = element_line(linewidth = 0.4)
  )

--------------------------------------------------------------------------------
##-----------Fig A — tSNE：PRFs 高亮为红色亮点
tsne_df <- as.data.frame(Embeddings(fib_prf, "tsne"))
tsne_df$celltype <- fib_prf$celltype_v3

p_3A <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2)) +
  # 1) Fibroblasts1/2 底色 —— 栅格化
  rasterise(
    geom_point(data = subset(tsne_df, celltype != "PRFs"),
               aes(color = celltype),
               size = 0.4, alpha = 0.35),
    dpi = 600) +
  # 2) PRFs 高亮点 —— 也栅格化(点数少, 但统一栅格避免混合渲染)
  rasterise(
    geom_point(data = subset(tsne_df, celltype == "PRFs"),
               aes(fill = celltype),
               shape = 21, size = 2.2, color = "black", stroke = 0.4),
    dpi = 600) +
  scale_color_manual(values = celltype_colors_v3, name = NULL) +
  scale_fill_manual(values = celltype_colors_v3, guide = "none") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "A. Three reparative cell populations on tSNE",
       subtitle = sprintf("Fib1: %d  |  Fib2: %d  |  PRFs: %d (highlighted)",
                          sum(tsne_df$celltype == "Fibroblasts1"),
                          sum(tsne_df$celltype == "Fibroblasts2"),
                          sum(tsne_df$celltype == "PRFs")),
       x = "tSNE_1", y = "tSNE_2") +
  theme_fig3 +
  theme(legend.position = "right",
        legend.text = element_text(size = 10, face = "bold"))


ggsave(paste0(out_dir, "Fig3A_v2.pdf"), p_3A,
       width = 6, height = 4.5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig3A_v2.png"), p_3A,
       width = 6, height = 4.5, dpi = 600)


--------------------------------------------------------------------------------
##---------Fig B — Reparative markers：6 个最判别 marker + y 轴重排

# ★ 精简 marker（12 → 6 个）
reparative_markers <- c("POSTN", "TNC", "COL1A1", "COL3A1", "SPARC", "ACTA2")
reparative_markers <- intersect(reparative_markers, rownames(fib_prf))

# ★ y 轴顺序：PRFs 在最上（视觉重点 = 最强表达）
fib_prf$celltype_v3 <- factor(fib_prf$celltype_v3,
                              levels = c("Fibroblasts1",
                                         "Fibroblasts2",
                                         "PRFs"))
Idents(fib_prf) <- fib_prf$celltype_v3

p_3B <- DotPlot(fib_prf,
                features  = reparative_markers,
                dot.scale = 10,
                col.min   = 0,
                col.max   = 2.5) +
  scale_color_gradientn(
    colours = c("#FFFFFF", "#FCD2C2", "#E59085",
                "#C04040", "#8B0000"),
    limits  = c(0, 2.5), oob = scales::squish,
    name    = "Avg. expression") +
  scale_size_continuous(range = c(0.5, 10), name = "% expressed") +
  labs(title = "B. Reparative markers — PRFs express highest levels",
       subtitle = "Core reparative + myofibroblast markers",
       x = NULL, y = NULL) +
  theme_fig3 +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,
                                   face = "italic", size = 11,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 12,
                                   color = "black"),
        legend.position = "right",
        legend.title    = element_text(size = 9),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

ggsave(paste0(out_dir, "Fig3B_v2.pdf"), p_3B,
       width = 7.5, height = 4, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig3B_v2.png"), p_3B,
       width = 7.5, height = 4, dpi = 600)


##---------Fig B — Proliferation markers：删除 vln，只留单行 dot plot
proliferation_markers <- c("MKI67", "TOP2A", "UBE2C", "BIRC5",
                           "PCNA", "CDK1", "CCNB1", "STMN1")
proliferation_markers <- intersect(proliferation_markers, rownames(fib_prf))

p_3C <- DotPlot(fib_prf,
                features  = proliferation_markers,
                dot.scale = 10,
                col.min   = 0,
                col.max   = 2.5) +
  scale_color_gradientn(
    colours = c("#FFFFFF", "#D6E4F2", "#9DC3E6",
                "#3A6FA0", "#0A3D6B"),     # ★ 蓝色系，与 B 区分
    limits  = c(0, 2.5), oob = scales::squish,
    name    = "Avg. expression") +
  scale_size_continuous(range = c(0.5, 10), name = "% expressed") +
  labs(title = "C. Proliferation markers — PRFs uniquely cycle",
       subtitle = "MKI67 mean = 0.689 in PRFs vs. < 0.01 in Fib1/2",
       x = NULL, y = NULL) +
  theme_fig3 +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,
                                   face = "italic", size = 11,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 12,
                                   color = "black"),
        legend.position = "right",
        legend.title    = element_text(size = 9),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

ggsave(paste0(out_dir, "Fig3C_v2.pdf"), p_3C,
       width = 8, height = 4, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig3C_v2.png"), p_3C,
       width = 8, height = 4, dpi = 600)

--------------------------------------------------------------------------------
##------Fig C dotplot -------------------
celltype_colors_v3 <- c("Fibroblasts1"="#D4A89AFF",
                        "Fibroblasts2"="#E2ED50FF",
                        "PRFs"="#524559FF")

## ── 两个 module score: 与正文 marker 完全一致 ──────────────────────
reparative_genes    <- c("POSTN","TNC","COL1A1","COL3A1","SPARC")
proliferation_genes <- c("MKI67","TOP2A","UBE2C","BIRC5")

set.seed(2026)
fib_prf <- AddModuleScore(fib_prf, features = list(intersect(reparative_genes, rownames(fib_prf))),
                          name = "Reparative", ctrl = 50)
fib_prf <- AddModuleScore(fib_prf, features = list(intersect(proliferation_genes, rownames(fib_prf))),
                          name = "Proliferation", ctrl = 50)

df <- data.frame(
  Reparative    = fib_prf$Reparative1,
  Proliferation = fib_prf$Proliferation1,
  celltype      = factor(fib_prf$celltype_v3,
                         levels = c("Fibroblasts1","Fibroblasts2","PRFs"))
)
## PRFs 后画, 确保叠在最上层不被遮
df <- df[order(df$celltype == "PRFs"), ]

## 象限分隔线: 用各 score 的中位/0 作参考 (这里用 0, 因 ModuleScore 以 0 为中心)
xint <- 0; yint <- 0

p_3B <- ggplot(df, aes(Reparative, Proliferation)) +
  ## 散点栅格化 (两万多细胞, 必须栅格化)
  rasterise(
    geom_point(aes(color = celltype), size = 0.5, alpha = 0.55),
    dpi = 600) +
  ## PRFs 再描一层大点突出主角
  rasterise(
    geom_point(data = subset(df, celltype == "PRFs"),
               aes(fill = celltype), shape = 21, size = 1.8,
               color = "black", stroke = 0.3),
    dpi = 600) +
  geom_hline(yintercept = yint, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = xint, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  scale_color_manual(values = celltype_colors_v3, name = NULL) +
  scale_fill_manual(values = celltype_colors_v3, guide = "none") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  annotate("text", x = max(df$Reparative)*0.95, y = max(df$Proliferation)*0.95,
           label = "Reparative +\nProliferating", hjust = 1, vjust = 1,
           fontface = "bold", size = 3.3, color = "grey25") +
  labs(title = "B. PRFs uniquely co-express reparative and proliferation programs",
       subtitle = "Each dot = one cell; PRFs occupy the dual-high quadrant",
       x = "Reparative module score",
       y = "Proliferation module score") +
  theme_fig3 +
  theme(legend.position = "right",
        legend.text = element_text(size = 10, face = "bold"))

ggsave(paste0(out_dir, "Fig3B_dualscore.pdf"), p_3B, width = 6, height = 5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig3B_dualscore.png"), p_3B, width = 6, height = 5, dpi = 600)
print(p_3B)

## ── 顺手验证: 给正文用的 proliferation score 分亚群均值 ──────────────
df %>% group_by(celltype) %>%
  summarise(rep_med  = round(median(Reparative), 3),
            prol_med = round(median(Proliferation), 3),
            n = n()) %>% print()


------------------------------------------------------------------------------
  
##--------Fig D — Pseudotime：删 celltype 子图，箱线图按 PRFs→Fib2→Fib1 排
  
  # Fig. 3D replacement: donor-level pseudotime analysis
  #
  # Purpose
  #   Replace cell-level inference in the original Fig. 3D with donor-level
  #   summaries. Each point is one donor's median pseudotime within a cell
  #   population; cells are not treated as independent biological replicates.
  #
  # Statistical plan
  #   1. Restrict inference to plaque samples (ASYM and SYM).
  #   2. Require >= 5 cells per donor x population pseudobulk.
  #   3. Global comparison: Friedman test among donors containing all three
  #      populations.
  #   4. Planned pairwise comparisons: paired Wilcoxon signed-rank tests, using
  #      all donors available for each pair; Holm correction across three tests.
  #   5. Effect size: paired rank-biserial correlation. Also report a bootstrap
  #      95% CI for the mean within-donor pseudotime difference.
  #
  # Important interpretation
  #   The trajectory root in the upstream analysis was anchored to an MKI67+
  #   PRF cell. Therefore, an early PRF pseudotime is consistent with the chosen
  #   trajectory orientation, but it is not independent evidence that PRFs are
  #   progenitors.
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

# -----------------------------------------------------------------------------
# 1. Paths and analysis settings
# -----------------------------------------------------------------------------

input_rds <- "G:/13sample/Fig3/Fig3_v3/fib_prf_v3_pseudotime.rds"
output_dir <- file.path(getwd(), "Fig3D_donor_level_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_donor_population <- 5L
bootstrap_repetitions <- 10000L
set.seed(20260831)

population_levels <- c("PRFs", "Fibroblasts1", "Fibroblasts2")
population_colours <- c(
  "PRFs" = "#524559",
  "Fibroblasts1" = "#D4A89A",
  "Fibroblasts2" = "#E2ED50"
)

# -----------------------------------------------------------------------------
# 2. Load data and validate metadata
# -----------------------------------------------------------------------------

stopifnot(file.exists(input_rds))
fib_prf <- readRDS(input_rds)

required_columns <- c(
  "orig.ident", "group", "celltype_v3", "pseudotime"
)
missing_columns <- setdiff(required_columns, colnames(fib_prf@meta.data))
if (length(missing_columns) > 0L) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

cell_level <- fib_prf@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  transmute(
    cell_id = cell_id,
    donor_id = as.character(orig.ident),
    group = as.character(group),
    population = as.character(celltype_v3),
    pseudotime = as.numeric(pseudotime)
  ) %>%
  filter(
    group %in% c("ASYM", "SYM"),
    population %in% population_levels,
    is.finite(pseudotime)
  )

if (nrow(cell_level) == 0L) {
  stop("No eligible plaque cells were found after filtering.")
}

# A donor must map to only one clinical group.
donor_group_check <- cell_level %>%
  distinct(donor_id, group) %>%
  count(donor_id, name = "n_groups") %>%
  filter(n_groups != 1L)

if (nrow(donor_group_check) > 0L) {
  stop(
    "At least one donor maps to more than one clinical group: ",
    paste(donor_group_check$donor_id, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# 3. One pseudotime summary per donor x population
# -----------------------------------------------------------------------------

donor_population_all <- cell_level %>%
  group_by(donor_id, group, population) %>%
  summarise(
    n_cells = n(),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    mean_pseudotime = mean(pseudotime, na.rm = TRUE),
    q1_pseudotime = unname(quantile(pseudotime, 0.25, na.rm = TRUE)),
    q3_pseudotime = unname(quantile(pseudotime, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

excluded_low_cell_pseudobulks <- donor_population_all %>%
  filter(n_cells < min_cells_per_donor_population)

donor_population <- donor_population_all %>%
  filter(n_cells >= min_cells_per_donor_population) %>%
  mutate(population = factor(population, levels = population_levels))

write.csv(
  donor_population,
  file.path(output_dir, "Fig3D_donor_population_median_pseudotime.csv"),
  row.names = FALSE
)
write.csv(
  excluded_low_cell_pseudobulks,
  file.path(output_dir, "Fig3D_excluded_low_cell_pseudobulks.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Global donor-matched comparison: Friedman test
# -----------------------------------------------------------------------------

complete_donors <- donor_population %>%
  select(donor_id, population, median_pseudotime) %>%
  pivot_wider(names_from = population, values_from = median_pseudotime) %>%
  drop_na(all_of(population_levels)) %>%
  arrange(donor_id)

n_complete <- nrow(complete_donors)
k_populations <- length(population_levels)

if (n_complete < 2L) {
  stop("Fewer than two donors contain all three populations.")
}

friedman_matrix <- as.matrix(
  complete_donors[, population_levels, drop = FALSE]
)
friedman_result <- friedman.test(friedman_matrix)
kendalls_w <- unname(friedman_result$statistic) /
  (n_complete * (k_populations - 1L))

global_statistics <- data.frame(
  test = "Friedman test",
  analysis_unit = "donor median pseudotime",
  complete_donors = n_complete,
  populations = k_populations,
  statistic_Q = unname(friedman_result$statistic),
  df = unname(friedman_result$parameter),
  p_value = friedman_result$p.value,
  kendalls_W = kendalls_w
)

write.csv(
  global_statistics,
  file.path(output_dir, "Fig3D_global_friedman_test.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Planned donor-matched pairwise comparisons
# -----------------------------------------------------------------------------

bootstrap_mean_difference <- function(differences, repetitions = 10000L) {
  differences <- differences[is.finite(differences)]
  n <- length(differences)
  if (n < 2L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  bootstrap_means <- replicate(
    repetitions,
    mean(sample(differences, size = n, replace = TRUE))
  )
  unname(quantile(bootstrap_means, c(0.025, 0.975), na.rm = TRUE))
}

paired_wilcoxon_summary <- function(data, lower_population, higher_population) {
  paired_data <- data %>%
    filter(population %in% c(lower_population, higher_population)) %>%
    select(donor_id, population, median_pseudotime) %>%
    pivot_wider(names_from = population, values_from = median_pseudotime) %>%
    drop_na(all_of(c(lower_population, higher_population))) %>%
    arrange(donor_id)
  
  lower_values <- paired_data[[lower_population]]
  higher_values <- paired_data[[higher_population]]
  differences <- higher_values - lower_values
  n_pairs <- length(differences)
  
  if (n_pairs < 2L) {
    stop(
      "Fewer than two matched donors for ", lower_population,
      " versus ", higher_population, "."
    )
  }
  
  # R's exact signed-rank calculation is valid only without zero differences
  # or tied absolute differences. Otherwise use the asymptotic calculation.
  exact_ok <- all(differences != 0) &&
    !anyDuplicated(abs(differences))
  
  wilcoxon_result <- wilcox.test(
    higher_values,
    lower_values,
    paired = TRUE,
    alternative = "two.sided",
    exact = exact_ok,
    correct = !exact_ok,
    conf.int = FALSE
  )
  
  rank_sum_total <- n_pairs * (n_pairs + 1) / 2
  positive_rank_sum <- unname(wilcoxon_result$statistic)
  paired_rank_biserial <- 2 * positive_rank_sum / rank_sum_total - 1
  bootstrap_ci <- bootstrap_mean_difference(
    differences,
    repetitions = bootstrap_repetitions
  )
  
  data.frame(
    comparison = paste(higher_population, "minus", lower_population),
    lower_population = lower_population,
    higher_population = higher_population,
    paired_donors = n_pairs,
    wilcoxon_V = positive_rank_sum,
    p_value = wilcoxon_result$p.value,
    exact_test = exact_ok,
    paired_rank_biserial = paired_rank_biserial,
    mean_paired_difference = mean(differences),
    mean_difference_ci95_lower = bootstrap_ci[1],
    mean_difference_ci95_upper = bootstrap_ci[2],
    median_paired_difference = median(differences),
    stringsAsFactors = FALSE
  )
}

pairwise_statistics <- bind_rows(
  paired_wilcoxon_summary(donor_population, "PRFs", "Fibroblasts1"),
  paired_wilcoxon_summary(donor_population, "PRFs", "Fibroblasts2"),
  paired_wilcoxon_summary(donor_population, "Fibroblasts1", "Fibroblasts2")
) %>%
  mutate(holm_q_value = p.adjust(p_value, method = "holm"))

write.csv(
  pairwise_statistics,
  file.path(output_dir, "Fig3D_pairwise_paired_wilcoxon_tests.csv"),
  row.names = FALSE
)

population_summary <- donor_population %>%
  group_by(population) %>%
  summarise(
    donors = n(),
    median_of_donor_medians = median(median_pseudotime),
    mean_of_donor_medians = mean(median_pseudotime),
    sd_of_donor_medians = sd(median_pseudotime),
    .groups = "drop"
  )

write.csv(
  population_summary,
  file.path(output_dir, "Fig3D_population_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Publication-ready replacement panel
# -----------------------------------------------------------------------------

format_p <- function(p_value) {
  if (is.na(p_value)) return("NA")
  if (p_value < 0.001) return("<0.001")
  formatC(p_value, format = "f", digits = 3)
}

n_by_population <- donor_population %>%
  count(population, name = "n_donors") %>%
  complete(
    population = factor(population_levels, levels = population_levels),
    fill = list(n_donors = 0L)
  ) %>%
  arrange(population)

x_axis_labels <- setNames(
  paste0(as.character(n_by_population$population), "\n(n=", n_by_population$n_donors, " donors)"),
  as.character(n_by_population$population)
)

plot_range <- range(donor_population$median_pseudotime, na.rm = TRUE)
plot_span <- diff(plot_range)
if (!is.finite(plot_span) || plot_span == 0) plot_span <- 1

bracket_step <- 0.105 * plot_span
bracket_tip <- 0.025 * plot_span
bracket_base <- plot_range[2] + 0.08 * plot_span

annotation_data <- pairwise_statistics %>%
  mutate(
    x_start = c(1, 1, 2),
    x_end = c(2, 3, 3),
    y = bracket_base + c(0, 2, 1) * bracket_step,
    label = paste0(
      "Holm q ", vapply(holm_q_value, format_p, character(1)),
      "; n=", paired_donors
    )
  )

subtitle_text <- paste0(
  "Friedman P ", format_p(global_statistics$p_value),
  "; Kendall's W=", sprintf("%.2f", global_statistics$kendalls_W),
  "; complete donors n=", global_statistics$complete_donors
)

fig3d <- ggplot(
  donor_population,
  aes(x = population, y = median_pseudotime)
) +
  geom_line(
    aes(group = donor_id),
    colour = "grey70",
    linewidth = 0.55,
    alpha = 0.75,
    na.rm = TRUE
  ) +
  geom_boxplot(
    aes(fill = population),
    width = 0.50,
    outlier.shape = NA,
    linewidth = 0.75,
    alpha = 0.85
  ) +
  geom_point(
    aes(fill = population),
    shape = 21,
    size = 3.0,
    stroke = 0.7,
    colour = "black",
    position = position_jitter(width = 0.035, height = 0),
    show.legend = FALSE
  ) +
  geom_segment(
    data = annotation_data,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.55
  ) +
  geom_segment(
    data = annotation_data,
    aes(x = x_start, xend = x_start, y = y, yend = y - bracket_tip),
    inherit.aes = FALSE,
    linewidth = 0.55
  ) +
  geom_segment(
    data = annotation_data,
    aes(x = x_end, xend = x_end, y = y, yend = y - bracket_tip),
    inherit.aes = FALSE,
    linewidth = 0.55
  ) +
  geom_text(
    data = annotation_data,
    aes(x = (x_start + x_end) / 2, y = y + 0.018 * plot_span, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 0
  ) +
  scale_fill_manual(values = population_colours, drop = FALSE) +
  scale_x_discrete(labels = x_axis_labels, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.25))
  ) +
  labs(
    title = "D. Donor-level pseudotime ordering",
    subtitle = subtitle_text,
    x = NULL,
    y = "Median pseudotime per donor",
    caption = paste0(
      "Each point is one donor; grey lines connect matched donor populations. ",
      "Only donor-population summaries with >=",
      min_cells_per_donor_population,
      " cells are included. Trajectory root was anchored to an MKI67+ PRF cell."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10.5, colour = "grey30"),
    plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0),
    axis.text.x = element_text(size = 10.5, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(size = 11.5),
    legend.position = "none",
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  filename = file.path(output_dir, "Fig3D_donor_median_pseudotime.pdf"),
  plot = fig3d,
  width = 7.2,
  height = 6.0,
  units = "in",
  device = cairo_pdf
)
ggsave(
  filename = file.path(output_dir, "Fig3D_donor_median_pseudotime.png"),
  plot = fig3d,
  width = 7.2,
  height = 6.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 7. Save a concise methods/statistics record and session information
# -----------------------------------------------------------------------------

methods_text <- c(
  "Fig. 3D donor-level pseudotime analysis",
  "",
  paste0(
    "Pseudotime was summarized as the median for each donor within each ",
    "fibroblast population. Donor-population summaries containing fewer than ",
    min_cells_per_donor_population,
    " cells were excluded."
  ),
  paste0(
    "A Friedman test compared the three populations among donors represented ",
    "in all three populations (n=", n_complete, ")."
  ),
  paste0(
    "Planned pairwise comparisons used two-sided paired Wilcoxon signed-rank ",
    "tests with Holm correction across the three comparisons."
  ),
  paste0(
    "The paired rank-biserial correlation and a nonparametric bootstrap 95% ",
    "confidence interval for the mean within-donor difference (",
    bootstrap_repetitions,
    " bootstrap samples) were reported as effect-size summaries."
  ),
  "Cells were not treated as independent biological replicates.",
  "The trajectory root was anchored upstream to an MKI67+ PRF cell; therefore, trajectory direction is conditional on that rooting choice."
)

writeLines(
  methods_text,
  con = file.path(output_dir, "Fig3D_methods_text.txt"),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "Fig3D_sessionInfo.txt")
)

print(fig3d)
message("Fig. 3D donor-level analysis completed. Results: ", output_dir)


------------------------------------------------------------------------------------
##----------------- FigE: donor-level pseudobulk analysis
#
# Replaces cell-level Wilcoxon/Kruskal-Wallis tests in the original Fig. 3E.
# Raw RNA counts are summed within donor x fibroblast population, normalized by
# edgeR TMM, converted to log2 CPM, and summarized as a prespecified healing
# signature score. Donors, not cells, are the independent observations.

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

input_rds <- "G:/13sample/Fig3/Fig3_v3/fib_prf_v3_pseudotime.rds"

out_dir <- file.path(
  "G:/13sample/A 颈动脉单细胞文章/TREM2??? macrophage–reparative fibroblast healing axis",
  "Communications_Biology_submission/Communications_Biology_submission",
  "Fig3E_donor_pseudobulk_results"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260831)
bootstrap_iterations <- 10000L
minimum_cells_per_pseudobulk <- 5L

healing_signature <- c(
  "POSTN", "TNC", "LOXL2", "COMP", "ASPN",
  "CTGF", "SPARC", "THBS1", "MMP2", "TIMP1",
  "COL1A1", "COL3A1", "DCN", "LUM"
)

celltype_levels <- c("Fibroblasts1", "Fibroblasts2", "PRFs")

celltype_colors <- c(
  "Fibroblasts1" = "#D4A89AFF",
  "Fibroblasts2" = "#E2ED50FF",
  "PRFs"         = "#524559FF"
)

group_colors <- c(
  "ASYM" = "#F9C80EFF",
  "SYM"  = "#C71000FF"
)

theme_fig3 <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey35"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4)
  )

# ==============================================================================
# 2. Load object and validate biological-replicate metadata
# ==============================================================================

fib_prf <- readRDS(input_rds)

required_metadata <- c("orig.ident", "sample", "group", "celltype_v3")
missing_metadata <- setdiff(required_metadata, colnames(fib_prf@meta.data))
if (length(missing_metadata) > 0L) {
  stop("Missing metadata column(s): ", paste(missing_metadata, collapse = ", "))
}

if (!"RNA" %in% Assays(fib_prf)) {
  stop("RNA assay is absent; raw-count pseudobulk cannot be generated.")
}

rna_counts <- GetAssayData(fib_prf, assay = "RNA", layer = "counts")
if (nrow(rna_counts) == 0L || ncol(rna_counts) == 0L) {
  stop("RNA counts are empty; raw-count pseudobulk cannot be generated.")
}

md <- fib_prf@meta.data %>%
  mutate(
    cell_id = rownames(fib_prf@meta.data),
    donor_id = as.character(orig.ident),
    group = as.character(group),
    celltype_v3 = as.character(celltype_v3)
  ) %>%
  filter(
    group %in% c("ASYM", "SYM"),
    celltype_v3 %in% celltype_levels,
    !is.na(donor_id)
  )

# One donor must belong to exactly one disease group.
donor_group_check <- md %>% distinct(donor_id, group) %>% count(donor_id)
if (any(donor_group_check$n != 1L)) {
  stop("At least one donor is assigned to more than one ASYM/SYM group.")
}

# ==============================================================================
# 3. Construct donor x population raw-count pseudobulks
# ==============================================================================

md <- md %>%
  mutate(
    pseudobulk_id = paste(donor_id, group, celltype_v3, sep = "___"),
    pseudobulk_id = factor(pseudobulk_id, levels = unique(pseudobulk_id))
  )

# Match count-matrix columns exactly to metadata rows before aggregation.
if (!all(md$cell_id %in% colnames(rna_counts))) {
  stop("Some selected metadata cell IDs are absent from the RNA count matrix.")
}
rna_counts <- rna_counts[, md$cell_id, drop = FALSE]

design_pb <- Matrix::sparse.model.matrix(~ 0 + pseudobulk_id, data = md)
colnames(design_pb) <- levels(md$pseudobulk_id)

# genes x cells multiplied by cells x pseudobulk samples
pseudobulk_counts <- rna_counts %*% design_pb

pseudobulk_meta <- md %>%
  count(donor_id, group, celltype_v3, pseudobulk_id, name = "n_cells") %>%
  mutate(pseudobulk_id = as.character(pseudobulk_id))

pseudobulk_meta <- pseudobulk_meta[
  match(colnames(pseudobulk_counts), pseudobulk_meta$pseudobulk_id),
  , drop = FALSE
]

if (any(is.na(pseudobulk_meta$pseudobulk_id))) {
  stop("Pseudobulk count columns could not be matched to pseudobulk metadata.")
}

# A score cannot be estimated for a population absent from a donor. Existing
# pseudobulks with fewer than five cells are excluded prospectively and listed.
excluded_low_cell <- pseudobulk_meta %>%
  filter(n_cells < minimum_cells_per_pseudobulk)

keep_pb <- pseudobulk_meta$n_cells >= minimum_cells_per_pseudobulk
pseudobulk_counts <- pseudobulk_counts[, keep_pb, drop = FALSE]
pseudobulk_meta <- pseudobulk_meta[keep_pb, , drop = FALSE]

if (nrow(excluded_low_cell) > 0L) {
  warning(
    "Excluded pseudobulk samples with fewer than ",
    minimum_cells_per_pseudobulk, " cells: ",
    paste(excluded_low_cell$pseudobulk_id, collapse = ", ")
  )
}

# ==============================================================================
# 4. edgeR TMM normalization and pseudobulk healing score
# ==============================================================================

# Remove only genes with zero counts in every pseudobulk sample. No DEG-driven
# feature selection is used for the prespecified signature score.
nonzero_gene <- Matrix::rowSums(pseudobulk_counts) > 0
pseudobulk_counts <- pseudobulk_counts[nonzero_gene, , drop = FALSE]

dge <- edgeR::DGEList(counts = pseudobulk_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

signature_present <- intersect(healing_signature, rownames(log_cpm))
signature_missing <- setdiff(healing_signature, signature_present)

if (length(signature_present) < 5L) {
  stop(
    "Fewer than five healing-signature genes are present in pseudobulk data: ",
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

# Each prespecified gene contributes equally: z transform each gene across
# pseudobulk samples, then average the gene z scores within each sample.
signature_gene_z <- t(scale(t(signature_log_cpm)))
healing_score <- colMeans(signature_gene_z, na.rm = TRUE)

pseudobulk_meta <- pseudobulk_meta %>%
  mutate(
    healing_pb_score = healing_score[pseudobulk_id],
    group = factor(group, levels = c("ASYM", "SYM")),
    celltype_v3 = factor(celltype_v3, levels = celltype_levels)
  )

if (any(!is.finite(pseudobulk_meta$healing_pb_score))) {
  stop("At least one pseudobulk healing score is non-finite.")
}

# ==============================================================================
# 5. Statistical helpers: exact donor-level tests, effects and bootstrap CIs
# ==============================================================================

bootstrap_mean_difference <- function(x1, x0, paired = FALSE,
                                      iterations = bootstrap_iterations) {
  if (paired) {
    stopifnot(length(x1) == length(x0))
    delta <- x1 - x0
    boot <- replicate(
      iterations,
      mean(sample(delta, length(delta), replace = TRUE))
    )
    estimate <- mean(delta)
  } else {
    boot <- replicate(
      iterations,
      mean(sample(x1, length(x1), replace = TRUE)) -
        mean(sample(x0, length(x0), replace = TRUE))
    )
    estimate <- mean(x1) - mean(x0)
  }
  
  ci <- unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
  c(estimate = estimate, ci_low = ci[1], ci_high = ci[2])
}

paired_population_test <- function(data, level0, level1) {
  wide <- data %>%
    filter(celltype_v3 %in% c(level0, level1)) %>%
    select(donor_id, celltype_v3, healing_pb_score) %>%
    pivot_wider(names_from = celltype_v3, values_from = healing_pb_score) %>%
    drop_na(all_of(c(level0, level1))) %>%
    arrange(donor_id)
  
  x0 <- wide[[level0]]
  x1 <- wide[[level1]]
  delta <- x1 - x0
  exact_ok <- all(delta != 0) && !anyDuplicated(abs(delta))
  
  fit <- wilcox.test(
    x1, x0, paired = TRUE, alternative = "two.sided",
    exact = exact_ok, correct = !exact_ok
  )
  
  nz <- delta != 0
  ranks <- rank(abs(delta[nz]), ties.method = "average")
  rank_biserial <- if (any(nz)) {
    (sum(ranks[delta[nz] > 0]) - sum(ranks[delta[nz] < 0])) / sum(ranks)
  } else {
    0
  }
  
  boot <- bootstrap_mean_difference(x1, x0, paired = TRUE)
  
  tibble(
    comparison = paste0(level1, " - ", level0),
    level0 = level0,
    level1 = level1,
    n_paired_donors = length(x0),
    test = ifelse(exact_ok, "Exact paired Wilcoxon signed-rank",
                  "Paired Wilcoxon signed-rank with continuity correction"),
    statistic = unname(fit$statistic),
    p_value = fit$p.value,
    rank_biserial = rank_biserial,
    mean_difference = unname(boot["estimate"]),
    mean_difference_ci_low = unname(boot["ci_low"]),
    mean_difference_ci_high = unname(boot["ci_high"])
  )
}

independent_group_test <- function(data, population) {
  d <- data %>% filter(celltype_v3 == population)
  x0 <- d$healing_pb_score[d$group == "ASYM"]
  x1 <- d$healing_pb_score[d$group == "SYM"]
  exact_ok <- !anyDuplicated(c(x0, x1))
  
  fit <- wilcox.test(
    x1, x0, paired = FALSE, alternative = "two.sided",
    exact = exact_ok, correct = !exact_ok
  )
  
  u <- unname(fit$statistic)
  rank_biserial <- 2 * u / (length(x1) * length(x0)) - 1
  boot <- bootstrap_mean_difference(x1, x0, paired = FALSE)
  
  tibble(
    celltype_v3 = population,
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

# ==============================================================================
# 6. FigE1: within-donor comparison across fibroblast populations
# ==============================================================================

complete_e1 <- pseudobulk_meta %>%
  select(donor_id, celltype_v3, healing_pb_score) %>%
  pivot_wider(names_from = celltype_v3, values_from = healing_pb_score) %>%
  drop_na(all_of(celltype_levels))

if (nrow(complete_e1) < 3L) {
  stop("Fewer than three donors have all three fibroblast populations.")
}

friedman_matrix <- as.matrix(complete_e1[, celltype_levels, drop = FALSE])
friedman_fit <- friedman.test(friedman_matrix)
kendalls_w <- unname(friedman_fit$statistic) /
  (nrow(friedman_matrix) * (ncol(friedman_matrix) - 1))

e1_pairs <- list(
  c("Fibroblasts1", "Fibroblasts2"),
  c("Fibroblasts1", "PRFs"),
  c("Fibroblasts2", "PRFs")
)

e1_pairwise <- bind_rows(lapply(
  e1_pairs,
  function(x) paired_population_test(pseudobulk_meta, x[1], x[2])
)) %>%
  mutate(p_holm_three_pairs = p.adjust(p_value, method = "holm"))

e1_global <- tibble(
  test = "Friedman test on donors with all three populations",
  n_complete_donors = nrow(friedman_matrix),
  statistic = unname(friedman_fit$statistic),
  df = unname(friedman_fit$parameter),
  p_value = friedman_fit$p.value,
  kendalls_w = kendalls_w
)

# ==============================================================================
# 7. FigE2: ASYM versus SYM within each population
# ==============================================================================

e2_tests <- bind_rows(lapply(
  celltype_levels,
  function(pop) independent_group_test(pseudobulk_meta, pop)
)) %>%
  mutate(p_holm_three_populations = p.adjust(p_value, method = "holm"))

descriptive_stats <- pseudobulk_meta %>%
  group_by(celltype_v3, group) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    total_cells = sum(n_cells),
    mean_score = mean(healing_pb_score),
    sd_score = sd(healing_pb_score),
    median_score = median(healing_pb_score),
    q1_score = unname(quantile(healing_pb_score, 0.25)),
    q3_score = unname(quantile(healing_pb_score, 0.75)),
    .groups = "drop"
  )

# ==============================================================================
# 8. Publication figure: all plotted observations are donor pseudobulks
# ==============================================================================

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
}

# E1: donor trajectories across populations. Donors without detectable PRFs are
# not imputed; their Fibroblasts1/Fibroblasts2 values remain visible.
e1_range <- range(pseudobulk_meta$healing_pb_score, na.rm = TRUE)
e1_span <- diff(e1_range)
if (!is.finite(e1_span) || e1_span == 0) e1_span <- 1

e1_brackets <- e1_pairwise %>%
  mutate(
    x_start = match(level0, celltype_levels),
    x_end = match(level1, celltype_levels),
    y = e1_range[2] + c(0.12, 0.30, 0.48) * e1_span,
    tick = 0.035 * e1_span,
    label = paste0("Holm q=", vapply(p_holm_three_pairs, format_p, character(1)))
  )

p_3E1 <- ggplot(
  pseudobulk_meta,
  aes(x = celltype_v3, y = healing_pb_score, group = donor_id)
) +
  geom_boxplot(
    aes(fill = celltype_v3, group = celltype_v3),
    width = 0.52, outlier.shape = NA, alpha = 0.42,
    linewidth = 0.5, color = "black"
  ) +
  geom_line(color = "grey60", linewidth = 0.45, alpha = 0.65) +
  geom_point(
    aes(fill = celltype_v3),
    shape = 21, size = 2.5, stroke = 0.45, color = "black", alpha = 0.95
  ) +
  geom_segment(
    data = e1_brackets,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.4
  ) +
  geom_segment(
    data = e1_brackets,
    aes(x = x_start, xend = x_start, y = y, yend = y - tick),
    inherit.aes = FALSE, linewidth = 0.4
  ) +
  geom_segment(
    data = e1_brackets,
    aes(x = x_end, xend = x_end, y = y, yend = y - tick),
    inherit.aes = FALSE, linewidth = 0.4
  ) +
  geom_text(
    data = e1_brackets,
    aes(x = (x_start + x_end) / 2, y = y + 0.035 * e1_span, label = label),
    inherit.aes = FALSE, size = 2.8
  ) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.55))) +
  labs(
    title = "E1. Donor-pseudobulk healing score across populations",
    subtitle = paste0(
      "Friedman P=", format_p(e1_global$p_value),
      "; Kendall's W=", formatC(e1_global$kendalls_w, digits = 2, format = "f"),
      "; complete donors n=", e1_global$n_complete_donors
    ),
    x = NULL,
    y = "Healing signature score\n(TMM logCPM, gene-wise z mean)"
  ) +
  theme_fig3 +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 10)
  )

# E2: independent donor comparisons, separately within each population.
e2_annotation <- pseudobulk_meta %>%
  group_by(celltype_v3) %>%
  summarise(
    y_max = max(healing_pb_score),
    y_min = min(healing_pb_score),
    .groups = "drop"
  ) %>%
  left_join(e2_tests, by = "celltype_v3") %>%
  mutate(
    span = pmax(y_max - y_min, 0.2),
    y = y_max + 0.18 * span,
    label = paste0(
      "n=", n_ASYM, " vs ", n_SYM,
      "\nP=", vapply(p_value, format_p, character(1)),
      "; q=", vapply(p_holm_three_populations, format_p, character(1))
    )
  )

p_3E2 <- ggplot(
  pseudobulk_meta,
  aes(x = group, y = healing_pb_score, fill = group)
) +
  geom_boxplot(
    width = 0.54, outlier.shape = NA,
    alpha = 0.48, linewidth = 0.5
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 20260831),
    shape = 21, size = 2.6, color = "black", stroke = 0.45
  ) +
  geom_text(
    data = e2_annotation,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE, size = 2.75, lineheight = 0.95
  ) +
  facet_wrap(~ celltype_v3, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.30))) +
  labs(
    title = "E2. Donor-pseudobulk healing score: ASYM versus SYM",
    subtitle = "Exact donor-level Wilcoxon tests; q = Holm adjustment across three populations",
    x = NULL,
    y = "Healing signature score"
  ) +
  theme_fig3 +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(face = "bold", size = 9),
    panel.spacing.x = unit(0.65, "lines")
  )

p_3E <- p_3E1 + p_3E2 + plot_layout(widths = c(1.0, 1.45))

# ==============================================================================
# 9. Export figure, source data, statistics and methods text
# ==============================================================================

ggsave(
  file.path(out_dir, "Fig3E_donor_pseudobulk.pdf"),
  p_3E, width = 12.2, height = 5.0, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)

ggsave(
  file.path(out_dir, "Fig3E_donor_pseudobulk.png"),
  p_3E, width = 12.2, height = 5.0, units = "in",
  dpi = 600, bg = "white"
)

write.csv(
  pseudobulk_meta,
  file.path(out_dir, "Fig3E_donor_population_pseudobulk_scores.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  as.data.frame(signature_log_cpm),
  file.path(out_dir, "Fig3E_signature_gene_logCPM.csv"),
  row.names = TRUE, fileEncoding = "UTF-8"
)

write.csv(
  e1_global,
  file.path(out_dir, "Fig3E1_Friedman_global_test.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  e1_pairwise,
  file.path(out_dir, "Fig3E1_paired_Wilcoxon_tests.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  e2_tests,
  file.path(out_dir, "Fig3E2_ASYM_vs_SYM_Wilcoxon_tests.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  descriptive_stats,
  file.path(out_dir, "Fig3E_descriptive_statistics.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  excluded_low_cell,
  file.path(out_dir, "Fig3E_excluded_low_cell_pseudobulks.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

methods_text <- c(
  "Fig. 3E donor-pseudobulk analysis",
  "",
  paste0(
    "Raw RNA counts were summed separately for each donor and fibroblast ",
    "population. Pseudobulk libraries containing fewer than five cells were ",
    "excluded without imputing scores for populations absent from a donor."
  ),
  paste0(
    "Pseudobulk counts were normalized using the trimmed mean of M-values ",
    "method in edgeR and transformed to log2 counts per million with a prior ",
    "count of 2."
  ),
  paste0(
    "For the prespecified healing signature, each gene was standardized across ",
    "pseudobulk samples and the mean gene z score was calculated per sample."
  ),
  paste0(
    "Population differences were tested using a Friedman test among donors ",
    "with all three populations, followed by planned paired Wilcoxon signed-rank ",
    "tests with Holm correction."
  ),
  paste0(
    "Within each population, ASYM and SYM donors were compared using two-sided ",
    "exact Wilcoxon rank-sum tests with Holm correction across the three ",
    "populations."
  ),
  paste0(
    "Rank-biserial correlations and bootstrap 95% confidence intervals for ",
    "mean donor-level differences were reported. No cell was treated as an ",
    "independent biological replicate."
  ),
  paste0("Signature genes present: ", paste(rownames(signature_log_cpm), collapse = ", ")),
  paste0("Signature genes absent: ", ifelse(length(signature_missing) == 0,
                                            "none", paste(signature_missing, collapse = ", ")))
)

writeLines(
  methods_text,
  file.path(out_dir, "Fig3E_methods_text.txt"),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "Fig3E_R_sessionInfo.txt")
)

message("Completed Fig. 3E donor-pseudobulk analysis.")
message("Results saved to: ", out_dir)
print(e1_global)
print(e1_pairwise)
print(e2_tests)


------------------------------------------------------------------------
##-----------FigF — 单 panel + per 1,000 cells 单位
# Panel F final — Bar plot, ASYM vs SYM 
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsignif)
  library(patchwork)
})

# ── 配色（与 D/E 一致）────────────────────────────────────────────
group_colors_2 <- c(
  "ASYM" = "#F9C80EFF",   # 黄
  "SYM"  = "#C71000FF"    # 红
)

theme_fig3 <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text     = element_text(color = "black"),
    axis.line     = element_line(linewidth = 0.4),
    axis.ticks    = element_line(linewidth = 0.4)
  )

# ── 数据准备 ─────────────────────────────────────────────────────────
sample_col <- if ("orig.ident" %in% colnames(fib_prf@meta.data)) "orig.ident" else "sample"

prop_long <- fib_prf@meta.data %>%
  group_by(!!sym(sample_col), group, celltype_v3) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(!!sym(sample_col), group) %>%
  mutate(per_1000 = n_cells / sum(n_cells) * 1000) %>%
  ungroup()
colnames(prop_long)[1] <- "sample"

# ★ 删除 Normal，只保留 ASYM/SYM
prop_long <- prop_long %>%
  filter(group %in% c("ASYM", "SYM")) %>%
  mutate(group       = factor(group, levels = c("ASYM", "SYM")),
         celltype_v3 = factor(celltype_v3,
                              levels = c("Fibroblasts1",
                                         "Fibroblasts2",
                                         "PRFs")))

# 均值 + SE（每个 celltype × group）
summary_df <- prop_long %>%
  group_by(celltype_v3, group) %>%
  summarise(
    mean_val = mean(per_1000),
    se_val   = sd(per_1000) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  )
print(summary_df)

# P 值（每个亚群 ASYM vs SYM）
p_vals_F <- prop_long %>%
  group_by(celltype_v3) %>%
  summarise(
    p_val = tryCatch(wilcox.test(per_1000 ~ group)$p.value,
                     error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    is.na(p_val)        ~ "n.s.",
    p_val < 0.001       ~ sprintf("P = %.1e ***", p_val),
    p_val < 0.01        ~ sprintf("P = %.3f **",  p_val),
    p_val < 0.05        ~ sprintf("P = %.3f *",   p_val),
    p_val < 0.10        ~ sprintf("P = %.3f (trend)", p_val),
    TRUE                ~ sprintf("P = %.3f n.s.", p_val)
  ))
print(p_vals_F)

# ── 绘图：柱状图 + 误差棒 + 个体点 + P 值 bracket ────────────────────
p_3F <- ggplot() +
  # 1) 柱：组均值（与第二张图同款）
  geom_col(data = summary_df,
           aes(x = group, y = mean_val, fill = group),
           width = 0.65, color = "black", linewidth = 0.4,
           alpha = 0.9) +
  # 2) 误差棒（mean ± SE）
  geom_errorbar(data = summary_df,
                aes(x = group,
                    ymin = pmax(0, mean_val - se_val),
                    ymax = mean_val + se_val),
                width = 0.18, linewidth = 0.5,
                color = "black") +
  # 3) 个体散点（每个样本一个点）
  geom_jitter(data = prop_long,
              aes(x = group, y = per_1000),
              width = 0.12, height = 0,
              size = 2.2, shape = 21,
              fill = "white", color = "grey20",
              stroke = 0.5, alpha = 0.95) +
  # 4) P 值 bracket（ASYM vs SYM）
  geom_signif(
    data = prop_long,
    aes(x = group, y = per_1000),
    comparisons      = list(c("ASYM", "SYM")),
    test             = "wilcox.test",
    map_signif_level = function(p) {
      if (is.na(p)) "n.s."
      else if (p < 0.001) sprintf("P = %.1e", p)
      else sprintf("P = %.3f", p)
    },
    textsize    = 3.0,
    tip_length  = 0.015,
    vjust       = -0.1,
    color       = "black",
    inherit.aes = FALSE) +
  # 5) facet by celltype
  facet_wrap(~ celltype_v3, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors_2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "F. Reparative compartment composition: ASYM vs SYM",
       subtitle = "Fibroblasts1 trends upward in SYM; Fibroblasts2 and PRFs remain stable",
       x = NULL, y = "Cells per 1,000 reparative cells") +
  theme_fig3 +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey95", color = "grey80",
                                    linewidth = 0.3),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(face = "bold", size = 11),
    panel.spacing.x  = unit(0.8, "lines")
  )

ggsave(paste0(out_dir, "Fig3F_final.pdf"), p_3F,
       width = 9, height = 4.2, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig3F_final.png"), p_3F,
       width = 9, height = 4.2, dpi = 600)

cat("\n??? Fig 3F final version saved\n")

--------------------------------------------------------------------------------
  
##------------Fig G PRFs absolute count — Bar plot, ASYM vs SYM 

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsignif)
})

# ── 配色（与 Panel F 一致）──────────────────────────────────────────
group_colors_2 <- c(
  "ASYM" = "#F9C80EFF",   # 黄
  "SYM"  = "#C71000FF"    # 红
)

theme_fig3 <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text     = element_text(color = "black"),
    axis.line     = element_line(linewidth = 0.4),
    axis.ticks    = element_line(linewidth = 0.4)
  )

# ── 数据准备 ─────────────────────────────────────────────────────────
sample_col <- if ("orig.ident" %in% colnames(fib_prf@meta.data)) "orig.ident" else "sample"

# 每个样本的 PRFs 绝对数
abs_count <- fib_prf@meta.data %>%
  filter(celltype_v3 == "PRFs") %>%
  group_by(!!sym(sample_col), group) %>%
  summarise(n_PRFs = n(), .groups = "drop")
colnames(abs_count)[1] <- "sample"

# 加入没有 PRFs 的样本（n_PRFs = 0）
all_samples <- fib_prf@meta.data %>%
  distinct(!!sym(sample_col), group)
colnames(all_samples)[1] <- "sample"

abs_count <- all_samples %>%
  left_join(abs_count, by = c("sample", "group")) %>%
  mutate(n_PRFs = ifelse(is.na(n_PRFs), 0, n_PRFs))

# ★ 删除 Normal，只保留 ASYM/SYM
abs_count <- abs_count %>%
  filter(group %in% c("ASYM", "SYM")) %>%
  mutate(group = factor(group, levels = c("ASYM", "SYM")))

print(abs_count)

# 均值 + SE
abs_summary <- abs_count %>%
  group_by(group) %>%
  summarise(
    mean_val = mean(n_PRFs),
    se_val   = sd(n_PRFs) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  )
print(abs_summary)

# P 值
p_abs <- tryCatch(
  wilcox.test(n_PRFs ~ group, data = abs_count)$p.value,
  error = function(e) NA
)
cat(sprintf("PRFs absolute count — Wilcoxon P = %.3f\n", p_abs))

# ── 绘图：柱状图 + 误差棒 + 个体点 + P 值 ────────────────────────────
p_abs_bar <- ggplot() +
  # 1) 柱：组均值
  geom_col(data = abs_summary,
           aes(x = group, y = mean_val, fill = group),
           width = 0.55, color = "black", linewidth = 0.4,
           alpha = 0.9) +
  # 2) 误差棒（mean ± SE）
  geom_errorbar(data = abs_summary,
                aes(x = group,
                    ymin = pmax(0, mean_val - se_val),
                    ymax = mean_val + se_val),
                width = 0.15, linewidth = 0.5,
                color = "black") +
  # 3) 个体散点（每个样本一个点）
  geom_jitter(data = abs_count,
              aes(x = group, y = n_PRFs),
              width = 0.10, height = 0,
              size = 2.5, shape = 21,
              fill = "white", color = "grey20",
              stroke = 0.5, alpha = 0.95) +
  # 4) P 值 bracket
  geom_signif(
    data = abs_count,
    aes(x = group, y = n_PRFs),
    comparisons      = list(c("ASYM", "SYM")),
    test             = "wilcox.test",
    map_signif_level = function(p) {
      if (is.na(p)) "n.s."
      else if (p < 0.001) sprintf("P = %.1e", p)
      else sprintf("P = %.3f", p)
    },
    textsize    = 3.0,
    tip_length  = 0.015,
    vjust       = -0.1,
    color       = "black",
    inherit.aes = FALSE) +
  scale_fill_manual(values = group_colors_2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "PRFs absolute count per sample",
       subtitle = "No numerical expansion in symptomatic disease",
       x = NULL, y = "Number of PRFs per sample") +
  theme_fig3 +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(face = "bold", size = 11)
  )

ggsave(paste0(out_dir, "Fig3F_PRFs_abscount_bar.pdf"), p_abs_bar,
       width = 4.5, height = 4.2, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig3F_PRFs_abscount_bar.png"), p_abs_bar,
       width = 4.5, height = 4.2, dpi = 600)

cat("\n??? PRFs absolute count bar plot saved\n")
