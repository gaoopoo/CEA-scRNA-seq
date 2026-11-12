# Project root: the folder containing SMC/, macro/, Fig3/, Fig5_CellChat/, etc.
# Every path in this script is built from this one setting, so the code no longer
# depends on a single hard-coded drive letter. Either start R/RStudio in the
# project root, or edit the one line below to point at it directly.
# Stored as an R option so it survives the rm(list = ls()) calls between sections.
options(cea.project_root = normalizePath(".", mustWork = FALSE))

#--Figure 7. The TREM2??? lipid-associated macrophage (LAM) signature is recoverable across independent bulk transcriptomic cohorts

# FigureA — Single-cell to clinical translation workflow

suppressMessages({
  library(ggplot2)
  library(tidyverse)
})

out_dir <- paste0(getOption("cea.project_root"), "/Fig6/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ── 主题 ──────────────────────────────────────────────────────────
theme_schematic <- theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 13,
                                     hjust = 0, margin = margin(b = 8)),
        plot.subtitle = element_text(size = 10, color = "grey40",
                                     hjust = 0, margin = margin(b = 12)),
        plot.margin   = margin(15, 15, 15, 15))

## ── 配色（与 Fig4-5 一致）─────────────────────────────────────────
cols <- list(
  sc_data    = "#5DC8B8",
  fibroblast = "#F0B5AD",
  trem2_lam  = "#8B0000",
  combined   = "#E59085",
  cohort1    = "#3A6FA0",   # GSE100927  presence
  cohort2    = "#5C8AB8",   # GSE28829   severity（新增）
  cohort3    = "#7FA5C7",   # GSE43292   within-patient
  arrow      = "grey40"
)

## ── box 坐标 / 文本 ───────────────────────────────────────────────
box_df <- tibble::tribble(
  ~id,      ~x,    ~y,    ~w,    ~h,   ~fill,            ~label,
  # Row1 数据源
  "sc",     5.5,   10,    4.6,   1.4,  cols$sc_data,     "Single-cell atlas\n(12 CEA: 6 ASYM + 6 SYM)",
  # Row2 三个 signature
  "fib",    2.0,   7.3,   2.7,   1.5,  cols$fibroblast,  "Fibroblast\nhealing\n10 genes",
  "trem2",  5.5,   7.3,   2.7,   1.5,  cols$trem2_lam,   "TREM2+ LAM\n21 genes",
  "comb",   9.0,   7.3,   2.7,   1.5,  cols$combined,    "Combined\nhealing\n28 genes",
  # Row3 三个 cohort（按严格度阶梯排序）
  "gse1",   2.0,   4.4,   3.0,   1.7,  cols$cohort1,     "GSE100927\n29 plaque vs 12 control\n(presence)",
  "gse2",   5.5,   4.4,   3.0,   1.7,  cols$cohort2,     "GSE28829\n16 advanced vs 13 early\n(severity)",
  "gse3",   9.0,   4.4,   3.0,   1.7,  cols$cohort3,     "GSE43292\n32 paired atheroma/intact\n(within-patient)",
  # Row4 输出
  "output", 5.5,   1.6,   5.6,   1.4,  "grey25",         "ssGSEA + bootstrap ROC\nAUC validation"
)

## ── 箭头 ──────────────────────────────────────────────────────────
arrow_df <- tibble::tribble(
  ~x,    ~y,    ~xend,  ~yend,
  # sc -> 3 signatures
  5.5,   9.3,   2.0,    8.05,
  5.5,   9.3,   5.5,    8.05,
  5.5,   9.3,   9.0,    8.05,
  # signatures(中心) -> 3 cohorts
  5.5,   6.55,  2.0,    5.25,
  5.5,   6.55,  5.5,    5.25,
  5.5,   6.55,  9.0,    5.25,
  # 3 cohorts -> output
  2.0,   3.55,  5.5,    2.3,
  5.5,   3.55,  5.5,    2.3,
  9.0,   3.55,  5.5,    2.3
)

## ── 绘图 ──────────────────────────────────────────────────────────
p_6A <- ggplot() +
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = cols$arrow, linewidth = 0.6, alpha = 0.7,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_tile(data = box_df,
            aes(x = x, y = y, width = w, height = h, fill = fill),
            color = "black", linewidth = 0.5, alpha = 0.92) +
  geom_text(data = box_df,
            aes(x = x, y = y, label = label),
            color = "white", fontface = "bold",
            size = 3.2, lineheight = 0.95) +
  annotate("text", x = 0.3, y = 10,  label = "Source",     angle = 90,
           fontface = "bold", size = 4, color = "grey30") +
  annotate("text", x = 0.3, y = 7.3, label = "Signatures", angle = 90,
           fontface = "bold", size = 4, color = "grey30") +
  annotate("text", x = 0.3, y = 4.4, label = "Cohorts",    angle = 90,
           fontface = "bold", size = 4, color = "grey30") +
  annotate("text", x = 0.3, y = 1.6, label = "Analysis",   angle = 90,
           fontface = "bold", size = 4, color = "grey30") +
  scale_fill_identity() +
  coord_fixed(ratio = 0.8, xlim = c(0, 11.5), ylim = c(0.5, 11)) +
  labs(title = "A. Single-cell to clinical translation workflow",
       subtitle = "Three healing signatures evaluated across three cohorts of increasing stringency") +
  theme_schematic

## ── 保存（默认 pdf 设备）─────────────────────────────────────────
ggsave(paste0(out_dir, "Fig6A.pdf"), p_6A, width = 7.5, height = 6)
ggsave(paste0(out_dir, "Fig6A.png"), p_6A, width = 7.5, height = 6, dpi = 600)
cat("\nFig 6A (V3) saved\n")

--------------------------------------------------------------------------------
# FigureB — ssGSEA scores 综合图 (3 sig × 2 cohort)
################################################################################

suppressMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggh4x)        # 用于 facet header 个性化配色
  library(patchwork)
})


out_dir <- paste0(getOption("cea.project_root"), "/Fig6/")

# ── 1. 读取已有 ssGSEA 数据 ────────────────────────────────────────
score_100927 <- read.csv(paste0(out_dir, "Fig6_GSE100927_scores.csv"),
                         stringsAsFactors = FALSE)
score_43292  <- read.csv(paste0(out_dir, "Fig6_GSE43292_scores.csv"),
                         stringsAsFactors = FALSE)

# ── 2. 整合长格式 ─────────────────────────────────────────────────
prepare_long <- function(df, cohort_name) {
  df %>%
    pivot_longer(cols = c("Fibroblast_Healing", "TREM2_LAM", "Combined_Healing"),
                 names_to = "Signature",
                 values_to = "Score") %>%
    mutate(Cohort = cohort_name)
}

score_long <- bind_rows(
  prepare_long(score_100927, "GSE100927"),
  prepare_long(score_43292,  "GSE43292")
)

# 统一 group 标签（control/intact = "Healthy"; plaque/atheroma = "Disease"）
score_long <- score_long %>%
  mutate(group_unified = case_when(
    group %in% c("control", "intact") ~ "Healthy",
    group %in% c("plaque", "atheroma") ~ "Disease",
    TRUE ~ as.character(group)
  ))

# 设 factor 顺序
score_long$Signature <- factor(score_long$Signature,
                               levels = c("Fibroblast_Healing",
                                          "TREM2_LAM",
                                          "Combined_Healing"))
score_long$Cohort <- factor(score_long$Cohort,
                            levels = c("GSE100927", "GSE43292"))
score_long$group_unified <- factor(score_long$group_unified,
                                   levels = c("Healthy", "Disease"))

# 重命名 facet header（更可读）
levels(score_long$Signature) <- c(
  "Fibroblast healing (10 genes)",
  "TREM2??? LAM (21 genes)",
  "Combined healing (28 genes)"
)
levels(score_long$Cohort) <- c(
  "GSE100927 (plaque vs control)",
  "GSE43292 (atheroma vs intact)"
)

# ── 3. 计算每个 facet 的 Wilcoxon P 值 ─────────────────────────────
p_values <- score_long %>%
  group_by(Signature, Cohort) %>%
  summarise(
    p_value = wilcox.test(Score ~ group_unified)$p.value,
    y_pos   = max(Score) * 1.04,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      p_value < 1e-9  ~ sprintf("P = %.1e", p_value),
      p_value < 0.001 ~ sprintf("P = %.1e", p_value),
      p_value < 0.05  ~ sprintf("P = %.3f", p_value),
      TRUE            ~ sprintf("P = %.2f (n.s.)", p_value)
    )
  )

print(p_values)

# ── 4. 配色 ────────────────────────────────────────────────────────
group_col <- c("Healthy" = "#F9C80EFF", "Disease" = "#C71000FF")

# Facet header 配色（强调 TREM2??? LAM 行）
strip_colors <- strip_themed(
  background_y = elem_list_rect(fill = c(
    "#F0B5AD",                  # Fibroblast: 浅红
    "#8B0000",                  # TREM2??? LAM: 深红（强调）
    "#E59085"                   # Combined: 中红
  )),
  text_y = elem_list_text(
    color = c("grey15", "white", "grey15"),
    face  = c("bold", "bold", "bold"),
    size  = c(11, 11.5, 11)
  ),
  background_x = elem_list_rect(fill = c("grey90", "grey90")),
  text_x = elem_list_text(face = "bold", size = 11, color = "grey15")
)

# ── 5. 绘图 ────────────────────────────────────────────────────────
p_6B <- ggplot(score_long,
               aes(x = group_unified, y = Score, fill = group_unified)) +
  # 1) 小提琴
  geom_violin(scale = "width", trim = TRUE,
              color = "black", linewidth = 0.3,
              alpha = 0.85) +
  # 2) 嵌套箱体
  geom_boxplot(width = 0.15, fill = "white",
               outlier.size = 0.3, linewidth = 0.3) +
  # 3) P 值标注
  geom_text(data = p_values,
            aes(x = 1.5, y = y_pos, label = p_label),
            inherit.aes = FALSE,
            size = 3.5, fontface = "bold",
            color = "grey15") +
  scale_fill_manual(values = group_col, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  # 4) Facet (用 ggh4x 的 facet_grid2 实现个性化 header)
  facet_grid2(Signature ~ Cohort,
              scales = "free_y",
              strip = strip_colors) +
  # 5) 标题
  labs(title = "B. ssGSEA scores across signatures and cohorts",
       subtitle = "Three signatures (rows) × Two independent cohorts (columns)",
       x = NULL,
       y = "ssGSEA score") +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    axis.text        = element_text(color = "black"),
    axis.text.x      = element_text(face = "bold", size = 11),
    axis.text.y      = element_text(size = 9),
    axis.line        = element_line(linewidth = 0.4),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.spacing    = unit(0.4, "cm")
  )

# ── 6. 保存 ────────────────────────────────────────────────────────
ggsave(paste0(out_dir, "Fig6B.pdf"), p_6B,
       width = 8.5, height = 8, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig6B.png"), p_6B,
       width = 8.5, height = 8, dpi = 600)

cat("\n??? Fig 6B saved\n")


--------------------------------------------------------------------------------
  # Figure C-D-E from the existing per-sample score files.
  
  required <- c("pROC", "ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install the missing packages first: ", paste(missing, collapse = ", "))
}
suppressPackageStartupMessages(library(ggplot2))

input_dir <- getOption("roc_fix.input_dir", paste0(getOption("cea.project_root"), "/Fig6"))
out_dir <- getOption("roc_fix.output_dir", file.path(input_dir, "ROC_fixed"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
n_boot <- 1000L
# A plotting-only repair should not replace previously reported bootstrap CIs.
# Set options(roc_fix.recompute_ci = TRUE) before sourcing to recompute them.
recompute_ci <- isTRUE(getOption("roc_fix.recompute_ci", FALSE))

cohorts <- list(
  GSE100927 = list(file = "Fig6_GSE100927_scores.csv",
                   lv = c("control", "plaque"), panel = "C"),
  GSE43292 = list(file = "Fig6_GSE43292_scores.csv",
                  lv = c("intact", "atheroma"), panel = "D"),
  GSE28829 = list(file = "Fig6_GSE28829_scores.csv",
                  lv = c("early", "advanced"), panel = "E")
)
# Preserve the original calculation and drawing order (LAM drawn last).
sigs <- c("Fibroblast_Healing", "Combined_Healing", "TREM2_LAM")
sig_col <- c(Fibroblast_Healing = "#9E9E9E",
             Combined_Healing = "#E59085", TREM2_LAM = "#8B0000")
sig_lw <- c(Fibroblast_Healing = 0.55,
            Combined_Healing = 0.85, TREM2_LAM = 1.25)
sig_lab <- c(Fibroblast_Healing = "Fibroblast healing",
             Combined_Healing = "Combined healing", TREM2_LAM = "TREM2+ LAM")
legend_order <- c("TREM2_LAM", "Fibroblast_Healing", "Combined_Healing")
ypos <- c(TREM2_LAM = 0.18, Fibroblast_Healing = 0.10, Combined_Healing = 0.02)

roc_list <- auc_list <- check_list <- object_list <- list()
tol <- 1e-10
old_auc_file <- file.path(input_dir, "Fig6_AUC_3cohorts.csv")
old_auc <- if (file.exists(old_auc_file)) read.csv(old_auc_file) else NULL
if (!recompute_ci && is.null(old_auc)) {
  stop("Existing AUC table is missing. To compute new CIs, set ",
       "options(roc_fix.recompute_ci = TRUE) before sourcing this script.")
}

for (co in names(cohorts)) {
  cfg <- cohorts[[co]]
  f <- file.path(input_dir, cfg$file)
  if (!file.exists(f)) stop("Missing score file: ", f)
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("group", sigs) %in% names(d))) stop("Missing columns in ", f)
  if (anyNA(d$group) || !all(d$group %in% cfg$lv)) {
    stop("Unexpected or missing group labels in ", f)
  }
  if ("sample" %in% names(d) && anyDuplicated(d$sample)) {
    stop("Duplicate sample identifiers in ", f)
  }
  d$group <- factor(d$group, levels = cfg$lv)
  
  for (s in sigs) {
    if (!is.numeric(d[[s]]) || any(is.infinite(d[[s]]))) {
      stop("Scores must be numeric and finite (NA is allowed): ", co, " / ", s)
    }
    # Match pROC's original default removal of missing scores, but report it.
    ok <- !is.na(d[[s]])
    if (!all(ok)) warning(co, " / ", s, ": excluded ", sum(!ok), " missing scores")
    if (any(table(d$group[ok]) == 0L)) stop("Both classes are required: ", co, " / ", s)
    r <- pROC::roc(d$group[ok], d[[s]][ok], levels = cfg$lv,
                   direction = "<", quiet = TRUE)
    actual_auc <- as.numeric(pROC::auc(r))
    prior <- if (!is.null(old_auc)) {
      old_auc[old_auc$Cohort == co & old_auc$Signature == s, , drop = FALSE]
    } else data.frame()
    if (recompute_ci) {
      ci <- pROC::ci.auc(r, method = "bootstrap", boot.n = n_boot,
                         boot.stratified = TRUE)
      ci_source <- "Recomputed: stratified bootstrap, 1000 replicates"
    } else {
      if (nrow(prior) != 1L || !all(c("AUC", "lo", "hi") %in% names(prior))) {
        stop("Existing AUC/CI entry is missing or duplicated: ", co, " / ", s)
      }
      vals <- unlist(prior[c("AUC", "lo", "hi")], use.names = FALSE)
      if (any(!is.finite(vals)) || any(vals < 0 | vals > 1) ||
          prior$lo > prior$hi || abs(prior$AUC - actual_auc) > tol) {
        stop("Existing AUC/CI table does not match this ROC: ", co, " / ", s,
             ". Verify inputs, or explicitly request CI recomputation.")
      }
      ci <- c(prior$lo, actual_auc, prior$hi)
      ci_source <- "Existing AUC table; point estimate checked against current scores"
    }
    key <- paste(co, s, sep = "__")
    object_list[[key]] <- r
    z <- data.frame(Cohort = co, Signature = s,
                    Source_index = seq_along(r$sensitivities),
                    Threshold = r$thresholds,
                    Specificity = r$specificities,
                    Sensitivity = r$sensitivities,
                    FPR = 1 - r$specificities)
    
    # Reproduce the original x-only sort for diagnosis (NOT for plotting).
    old_order <- z[order(z$FPR, method = "radix"), ]
    old_drops <- sum(diff(old_order$Sensitivity) < -tol)
    
    # KEY FIX: sort complete coordinate rows by BOTH x and y within each ROC.
    # Never sort the x and y vectors independently.
    z <- z[order(z$FPR, z$Sensitivity, method = "radix"), ]
    rownames(z) <- NULL
    monotone <- all(diff(z$FPR) >= -tol) && all(diff(z$Sensitivity) >= -tol)
    endpoints <- all(abs(c(z$FPR[1], z$Sensitivity[1],
                           tail(z$FPR, 1) - 1, tail(z$Sensitivity, 1) - 1)) < tol)
    # Confirm the plotted path's trapezoidal area equals the pROC AUC.
    path_auc <- sum(diff(z$FPR) *
                      (head(z$Sensitivity, -1) + tail(z$Sensitivity, -1)) / 2)
    if (!monotone || !endpoints || abs(path_auc - actual_auc) > tol) {
      stop("ROC integrity check failed: ", co, " / ", s)
    }
    roc_list[[key]] <- z
    auc_list[[key]] <- data.frame(Cohort = co, Signature = s,
                                  AUC = actual_auc, lo = as.numeric(ci[1]),
                                  hi = as.numeric(ci[3]))
    check_list[[key]] <- data.frame(
      Cohort = co, Signature = s, N_control = sum(d$group[ok] == cfg$lv[1]),
      N_case = sum(d$group[ok] == cfg$lv[2]), N_excluded = sum(!ok),
      Old_x_only_downward_segments = old_drops,
      Fixed_downward_segments = sum(diff(z$Sensitivity) < -tol),
      Monotone = monotone, Endpoints_valid = endpoints,
      Path_AUC = path_auc, pROC_AUC = actual_auc, CI_source = ci_source,
      AUC_difference_from_existing = if (nrow(prior) == 1L) actual_auc - prior$AUC else NA_real_,
      CI_low_difference_from_existing = if (nrow(prior) == 1L) as.numeric(ci[1]) - prior$lo else NA_real_,
      CI_high_difference_from_existing = if (nrow(prior) == 1L) as.numeric(ci[3]) - prior$hi else NA_real_
    )
  }
}

roc_df <- do.call(rbind, roc_list)
auc_df <- do.call(rbind, auc_list)
checks <- do.call(rbind, check_list)
rownames(roc_df) <- rownames(auc_df) <- rownames(checks) <- NULL
roc_df$Signature <- factor(roc_df$Signature, levels = sigs)
auc_df$Signature <- factor(auc_df$Signature, levels = sigs)

# KEY FIX: every annotation is generated from the calculated AUC and CI.
# The original hard-coded GSE100927 LAM label is deliberately removed.
auc_df$disp <- sprintf("AUC = %.2f (%.2f-%.2f)", auc_df$AUC, auc_df$lo, auc_df$hi)
auc_df$y_pos <- unname(ypos[as.character(auc_df$Signature)])

make_panel <- function(co, show_y = TRUE) {
  z <- roc_df[roc_df$Cohort == co, ]
  a <- auc_df[auc_df$Cohort == co, ]
  p <- ggplot(z, aes(FPR, Sensitivity, color = Signature,
                     linewidth = Signature, group = Signature)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey60", linewidth = 0.4) +
    # KEY FIX: geom_path preserves the explicitly validated row order.
    geom_path(linejoin = "round", lineend = "round") +
    geom_text(data = a, aes(x = 0.98, y = y_pos, label = disp, color = Signature),
              inherit.aes = FALSE, hjust = 1, size = 2.75,
              fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = sig_col, labels = sig_lab,
                       breaks = legend_order, name = NULL) +
    scale_linewidth_manual(values = sig_lw, guide = "none") +
    scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                       limits = c(0, 1), expand = expansion(mult = 0.04)) +
    scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                       limits = c(0, 1), expand = expansion(mult = 0.04)) +
    facet_wrap(~ Cohort) + coord_fixed() +
    labs(x = "1 - Specificity (False positive rate)",
         y = if (show_y) "Sensitivity (True positive rate)" else NULL,
         tag = cohorts[[co]]$panel) +
    theme_bw(base_size = 10) +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          strip.background = element_rect(fill = "grey92", color = "grey80", linewidth = 0.3),
          strip.text = element_text(size = 10),
          panel.grid = element_blank(),
          plot.tag = element_text(face = "bold", size = 19),
          plot.margin = margin(7, 7, 7, 7),
          legend.text = element_text(size = 9)) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.1)))
  p
}

plots <- lapply(names(cohorts), function(co) make_panel(co, co == "GSE100927"))
names(plots) <- names(cohorts)
p_CDE <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
  theme(legend.position = "right")
p_CD <- patchwork::wrap_plots(plots[c("GSE100927", "GSE43292")],
                              nrow = 1, guides = "collect") &
  theme(legend.position = "right")

save_pair <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, units = "in", bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), p,
         width = width, height = height, units = "in", dpi = 600, bg = "white")
}
save_pair(p_CDE, "ROC_CDE_fixed", 12, 4.1)
save_pair(p_CD, "ROC_CD_fixed", 9, 4.1)
for (co in names(cohorts)) {
  save_pair(make_panel(co) + theme(legend.position = "bottom"),
            paste0("ROC_", cohorts[[co]]$panel, "_", co, "_fixed"), 4.6, 4.9)
}
write.csv(auc_df[c("Cohort", "Signature", "AUC", "lo", "hi")],
          file.path(out_dir, "ROC_AUC_3cohorts_fixed.csv"), row.names = FALSE)
write.csv(roc_df, file.path(out_dir, "ROC_coordinates_fixed.csv"), row.names = FALSE)
write.csv(checks, file.path(out_dir, "ROC_integrity_checks.csv"), row.names = FALSE)
saveRDS(object_list, file.path(out_dir, "ROC_objects.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
input_files <- vapply(cohorts, function(x) file.path(input_dir, x$file), character(1))
if (!recompute_ci) input_files <- c(input_files, old_auc_file)
write.csv(data.frame(File = input_files, MD5 = unname(tools::md5sum(input_files))),
          file.path(out_dir, "input_file_checksums.csv"), row.names = FALSE)
print(checks[c("Cohort", "Signature", "Old_x_only_downward_segments",
               "Fixed_downward_segments", "pROC_AUC")], row.names = FALSE)
print(auc_df[c("Cohort", "Signature", "disp")], row.names = FALSE)
if (interactive()) print(p_CDE)
message("Saved corrected figures and audit tables to: ", normalizePath(out_dir, winslash = "/"))

--------------------------------------------------------------------------------
  # Figure F — 单图版 (仅配对连线图) 
suppressMessages({
  library(tidyverse); library(ggplot2)
})

out_dir <- paste0(getOption("cea.project_root"), "/Fig6/")

# ── 数据准备 ─────────────────────────────────────────────────────
score_43292 <- read.csv(paste0(out_dir, "Fig6_GSE43292_scores.csv"), stringsAsFactors = FALSE)
paired_info <- read.csv(paste0(out_dir, "GSE43292_paired_info.csv"), stringsAsFactors = FALSE)

df_paired <- score_43292 %>%
  left_join(paired_info %>% dplyr::select(sample, patient_id, tissue), by = "sample") %>%
  mutate(tissue_short = ifelse(grepl("Atheroma", tissue), "Atheroma", "Intact"),
         tissue_short = factor(tissue_short, levels = c("Intact", "Atheroma")))

diff_summary <- df_paired %>%
  dplyr::select(patient_id, tissue_short, TREM2_LAM) %>%
  pivot_wider(names_from = tissue_short, values_from = TREM2_LAM) %>%
  mutate(diff = Atheroma - Intact,
         direction = case_when(
           diff > 0.05  ~ "Atheroma > Intact",
           diff < -0.05 ~ "Atheroma < Intact",
           TRUE         ~ "Comparable"))

n_total      <- nrow(diff_summary)
n_higher     <- sum(diff_summary$direction == "Atheroma > Intact")
n_lower      <- sum(diff_summary$direction == "Atheroma < Intact")
n_comparable <- sum(diff_summary$direction == "Comparable")
p_paired <- wilcox.test(diff_summary$Atheroma, diff_summary$Intact, paired = TRUE)$p.value

# ── 配色逻辑优化: 点(组织)与线(方向)用两套不冲突的色系 ──────────────
# 点 = 组织类型: Intact 冷灰蓝 / Atheroma 暖红 (冷暖对比)
tissue_col <- c("Intact" = "#5B8BA6", "Atheroma" = "#D1495B")
# 线 = 变化方向: 升高 橙 / 持平 浅灰 / 反向 深蓝紫 (与点色错开, 避免混淆)
direction_col <- c("Atheroma > Intact" = "#E8772E",   # 橙 (升高)
                   "Comparable"        = "#C9C9C9",   # 浅灰 (持平)
                   "Atheroma < Intact" = "#3B5BA5")   # 蓝紫 (反向)

df_plot <- df_paired %>%
  left_join(diff_summary %>% dplyr::select(patient_id, direction), by = "patient_id")
df_plot$direction <- factor(df_plot$direction,
                            levels = c("Atheroma > Intact","Comparable","Atheroma < Intact"))

# ── 配对连线图 ───────────────────────────────────────────────────
p_D <- ggplot(df_plot, aes(x = tissue_short, y = TREM2_LAM, group = patient_id)) +
  geom_line(aes(color = direction, linewidth = direction, alpha = direction)) +
  geom_point(aes(fill = tissue_short), shape = 21, size = 3.0,
             stroke = 0.4, color = "black", alpha = 0.95) +
  annotate("text", x = 1.5, y = max(df_paired$TREM2_LAM) * 1.07,
           label = sprintf("Paired Wilcoxon\nP = %.2e", p_paired),
           size = 4, fontface = "bold", color = "grey15", lineheight = 0.95) +
  scale_color_manual(values = direction_col, name = "Within-patient direction") +
  scale_fill_manual(values = tissue_col, name = "Tissue") +
  scale_linewidth_manual(values = c("Atheroma > Intact"=0.8,"Comparable"=0.4,"Atheroma < Intact"=0.7),
                         guide = "none") +
  scale_alpha_manual(values = c("Atheroma > Intact"=0.85,"Comparable"=0.4,"Atheroma < Intact"=0.8),
                     guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(title = "GSE43292 paired analysis (n = 32 patients)",
       subtitle = sprintf("%d atheroma > intact (%.0f%%) | %d comparable (%.0f%%) | %d reverse (%.0f%%)",
                          n_higher, n_higher/n_total*100,
                          n_comparable, n_comparable/n_total*100,
                          n_lower, n_lower/n_total*100),
       x = NULL, y = "TREM2\u207A LAM ssGSEA score") +
  guides(fill = guide_legend(override.aes = list(size = 3.5), order = 1),
         color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1), order = 2)) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        axis.text     = element_text(color = "black"),
        axis.text.x   = element_text(face = "bold", size = 12),
        axis.line     = element_line(linewidth = 0.4),
        axis.ticks    = element_line(linewidth = 0.4),
        legend.position = "right",
        legend.title  = element_text(size = 10, face = "bold"),
        legend.text   = element_text(size = 9),
        panel.grid.major.y = element_line(color = "grey94", linewidth = 0.3),
        plot.margin   = margin(12, 12, 12, 12))

ggsave(paste0(out_dir, "Fig6D_single.pdf"), p_D, width = 7, height = 5.5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig6D_single.png"), p_D, width = 7, height = 5.5, dpi = 600)
cat("\nFig 6D single saved\n")

--------------------------------------------------------------------------------
  # FigG — GSE43292 within-patient concordance（单图：横向堆叠柱）
  # 从左到右：红 Atheroma>Intact / 浅灰 Comparable / 蓝 Atheroma<Intact
suppressMessages({ library(tidyverse); library(ggplot2) })
out_dir <- paste0(getOption("cea.project_root"), "/Fig6/")
ver     <- "v2"
## ── 数据准备 ─────────────────────────────────────────────────────
score_43292 <- read.csv(paste0(out_dir, "Fig6_GSE43292_scores.csv"), stringsAsFactors = FALSE)
paired_info <- read.csv(paste0(out_dir, "GSE43292_paired_info.csv"), stringsAsFactors = FALSE)

df_paired <- score_43292 %>%
  left_join(paired_info %>% dplyr::select(sample, patient_id, tissue), by = "sample") %>%
  mutate(tissue_short = ifelse(grepl("Atheroma", tissue), "Atheroma", "Intact"))

diff_summary <- df_paired %>%
  dplyr::select(patient_id, tissue_short, TREM2_LAM) %>%
  pivot_wider(names_from = tissue_short, values_from = TREM2_LAM) %>%
  mutate(diff = Atheroma - Intact,
         direction = case_when(diff > 0.05  ~ "Atheroma > Intact",
                               diff < -0.05 ~ "Atheroma < Intact",
                               TRUE         ~ "Comparable"))

n_total      <- nrow(diff_summary)
n_higher     <- sum(diff_summary$direction == "Atheroma > Intact")
n_lower      <- sum(diff_summary$direction == "Atheroma < Intact")
n_comparable <- sum(diff_summary$direction == "Comparable")

## ── 统计数据框 ───────────────────────────────────────────────────
## 横向（coord_flip 后）：factor 的第一个 level 排在最左。
## 要"左红右蓝"，factor 顺序设为 红、灰、蓝。
concordance_df <- data.frame(
  category = c("Atheroma > Intact", "Comparable", "Atheroma < Intact"),
  count    = c(n_higher, n_comparable, n_lower)
) %>%
  mutate(label = paste0(count, "/", n_total, "\n(", round(count/n_total*100, 0), "%)"),
         category = factor(category,
                           levels = c("Atheroma > Intact", "Comparable", "Atheroma < Intact")))

direction_col <- c("Atheroma > Intact" = "#8B0000",   # 红
                   "Comparable"        = "#BBBBBB",   # 浅灰
                   "Atheroma < Intact" = "#3A6FA0")   # 蓝

## ── 横向堆叠柱 ───────────────────────────────────────────────────
p <- ggplot(concordance_df, aes(x = "", y = count, fill = category)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.5) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.6, lineheight = 0.95) +
  scale_fill_manual(values = direction_col, name = "Within-patient direction",
                    breaks = c("Atheroma > Intact", "Comparable", "Atheroma < Intact")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, n_total)) +
  coord_flip() +
  labs(title = "Within-patient concordance",
       subtitle = sprintf("%d/%d (%.0f%%) patients do not show Intact > Atheroma",
                          n_higher + n_comparable, n_total,
                          (n_higher + n_comparable)/n_total*100),
       x = NULL, y = "Number of patients") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "grey40"),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.line.y   = element_blank(),
        axis.text.x   = element_text(size = 10, color = "black"),
        axis.line.x   = element_line(linewidth = 0.4),
        axis.ticks.x  = element_line(linewidth = 0.4),
        legend.position = "bottom",
        legend.title  = element_text(size = 9, face = "bold"),
        legend.text   = element_text(size = 8.5),
        legend.key.size = unit(0.4, "cm"),
        plot.margin   = margin(12, 16, 12, 16)) +
  guides(fill = guide_legend(nrow = 1))

ggsave(paste0(out_dir, "Fig6_concordance_", ver, ".pdf"), p, width = 6.5, height = 2.6)
ggsave(paste0(out_dir, "Fig6_concordance_", ver, ".png"), p, width = 6.5, height = 2.6, dpi = 600)
cat("\nsaved Fig6_concordance_", ver, " (pdf/png)\n", sep = "")

--------------------------------------------------------------------------------
  #FigH--GSE28829 — TREM2+ LAM score 的"进展性"验证
# 目的：检验 LAM ssGSEA score 是否随颈动脉斑块进展升高
#       （early / pathological intimal thickening  ->  advanced / fibrous cap atheroma）
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("GEOquery", "GSVA"))
# install.packages(c("ggplot2", "pROC"))
rm(list = ls())
gc()
library(GEOquery)
library(GSVA)
library(ggplot2)

outdir <- paste0(getOption("cea.project_root"), "/Fig6")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


## ---- 1. 拉取 ------------------------------------------------------------
gse  <- getGEO("GSE28829", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
expr <- Biobase::exprs(gse)        # 探针级矩阵（series matrix 已归一化）
fd   <- Biobase::fData(gse)
pd   <- Biobase::pData(gse)

## ---- 2. 先看分组字段（务必先跑这段，确认 early/advanced 怎么编码）-------
cat("\n--- pData 候选分组列 ---\n")
cand <- grep("characteristics|title|source", colnames(pd), ignore.case = TRUE)
print(head(pd[, cand, drop = FALSE], 4))


## ---- 3. 分组（按需调整 pattern / 列名）----------------------------------
lab <- tolower(paste(pd$title, pd$source_name_ch1, pd$characteristics_ch1))
grp <- ifelse(grepl("advanced|fibrous", lab), "advanced",
              ifelse(grepl("early|intimal",   lab), "early", NA))
grp <- factor(grp, levels = c("early", "advanced"))

cat("\n--- 分组计数（应为 early=13, advanced=16）---\n")
print(table(grp, useNA = "ifany"))
stopifnot(sum(!is.na(grp)) == ncol(expr))   # 没有漏分；若报错说明 pattern 要改

## ---- 4. 探针 -> 基因符号（取表达均值最高的探针）-------------------------
sym_col <- grep("^gene symbol$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
sym <- as.character(fd[[sym_col]])
keep <- !is.na(sym) & sym != "" & !grepl("///", sym)   # 去掉空/多映射探针
expr <- expr[keep, ]; sym <- sym[keep]
rownames(expr) <- sym
expr <- expr[order(rowMeans(expr), decreasing = TRUE), ]
expr <- expr[!duplicated(rownames(expr)), ]            # 每基因保留最高表达探针
cat(sprintf("\n基因数（去重后）：%d\n", nrow(expr)))

## ---- 5. 21 基因 LAM 清单 ----------
lam_genes <- c(
  "TREM2", "APOE", "SPP1", "GPNMB", "LIPA", "FABP5",
  "CTSD", "CTSB", "CD9", "CD63", "LGALS3", "FABP4",
  "PLIN2", "CD36", "MMP9", "ACP5", "LPL",
  "SLC11A1", "CTSL", "CTSK", "CTSS"
)
stopifnot(length(lam_genes) >= 1)   # 提醒：没填会在此停下
present <- intersect(lam_genes, rownames(expr))
cat(sprintf("LAM 基因映射上：%d / %d\n", length(present), length(lam_genes)))
if (length(present) < length(lam_genes))
  cat("未映射：", setdiff(lam_genes, rownames(expr)), "\n")

## ---- 6. ssGSEA 打分）-------------------------------------
score <- as.numeric(gsva(expr, list(LAM = present), method = "ssgsea"))
## --- 若你的 GSVA 是旧版，注释上面两行，改用： ---
# score <- as.numeric(gsva(expr, list(LAM = present), method = "ssgsea"))

df <- data.frame(sample = colnames(expr), LAM = score, group = grp)

## ---- 7. 统计 + AUC）----------------------------------
w <- wilcox.test(LAM ~ group, data = df)            # early vs advanced
cat("\n=== Wilcoxon (early vs advanced) ===\n"); print(w)

# 方向核对：中位数应当 advanced > early
print(tapply(df$LAM, df$group, median))

# 与 Fig6 一致的 AUC（advanced 为阳性类）
if (requireNamespace("pROC", quietly = TRUE)) {
  roc <- pROC::roc(df$group, df$LAM, levels = c("early", "advanced"),
                   direction = "<", quiet = TRUE)
  ci  <- pROC::ci.auc(roc, method = "bootstrap", boot.n = 1000)
  cat(sprintf("AUC = %.2f (95%% CI %.2f-%.2f)\n", as.numeric(roc$auc), ci[1], ci[3]))
}

## ---- 8. 箱线图-------------------
p <- ggplot(df, aes(group, LAM, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.85) +
  scale_fill_manual(values = c(early = "#F2C200", advanced = "#D7261E")) +
  labs(x = NULL, y = "TREM2+ LAM ssGSEA score",
       title = sprintf("GSE28829  (Wilcoxon P = %.2g)", w$p.value)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
print(p)

ggsave(file.path(outdir, "GSE28829_LAM_early_vs_advanced.pdf"),
       p, width = 3.2, height = 3.6)            
ggsave(file.path(outdir, "GSE28829_LAM_early_vs_advanced.png"),
       p, width = 3.2, height = 3.6, dpi = 600)


