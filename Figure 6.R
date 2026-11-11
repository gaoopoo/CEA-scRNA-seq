#Figure 6:A three-layered cell–cell communication architecture characterizes a repair-associated program

# FigA — ASYM vs SYM dual-axis line plot
#
# Number of interactions: solid line + circles (left y-axis)
# Interaction strength:   dashed line + triangles (right y-axis)
# ASYM nodes: yellow; SYM nodes: red
# Raw values and SYM/ASYM fold changes are both displayed.

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(dplyr)
})

setwd("G:/13sample")
out_dir <- "G:/13sample/Fig5_v2/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(
  "ASYM" = "#FFE705FF",
  "SYM"  = "#C71000FF"
)

metric_colors <- c(
  "Number of interactions" = "#4D4D4D",
  "Interaction strength"   = "#8C8C8C"
)

# ============================================================================
# Load only ASYM and SYM CellChat objects — NP/Normal is completely removed
# ============================================================================
cellchat_asym <- readRDS("cellchat/cellchat_asym.rds")
cellchat_sym  <- readRDS("cellchat/cellchat_sym.rds")

cellchat_list <- list(
  ASYM = cellchat_asym,
  SYM  = cellchat_sym
)

# Extract values directly from the CellChat objects; no hard-coded values.
df_5A <- data.frame(
  group = factor(names(cellchat_list), levels = c("ASYM", "SYM")),
  x = c(1, 2),
  count = vapply(
    cellchat_list,
    function(x) sum(x@net$count, na.rm = TRUE),
    numeric(1)
  ),
  strength = vapply(
    cellchat_list,
    function(x) sum(x@net$weight, na.rm = TRUE),
    numeric(1)
  ),
  row.names = NULL
)

if (any(!is.finite(df_5A$count)) || any(!is.finite(df_5A$strength))) {
  stop("The extracted interaction number or strength contains non-finite values.")
}

if (any(df_5A$count <= 0) || any(df_5A$strength <= 0)) {
  stop("This dual-axis/fold-change plot requires positive count and strength values.")
}

print(df_5A)

# ============================================================================
# Dual-axis transformation
#
# Align the two axes at ASYM. This is a fixed, zero-preserving linear mapping:
#   displayed strength = raw strength * scale_factor
#   right-axis strength = displayed value / scale_factor
# Therefore, separation at SYM reflects the different SYM/ASYM fold changes.
# ============================================================================
asym_count    <- df_5A$count[df_5A$group == "ASYM"]
asym_strength <- df_5A$strength[df_5A$group == "ASYM"]

scale_factor <- asym_count / asym_strength

count_fold <-
  df_5A$count[df_5A$group == "SYM"] /
  df_5A$count[df_5A$group == "ASYM"]

strength_fold <-
  df_5A$strength[df_5A$group == "SYM"] /
  df_5A$strength[df_5A$group == "ASYM"]

message(sprintf("Number of interactions, SYM/ASYM: %.3fx", count_fold))
message(sprintf("Interaction strength, SYM/ASYM: %.3fx", strength_fold))

# Long-format plotting data.
metric_levels <- c("Number of interactions", "Interaction strength")

df_plot <- bind_rows(
  df_5A %>%
    transmute(
      group,
      x,
      metric = "Number of interactions",
      raw_value = count,
      display_y = count,
      raw_label = formatC(count, format = "f", digits = 0, big.mark = ",")
    ),
  df_5A %>%
    transmute(
      group,
      x,
      metric = "Interaction strength",
      raw_value = strength,
      display_y = strength * scale_factor,
      raw_label = formatC(strength, format = "f", digits = 2)
    )
) %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    # Tiny horizontal separation keeps the circle and triangle visible where
    # the two metrics coincide at the ASYM baseline.
    plot_x = x + if_else(metric == "Number of interactions", -0.018, 0.018)
  )

# Positions for direct fold-change labels between ASYM and SYM.
display_span <- diff(range(df_plot$display_y))
if (display_span == 0) display_span <- max(df_plot$display_y) * 0.1

fold_labels <- data.frame(
  metric = factor(metric_levels, levels = metric_levels),
  x = c(1.50, 1.50),
  y = c(
    mean(df_5A$count) - 0.08 * display_span,
    mean(df_5A$strength * scale_factor) + 0.08 * display_span
  ),
  label = c(
    sprintf("Number: %.2f\u00d7", count_fold),
    sprintf("Strength: %.2f\u00d7", strength_fold)
  )
)

# Honest viewing range for a point/line plot. The transformation remains
# anchored at zero even though the displayed window focuses on observed data.
y_range <- range(df_plot$display_y)
y_pad <- max(diff(y_range) * 0.24, max(y_range) * 0.04)
y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)

# ============================================================================
# Plot
# ============================================================================
p_5A <- ggplot(
  df_plot,
  aes(
    x = plot_x,
    y = display_y,
    group = metric,
    color = metric,
    linetype = metric,
    shape = metric,
    fill = group
  )
) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 5.8, stroke = 0.85) +
  
  # Raw values: number above its circles; strength below its triangles.
  geom_text(
    data = filter(df_plot, metric == "Number of interactions"),
    aes(label = raw_label),
    color = "black",
    vjust = -1.25,
    size = 3.8,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(df_plot, metric == "Interaction strength"),
    aes(label = raw_label),
    color = "black",
    vjust = 2.15,
    size = 3.6,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  # Direct SYM/ASYM fold-change labels.
  geom_label(
    data = fold_labels,
    aes(x = x, y = y, label = label, color = metric),
    inherit.aes = FALSE,
    size = 3.3,
    fontface = "bold",
    fill = scales::alpha("white", 0.88),
    label.size = 0,
    show.legend = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("ASYM", "SYM"),
    limits = c(0.72, 2.28),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    name = "Number of interactions",
    limits = y_limits,
    expand = expansion(mult = c(0.01, 0.01)),
    labels = function(x) formatC(x, format = "f", digits = 0, big.mark = ","),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Interaction strength",
      labels = function(x) formatC(x, format = "f", digits = 0)
    )
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_manual(
    values = metric_colors,
    breaks = metric_levels,
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "Number of interactions" = "solid",
      "Interaction strength" = "dashed"
    ),
    breaks = metric_levels,
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Number of interactions" = 21,
      "Interaction strength" = 24
    ),
    breaks = metric_levels,
    name = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(fill = "white", size = 4.2, linewidth = 1.1)
    ),
    linetype = guide_legend(
      override.aes = list(fill = "white", size = 4.2, linewidth = 1.1)
    ),
    shape = guide_legend(
      override.aes = list(fill = "white", size = 4.2, linewidth = 1.1)
    )
  ) +
  labs(
    x = NULL,
    title = "Intercellular communication increases from ASYM to SYM",
    subtitle = paste0(
      "Axes are aligned at ASYM; solid circles use the left axis and ",
      "dashed triangles use the right axis"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold", size = 13.5, hjust = 0.5,
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 9.8, hjust = 0.5, color = "grey30",
      margin = margin(b = 10)
    ),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(size = 11.5, margin = margin(r = 7)),
    axis.title.y.right = element_text(size = 11.5, margin = margin(l = 7)),
    axis.line = element_line(linewidth = 0.55, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10.5),
    legend.key.width = grid::unit(1.3, "cm"),
    legend.margin = margin(b = 5),
    plot.margin = margin(8, 12, 8, 10)
  )

# ============================================================================
# Export
# ============================================================================
pdf_file <- file.path(out_dir, "Fig5A_ASYM_SYM_dualaxis_lines.pdf")
png_file <- file.path(out_dir, "Fig5A_ASYM_SYM_dualaxis_lines.png")

ggsave(
  filename = pdf_file,
  plot = p_5A,
  width = 7.2,
  height = 5.2,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = p_5A,
  width = 7.2,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

print(p_5A)
message("Saved PDF: ", pdf_file)
message("Saved PNG: ", png_file)


--------------------------------------------------------------------------------
  # Figure B — Top 15 differential pathways (SYM vs ASYM)
  # 0. 环境
  # ============================================================================
rm(list = ls())
gc()

suppressMessages({
  library(CellChat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

setwd("G:/13sample")
out_dir <- "G:/13sample/Fig5_v2/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. 加载已有 CellChat 对象 ───────────────────────────────────────
cellchat_asym <- readRDS("cellchat/cellchat_asym.rds")
cellchat_sym  <- readRDS("cellchat/cellchat_sym.rds")

cellchat_2grp <- mergeCellChat(
  list(ASYM = cellchat_asym, SYM = cellchat_sym),
  add.names = c("ASYM", "SYM"),
  cell.prefix = TRUE
)

# ── 2. 跑 rankNet 并提取数据 ───────────────────────────────────────
# rankNet 返回 ggplot 对象，data 在 $data
rank_obj <- rankNet(cellchat_2grp, mode = "comparison",
                    stacked = FALSE, do.stat = TRUE,
                    return.data = TRUE)

# 提取数据 (rankNet 的 return.data 行为)
if (is.list(rank_obj) && "signaling.contribution" %in% names(rank_obj)) {
  rank_df <- rank_obj$signaling.contribution
} else if (inherits(rank_obj, "ggplot")) {
  rank_df <- rank_obj$data
} else {
  rank_df <- rank_obj
}

cat("rankNet data structure:\n")
print(head(rank_df))
print(colnames(rank_df))

# ── 3. 计算 SYM vs ASYM 差异 ───────────────────────────────────────
# 假设列名：name (pathway), group (ASYM/SYM), contribution / contribution.relative
# 转成宽表
rank_wide <- rank_df %>%
  dplyr::select(name, group, contribution) %>%
  pivot_wider(names_from = group, values_from = contribution,
              values_fill = 0) %>%
  mutate(
    diff_strength = SYM - ASYM,
    direction = case_when(
      diff_strength >  0 ~ "SYM_up",
      diff_strength <  0 ~ "ASYM_up",
      TRUE ~ "ns"
    )
  )

cat("\nTotal pathways:", nrow(rank_wide), "\n")

# 取 SYM 上调 top 8 + ASYM 上调 top 7
top_sym <- rank_wide %>%
  filter(direction == "SYM_up") %>%
  arrange(desc(diff_strength)) %>%
  head(8)

top_asym <- rank_wide %>%
  filter(direction == "ASYM_up") %>%
  arrange(diff_strength) %>%
  head(7)

top_combined <- bind_rows(top_sym, top_asym) %>%
  arrange(diff_strength) %>%
  mutate(name = factor(name, levels = name))

cat("\n=== Top 15 differential pathways ===\n")
print(top_combined)

write.csv(top_combined, paste0(out_dir, "Fig5B_top15_pathways.csv"),
          row.names = FALSE)

# ── 4. 画 bar plot ─────────────────────────────────────────────────
p_5B <- ggplot(top_combined,
               aes(x = diff_strength, y = name, fill = direction)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "SYM_up"   = "#C71000",
    "ASYM_up"  = "#3A6FA0"),
    labels = c(
      "SYM_up"   = "Enriched in SYM (repair / matrix)",
      "ASYM_up"  = "Enriched in ASYM (inflammation)"),
    name = NULL) +
  # 通路标签放在柱子的另一侧
  geom_text(aes(label = name,
                x = ifelse(diff_strength > 0, -0.5, 0.5),
                hjust = ifelse(diff_strength > 0, 1, 0)),
            size = 3.3, fontface = "bold", color = "grey15") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(title = "B. Top 15 differential signaling pathways",
       subtitle = "SYM vs ASYM (positive = enriched in SYM)",
       x = "Information flow difference (SYM ??? ASYM)",
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.text.x   = element_text(color = "black"),
    axis.line     = element_line(linewidth = 0.4),
    legend.position = "top",
    legend.text     = element_text(size = 9, face = "bold")
  )

ggsave(paste0(out_dir, "Fig5B_v2.pdf"), p_5B,
       width = 8, height = 6, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig5B_v2.png"), p_5B,
       width = 8, height = 6, dpi = 600)

cat("\n??? Fig 5B v2 saved\n")

--------------------------------------------------------------------------------
# FigureC — Three-panel layered architecture
# 三栏独立子图，每栏一个 layer，sender-receiver 关系清晰

suppressMessages({
  library(CellChat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

setwd("G:/13sample")
out_dir <- "G:/13sample/Fig5_v2/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cellchat_sym <- readRDS("cellchat/cellchat_sym.rds")

# ── 1. 数据准备（沿用之前的 top_pairs_all）─────────────────────────
layer_pathways <- list(
  "Layer 1: Initiator"   = c("PERIOSTIN", "TENASCIN"),
  "Layer 2: Execute"     = c("COLLAGEN", "FN1", "LAMININ", "SPP1"),
  "Layer 3: Integrate"   = c("MHC-II", "APP", "CXCL")
)

core_celltypes <- c(
  "Fibroblasts1", "Fibroblasts2", "Proliferating cells2",
  "SMCs1", "SMCs3", "Macrophages2", "ECs", "T cells1"
)

celltype_short <- c(
  "Fibroblasts1" = "Fibro1", "Fibroblasts2" = "Fibro2",
  "Proliferating cells2" = "PRFs",
  "SMCs1" = "SMC1", "SMCs3" = "SMC3",
  "Macrophages2" = "mLAM", "ECs" = "ECs", "T cells1" = "T cells"
)

# 提取每层 top pairs
get_top_pairs <- function(obj, pathways, cells, top_n = 5) {
  pairs <- data.frame()
  for (pwy in intersect(pathways, obj@netP$pathways)) {
    pwy_idx <- which(obj@netP$pathways == pwy)
    mat <- obj@netP$prob[, , pwy_idx]
    common_r <- intersect(rownames(mat), cells)
    common_c <- intersect(colnames(mat), cells)
    if (length(common_r) == 0 || length(common_c) == 0) next
    sub <- mat[common_r, common_c, drop = FALSE]
    df_p <- as.data.frame(as.table(sub))
    colnames(df_p) <- c("from", "to", "strength")
    df_p$pathway <- pwy
    df_p <- df_p[df_p$strength > 0, ]
    pairs <- bind_rows(pairs, df_p)
  }
  if (nrow(pairs) == 0) return(NULL)
  pairs %>% arrange(desc(strength)) %>% head(top_n)
}

# ── 2. 通用绘图函数（每个 layer 一个 panel）────────────────────────
# ── 2. 优化后的通用绘图函数（增强箭头、防文字重叠、提升高级感） ────────────────────────
# ── 2. 终极优化版绘图函数：引入高级曲线、扩充安全边距防止文字被裁 ────────────────────────
plot_layer <- function(layer_idx, layer_color, sub_title) {
  layer_name <- names(layer_pathways)[layer_idx]
  pwys <- layer_pathways[[layer_idx]]
  
  pairs <- get_top_pairs(cellchat_sym, pwys, core_celltypes, top_n = 5)
  if (is.null(pairs)) return(NULL)
  
  involved <- unique(c(as.character(pairs$from), as.character(pairs$to)))
  involved <- intersect(core_celltypes, involved)
  
  cell_roles <- data.frame(cell = involved) %>%
    rowwise() %>%
    mutate(
      out_str = sum(pairs$strength[as.character(pairs$from) == cell & as.character(pairs$to) != cell]),
      in_str  = sum(pairs$strength[as.character(pairs$to)   == cell & as.character(pairs$from) != cell]),
      auto_str = sum(pairs$strength[as.character(pairs$from) == cell & as.character(pairs$to)   == cell]),
      total = out_str + in_str + auto_str
    ) %>%
    ungroup() %>%
    mutate(role = case_when(
      out_str > in_str & out_str > 0  ~ "Sender",
      in_str  > out_str & in_str  > 0 ~ "Receiver",
      auto_str > 0                    ~ "Autocrine",
      TRUE                            ~ "Both"
    ))
  
  senders   <- cell_roles$cell[cell_roles$out_str > 0]
  receivers <- cell_roles$cell[cell_roles$in_str  > 0 & !cell_roles$cell %in% senders]
  
  sender_only   <- setdiff(senders, cell_roles$cell[cell_roles$in_str > 0])
  receiver_only <- setdiff(receivers, senders)
  both_role     <- intersect(senders, cell_roles$cell[cell_roles$in_str > 0])
  
  cell_pos <- data.frame()
  if (length(sender_only) > 0) {
    cell_pos <- rbind(cell_pos, data.frame(cell = sender_only, x = 0,
                                           y = seq(1, length(sender_only), by = 1)))
  }
  if (length(both_role) > 0) {
    cell_pos <- rbind(cell_pos, data.frame(cell = both_role, x = 1,
                                           y = seq(1, length(both_role), by = 1) + 0.3))
  }
  if (length(receiver_only) > 0) {
    cell_pos <- rbind(cell_pos, data.frame(cell = receiver_only, x = 2,
                                           y = seq(1, length(receiver_only), by = 1)))
  }
  
  cell_pos$y <- cell_pos$y - mean(cell_pos$y)
  cell_pos <- cell_pos %>% left_join(cell_roles, by = "cell") %>%
    mutate(short = celltype_short[cell])
  
  # ★ 关键诊断：打印 cell_pos 检查是否有异常
  cat("\n--- Layer", layer_idx, "cell_pos ---\n")
  print(cell_pos)
  
  pairs_plot <- pairs %>%
    mutate(from = as.character(from), to = as.character(to)) %>%
    left_join(cell_pos %>% dplyr::select(cell, x, y) %>%
                dplyr::rename(x_from = x, y_from = y), by = c("from" = "cell")) %>%
    left_join(cell_pos %>% dplyr::select(cell, x, y) %>%
                dplyr::rename(x_to = x, y_to = y), by = c("to" = "cell"))
  
  auto_pairs <- pairs_plot %>% filter(from == to)
  para_pairs <- pairs_plot %>% filter(from != to)
  
  # ★ 修复 1：给每条 paracrine 边分配不同的 curvature 让它们错开
  if (nrow(para_pairs) > 0) {
    para_pairs <- para_pairs %>%
      mutate(
        dx = x_to - x_from,
        dy = y_to - y_from,
        len = sqrt(dx^2 + dy^2),
        x_from_new = x_from + (dx / len) * 0.14,
        y_from_new = y_from + (dy / len) * 0.14,
        x_to_new   = x_to - (dx / len) * 0.16,
        y_to_new   = y_to - (dy / len) * 0.16,
        pair_id = row_number()
      )
  }
  
  y_min <- min(cell_pos$y) - 0.55
  y_max <- max(cell_pos$y) + 0.55
  
  p <- ggplot()
  
  # ★ 修复 2：用循环画曲线，每条用不同 curvature
  # curvature 在 -0.3 到 0.3 之间分布，确保多条线之间分开
  if (nrow(para_pairs) > 0) {
    n_edges <- nrow(para_pairs)
    # 生成不同的弧度值
    curv_values <- seq(-0.25, 0.25, length.out = n_edges)
    if (n_edges == 1) curv_values <- 0.08
    
    for (i in seq_len(n_edges)) {
      edge <- para_pairs[i, ]
      p <- p + geom_curve(data = edge,
                          aes(x = x_from_new, y = y_from_new,
                              xend = x_to_new, yend = y_to_new,
                              linewidth = strength, alpha = strength),
                          color = layer_color,
                          curvature = curv_values[i],
                          arrow = arrow(length = unit(0.18, "cm"),
                                        type = "closed"),
                          inherit.aes = FALSE)
    }
  }
  
  # 自分泌环
  if (nrow(auto_pairs) > 0) {
    auto_pos <- auto_pairs %>%
      mutate(from = as.character(from)) %>%
      left_join(cell_pos %>% dplyr::select(cell, x, y), by = c("from" = "cell"))
    
    p <- p + 
      geom_point(data = auto_pos, aes(x = x + 0.16, y = y + 0.16),
                 shape = 21, size = 5.5, color = layer_color,
                 fill = "white", stroke = 0.9, inherit.aes = FALSE) +
      geom_text(data = auto_pos, aes(x = x + 0.16, y = y + 0.16),
                label = "???", color = layer_color, size = 3.5,
                fontface = "bold", inherit.aes = FALSE)
  }
  
  # 节点
  p <- p +
    geom_point(data = cell_pos, aes(x = x, y = y, size = total),
               shape = 21, fill = layer_color, color = "black",
               stroke = 0.6, alpha = 0.95, inherit.aes = FALSE)
  
  # 标签
  p <- p +
    geom_text(data = cell_pos,
              aes(x = x, y = y - 0.28, label = short),
              fontface = "bold", color = "black", size = 3.8,
              inherit.aes = FALSE)
  
  # ★ 修复 3：缩小 plot.margin + 显式 panel.background
  p <- p +
    scale_linewidth_continuous(range = c(0.6, 3.0), guide = "none") + 
    scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
    scale_size_continuous(range = c(6, 13), guide = "none") +
    coord_cartesian(xlim = c(-0.5, 2.5),
                    ylim = c(y_min, y_max),
                    clip = "off") +
    labs(title = layer_name, subtitle = sub_title) +
    theme_void(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 12,
                                       color = layer_color),
          plot.subtitle = element_text(size = 9.5, color = "grey30",
                                       lineheight = 1.1),
          plot.margin   = margin(8, 12, 8, 12),   # ★ 减小 margin
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA))
  
  return(p)
}

# ── 3. 生成 3 个 layer panel ───────────────────────────────────────
p_L1 <- plot_layer(1, "#8B0000",
                   "PERIOSTIN + TENASCIN\nNiche-initiating signals")

p_L2 <- plot_layer(2, "#C04040",
                   "COLLAGEN, FN1, LAMININ, SPP1\nPRFs-centric ECM hub")

p_L3 <- plot_layer(3, "#3A6FA0",
                   "MHC-II, APP, CXCL\nImmune integration")

p_5C_v5 <- p_L1 | p_L2 | p_L3
p_5C_v5 <- p_5C_v5 + plot_annotation(
  title = "C. Three-layered communication architecture",
  subtitle = "Each layer shows top L-R signaling pairs in SYM | Arrow = paracrine | ??? = autocrine | Node size = total strength",
  theme = theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )
)

ggsave(paste0(out_dir, "Fig5C_v5.pdf"), p_5C_v5,
       width = 14, height = 5, device = cairo_pdf)   # ★ height: 6 → 5
ggsave(paste0(out_dir, "Fig5C_v5.png"), p_5C_v5,
       width = 14, height = 5, dpi = 600)

cat("\n??? Fig 5C v5 saved\n")


-------------------------------------------------------------------------------
# FigureD and E — Pathway strength: ASYM vs SYM (6 pathways)

suppressMessages({
  library(CellChat)
  library(tidyverse)
  library(ggplot2)
  library(ggsignif)
  library(patchwork)
})

setwd("G:/13sample")
out_dir <- "G:/13sample/Fig5_v2/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cellchat_asym <- readRDS("cellchat/cellchat_asym.rds")
cellchat_sym  <- readRDS("cellchat/cellchat_sym.rds")

# Step 1: 定义 9 条核心通路
# ════════════════════════════════════════════════════════════════════
pathway_layers <- data.frame(
  pathway = c("PERIOSTIN", "TENASCIN",
              "COLLAGEN", "FN1", "LAMININ", "SPP1",
              "MHC-II",    "APP",  "CXCL"),
  layer   = c("Layer 1", "Layer 1",
              "Layer 2", "Layer 2", "Layer 2", "Layer 2",
              "Layer 3", "Layer 3", "Layer 3"),
  stringsAsFactors = FALSE
)


# ════════════════════════════════════════════════════════════════════
# Step 2: 计算每条通路的总通讯强度（ASYM 和 SYM）
# ════════════════════════════════════════════════════════════════════
calc_strength <- function(obj, pathway) {
  if (!(pathway %in% obj@netP$pathways)) return(0)
  pwy_idx <- which(obj@netP$pathways == pathway)
  mat <- obj@netP$prob[, , pwy_idx]
  return(sum(mat, na.rm = TRUE))
}

# 构建原始数据框
df_raw <- pathway_layers %>%
  rowwise() %>%
  mutate(
    ASYM = calc_strength(cellchat_asym, pathway),
    SYM  = calc_strength(cellchat_sym,  pathway)
  ) %>%
  ungroup() %>%
  mutate(fc = SYM / ASYM)

cat("=== Raw pathway strength ===\n")
print(df_raw)


# ════════════════════════════════════════════════════════════════════
# Step 3: 转长格式 + 归一化（每条通路 ASYM = 1.0）
# ════════════════════════════════════════════════════════════════════
df_long <- df_raw %>%
  pivot_longer(cols = c(ASYM, SYM),
               names_to = "group",
               values_to = "strength") %>%
  group_by(pathway) %>%
  mutate(rel_strength = strength / strength[group == "ASYM"]) %>%
  ungroup() %>%
  mutate(
    group = factor(group, levels = c("ASYM", "SYM")),
    # 高幅 vs 中等幅
    panel_side = ifelse(fc > 1.4, "Left", "Right")
  )

cat("\n=== Normalized data ===\n")
print(df_long)


# ════════════════════════════════════════════════════════════════════
# Step 4: 配色（9 条通路，按 layer 分色系）
# ════════════════════════════════════════════════════════════════════
pathway_col <- c(
  # Layer 1 - 深红系
  "PERIOSTIN" = "#8B0000",
  "TENASCIN"  = "#C04040",
  # Layer 2 - 红橙系
  "COLLAGEN"  = "#E59085",
  "FN1"       = "#F0B5AD",
  "LAMININ"   = "#F8C8B0",
  "SPP1"      = "#D4906E",
  # Layer 3 - 蓝系
  "MHC-II"    = "#3A6FA0",
  "APP"       = "#7FA5C7",
  "CXCL"      = "#A5C2D9"
)


# ════════════════════════════════════════════════════════════════════
# Step 5: 通用主题
# ════════════════════════════════════════════════════════════════════
theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9.5, color = "grey40"),
    axis.text          = element_text(color = "black"),
    axis.text.x        = element_text(face = "bold", size = 12),
    axis.text.y        = element_text(size = 10),
    axis.line          = element_line(linewidth = 0.4),
    axis.ticks         = element_line(linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    legend.position    = "none",
    plot.margin        = margin(10, 10, 10, 10)
  )


# ════════════════════════════════════════════════════════════════════
# Step 6: 左 panel — 高幅通路 (PERIOSTIN, CXCL, LAMININ, SPP1)
# ════════════════════════════════════════════════════════════════════
df_L <- df_long %>% filter(panel_side == "Left")
labels_L <- df_L %>% filter(group == "SYM")

p_L <- ggplot(df_L,
              aes(x = group, y = rel_strength,
                  group = pathway, color = pathway)) +
  # 基准线
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  # 线
  geom_line(linewidth = 1.4, alpha = 0.92) +
  # 点
  geom_point(size = 4.2, shape = 21, fill = "white", stroke = 1.4) +
  # 标签（不用 ggrepel，手动定位避免幽灵点）
  geom_text(data = labels_L,
            aes(label = sprintf("%s ×%.2f", pathway, rel_strength)),
            x = 2.18, hjust = 0,
            fontface = "bold", size = 3.8) +
  scale_color_manual(values = pathway_col) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.85))) +
  scale_y_continuous(breaks = c(1, 1.5, 2, 2.5),
                     labels = c("1×", "1.5×", "2×", "2.5×"),
                     limits = c(0.95, 2.5)) +
  labs(title = "High-amplitude pathways",
       subtitle = "Large relative gain (fc > 1.4)",
       x = NULL,
       y = "Relative strength (vs ASYM)") +
  theme_panel


# ════════════════════════════════════════════════════════════════════
# Step 7: 右 panel — 中等幅通路 (COLLAGEN, FN1, TENASCIN, MHC-II, APP)
# ════════════════════════════════════════════════════════════════════
df_R <- df_long %>% filter(panel_side == "Right")
labels_R <- df_R %>% filter(group == "SYM") %>%
  arrange(desc(rel_strength)) %>%
  # 给每条线的标签做 y 微调（避免重叠）
  mutate(label_y = rel_strength + c(0.005, 0.000, 0.000, -0.005, -0.005))

p_R <- ggplot(df_R,
              aes(x = group, y = rel_strength,
                  group = pathway, color = pathway)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 1.4, alpha = 0.92) +
  geom_point(size = 4.2, shape = 21, fill = "white", stroke = 1.4) +
  geom_text(data = labels_R,
            aes(label = sprintf("%s ×%.2f", pathway, rel_strength),
                y = label_y),
            x = 2.18, hjust = 0,
            fontface = "bold", size = 3.8) +
  scale_color_manual(values = pathway_col) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.85))) +
  scale_y_continuous(breaks = seq(1, 1.20, by = 0.05),
                     labels = c("1×", "1.05×", "1.10×", "1.15×", "1.20×"),
                     limits = c(0.99, 1.22)) +
  labs(title = "Moderate-amplitude pathways",
       subtitle = "Concordant system-wide upregulation",
       x = NULL,
       y = NULL) +
  theme_panel


# ════════════════════════════════════════════════════════════════════
# Step 8: 拼接 + 保存
# ════════════════════════════════════════════════════════════════════
p_5D <- p_L + p_R + plot_layout(widths = c(1, 1.5))

p_5D <- p_5D + plot_annotation(
  title = "D. Pathway-level relative change from ASYM to SYM",
  subtitle = "Nine core pathways | Each normalized to its own ASYM baseline (×1.00)",
  theme = theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )
)

ggsave(paste0(out_dir, "Fig5Dv5.pdf"), p_5D,
       width = 12, height = 5.5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig5Dv5.png"), p_5D,
       width = 12, height = 5.5, dpi = 600)

cat("\n??? Figure 5D saved to:", out_dir, "\n")

--------------------------------------------------------------------------------
# FigureF — 中心性 ASYM vs SYM 对比 (移动轨迹散点)
# 展示从 ASYM → SYM, 哪些细胞向"高发送+高接收"(hub地位)移动
suppressMessages({
  library(CellChat); library(tidyverse); library(ggplot2); library(ggrepel)
})

setwd("G:/13sample")
out_dir <- "G:/13sample/Fig5_centrality/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 读取两组 cellchat ────────────────────────────────────────────
cc_asym <- readRDS("cellchat/cellchat_asym.rds")
cc_sym  <- readRDS("cellchat/cellchat_sym.rds")
cc_asym <- netAnalysis_computeCentrality(cc_asym, slot.name = "netP")
cc_sym  <- netAnalysis_computeCentrality(cc_sym,  slot.name = "netP")

core_celltypes <- c("Fibroblasts1","Fibroblasts2","Proliferating cells2",
                    "SMCs1","SMCs3","Macrophages2","ECs","T cells1")
celltype_short <- c("Fibroblasts1"="Fibro1","Fibroblasts2"="Fibro2",
                    "Proliferating cells2"="PRFs","SMCs1"="SMC1","SMCs3"="SMC3",
                    "Macrophages2"="mLAM","ECs"="ECs","T cells1"="T cells")

# ── 通用: 提取某 cellchat 的 Sender/Receiver ─────────────────────
get_SR <- function(cc, grp) {
  centr <- cc@netP$centr
  all_cells <- levels(cc@idents); n <- length(all_cells)
  agg <- function(role) {
    mat <- sapply(centr, function(x){
      v <- x[[role]]; out <- setNames(rep(0,n), all_cells); out[names(v)] <- v; out
    }); rowSums(mat)
  }
  data.frame(cell=all_cells, Sender=agg("outdeg"), Receiver=agg("indeg"), group=grp) %>%
    filter(cell %in% core_celltypes) %>% mutate(short=celltype_short[cell])
}

sr_asym <- get_SR(cc_asym, "ASYM")
sr_sym  <- get_SR(cc_sym,  "SYM")
sr_all  <- bind_rows(sr_asym, sr_sym)

# ── 核对表 ───────────────────────────────────────────────────────
compare_tbl <- sr_asym %>% dplyr::select(short, Sender_ASYM=Sender, Receiver_ASYM=Receiver) %>%
  left_join(sr_sym %>% dplyr::select(short, Sender_SYM=Sender, Receiver_SYM=Receiver), by="short") %>%
  mutate(dSender = Sender_SYM - Sender_ASYM,
         dReceiver = Receiver_SYM - Receiver_ASYM,
         total_shift = sqrt(dSender^2 + dReceiver^2)) %>%
  arrange(desc(total_shift))
cat("\n=== ASYM → SYM 中心性变化 (按移动距离排序) ===\n")
print(compare_tbl %>% mutate(across(where(is.numeric), ~round(.,2))))
write.csv(compare_tbl %>% mutate(across(where(is.numeric),~round(.,4))),
          paste0(out_dir,"centrality_ASYM_vs_SYM.csv"), row.names=FALSE)

# ── 准备轨迹数据 (每个细胞 ASYM点→SYM点) ─────────────────────────
traj <- sr_asym %>% dplyr::select(short, x0=Sender, y0=Receiver) %>%
  left_join(sr_sym %>% dplyr::select(short, x1=Sender, y1=Receiver), by="short")
# 高亮 hub
hub_cols <- c("PRFs"="#00A878","mLAM"="#3A6FA0","Other"="#C0C0C0")
traj$col <- ifelse(traj$short=="PRFs","PRFs", ifelse(traj$short=="mLAM","mLAM","Other"))
sr_all$col <- ifelse(sr_all$short=="PRFs","PRFs", ifelse(sr_all$short=="mLAM","mLAM","Other"))
sr_sym$col <- ifelse(sr_sym$short=="PRFs","PRFs", ifelse(sr_sym$short=="mLAM","mLAM","Other"))

# ── 移动轨迹散点 ─────────────────────────────────────────────────
p_traj <- ggplot() +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey85") +
  # 箭头: ASYM -> SYM
  geom_segment(data=traj,
               aes(x=x0, y=y0, xend=x1, yend=y1, color=col),
               arrow=arrow(length=unit(0.25,"cm"), type="closed"),
               linewidth=0.9, alpha=0.8) +
  # ASYM 点 (空心)
  geom_point(data=sr_all %>% filter(group=="ASYM"),
             aes(Sender, Receiver, color=col), shape=1, size=3.5, stroke=1.2) +
  # SYM 点 (实心)
  geom_point(data=sr_all %>% filter(group=="SYM"),
             aes(Sender, Receiver, fill=col), shape=21, size=5, color="white", stroke=0.8) +
  # 标签 (标在SYM终点)
  geom_text_repel(data=sr_sym, aes(Sender, Receiver, label=short, color=col),
                  fontface="bold", size=4, box.padding=0.6, seed=1, show.legend=FALSE) +
  scale_color_manual(values=hub_cols, guide="none") +
  scale_fill_manual(values=hub_cols, guide="none") +
  labs(title="Shift in network centrality from ASYM to SYM",
       subtitle="Open circle = ASYM, filled circle = SYM; arrow shows the trajectory. PRFs and mLAM move toward hub positions.",
       x="Outgoing strength (Sender)", y="Incoming strength (Receiver)") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=9,color="grey40"),
        axis.text=element_text(color="black"),
        axis.line=element_line(linewidth=0.4))

ggsave(paste0(out_dir,"centrality_trajectory.pdf"), p_traj, width=7.5, height=6.5, device=cairo_pdf)
ggsave(paste0(out_dir,"centrality_trajectory.png"), p_traj, width=7.5, height=6.5, dpi=600)

# ── 补充: 分组条形图 (每细胞 Sender/Receiver 在两组对比) ──────────
bar_df <- sr_all %>%
  pivot_longer(c(Sender,Receiver), names_to="role", values_to="score")
bar_df$group <- factor(bar_df$group, levels=c("ASYM","SYM"))
# 细胞按 SYM receiver 排序
ord <- sr_sym %>% arrange(desc(Receiver)) %>% pull(short)
bar_df$short <- factor(bar_df$short, levels=ord)

p_bar <- ggplot(bar_df, aes(short, score, fill=group)) +
  geom_col(position=position_dodge(0.75), width=0.7, color="white", linewidth=0.3) +
  facet_wrap(~role, ncol=1, scales="free_y") +
  scale_fill_manual(values=c("ASYM"="#E8A33D","SYM"="#C0392B"), name=NULL) +
  labs(title="Sender / receiver strength: ASYM vs SYM",
       x=NULL, y="Interaction strength") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold",size=12),
        axis.text.x=element_text(angle=35, hjust=1, face="bold", color="black"),
        axis.text.y=element_text(color="black"),
        strip.text=element_text(face="bold", size=11),
        strip.background=element_rect(fill="grey92", color=NA),
        legend.position="top")

ggsave(paste0(out_dir,"centrality_barplot.pdf"), p_bar, width=7, height=6.5, device=cairo_pdf)
ggsave(paste0(out_dir,"centrality_barplot.png"), p_bar, width=7, height=6.5, dpi=600)

cat("\n完成. 输出:\n  centrality_trajectory (移动轨迹散点-主图)\n  centrality_barplot (分组条形图-备选)\n  centrality_ASYM_vs_SYM.csv (核对表)\n")
print(p_traj); print(p_bar)
  
--------------------------------------------------------------------------------
##---FigG and H_candidate_donor_balanced_lineage_curves
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

  
  
  
  
