# Project root: the folder containing SMC/, macro/, Fig3/, Fig5_CellChat/, etc.
# Every path in this script is built from this one setting, so the code no longer
# depends on a single hard-coded drive letter. Either start R/RStudio in the
# project root, or edit the one line below to point at it directly.
# Stored as an R option so it survives the rm(list = ls()) calls between sections.
options(cea.project_root = normalizePath(".", mustWork = FALSE))

#Supplemental Figure 2

## FigureA 
##------------细胞群之间的相关性热图--------------
av <-AverageExpression(combined131,group.by = "celltype", assays = "RNA")
av<-av[[1]]

cg<-names(tail(sort(apply(av, 1, sd)),1000))#选出标准差最大的1000个基因
# 查看这1000个基因在各细胞群中的表达矩阵
expression_matrix <- av[cg,]

# 确保 expression_matrix 是矩阵
if (!is.matrix(expression_matrix)) {
  expression_matrix <- as.matrix(expression_matrix)
}

# 计算细胞群的相关性矩阵（使用Spearman相关性）
cor_matrix <- cor(expression_matrix, method = 'spearman')

#pheatmap绘制热图
library(pheatmap)

# 生成热图
pdf("combined131/celltype_correlation_heatmap.pdf", width = 10, height = 8)
pheatmap(cor_matrix, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         clustering_method = "complete", 
         display_numbers = F,
         main = "Cell Type Correlation Heatmap (Spearman)")
dev.off()


# 生成调色板，类似于Nature风格
library(RColorBrewer)
nature_colors <- colorRampPalette(c("blue", "white", "red"))(100)
# 生成热图并保存为PDF
pdf("combined131/celltype_correlation_heatmap_nature_style.pdf", width = 10, height = 8)
pheatmap(cor_matrix, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         clustering_method = "complete", 
         display_numbers = F,
         color = nature_colors,
         main = "Cell Type Correlation Heatmap (Spearman, Nature Style)")
dev.off()


--------------------------------------------------------------------------------
#----------- Fig B UMAP  ggrastr 栅格化
#

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrastr)        # ★ 关键
})

out_dir <- paste0(getOption("cea.project_root"), "/Alsaigh_integration/")

# 假设 combined 对象已加载，含 dataset / your_celltype / seurat_clusters / umap_1 / umap_2

# ── 提取 UMAP 坐标 ────────────────────────────────────────────────
umap_coords <- as.data.frame(Embeddings(combined, "umap"))
umap_coords$dataset <- combined$dataset
umap_coords$your_celltype <- combined$your_celltype
umap_coords$seurat_clusters <- combined$seurat_clusters

umap_coords$celltype_simple <- case_when(
  grepl("^SMC", umap_coords$your_celltype) ~ "SMCs",
  grepl("^Fibroblast", umap_coords$your_celltype) ~ "Fibroblasts",
  grepl("^Macrophage", umap_coords$your_celltype) ~ "Macrophages",
  grepl("^T cell", umap_coords$your_celltype) ~ "T cells",
  grepl("^B cell", umap_coords$your_celltype) ~ "B cells",
  grepl("^Proliferating", umap_coords$your_celltype) ~ "Proliferating",
  TRUE ~ umap_coords$your_celltype
)

simple_colors <- c(
  "SMCs"="#3C999CFF", "Fibroblasts"="#B46356FF", "Macrophages"="#CA480DFF",
  "Monocytes"="#F9AB3CFF", "DCs"="#13A4EBFF", "T cells"="#944D40FF",
  "B cells"="#A1B3E2FF", "NK cells"="#2F7270FF", "Neutrophils"="#FFF7A3FF",
  "Mast cells"="#8BEAFFFF", "ECs"="#FFE705FF", "Proliferating"="#524559FF"
)

# ── Panel 1A: By dataset (栅格化) ─────────────────────────────────
p_A <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = dataset)) +
  ggrastr::rasterise(                                                # ★ 关键
    geom_point(size = 0.25, alpha = 0.75),
    dpi = 600
  ) +
  scale_color_manual(values = c("Our_study" = "#C71000FF",
                                "Alsaigh_2022" = "#3A6FA0FF")) +
  labs(title = "Dataset (Harmony-integrated)",
       x = "umap_1", y = "umap_2") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# ── Panel 1B: Our cell types overlaid (栅格化) ────────────────────
p_B <- ggplot() +
  # 底层：Alsaigh 灰色背景（栅格化）
  ggrastr::rasterise(
    geom_point(data = umap_coords[umap_coords$dataset == "Alsaigh_2022", ],
               aes(x = umap_1, y = umap_2),
               color = "grey85", size = 0.15, alpha = 0.6),
    dpi = 600
  ) +
  # 顶层：Our 有色（栅格化）
  ggrastr::rasterise(
    geom_point(data = umap_coords[umap_coords$dataset == "Our_study" &
                                    !is.na(umap_coords$celltype_simple), ],
               aes(x = umap_1, y = umap_2, color = celltype_simple),
               size = 0.4, alpha = 0.85),
    dpi = 600
  ) +
  scale_color_manual(values = simple_colors, name = NULL, na.translate = FALSE) +
  labs(title = "Our cell types overlaid on integrated space",
       subtitle = "Grey background = Alsaigh cells",
       x = "umap_1", y = "umap_2") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0.5),
        legend.position = "right",
        legend.text   = element_text(size = 9, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

# ── Panel 1C: Integrated clusters (栅格化) ────────────────────────
cluster_centers <- umap_coords %>%
  group_by(seurat_clusters) %>%
  summarise(x = median(umap_1), y = median(umap_2), n = n(), .groups = "drop") %>%
  filter(n > 50)

p_C <- ggplot(umap_coords, aes(x = umap_1, y = umap_2,
                               color = factor(seurat_clusters))) +
  ggrastr::rasterise(                                                # ★ 关键
    geom_point(size = 0.25, alpha = 0.85),
    dpi = 600
  ) +
  geom_label(data = cluster_centers,                                 # 矢量保留
             aes(x = x, y = y, label = seurat_clusters),
             fill = "white", color = "black",
             size = 3, fontface = "bold",
             label.padding = unit(0.15, "lines"),
             alpha = 0.85,
             inherit.aes = FALSE) +
  labs(title = "Integrated clusters (Seurat resolution 0.5)",
       x = "umap_1", y = "umap_2") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "none")

# ── 拼图 ──────────────────────────────────────────────────────────
fig1F_umap <- p_A | p_B | p_C
fig1F_umap <- fig1F_umap + plot_annotation(
  title = "F. Integration with Alsaigh 2022 (GSE159677) confirms cell type taxonomy",
  subtitle = sprintf("Harmony integration of %s cells (Our: %s + Alsaigh: %s) from 15 patients",
                     format(ncol(combined), big.mark = ","),
                     format(sum(combined$dataset == "Our_study"), big.mark = ","),
                     format(sum(combined$dataset == "Alsaigh_2022"), big.mark = ",")),
  theme = theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10.5, color = "grey40")
  )
)

ggsave(paste0(out_dir, "Fig1F_v2_UMAP_rasterized.pdf"), fig1F_umap,
       width = 17, height = 5.5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig1F_v2_UMAP_rasterized.png"), fig1F_umap,
       width = 17, height = 5.5, dpi = 600)

cat("\n??? Fig 1F UMAP (rasterized) saved\n") 
  
--------------------------------------------------------------------------------
  
  
# Fig C-D — 组间差异可视化
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggsignif)
  library(ggdist)         # raincloud / half-violin (现代风)
  library(ggbeeswarm)     # sina plot
  library(gghalves)       # half-violin 拼接
})

out_dir <- paste0(getOption("cea.project_root"), "/Alsaigh_integration/")

# 配色
group_colors <- c(
  "ASYM"    = "#FFE705FF",   # 黄
  "SYM"     = "#C71000FF"    # 红
)
tissue_colors <- c(
  "control" = "#A8A8A8",
  "plaque"  = "#C71000FF"
)

# ============================================================================
# 数据准备（接续您的整合脚本）
# ============================================================================
if (!exists("combined")) {
  combined <- readRDS(paste0(out_dir, "integrated_obj_with_Alsaigh.rds"))
}

fibroblast_smc_clusters <- c("3","5","7","9","11","12","16","20")
macrophage_clusters     <- c("4","6","13","14")

combined$broad_type <- "Other"
combined$broad_type[combined$seurat_clusters %in% fibroblast_smc_clusters] <- "Fibroblast_SMC"
combined$broad_type[combined$seurat_clusters %in% macrophage_clusters]     <- "Macrophage"


#方案 1：RAINCLOUD PLOT（云雨图）—— 最推荐

make_raincloud <- function(df, x_var, y_var, fill_colors,
                           comparisons_list, title, y_label) {
  
  ggplot(df, aes_string(x = x_var, y = y_var, fill = x_var, color = x_var)) +
    
    # 雨：散点（左下方）
    ggdist::stat_dots(
      side = "left",
      justification = 1.12,
      binwidth = NA,
      dotsize = 0.8,
      alpha = 0.35
    ) +
    
    # 云：半小提琴（右上方密度分布）
    ggdist::stat_halfeye(
      adjust = 0.8,
      width = 0.6,
      .width = 0,           # 不画中位数线
      justification = -0.15,
      point_colour = NA,
      alpha = 0.7
    ) +
    
    # 箱：迷你箱线图（中间）
    geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.5,
      color = "black",
      size = 0.4
    ) +
    
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    
    # P 值标注
    geom_signif(
      comparisons = comparisons_list,
      map_signif_level = function(p) {
        if (p < 0.001) "***"
        else if (p < 0.01) "**"
        else if (p < 0.05) "*"
        else "ns"
      },
      test = "wilcox.test",
      textsize = 5,
      tip_length = 0.01,
      step_increase = 0.08,
      color = "black"
    ) +
    
    labs(title = title, y = y_label, x = "") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 11, color = "black"),
      axis.title.y = element_text(size = 11),
      axis.line = element_line(size = 0.6)
    )
}

# 4 个 Panel
fibro_alsaigh <- combined@meta.data %>%
  filter(broad_type == "Fibroblast_SMC", dataset == "Alsaigh_2022") %>%
  mutate(tissue_type = factor(tissue_type, levels = c("control","plaque")))

fibro_ours <- combined@meta.data %>%
  filter(broad_type == "Fibroblast_SMC", dataset == "Our_study", !is.na(group)) %>%
  mutate(group = factor(group, levels = c("ASYM","SYM")))

mac_alsaigh <- combined@meta.data %>%
  filter(broad_type == "Macrophage", dataset == "Alsaigh_2022") %>%
  mutate(tissue_type = factor(tissue_type, levels = c("control","plaque")))

mac_ours <- combined@meta.data %>%
  filter(broad_type == "Macrophage", dataset == "Our_study", !is.na(group)) %>%
  mutate(group = factor(group, levels = c("ASYM","SYM")))

p1 <- make_raincloud(fibro_alsaigh, "tissue_type", "FibroHealing1", tissue_colors,
                     list(c("control","plaque")),
                     "Alsaigh 2022 (Tissue-type)", "Fibroblast healing score")

p2 <- make_raincloud(fibro_ours, "group", "FibroHealing1", group_colors,
                     list( c("ASYM","SYM")),
                     "Our study (Pathology-anchored)", "Fibroblast healing score")

p3 <- make_raincloud(mac_alsaigh, "tissue_type", "TREM2_LAM2", tissue_colors,
                     list(c("control","plaque")),
                     "Alsaigh 2022 (Tissue-type)", "TREM2+ LAM score")

p4 <- make_raincloud(mac_ours, "group", "TREM2_LAM2", group_colors,
                     list( c("ASYM","SYM")),
                     "Our study (Pathology-anchored)", "TREM2+ LAM score")

fig1F_raincloud <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Fig 1F. Pathology-anchored grouping reveals healing program",
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  )

ggsave(paste0(out_dir, "Fig1F_RAINCLOUD.pdf"), fig1F_raincloud,
       width = 13, height = 11, device = "pdf")
ggsave(paste0(out_dir, "Fig1F_RAINCLOUD.png"), fig1F_raincloud,
       width = 13, height = 11, dpi = 300)


##方案 2：Half-Violin + Dot + Boxplot（三合一精简版）

# 1. 重新定义函数，使用 ggdist 替代 gghalves
make_halfviolin <- function(df, x_var, y_var, fill_colors,
                            comparisons_list, title, y_label) {
  
  # 使用官方推荐的 .data[[]] 替代即将被淘汰的 aes_string
  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]], color = .data[[x_var]])) +
    
    # 左半小提琴 (ggdist 语法)
    ggdist::stat_halfeye(
      adjust = 0.5,
      width = 0.5,
      .width = 0,             # 不显示底部的区间条
      justification = 1.1,    # 向左偏移
      point_colour = NA,      # 不显示中位数点
      alpha = 0.7
    ) +
    
    # 右半散点 (ggdist 语法)
    ggdist::stat_dots(
      side = "right",
      justification = -0.15,  # 向右偏移，给箱线图留出空间
      binwidth = NA,
      dotsize = 0.6,
      alpha = 0.4
    ) +
    
    # 中间小箱线
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      size = 0.4
    ) +
    
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    
    # 显著性检验
    geom_signif(
      comparisons = comparisons_list,
      map_signif_level = function(p) {
        if (p < 0.001) "***"
        else if (p < 0.01) "**"
        else if (p < 0.05) "*"
        else "ns"
      },
      test = "wilcox.test",
      textsize = 5,
      tip_length = 0.01,
      step_increase = 0.08,
      color = "black"
    ) +
    
    labs(title = title, y = y_label, x = "") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 11, color = "black"),
      axis.title.y = element_text(size = 11),
      axis.line = element_line(size = 0.6)
    )
}

p1b <- make_halfviolin(fibro_alsaigh, "tissue_type", "FibroHealing1", tissue_colors,
                       list(c("control","plaque")),
                       "Alsaigh 2022", "Fibroblast healing score")

p2b <- make_halfviolin(fibro_ours, "group", "FibroHealing1", group_colors,
                       list( c("ASYM","SYM")),
                       "Our study", "Fibroblast healing score")

p3b <- make_halfviolin(mac_alsaigh, "tissue_type", "TREM2_LAM2", tissue_colors,
                       list(c("control","plaque")),
                       "Alsaigh 2022", "TREM2+ LAM score")

p4b <- make_halfviolin(mac_ours, "group", "TREM2_LAM2", group_colors,
                       list( c("ASYM","SYM")),
                       "Our study", "TREM2+ LAM score")

# 3. 拼图并添加主标题
fig1F_halfviolin <- (p1b | p2b) / (p3b | p4b) +
  plot_annotation(
    title = "Fig 1F. Half-violin + dot + boxplot",
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  )

# 4. 完美保存（此时不再报错）
ggsave(paste0(out_dir, "Fig1F_HALFVIOLIN.pdf"), fig1F_halfviolin,
       width = 13, height = 11, device = "pdf")

ggsave(paste0(out_dir, "Fig1F_HALFVIOLIN.png"), fig1F_halfviolin,
       width = 13, height = 11, dpi = 300)


# ----方案 3：Sina Plot（蜂群分布图）
#Sina Plot — 修长 + 紧凑优化版

make_sinaplot <- function(df, x_var, y_var, fill_colors,
                          comparisons_list, title, y_label) {
  
  ggplot(df, aes_string(x = x_var, y = y_var, color = x_var)) +
    
    # 背景小提琴（更窄）
    geom_violin(aes_string(fill = x_var), 
                alpha = 0.15, color = NA, trim = FALSE,
                width = 0.6) +              # ★ 从默认0.9改为0.6
    
    # Sina：按密度分布的散点（更窄）
    ggbeeswarm::geom_quasirandom(
      width = 0.18,                          # ★ 从0.3改为0.18
      size = 0.4,                            # 略小
      alpha = 0.5
    ) +
    
    # 中位数 + 四分位线（更窄）
    stat_summary(
      fun = median,
      fun.min = function(x) quantile(x, 0.25),
      fun.max = function(x) quantile(x, 0.75),
      geom = "crossbar",
      width = 0.28,                          # ★ 从0.4改为0.28
      color = "black",
      fill = NA,
      size = 0.35
    ) +
    
    # 均值标记
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 2.5,                            # ★ 略小
      color = "black"
    ) +
    
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    
    # P 值标注
    geom_signif(
      comparisons = comparisons_list,
      map_signif_level = function(p) {
        if (p < 0.001) "***"
        else if (p < 0.01) "**"
        else if (p < 0.05) "*"
        else "ns"
      },
      test = "wilcox.test",
      textsize = 4.5,
      tip_length = 0.008,                    # ★ 连接线更短
      step_increase = 0.07,                  # ★ 多组间隔更紧
      color = "black",
      size = 0.35
    ) +
    
    # x 轴留白（防止散点贴边）
    scale_x_discrete(expand = expansion(add = 0.6)) +   # ★ 关键：左右留白
    
    labs(title = title, y = y_label, x = "") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5,
                                margin = margin(b = 4)),
      axis.text.x = element_text(face = "bold", size = 10, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 10, margin = margin(r = 4)),
      axis.line = element_line(size = 0.5),
      axis.ticks = element_line(size = 0.4),
      
      # ★★★ 关键三个参数 ★★★
      aspect.ratio = 1.5,                    # ★ 强制高宽比，让图更修长
      plot.margin = margin(5, 3, 5, 3),      # ★ 单 panel 内边距收紧
      panel.spacing = unit(0.2, "lines")     # ★ panel 之间间距
    )
}


# ============================================================================
# 生成 4 个 panel（数据准备部分保持不变）
# ============================================================================
p1c <- make_sinaplot(fibro_alsaigh, "tissue_type", "FibroHealing1", tissue_colors,
                     list(c("control","plaque")),
                     "Alsaigh 2022", "Fibroblast healing score")

p2c <- make_sinaplot(fibro_ours, "group", "FibroHealing1", group_colors,
                     list( c("ASYM","SYM")),
                     "Our study", "Fibroblast healing score")

p3c <- make_sinaplot(mac_alsaigh, "tissue_type", "TREM2_LAM2", tissue_colors,
                     list(c("control","plaque")),
                     "Alsaigh 2022", "TREM2+ LAM score")

p4c <- make_sinaplot(mac_ours, "group", "TREM2_LAM2", group_colors,
                     list( c("ASYM","SYM")),
                     "Our study", "TREM2+ LAM score")


# ============================================================================
# 拼图：4 panel 紧凑布局 (2×2)
# ============================================================================
fig1F_sina <- (p1c | p2c) / (p3c | p4c) +
  plot_layout(guides = "collect") &        # 统一图例（虽然无图例）
  theme(
    plot.margin = margin(2, 2, 2, 2)       # ★ 每个子图外边距
  )

fig1F_sina <- fig1F_sina +
  plot_annotation(
    title = "Fig 1F. Pathology-anchored grouping reveals healing program",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5,
                                margin = margin(b = 6))
    )
  )


# ============================================================================
# 保存：宽度大幅缩减 + 高度调整
# ============================================================================
ggsave(paste0(out_dir, "Fig1F_SINA_compact.pdf"), fig1F_sina,
       width = 8.5,                          # ★ 从 13 缩到 8.5 寸
       height = 10,                          # ★ 高度保持 → 修长效果
       device = "pdf")

ggsave(paste0(out_dir, "Fig1F_SINA_compact.png"), fig1F_sina,
       width = 8.5, 
       height = 10, 
       dpi = 300)

