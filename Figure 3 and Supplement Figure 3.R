# Project root: the folder containing SMC/, macro/, Fig3/, Fig5_CellChat/, etc.
# Every path in this script is built from this one setting, so the code no longer
# depends on a single hard-coded drive letter. Either start R/RStudio in the
# project root, or edit the one line below to point at it directly.
# Stored as an R option so it survives the rm(list = ls()) calls between sections.
options(cea.project_root = normalizePath(".", mustWork = FALSE))


#---Figure 3. VSMC and fibroblast states define a contractile-to-reparative transcriptional continuum---

rm(list = ls())
setwd(paste0(getOption("cea.project_root"), "/SMC"))
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidydr)
plaque<- readRDS(file="combined131.rds")
##VSMC亚群分析
smc = plaque[, Idents(plaque) %in% c("Fibroblasts1" , "Fibroblasts2" ,
                                     "SMCs1","SMCs2","SMCs3")]

smc<- readRDS(file="smc.rds")


--------------------------------------------------------------------------------
## Figure A
##---------------绘制celltype tsne图
smc_col=c("#B46356FF","#7A4B42FF","#E2ED50FF","#66FFFFFF", "#3C999CFF")

tsne = smc@reductions$tsne @cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(celltype = smc@meta.data$celltype) # 注释后的label信息 ，改为celltype

p <- ggplot(tsne,aes(x= tSNE_1 , y = tSNE_2 ,color = celltype)) +
  geom_point(size = 1 , alpha =1 )+scale_color_manual(values = smc_col)

#三 ggplot2-umap-调整
#3.1 调整umap图 - theme
p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 轴刻度文本
        axis.ticks = element_blank(),# 轴刻度
        axis.line= element_blank(),#坐标轴
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
##3.2 调整umap图 - legend
p3 <- p2 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=3))) #设置legend中 点的大小 

#3.3 调整umap图 - annotation
#坐标轴放到左下角可以通过ggplot2添加箭头和文本实现。
p4 <- p3 + 
  geom_segment(aes(x = min(tsne$tSNE_1) , y = min(tsne$tSNE_2) ,
                   xend = min(tsne$tSNE_1) +10, yend = min(tsne$tSNE_2) ),#X轴线长，+10起点正向＋3
               colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm")))+ #X轴线宽、箭头
  geom_segment(aes(x = min(tsne$tSNE_1)  , y = min(tsne$tSNE_2)  ,
                   xend = min(tsne$tSNE_1) , yend = min(tsne$tSNE_2) + 10),#y轴线长
               colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(tsne$tSNE_1) +5, y = min(tsne$tSNE_2) -2, label = "tSNE_1",
           color="black",size = 3.5, fontface="plain" ) + #X轴注释tSNE_1,+/-数字调整距起点位置
  annotate("text", x = min(tsne$tSNE_1) -2, y = min(tsne$tSNE_2) + 5, label = "tSNE_2",
           color="black",size = 3.5, fontface="plain" ,angle=90) #Y轴注释tSNE_1
p4

ggsave("SMCtsne.PDF",width=8,height=6,device = "pdf")

--------------------------------------------------------------------------------
## FigureB
## ------------smc组间细胞比例柱状图--------------
# 计算每组中每个细胞类型的个数
cell_counts <- table(Idents(smc), smc$group)

# 计算每组中各细胞类型的比例
cell_prop <- prop.table(cell_counts, margin = 2)

# 转换为数据框格式
cell_prop_df <- as.data.frame(cell_prop)
colnames(cell_prop_df) <- c("CellType", "Group", "Proportion")

# 调整组别的顺序，使其从上到下为ASYM, SYM
cell_prop_df$Group <- factor(cell_prop_df$Group, levels = c( "SYM", "ASYM"))

# 绘制横向堆叠条形图
library(ggplot2)

p <- ggplot(cell_prop_df, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Group", y = "Proportion", title = "") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.5, 'cm')) +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = use_col) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  scale_y_reverse()

# 查看绘图
print(p)
# 保存图像为PDF文件
ggsave("smcproportion.pdf", plot = p, width = 8, height = 2, device = "pdf")

#只显示Var2中，ASYM和SYM两组的比例，使用filter函数
cell.prop <- filter(cell.prop, Var2 %in% c("ASYM", "SYM"))


-------------------------------------------------------------------------------
# Figure C — SMC/Fibroblast subset Dot Plot
rm(list = ls())
gc()
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(patchwork)
  library(scales)
})

smc <- readRDS(paste0(getOption("cea.project_root"), "/SMC/smc.rds"))

out_dir <- paste0(getOption("cea.project_root"), "/SMC/Fig2_dotplot/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. 细胞顺序（y 轴自上而下：Fibro1 → Fibro2 → SMC1 → SMC2 → SMC3）──
smc$celltype <- factor(smc$celltype,
                       levels = c("SMCs3","SMCs2","SMCs1",
                                  "Fibroblasts2","Fibroblasts1"))

# ── 2. 精简到 9 个最判别的 marker（删 SFRP2/FBLN1/MYLK/MYL9）──────
markers_smc_fibro <- c(
  # Reparative / Fibroblast (5)
  "POSTN", "TNC", "COL1A1", "COL3A1", "LOXL2",
  # Mesenchymal canonical (2)
  "DCN", "LUM",
  # Contractile / SMC (4)
  "ACTA2", "MYH11", "CNN1", "TAGLN"
)
markers_smc_fibro <- intersect(markers_smc_fibro, rownames(smc))

# ── 3. Dot plot（关键参数全部调整）─────────────────────────────────
p_dot <- DotPlot(
  smc,
  features  = markers_smc_fibro,
  group.by  = "celltype",
  dot.scale = 10,
  col.min   = 0,      # ★ 不显示负 z-score（白色 = 0 起点）
  col.max   = 2.5     # ★ 高表达上限截断，避免极端值压缩中段
) +
  # ★ 多档蓝色梯度，让强弱差异更鲜明
  scale_color_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6",
                "#2171B5", "#08306B"),
    limits  = c(0, 2.5),
    oob     = scales::squish,
    name    = "Avg. expression"
  ) +
  # ★ 点尺寸的动态范围（最小 → 最大几乎 20×）
  scale_size_continuous(range = c(0.5, 9), name = "% expressed") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    # ★ 去除网格线
    panel.grid       = element_blank(),
    # 保留横纵坐标线与刻度
    axis.line        = element_line(linewidth = 0.4, color = "black"),
    axis.ticks       = element_line(linewidth = 0.4, color = "black"),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                    size = 11, face = "italic",
                                    color = "black"),
    axis.text.y      = element_text(size = 12, face = "bold",
                                    color = "black"),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(8, 12, 8, 12)
  )

# ── 4. 顶部 marker 分组色带 ────────────────────────────────────────
n_fibro <- 7    # POSTN/TNC/COL1A1/COL3A1/LOXL2/DCN/LUM
n_smc   <- 4    # ACTA2/MYH11/CNN1/TAGLN
total   <- n_fibro + n_smc

top_band <- data.frame(
  xmin  = c(0.5, n_fibro + 0.5),
  xmax  = c(n_fibro + 0.5, total + 0.5),
  group = c("Reparative / Fibroblast", "Contractile / SMC"),
  fill  = c("#C04A4A", "#3A6FA0")
)

top_strip <- ggplot(top_band,
                    aes(xmin = xmin, xmax = xmax,
                        ymin = 0,   ymax = 1, fill = fill)) +
  geom_rect(color = "white", linewidth = 0.6) +
  geom_text(aes(x = (xmin + xmax)/2, y = 0.5, label = group),
            color = "white", fontface = "bold", size = 4) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.5, total + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 12, 0, 12))

# ── 5. 拼接 + 输出 ─────────────────────────────────────────────────
p_dot_enhanced <- top_strip / p_dot +
  plot_layout(heights = c(0.5, 10))

ggsave(paste0(out_dir, "Fig2BC_dotplot_v2.pdf"),
       p_dot_enhanced, width = 8.5, height = 5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig2BC_dotplot_v2.png"),
       p_dot_enhanced, width = 8.5, height = 5, dpi = 600)

cat("\n??? Saved Fig2BC_dotplot_v2.pdf\n")



--------------------------------------------------------------------------------
# Fig 2D -

suppressMessages({
  library(monocle3); library(Seurat); library(tidyverse)
  library(igraph); library(ggrastr); library(viridis)
})

out_dir <- paste0(getOption("cea.project_root"), "/SMC/Fig2_v2/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(11)

## ── 0. 读取 5 亚群对象 ───────────────────────────────────────────
smc_all <- readRDS(paste0(getOption("cea.project_root"), "/SMC/smc.rds"))
smc_all <- subset(smc_all, celltype %in% c("SMCs1","SMCs2","SMCs3",
                                           "Fibroblasts1","Fibroblasts2"))
smc_all$celltype <- droplevels(factor(smc_all$celltype))
if (!"data" %in% Layers(smc_all[["RNA"]])) smc_all <- NormalizeData(smc_all)

## ── 1. 取主流四群学轨迹(排除 SMCs2) ─────────────────────────────
traj_types <- c("SMCs1","SMCs3","Fibroblasts1","Fibroblasts2")
smc <- subset(smc_all, celltype %in% traj_types)
smc$celltype <- droplevels(factor(smc$celltype))
if (inherits(smc[["RNA"]], "Assay5") &&
    length(grep("^counts", Layers(smc[["RNA"]]))) > 1) smc[["RNA"]] <- JoinLayers(smc[["RNA"]])
expr <- GetAssayData(smc, assay = "RNA", layer = "counts")
gene_md <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr))
cds <- new_cell_data_set(expr, cell_metadata = smc@meta.data, gene_metadata = gene_md)

cds <- preprocess_cds(cds, num_dim = 50)
tsne4 <- Embeddings(smc, reduction = "tsne")[colnames(cds), ]; colnames(tsne4) <- c("tSNE_1","tSNE_2")
reducedDims(cds)$UMAP <- tsne4
cds <- cluster_cells(cds, reduction_method = "UMAP")

## ★★★ 唯一实质改动：调参让分支延伸到 SMCs3 ★★★
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE,
                   learn_graph_control = list(minimal_branch_len = 6,   # 15 -> 6
                                              ncenter = 200))            # 120 -> 200

get_root_nodes <- function(cds, root_cells) {
  cv <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds), ])
  paste0("Y_", names(table(cv[root_cells, ]))[which.max(table(cv[root_cells, ]))])
}
root_cells <- colnames(smc)[smc$celltype %in% c("SMCs1","SMCs3")]
cds <- order_cells(cds, root_pr_nodes = get_root_nodes(cds, root_cells))
smc$pseudotime <- pseudotime(cds)[colnames(smc)]

cat("\n四群中位 pseudotime:\n")
print(smc@meta.data %>% group_by(celltype) %>%
        summarise(med_pt = median(pseudotime, na.rm = TRUE), n = n()) %>% arrange(med_pt))

## ── 2. 轨迹骨架 ──────────────────────────────────────────────────
ica <- as.data.frame(t(cds@principal_graph_aux[["UMAP"]]$dp_mst))
colnames(ica) <- c("prin_x","prin_y"); ica$node <- rownames(ica)
seg <- igraph::as_data_frame(principal_graph(cds)[["UMAP"]], what = "edges")
seg$x    <- ica$prin_x[match(seg$from, ica$node)]; seg$y    <- ica$prin_y[match(seg$from, ica$node)]
seg$xend <- ica$prin_x[match(seg$to,   ica$node)]; seg$yend <- ica$prin_y[match(seg$to,   ica$node)]

## ── 3. 绘图数据 ──────────────────────────────────────────────────
tsne_all <- Embeddings(smc_all, reduction = "tsne"); colnames(tsne_all) <- c("tSNE_1","tSNE_2")
df4 <- data.frame(tSNE_1 = tsne_all[colnames(smc),1], tSNE_2 = tsne_all[colnames(smc),2],
                  pseudotime = smc$pseudotime, celltype = as.character(smc$celltype))
df4 <- df4[df4$tSNE_2 > -15 & is.finite(df4$pseudotime), ]
smc2_cells <- colnames(smc_all)[smc_all$celltype == "SMCs2"]
df2 <- data.frame(tSNE_1 = tsne_all[smc2_cells,1], tSNE_2 = tsne_all[smc2_cells,2])
df2 <- df2[df2$tSNE_2 > -15, ]

lab4 <- df4 %>% group_by(celltype) %>%
  summarise(x = median(tSNE_1), y = median(tSNE_2),
            med_pt = median(pseudotime, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sprintf("%s\nPT = %.2f", celltype, med_pt))
lab2 <- data.frame(x = median(df2$tSNE_1), y = median(df2$tSNE_2),
                   label = "SMCs2\n(off-trajectory)")

## ── 4. 绘图（同 v6 风格）─────────────────────────────────────────
p_D <- ggplot() +
  ggrastr::rasterise(geom_point(data = df2, aes(tSNE_1, tSNE_2),
                                color = "grey75", size = 0.6, alpha = 0.9), dpi = 600) +
  ggrastr::rasterise(geom_point(data = df4[order(df4$pseudotime), ],
                                aes(tSNE_1, tSNE_2, color = pseudotime),
                                size = 0.6, alpha = 0.9), dpi = 600) +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend),
               linewidth = 0.8, color = "black", inherit.aes = FALSE) +
  geom_point(data = ica, aes(prin_x, prin_y), size = 1.2, color = "black",
             fill = "white", shape = 21, stroke = 0.5, inherit.aes = FALSE) +
  geom_label(data = lab4, aes(x, y, label = label), fill = "white", color = "black",
             alpha = 0.85, size = 3.2, fontface = "bold", label.padding = unit(0.2, "lines")) +
  geom_label(data = lab2, aes(x, y, label = label), fill = "grey90", color = "grey25",
             alpha = 0.9, size = 3.2, fontface = "bold", label.padding = unit(0.2, "lines")) +
  scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
  labs(title = "Pseudotime trajectory (root: contractile SMCs1/3; SMCs2 off-trajectory)",
       x = "tSNE_1", y = "tSNE_2") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold"), legend.position = "right")

ggsave(paste0(out_dir, "Fig2D_v8_trajectory.pdf"), p_D, width = 8, height = 6.5)
ggsave(paste0(out_dir, "Fig2D_v8_trajectory.png"), p_D, width = 8, height = 6.5, dpi = 600)
cat("\nFig2D v8 saved (trajectory tuned to reach SMCs3; SMCs2 off-trajectory)\n")


--------------------------------------------------------------------------------

##-------------Figure E — 14-gene pseudotime 动态热图----------
# 14 个 marker（contractile ↓ + reparative ↑ + KLF4 transition）
pt_markers <- c(
  # Contractile (top, declining)
  "ACTA2", "MYH11", "CNN1", "TAGLN", "MYLK", "MYL9",
  # Transition factor
  "KLF4",
  # Reparative (bottom, rising)
  "POSTN", "TNC", "COL1A1", "COL3A1", "COL1A2", "LOXL2", "FN1"
)
pt_markers <- intersect(pt_markers, rownames(smc))

# 计算每个细胞每个 gene 的表达 + pseudotime
expr_long <- FetchData(smc, vars = c(pt_markers, "pseudotime")) %>%
  filter(!is.infinite(pseudotime), !is.na(pseudotime)) %>%
  pivot_longer(-pseudotime, names_to = "gene", values_to = "expr")

# 沿 pseudotime 做 50 个 bin 平滑
expr_smooth <- expr_long %>%
  mutate(pt_bin = cut(pseudotime,
                      breaks = quantile(pseudotime,
                                        probs = seq(0, 1, length.out = 51),
                                        na.rm = TRUE),
                      include.lowest = TRUE, labels = FALSE)) %>%
  group_by(gene, pt_bin) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(scaled_expr = scale(mean_expr)[,1]) %>%
  ungroup()

# 固定基因顺序：收缩在上，KLF4 中间，修复在下
gene_order <- c("ACTA2","MYH11","CNN1","TAGLN","MYLK","MYL9",
                "KLF4",
                "POSTN","TNC","COL1A1","COL3A1","COL1A2","LOXL2","FN1")
expr_smooth$gene <- factor(expr_smooth$gene,
                           levels = rev(intersect(gene_order, unique(expr_smooth$gene))))

p_E <- ggplot(expr_smooth,
              aes(x = pt_bin, y = gene, fill = scaled_expr)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("#2171B5", "#F7F7F7", "#C04A4A"),
    limits  = c(-2, 2), oob = scales::squish,
    name    = "Scaled\nexpression") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(1, 25, 50),
                     labels = c("Start", "Mid", "End")) +
  labs(x = "Pseudotime", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_text(color = "black"),
        axis.text.y  = element_text(face = "italic", color = "black",
                                    size = 9),
        axis.line    = element_line(linewidth = 0.4),
        axis.ticks   = element_line(linewidth = 0.4))

ggsave(paste0(out_dir, "Fig2E_14gene_heatmap.pdf"), p_E,
       width = 6, height = 5, device = cairo_pdf)


--------------------------------------------------------------------------------
#--------------Figure F — Contractile vs Reparative 模块分数交叉曲线（headline）-----

# 用 AddModuleScore 算两个程序的得分
contractile_genes <- c("ACTA2","MYH11","CNN1","TAGLN","MYLK","MYL9")
reparative_genes  <- c("POSTN","TNC","COL1A1","COL3A1","COL1A2","LOXL2","FN1","DCN","LUM")

smc <- AddModuleScore(smc, features = list(contractile_genes),
                      name = "Contractile_score", ctrl = 5,seed = 42)
smc <- AddModuleScore(smc, features = list(reparative_genes),
                      name = "Reparative_score", ctrl = 5,seed = 42)

# 沿 pseudotime 做 50 bin 平滑
score_df <- FetchData(smc, vars = c("pseudotime",
                                    "Contractile_score1",
                                    "Reparative_score1")) %>%
  filter(!is.infinite(pseudotime), !is.na(pseudotime)) %>%
  mutate(pt_bin = cut(pseudotime,
                      breaks = quantile(pseudotime,
                                        probs = seq(0, 1, length.out = 51),
                                        na.rm = TRUE),
                      include.lowest = TRUE, labels = FALSE)) %>%
  group_by(pt_bin) %>%
  summarise(pt_mid       = median(pseudotime),
            Contractile  = mean(Contractile_score1),
            Reparative   = mean(Reparative_score1),
            .groups = "drop") %>%
  pivot_longer(c(Contractile, Reparative),
               names_to = "program", values_to = "score")


--------------------------------------------------------------------------------
  
## FigureG  —  Phenotypic Switching Along Pseudotime
## VSMC contractile program  vs  reparative fibroblast-like program

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

## ---- 通用设置 -------------------------------------------------------------
WORKDIR <- paste0(getOption("cea.project_root"), "/SMC/Fig2_v2")      # 按需修改
setwd(WORKDIR)

## Program 配色 (与原图一致: contractile=深棕, repair=赭红)
prog_cols <- c(Contractile = "#5C3A30", Repair = "#C16B54")

## Program 基因集 (须与最初 AddModuleScore 用的完全一致)
contractile_genes <- c("ACTA2", "MYH11", "CNN1", "TAGLN", "MYL9")
repair_genes      <- c("POSTN", "TNC", "COL1A1", "COL3A1", "LOXL2", "SPARC")

## PART 1  —  从汇总 CSV 复现 Fig 2H (保证与之前输出完全一致)
## 直接读已固定的分箱均值 (20 个 pseudotime 等宽箱), 连点成线 —— 这一步不涉及
## AddModuleScore 的随机性, 因此结果与原图逐点重合.

summ <- read.csv(paste0(getOption("cea.project_root"), "/SMC/Pseudotime_ModuleScore/Pseudotime_ModuleScore_Summary.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
## 列: pseudotime_group, pseudotime, Contractile, Repair

## 转长表
summ_long <- summ %>%
  select(pseudotime, Contractile, Repair) %>%
  pivot_longer(c(Contractile, Repair),
               names_to = "Program", values_to = "ModuleScore") %>%
  mutate(Program = factor(Program, levels = c("Contractile", "Repair")))

## 交叉点 (见 PART 2 函数; 这里先调用得到 x_cross 以便画参考线)
find_crossover <- function(x, y1, y2) {
  d <- y1 - y2
  idx <- which(d[-1] * d[-length(d)] < 0)          # 符号翻转处
  out <- lapply(idx, function(i) {
    frac <- d[i] / (d[i] - d[i + 1])
    data.frame(pseudotime = x[i] + frac * (x[i + 1] - x[i]),
               score      = y1[i] + frac * (y1[i + 1] - y1[i]))
  })
  do.call(rbind, out)
}
cross <- find_crossover(summ$pseudotime, summ$Contractile, summ$Repair)
x_cross <- cross$pseudotime[1]; y_cross <- cross$score[1]
message(sprintf("Crossover: pseudotime = %.3f, score = %.3f", x_cross, y_cross))

p <- ggplot(summ_long, aes(pseudotime, ModuleScore, color = Program)) +
  geom_vline(xintercept = x_cross, linetype = "dotted",
             color = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 1.6, lineend = "round") +
  scale_color_manual(values = prog_cols) +
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10)) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1.0, 1.5)) +
  labs(title = "Phenotypic Switching Along Pseudotime",
       x = "Pseudotime", y = "Module Score", color = "Program") +
  theme_bw(base_size = 18) +
  theme(
    plot.title       = element_text(face = "bold", size = 22, hjust = 0),
    axis.title       = element_text(size = 20),
    axis.text        = element_text(size = 16, color = "black"),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.border     = element_rect(color = "black", linewidth = 0.8)
  )

## 线条本身是矢量(无需栅格化); 直接存 PDF/SVG 进 AI 即可
ggsave("Fig2H_phenotypic_switching.pdf", p, width = 8, height = 5.5)
ggsave("Fig2H_phenotypic_switching.png", p, width = 8, height = 5.5, dpi = 600)
print(p)


--------------------------------------------------------------------------------
#------------------FigH_pseudotime_density---------
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
PT_COL <- "pseudotime"
grp_cols <- c(ASYM = "#F9C80E", SYM = "#C71000")

df <- smc@meta.data %>%
  filter(group %in% c("ASYM","SYM"), is.finite(.data[[PT_COL]])) %>%
  transmute(pt = .data[[PT_COL]], group = factor(group, levels=c("ASYM","SYM")))

## 在统一网格上估计两组密度
rng <- range(df$pt); grid <- seq(rng[1], rng[2], length.out = 512)
dens <- df %>% group_by(group) %>%
  summarise(d = list(approx(density(pt, from=rng[1], to=rng[2], n=512), xout=grid)$y),
            .groups="drop") %>%
  mutate(grid = list(grid)) %>% unnest(c(d, grid))
diffdf <- dens %>% pivot_wider(names_from=group, values_from=d) %>%
  mutate(delta = SYM - ASYM)   # >0: SYM 更密

p_diff <- ggplot(diffdf, aes(grid, delta)) +
  geom_area(aes(fill = delta > 0), alpha = 0.7) +
  geom_hline(yintercept = 0, color = "grey30") +
  scale_fill_manual(values = c(`TRUE`="#C71000", `FALSE`="#F9C80E"),
                    labels = c("Enriched in ASYM (contractile end)",
                               "Enriched in SYM (reparative end)"),
                    name = NULL) +
  labs(x = "Pseudotime", y = "Density difference (SYM ??? ASYM)",
       title = "Ruptured plaques shift toward the reparative trajectory end") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=12), legend.position="bottom")
ggsave("Fig2_pseudotime_density_diff.pdf", p_diff, width = 7, height = 4.5)
print(p_diff)

#--------------"Fig2_contractile_fraction_ASYM_SYM-------------------

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

PT_COL <- "pseudotime"; GRP_COL <- "group"; SAMP_COL <- "sample"  # 按对象实际列名
grp_cols <- c(ASYM = "#F9C80E", SYM = "#C71000")

## ---- 1.按"样本"计算各亚群占比 ----
md <- smc@meta.data %>%
  filter(.data[[GRP_COL]] %in% c("ASYM","SYM")) %>%
  rename(group = !!GRP_COL, sample = !!SAMP_COL, celltype = celltype)

## 每个样本内各亚群比例
samp_comp <- md %>%
  count(group, sample, celltype) %>%
  group_by(sample) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

## 收缩型 compartment 占比 (SMCs1+SMCs3) —— 每个样本一个值
samp_contr <- md %>%
  mutate(comp = ifelse(celltype %in% c("SMCs1","SMCs3"),
                       "Contractile", "Reparative")) %>%
  count(group, sample, comp) %>%
  group_by(sample) %>% mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
  filter(comp == "Contractile")

print(samp_contr)   # 看每个样本的收缩型占比

## ---- 2. 统计: 样本水平 Wilcoxon (n=6 vs 6) ----
w <- wilcox.test(pct ~ group, data = samp_contr)
cat(sprintf("\n收缩型占比 ASYM vs SYM: Wilcoxon p = %.4f\n", w$p.value))
samp_contr %>% group_by(group) %>%
  summarise(median = median(pct), mean = mean(pct), n = n()) %>% print()


p_box <- ggplot(samp_contr, aes(group, pct, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 3, shape = 21, color = "black",
              aes(fill = group), show.legend = FALSE) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  annotate("text", x = 1.5, y = max(samp_contr$pct) * 1.08,
           label = sprintf("Wilcoxon p = %.3f", w$p.value), size = 4.5) +
  labs(x = NULL, y = "Contractile SMC fraction (%)",
       title = "Contractile SMCs are depleted in ruptured (SYM) plaques") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.major.x = element_blank())
ggsave("Fig2_contractile_fraction_ASYM_SYM.pdf", p_box, width = 4.5, height = 5)
print(p_box)


--------------------------------------------------------------------------------

# Supplementary Figure3 
# Panel A: MYOCD vs KLF4 同图（z-score 标准化，同一尺度看相反趋势）
#          MYOCD = 收缩主控（↓）, KLF4 = 去分化主控开关（↑）
# Panel B: ASYM vs SYM pseudotime 密度

rm(list = ls()); gc()
suppressMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggplot2); library(mgcv)
})
out_dir <- paste0(getOption("cea.project_root"), "/SMC/Fig2_v2/")
ver     <- "v6"

smc <- readRDS(paste0(getOption("cea.project_root"), "/SMC/Fig2_v2/smc_with_SMCs2_pseudotime.rds"))
smc <- subset(smc, celltype %in% c("SMCs1","SMCs3","Fibroblasts1","Fibroblasts2"))
smc <- subset(smc, cells = colnames(smc)[is.finite(smc$pseudotime)])

## ── Panel A：MYOCD + KLF4，各自 z-score 到同一尺度 ───────────────
genes <- c("MYOCD","KLF4"); genes <- genes[genes %in% rownames(smc)]
stopifnot(length(genes) == 2)
e <- FetchData(smc, vars = genes, layer = "data")
e <- as.data.frame(scale(e))                       # 每个基因 z-score（均值0 SD1）
e$pseudotime <- smc$pseudotime
longA <- pivot_longer(e, all_of(genes), names_to = "gene", values_to = "z")
longA$gene <- factor(longA$gene, levels = c("MYOCD","KLF4"))

tf_col <- c(MYOCD = "#3A6FA0", KLF4 = "#C0392B")   # 收缩蓝 / 去分化红

pA <- ggplot(longA, aes(pseudotime, z, color = gene)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6), se = TRUE,
              linewidth = 1.3, alpha = 0.12) +
  scale_color_manual(values = tf_col, name = NULL,
                     labels = c("MYOCD (contractile master regulator)",
                                "KLF4 (de-differentiation switch)")) +
  labs(x = "Pseudotime", y = "Scaled expression (z-score)",
       title = "Opposing dynamics of master regulators along the VSMC-to-fibroblast trajectory") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black"),
        legend.position = "top", legend.text = element_text(face = "bold"))+
  coord_cartesian(xlim = c(0, 18))      # 裁掉末端稀疏区 + 缩短横坐标

## ── Panel B：ASYM vs SYM 密度 ───────────────────────────────────
dfB <- data.frame(pseudotime = smc$pseudotime,
                  grp = as.character(smc[["group"]][,1])) %>% filter(grp %in% c("ASYM","SYM"))
dfB$grp <- factor(dfB$grp, levels = c("ASYM","SYM"))
med_line <- dfB %>% group_by(grp) %>% summarise(med = median(pseudotime), .groups="drop")

grp_cols <- c(ASYM = "#FFE705FF", SYM = "#C71000FF")      # 黄/红（如需切换取消注释）

pB <- ggplot(dfB, aes(pseudotime, fill = grp, color = grp)) +
  geom_density(alpha = 0.40, linewidth = 0.9) +
  geom_vline(data = med_line, aes(xintercept = med, color = grp),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_color_manual(values = grp_cols, name = NULL) +
  labs(x = "Pseudotime", y = "Cell density",
       title = "Symptomatic plaques shift toward the reparative end (dashed = group median)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black"),
        legend.position = "top", legend.text = element_text(face = "bold"))+
  coord_cartesian(xlim = c(0, 18))

p <- pA / pB + plot_layout(heights = c(1, 1))
ggsave(paste0(out_dir, "SuppFig_Fig2_MYOCD_KLF4_", ver, ".pdf"), p, width = 8, height = 7)
ggsave(paste0(out_dir, "SuppFig_Fig2_MYOCD_KLF4_", ver, ".png"), p, width = 8, height = 7, dpi = 600)
cat("\nsaved v6\n"); print(med_line)
################################################################################




################################################################################
# Supplementary Figure (Fig2 支撑) v5
# Panel A: 5 个【明星】VSMC 表型转换 TF，左右两张独立子图（各自 y 轴）
#   收缩维持↓: MYOCD, GATA6   |   去分化驱动↑: KLF4, TCF21, SOX9
# Panel B: ASYM vs SYM pseudotime 密度
################################################################################

rm(list = ls()); gc()
suppressMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggplot2); library(mgcv)
})
out_dir <- paste0(getOption("cea.project_root"), "/SMC/Fig2_v2/")
ver     <- "v5"

smc <- readRDS(paste0(getOption("cea.project_root"), "/SMC/Fig2_v2/smc_with_SMCs2_pseudotime.rds"))
smc <- subset(smc, celltype %in% c("SMCs1","SMCs3","Fibroblasts1","Fibroblasts2"))
smc <- subset(smc, cells = colnames(smc)[is.finite(smc$pseudotime)])

## ── 明星 TF（写死，审稿人一眼认得）──────────────────────────────
tf_contractile <- c("MYOCD","GATA6")            # 收缩维持，应↓
tf_modulation  <- c("KLF4","TCF21","SOX9")      # 去分化驱动，应↑
tf_contractile <- tf_contractile[tf_contractile %in% rownames(smc)]
tf_modulation  <- tf_modulation [tf_modulation  %in% rownames(smc)]

get_long <- function(genes) {
  e <- FetchData(smc, vars = genes, layer = "data"); e$pseudotime <- smc$pseudotime
  pivot_longer(e, all_of(genes), names_to = "gene", values_to = "expr")
}
longC <- get_long(tf_contractile); longC$gene <- factor(longC$gene, levels = tf_contractile)
longM <- get_long(tf_modulation);  longM$gene <- factor(longM$gene, levels = tf_modulation)

col_C <- setNames(c("#1F3F66","#3A6FA0")[seq_along(tf_contractile)], tf_contractile)
col_M <- setNames(c("#C0392B","#E8772E","#F2B705")[seq_along(tf_modulation)], tf_modulation)

base_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black"),
        legend.position = "top", legend.text = element_text(face = "bold", size = 9),
        legend.title = element_blank())

## 左：收缩维持（独立 y 轴）
pA1 <- ggplot(longC, aes(pseudotime, expr, color = gene)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6), se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = col_C) +
  labs(x = "Pseudotime", y = "Expression (smoothed)",
       title = "Contractile-maintenance TFs (decline)") + base_theme

## 右：去分化驱动（独立 y 轴）
pA2 <- ggplot(longM, aes(pseudotime, expr, color = gene)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6), se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = col_M) +
  labs(x = "Pseudotime", y = "Expression (smoothed)",
       title = "De-differentiation driver TFs (rise)") + base_theme

pA <- pA1 + pA2 + plot_layout(widths = c(1, 1))

## ── Panel B：ASYM vs SYM 密度 ───────────────────────────────────
dfB <- data.frame(pseudotime = smc$pseudotime,
                  grp = as.character(smc[["group"]][,1])) %>% filter(grp %in% c("ASYM","SYM"))
dfB$grp <- factor(dfB$grp, levels = c("ASYM","SYM"))
med_line <- dfB %>% group_by(grp) %>% summarise(med = median(pseudotime), .groups="drop")

grp_cols <- c(ASYM = "#FFE705FF", SYM = "#C71000FF")       # 论文标准红/蓝
# grp_cols <- c(ASYM = "#FFE705FF", SYM = "#C71000FF")      # 黄/红（如需切换）

pB <- ggplot(dfB, aes(pseudotime, fill = grp, color = grp)) +
  geom_density(alpha = 0.40, linewidth = 0.9) +
  geom_vline(data = med_line, aes(xintercept = med, color = grp),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_color_manual(values = grp_cols, name = NULL) +
  labs(x = "Pseudotime", y = "Cell density",
       title = "Symptomatic plaques shift toward the reparative end (dashed = group median)") +
  base_theme

p <- pA / pB + plot_layout(heights = c(1, 1))
ggsave(paste0(out_dir, "SuppFig_Fig2_TF_progression_", ver, ".pdf"), p, width = 9, height = 7.5)
ggsave(paste0(out_dir, "SuppFig_Fig2_TF_progression_", ver, ".png"), p, width = 9, height = 7.5, dpi = 600)
cat("\nsaved v5\n收缩维持 TF:", paste(tf_contractile, collapse=", "),
    "\n去分化驱动 TF:", paste(tf_modulation, collapse=", "), "\n"); print(med_line)
################################################################################


# Supplementary Figure Fig2 

rm(list = ls()); gc()
suppressMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggplot2); library(mgcv)
})
out_dir <- paste0(getOption("cea.project_root"), "/SMC/Fig2_v2/")
ver     <- "v4"

smc <- readRDS(paste0(getOption("cea.project_root"), "/SMC/Fig2_v2/smc_with_SMCs2_pseudotime.rds"))
smc <- subset(smc, celltype %in% c("SMCs1","SMCs3","Fibroblasts1","Fibroblasts2"))
smc <- subset(smc, cells = colnames(smc)[is.finite(smc$pseudotime)])

## ── VSMC 表型转换【转录因子】候选池（均为 TF，不含结构/marker 基因）──
# 收缩维持 TF（应↓）
tf_contractile <- c("MYOCD","SRF","TEAD1","GATA6","NFATC1","NFATC2","PRDM6")
# 去分化/modulation 驱动 TF（应↑）
tf_modulation  <- c("KLF4","KLF5","TCF21","SOX9","RUNX2","ELK1","ETS1",
                    "NR4A1","NR4A2","NR4A3","FOXO1","FOXO3","SMAD3","YAP1","FOSL1","JUN")

pick_present <- function(cands, obj, min_pct = 0.03) {
  cands <- cands[cands %in% rownames(obj)]
  if (length(cands) == 0) return(character(0))
  pct <- rowMeans(GetAssayData(obj, layer = "data")[cands, , drop = FALSE] > 0)
  names(pct)[pct > min_pct]                       # 全部满足表达比例的都保留
}

contractile_reg <- pick_present(tf_contractile, smc)
modulation_reg  <- pick_present(tf_modulation,  smc)
regs <- c(contractile_reg, modulation_reg)
cat("收缩维持 TF 选用:", paste(contractile_reg, collapse=", "), "\n")
cat("去分化驱动 TF 选用:", paste(modulation_reg,  collapse=", "), "\n")
stopifnot(length(regs) >= 2)

exprA <- FetchData(smc, vars = regs, layer = "data"); exprA$pseudotime <- smc$pseudotime
longA <- exprA %>% pivot_longer(all_of(regs), names_to = "gene", values_to = "expr") %>%
  mutate(group = ifelse(gene %in% contractile_reg,
                        "Contractile-maintenance TFs", "De-differentiation TFs"))
longA$gene  <- factor(longA$gene, levels = regs)
longA$group <- factor(longA$group, levels = c("Contractile-maintenance TFs","De-differentiation TFs"))

## 配色：收缩维持=蓝系渐变，去分化=红橙系渐变（按实际数量自动取色）
n_c <- length(contractile_reg); n_m <- length(modulation_reg)
blues <- colorRampPalette(c("#1F3F66","#3A6FA0","#7FA5C7"))(max(n_c,1))
reds  <- colorRampPalette(c("#7A1B12","#C0392B","#E8772E","#F2B705"))(max(n_m,1))
reg_col <- setNames(c(blues[seq_len(n_c)], reds[seq_len(n_m)]), regs)

pA <- ggplot(longA, aes(pseudotime, expr, color = gene)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6), se = FALSE, linewidth = 1.0) +
  scale_color_manual(values = reg_col, name = NULL) +
  facet_wrap(~ group, nrow = 1) +                 # 两组分面，方向对比更清晰
  labs(x = "Pseudotime", y = "Expression (smoothed)",
       title = "VSMC phenotype-switching transcription factors along pseudotime") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold", size = 10),
        axis.text = element_text(color = "black"),
        legend.position = "right", legend.text = element_text(size = 8))

## ── Panel B：ASYM vs SYM 密度 ───────────────────────────────────
dfB <- data.frame(pseudotime = smc$pseudotime,
                  grp = as.character(smc[["group"]][,1])) %>%
  filter(grp %in% c("ASYM","SYM"))
dfB$grp <- factor(dfB$grp, levels = c("ASYM","SYM"))
med_line <- dfB %>% group_by(grp) %>% summarise(med = median(pseudotime), .groups="drop")

## ★ 配色二选一（默认论文标准 红/蓝）：
grp_cols <- c(ASYM = "#FFE705FF", SYM = "#C71000FF")        # 黄/红

pB <- ggplot(dfB, aes(pseudotime, fill = grp, color = grp)) +
  geom_density(alpha = 0.40, linewidth = 0.9) +
  geom_vline(data = med_line, aes(xintercept = med, color = grp),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_color_manual(values = grp_cols, name = NULL) +
  labs(x = "Pseudotime", y = "Cell density",
       title = "Symptomatic plaques shift toward the reparative end of the trajectory (dashed = group median)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(color = "black"),
        legend.position = "top", legend.text = element_text(face = "bold"))

p <- pA / pB + plot_layout(heights = c(1, 1))
ggsave(paste0(out_dir, "SuppFig_Fig2_TF_progression_", ver, ".pdf"), p, width = 9, height = 7.5)
ggsave(paste0(out_dir, "SuppFig_Fig2_TF_progression_", ver, ".png"), p, width = 9, height = 7.5, dpi = 600)
cat("\nsaved v4 ; 组中位 pseudotime:\n"); print(med_line)

