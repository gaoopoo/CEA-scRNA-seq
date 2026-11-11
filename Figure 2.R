#Figure 2. Pathology-stratified single-cell atlas identifies repair-associated programs in human carotid plaques
## Figure A-C
#-------------基因表达可视化----------------
  #加载包
  rm(list = ls())
gc()
setwd("G:/13sample/")
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidydr)
combined131<- readRDS(file="combined131.rds")

use_col<- c(Fibroblasts1="#B46356FF",
            "T cells1"="#944D40FF", 
            SMCs1="#7A4B42FF", 
            Macrophages1= "#CA480DFF", 
            "T cells2"= "#9B5928FF", 
            Macrophages2= "#C86C18FF",
            "T cells3"= "#D78017FF",
            Monocytes= "#F9AB3CFF", 
            ECs= "#FFE705FF",
            Neutrophils="#FFF7A3FF",
            Fibroblasts2="#E2ED50FF",
            "NK cells"="#2F7270FF", 
            SMCs2="#66FFFFFF", 
            SMCs3= "#3C999CFF",  
            "Mast cells"= "#8BEAFFFF",
            DCs= "#13A4EBFF",
            "Proliferating cells1"= "#2E7BC5FF", 
            "B cells1"= "#A1B3E2FF", 
            "B cells2"= "#BFCAFEFF",
            Macrophages3="#8B799CFF",
            "Proliferating cells2"="#524559FF",
            "B cell"="#FF6F00FF",#maincelltype_col
            DC="#C71000FF",
            EC="#008EA0FF",
            Fibro="#8A4198FF",
            Macro="#5A9599FF",
            Mast="#FF6348FF",
            Mono="#84D7E1FF",
            Neutro="#FF95A8FF",
            NK="#3D3B25FF",
            Prolifing="#ADE2D0FF",
            SMC="#1A5354FF",
            "T cell"="#3F4041FF",
            ASYM="#C71000FF",#group_col
            SYM="#8BEAFFFF")

#主要细胞亚群dotplot
combined131$maincelltype <- fct_relevel(combined131$maincelltype, "Prolifing","Bcell","DC","Mast",
                                        "NK","Neutro","EC","Mono",
                                        "Macro","Tcell", "SMC","Fibro")

markers1 <- c( "LUM", "DCN", "SFRP2", "FBLN1",#Fibro
               "ACTA2", "MYH11", "TAGLN", "MYLK", #SMC
               "TRAC", "GZMA", "IL7R","CD8A", #Tcell
               "C1QC", "C1QB","C1QA","SELENOP", #Macro
               "EREG","IL1A","SERPINB2","KYNU",#Mono
               "MMRN1","CD34","VWF","TIE1",#EC
               "MMP9","CD36","LPL","MMP12",#Neutro
               "XCL1","KLRD1","NKG7","TRDC",#NK
               "TPSB2","ADCYAP1","CPA3","HDC",#Mast
               "CCR7","FCER1A","NRARP","GLA",#DC cell
               "JCHAIN","IGKC","MS4A1","LY9",#B cell
               "TOP2A","MKI67","TK1","TYMS")#profilng 
DotPlot(combined131, features = markers1, cols = c("#FFFFFF", "#00468BFF"),
        dot.scale = 8, group.by = "maincelltype") +
  theme(axis.text.x = element_text(angle = 90))           

ggsave("combined131/dotplot.pdf",width=14,height=6,device = "pdf") 


###celltype;TSNE坐标轴调整
library(SingleR)
library(Seurat)
library(tidyverse)
library(ggrepel)
#
DimPlot(combined131, group.by = "celltype",reduction = "tsne",raster=FALSE)
#二 ggplot2绘制tsne图
#2.1 查看绘制tsne的数据
tsne = combined131@reductions$tsne @cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(celltype = combined131@meta.data$celltype) #celltype 或者者group

p <- ggplot(tsne,aes(x= tSNE_1 , y = tSNE_2 ,color = celltype)) +
  geom_point(size = 0.1 , alpha =1 )+scale_color_manual(values = use_col)

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

ggsave("combined131/celltypetsne1.pdf",width=8,height=6,device = "pdf")


##------------maincelltype;TSNE坐标轴调整
tsne = combined131@reductions$tsne @cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(maincelltype = combined131@meta.data$maincelltype) #celltype 或者者group
head(tsne )
p <- ggplot(tsne,aes(x= tSNE_1 , y = tSNE_2 ,color = maincelltype)) +
  geom_point(size = 0.1, alpha =1)+scale_color_manual(values = use_col)

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

ggsave("combined131/maincelltypetsne1.pdf",width=8,height=6,device = "pdf")


##------------group;TSNE坐标轴调整
tsne = combined131@reductions$tsne @cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(group = combined131@meta.data$group ) #celltype 或者者group
head(tsne )
p <- ggplot(tsne,aes(x= tSNE_1 , y = tSNE_2 ,color = group )) +
  geom_point(size = 0.1, alpha =1)+scale_color_manual(values = use_col)

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

ggsave("combined131/grouptsne.pdf",width=8,height=6,device = "pdf")


#细胞亚群合并maincell
maincell<- c("0"="VSMC","1"="LYM", "2"="VSMC", "3"= "Macro", 
             "4"= "LYM", "5"= "Macro","6"= "LYM", "7"= "Macro", 
             "8"= "EC","9"="MYE","10"="VSMC","11"="LYM", 
             "12"="VSMC", "13"= "VSMC", "14"= "MYE","15"= "MYE",
             "16"= "Prolifing", "17"= "LYM", "18"= "LYM",
             "19"="Macro","20"="Prolifing")
combined131[['maincell']] = unname(maincell[combined131@meta.data$seurat_clusters])

table(combined131$maincell)#查看maincell细胞数,输入下面

#----------------maincell 比例
mydata<-data.frame(occupation=c(rep("Macro",27144),rep("VSMC",43491),
                                rep("EC",5378),rep("LYM",33696),
                                rep("MYE",8595),rep("Prolifing",1163)))
#绘图，但是按照字母排列

pdf("combined131/maincellproportion饼图.pdf", width = 8, height = 6)
PieChart(occupation, data = mydata,
         hole =0,fill="viridis",#内部圈为0.5颜色
         values_color="red",values_size=1.5,#数字标签的颜色，字号
         main="Cell Proportion",#标题
         main_cex =1.3) 
dev.off()

#ASYM组，maincell 比例
mydata1<-data.frame(occupation1=c(rep("Macro",12561),rep("VSMC",23241),
                                  rep("EC",2687),rep("LYM",20250),
                                  rep("MYE",3485),rep("Prolifing",396)))

pdf("combined131/ASYM比例饼图.pdf", width = 8, height = 6)
PieChart(occupation1, data = mydata1,
         hole =0.6,fill="viridis",#内部圈为0.5颜色
         values_color="red",values_size=1.5,#数字标签的颜色，字号
         main="ASYM Cell Proportion",#标题
         main_cex =1.3) 
dev.off()

#SYM组，maincelltype 比例
mydata2<-data.frame(occupation2=c(rep("Macro",13164),rep("VSMC",15174),
                                  rep("EC",1700),rep("LYM",10776),
                                  rep("MYE",4558),rep("Prolifing",686)))

pdf("combined131/SYM比例饼图.pdf", width = 8, height = 6)
PieChart(occupation2, data = mydata2,
         hole =0.6,fill="viridis",#内部圈为0.5颜色
         values_color="red",values_size=1.5,#数字标签的颜色，字号
         main="SYM Cell Proportion",#标题
         main_cex =1.3) 
dev.off()


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
# FigureD 

suppressMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

setwd("G:/13sample/")
out_dir <- "G:/13sample/combined131/"

# ── 1. 数据准备 ───────────────────────────────────────────────────
prop_data <- data.frame(
  group = c( rep("ASYM", 6), rep("SYM", 6)),
  celltype = rep(c("VSMC","LYM","Macro","MYE","EC","Proliferating"), 3),
  count = c(
    23241, 20250, 12561, 3485, 2687, 396,     # ASYM
    15174, 10776, 13164, 4558, 1700, 686      # SYM
  )
)

prop_data <- prop_data %>%
  group_by(group) %>%
  mutate(pct = count / sum(count) * 100,
         pct_label = sprintf("%.0f%%", pct)) %>%
  ungroup()

# ── 2. 因子顺序（决定堆叠顺序）────────────────────────────────────
# ggplot 堆叠顺序：factor levels 第一个在底层
celltype_stack_order <- c("VSMC", "LYM", "Macro", "MYE", "EC", "Proliferating")
prop_data$celltype <- factor(prop_data$celltype, levels = celltype_stack_order)
prop_data$group <- factor(prop_data$group, levels = c( "ASYM", "SYM"))

# ── 3. 配色 ───────────────────────────────────────────────────────
maincell_col <- c(
  "VSMC"          = "#FFE705FF",
  "LYM"           = "#3A6FA0FF",
  "Macro"         = "#CA480DFF",
  "MYE"           = "#F9AB3CFF",
  "EC"            = "#8A4198FF",
  "Proliferating" = "#524559FF"
)

# ── 4. 分两类：内部标签 vs 外部标签 ─────────────────────────────
# 阈值：≥3% 在柱内标，<3% 用引线引到外侧
labels_inside  <- prop_data %>% filter(pct >= 3)
labels_outside <- prop_data %>% filter(pct < 3)

cat("Labels inside (>=3%):\n")
print(labels_inside %>% dplyr::select(group, celltype, pct))
cat("\nLabels outside (<3%):\n")
print(labels_outside %>% dplyr::select(group, celltype, pct))

# ── 5. 给外部标签计算 y 坐标 ─────────────────────────────────────
# 关键：必须按 ggplot 的堆叠规则算
# ggplot 用 position_stack 时，是按 factor levels 顺序累加
# factor levels: VSMC, LYM, Macro, MYE, EC, Proliferating
# 所以 VSMC 在最底层 (0 - pct_VSMC)
# LYM 在 VSMC 上方 (pct_VSMC - pct_VSMC+pct_LYM)
# ...

prop_data <- prop_data %>%
  arrange(group, celltype) %>%   # 关键：按 factor 顺序排序
  group_by(group) %>%
  mutate(
    y_cum_end = cumsum(pct),               # 该色块顶部
    y_cum_start = y_cum_end - pct,         # 该色块底部
    y_center = (y_cum_start + y_cum_end) / 2  # 该色块中心
  ) %>%
  ungroup()

# 重新分类（带 y_center 信息）
labels_inside  <- prop_data %>% filter(pct >= 3)
labels_outside <- prop_data %>% filter(pct < 3)

cat("\nOutside labels with y_center:\n")
print(labels_outside %>% dplyr::select(group, celltype, pct, y_center))

# ── 6. 绘图 ───────────────────────────────────────────────────────
p_D <- ggplot(prop_data,
              aes(x = group, y = pct, fill = celltype)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.4) +
  
  # ★ 内部标签：用手动算的 y_center，统一黑色
  geom_text(data = labels_inside,
            aes(x = group, y = y_center, label = pct_label),
            color = "black", fontface = "bold", size = 3.5,
            inherit.aes = FALSE) +
  
  # ★ 外部标签：只显示比例（不带细胞名），用 ggrepel 引到右侧
  geom_text_repel(
    data = labels_outside,
    aes(x = group, y = y_center, label = pct_label),
    nudge_x = 0.4,                  # 推到柱右侧
    direction = "y",                # 只在 y 方向调整避免重叠
    hjust = 0,                      # 左对齐
    segment.size = 0.4,
    segment.color = "grey40",
    fontface = "bold",
    size = 3.2,
    color = "black",
    box.padding = 0.3,
    point.padding = 0.1,
    min.segment.length = 0,         # 即使近也画引线
    inherit.aes = FALSE,
    seed = 42
  ) +
  
  scale_fill_manual(values = maincell_col, name = NULL,
                    breaks = celltype_stack_order) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     limits = c(0, 101)) +
  coord_cartesian(clip = "off") +     # ★ 允许标签超出绘图区
  
  labs(title = "D. Cell composition",
       x = NULL,
       y = "Proportion of cells") +
  
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(face = "bold", size = 11, color = "black"),
    axis.text.y      = element_text(size = 9, color = "black"),
    axis.title.y     = element_text(size = 9.5, color = "grey30"),
    axis.line.y      = element_line(linewidth = 0.4),
    axis.line.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(linewidth = 0.4),
    legend.position  = "right",
    legend.text      = element_text(size = 9, face = "bold"),
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(15, 30, 10, 10)   # ★ 右边留出空间给外部标签
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(paste0(out_dir, "Fig1D_v7_composition.pdf"), p_D,
       width = 6.5, height = 5, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig1D_v7_composition.png"), p_D,
       width = 6.5, height = 5, dpi = 600)

--------------------------------------------------------------------------------
#Fig E-F  (GSE159677)

# 0. 环境准备
# ============================================================================
rm(list = ls())
gc()  # 清内存

suppressMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(GEOquery)
  library(Matrix)
})

# 创建工作目录
out_dir <- "G:/13sample/Alsaigh_integration/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 用户指定颜色
celltype_colors <- c(
  "Fibroblasts1" = "#B46356FF",
  "T cells1"     = "#944D40FF",
  "SMCs1"        = "#7A4B42FF",
  "Macrophages1" = "#CA480DFF",
  "T cells2"     = "#9B5928FF",
  "Macrophages2" = "#C86C18FF",
  "T cells3"     = "#D78017FF",
  "Monocytes"    = "#F9AB3CFF",
  "ECs"          = "#FFE705FF",
  "Neutrophils"  = "#FFF7A3FF",
  "Fibroblasts2" = "#E2ED50FF",
  "NK cells"     = "#2F7270FF",
  "SMCs2"        = "#66FFFFFF",
  "SMCs3"        = "#3C999CFF",
  "Mast cells"   = "#8BEAFFFF",
  "DCs"          = "#13A4EBFF",
  "Proliferating cells1" = "#2E7BC5FF",
  "B cells1"     = "#A1B3E2FF",
  "B cells2"     = "#BFCAFEFF",
  "Macrophages3" = "#8B799CFF",
  "Proliferating cells2" = "#524559FF"
)

# ============================================================================
# 1. 下载 Alsaigh 2022 数据 (GSE159677)
# ============================================================================
cat("\n========== Step 1: Downloading Alsaigh 2022 data ==========\n")

# GEO 数据下载（首次跑会下载，~500 MB-1 GB）
download_dir <- paste0(out_dir, "GSE159677_raw/")
dir.create(download_dir, showWarnings = FALSE)

# 方法 1：用 GEOquery 自动下载
# (如果失败，跳到方法 2 手动下载)
tryCatch({
  options(timeout = 600)  # 延长超时时间
  
  if (!file.exists(paste0(download_dir, "GSE159677_RAW.tar"))) {
    cat("Downloading GSE159677_RAW.tar (~500MB)...\n")
    cat("If download is slow, manually download from:\n")
    cat("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159677&format=file\n")
    cat("Save to:", download_dir, "\n\n")
    
    download.file(
      url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159677&format=file",
      destfile = paste0(download_dir, "GSE159677_RAW.tar"),
      mode = "wb")
  }
  
  # 解压
  if (!dir.exists(paste0(download_dir, "GSE159677_extracted"))) {
    untar(paste0(download_dir, "GSE159677_RAW.tar"),
          exdir = paste0(download_dir, "GSE159677_extracted"))
  }
  
}, error = function(e) {
  cat("\n!!! 自动下载失败 !!!\n")
  cat("请手动下载：\n")
  cat("1. 访问 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677\n")
  cat("2. 点击 'http' 下方的 GSE159677_RAW.tar\n")
  cat("3. 下载到：", download_dir, "\n")
  cat("4. 解压到：", paste0(download_dir, "GSE159677_extracted"), "\n")
  cat("5. 重新运行此脚本\n")
  stop("Manual download required")
})

# 列出解压后的文件
extracted_files <- list.files(paste0(download_dir, "GSE159677_extracted"), 
                              pattern = "\\.h5$|\\.gz$", full.names = TRUE)
cat("\nExtracted files:\n")
print(extracted_files)


# ============================================================================
# 2. 加载 Alsaigh 2022 — 6 个样本 (3 患者 × AC/PA)
# ============================================================================
cat("\n========== Step 2: Extracting and Loading Alsaigh samples ==========\n")

alsaigh_samples <- list()
ext_dir <- "G:/13sample/Alsaigh_integration/GSE159677_raw/GSE159677_extracted/"

# 1. 找到所有真正包含表达矩阵的压缩包 (排除 moleculeinfo)
tar_files <- list.files(ext_dir, pattern = "featurebcmatrixfiltered\\.tar\\.gz$", full.names = TRUE)
cat("Found", length(tar_files), "matrix archives.\n")

for (t in tar_files) {
  # 提取 GSM 编号
  gsm_m <- regexpr("GSM[0-9]+", basename(t))
  gsm_id <- regmatches(basename(t), gsm_m)
  
  if (length(gsm_id) == 0) next
  
  sample_dir <- file.path(ext_dir, gsm_id)
  
  # 2. 解压压缩包
  if (!dir.exists(sample_dir) || length(list.files(sample_dir)) == 0) {
    dir.create(sample_dir, showWarnings = FALSE)
    cat("Extracting", basename(t), "...\n")
    untar(t, exdir = sample_dir)
  }
  
  # 3. 深度寻找矩阵文件 (无视它解压后可能嵌套了几层文件夹)
  matrix_files <- list.files(sample_dir, pattern = "matrix\\.mtx", recursive = TRUE, full.names = TRUE)
  
  if (length(matrix_files) > 0) {
    # 真正的读取目录是 matrix.mtx 所在的那个文件夹
    true_data_dir <- dirname(matrix_files[1])
    
    # 4. 暴力修复文件名：去除各种奇奇怪怪的前缀，强行改成 Seurat 认识的名字
    all_files <- list.files(true_data_dir, full.names = TRUE)
    for (f in all_files) {
      base_name <- basename(f)
      # 如果名字里包含 barcodes/features/genes/matrix，就把前面的乱七八糟前缀都删掉
      if (grepl("barcodes|features|genes|matrix", base_name, ignore.case = TRUE)) {
        new_name <- gsub(".*(barcodes|features|genes|matrix)", "\\1", base_name, ignore.case = TRUE)
        new_name <- gsub("genes", "features", new_name, ignore.case = TRUE) # 兼容 Seurat V3+ 要求
        
        # 执行重命名
        if (base_name != new_name) {
          file.rename(f, file.path(true_data_dir, new_name))
        }
      }
    }
    
    # 5. 映射临床标签
    if (gsm_id == "GSM4837523") { meta_name <- "AC1"; tissue <- "plaque" }
    else if (gsm_id == "GSM4837524") { meta_name <- "PA1"; tissue <- "control" }
    else if (gsm_id == "GSM4837525") { meta_name <- "AC2"; tissue <- "plaque" }
    else if (gsm_id == "GSM4837526") { meta_name <- "PA2"; tissue <- "control" }
    else if (gsm_id == "GSM4837527") { meta_name <- "AC3"; tissue <- "plaque" }
    else if (gsm_id == "GSM4837528") { meta_name <- "PA3"; tissue <- "control" }
    else { meta_name <- gsm_id; tissue <- "unknown" }
    
    cat("Loading", meta_name, "from", true_data_dir, "...\n")
    
    # 6. 读取并创建 Seurat 对象
    tryCatch({
      counts <- Read10X(data.dir = true_data_dir)
      obj <- CreateSeuratObject(counts = counts, project = meta_name,
                                min.cells = 3, min.features = 200)
      obj$orig.ident <- meta_name
      obj$dataset <- "Alsaigh_2022"
      obj$tissue_type <- tissue
      
      alsaigh_samples[[meta_name]] <- obj
      cat("  -> Success! Extracted", ncol(obj), "cells.\n")
    }, error = function(e) {
      cat("  -> FAILED:", e$message, "\n")
    })
    
  } else {
    cat("  -> ERROR: Could not find matrix.mtx inside", sample_dir, "\n")
  }
}

cat("\n========== Loading Complete ==========\n")
cat("Successfully loaded", length(alsaigh_samples), "Alsaigh samples.\n")

if (length(alsaigh_samples) > 0) {
  cat("Total cells (Alsaigh):", sum(sapply(alsaigh_samples, ncol)), "\n\n")
  
  # 合并 Alsaigh 6 个样本
  alsaigh_merged <- merge(alsaigh_samples[[1]], 
                          y = alsaigh_samples[-1],
                          add.cell.ids = names(alsaigh_samples))
  cat("Merged Alsaigh object:", ncol(alsaigh_merged), "cells\n")
  
  # QC
  alsaigh_merged[["percent.mt"]] <- PercentageFeatureSet(alsaigh_merged, pattern = "^MT-")
  alsaigh_merged <- subset(alsaigh_merged, 
                           subset = nFeature_RNA > 200 & 
                             nFeature_RNA < 5000 & 
                             percent.mt < 15)
  cat("After QC:", ncol(alsaigh_merged), "cells\n")
} else {
  stop("FATAL ERROR: No samples were loaded. Please check the error messages above.")
}


## ============================================================================
# 3.5 内存优化：对象瘦身与均衡降采样 (在合并前执行)
# ============================================================================
cat("\n========== Step 3.5: Diet and Downsampling for memory limit ==========\n")

# 1. 极限瘦身 (DietSeurat)：
# 您的原对象可能自带了以前跑过的 SCT/integrated assays，以及庞大的 PCA/UMAP 数据。
# 这些在重新整合时都不需要了，果断剔除，只保留最干净的 RNA counts 矩阵。
your_obj <- DietSeurat(your_obj, assays = "RNA", dimreducs = NULL)
alsaigh_merged <- DietSeurat(alsaigh_merged, assays = "RNA", dimreducs = NULL)
gc() # 强制回收内存

# 2. 均衡降采样 (Stratified Downsampling)：
# 根据您的细胞注释，每个细胞亚群最多随机抽取 2500 个细胞保留。
# （细胞数不到 2500 的亚群，如少量的 ECs 或特定免疫细胞，会被 100% 完整保留）
Idents(your_obj) <- "celltype"
your_obj <- subset(your_obj, downsample = 2500)

# 对 Alsaigh 数据也按样本进行一定程度的抽样（每个样本最多保留 5000 个细胞，6个样本最多3万）
Idents(alsaigh_merged) <- "orig.ident"
alsaigh_merged <- subset(alsaigh_merged, downsample = 5000)

cat("\n--- After Downsampling ---\n")
cat("Your dataset reduced to:", ncol(your_obj), "cells\n")
cat("Alsaigh dataset reduced to:", ncol(alsaigh_merged), "cells\n")
gc()

# ============================================================================
# 4. 合并两个数据集

cat("\n========== Step 4: Merging two datasets ==========\n")

# 转换为 V4 Assay，彻底规避 V5 JoinLayers 的内存黑洞
your_obj[["RNA"]] <- as(your_obj[["RNA"]], "Assay")
alsaigh_merged[["RNA"]] <- as(alsaigh_merged[["RNA"]], "Assay")
gc()

cat("Merging matrices...\n")
combined <- merge(your_obj, y = alsaigh_merged, 
                  add.cell.ids = c("our", "alsaigh"))

# 合并后立即销毁原对象，释放巨大空间给后续的 ScaleData 和 PCA
rm(your_obj, alsaigh_merged)
gc()

cat("Merged total cells:", ncol(combined), "\n")
print(table(combined$dataset))

# 标准化 + 高变基因
cat("\nNormalizing + finding variable features...\n")
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, nfeatures = 3000, verbose = FALSE)

# 再次内存保护：ScaleData 只对高变基因进行 scale，不 scale 所有基因（大大节省内存）
combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

# 5. ★★★ Harmony 整合（
# ============================================================================
cat("\n========== Step 5: Harmony integration ==========\n")

# 
library(harmony)
library(patchwork)

# 运行 Harmony 消除数据集之间的批次效应
combined <- RunHarmony(combined, 
                       group.by.vars = "dataset",
                       reduction.use = "pca",
                       theta = 2,              # theta=2 适度对齐，防止过度矫正抹杀生物学差异
                       max_iter = 20,          # 【修改点】：将 max.iter.harmony 改为 max_iter
                       verbose = FALSE)

# 使用 Harmony 降维结果计算 UMAP 和聚类
cat("Running UMAP and Clustering...\n")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)


# 6. ★★★ Fig 1F 主MAP 三张并排
# ============================================================================
cat("\n========== Step 6: Generating Fig 1F UMAP ==========\n")

# Panel 1: 按数据集
p1 <- DimPlot(combined, reduction = "umap",
              group.by = "dataset",
              cols = c("Our_study" = "#C71000FF",  # 红
                       "Alsaigh_2022" = "#008EA0FF"),  # 青
              pt.size = 0.3, raster = FALSE) +
  ggtitle("By dataset (after Harmony)") +
  theme(plot.title = element_text(face = "bold", size = 14))

# Panel 2: 按您的 cell type 注释（仅您的细胞）
your_celltype <- combined$celltype
your_celltype[combined$dataset == "Alsaigh_2022"] <- NA

combined$your_celltype <- your_celltype

p2 <- DimPlot(combined, reduction = "umap",
              group.by = "your_celltype",
              cols = celltype_colors,
              pt.size = 0.3, raster = FALSE,
              na.value = "lightgrey") +
  ggtitle("Our cell types\n(grey = Alsaigh cells)") +
  theme(plot.title = element_text(face = "bold", size = 14))

# Panel 3: 按整合后的新 cluster
p3 <- DimPlot(combined, reduction = "umap",
              group.by = "seurat_clusters",
              label = TRUE, label.size = 4,
              pt.size = 0.3, raster = FALSE) +
  ggtitle("Integrated clusters") +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none")

# 拼接主图
fig1F_main <- p1 | p2 | p3 +
  plot_annotation(
    title = "Fig 1F. Integration with Alsaigh 2022 (GSE159677)",
    subtitle = sprintf("Total: %d cells from %d patients (12 ours + 3 Alsaigh)",
                       ncol(combined), 12 + 3),
    theme = theme(plot.title = element_text(face = "bold", size = 16)))

ggsave(paste0(out_dir, "Fig1F_integrated_UMAP.pdf"), fig1F_main,
       width = 18, height = 6, device = "pdf")
ggsave(paste0(out_dir, "Fig1F_integrated_UMAP.png"), fig1F_main,
       width = 18, height = 6, dpi = 300)

cat("??? Saved Fig1F_integrated_UMAP.pdf\n")


# ============================================================================
# 7. Cell type concordance plot
# ============================================================================
cat("\n========== Step 7: Cell type concordance ==========\n")

# 把您每个 cell type 的细胞分到整合后的 cluster 里，看分布
concordance_df <- combined@meta.data %>%
  filter(!is.na(your_celltype)) %>%
  group_by(your_celltype, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(your_celltype) %>%
  mutate(pct = n / sum(n) * 100)

p_concordance <- ggplot(concordance_df, 
                        aes(x = seurat_clusters, y = your_celltype, 
                            fill = pct)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#C71000FF",
                      name = "% of cells") +
  labs(title = "Cell type concordance: Our annotation vs Integrated clusters",
       x = "Integrated cluster", y = "Our cell type") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(out_dir, "Fig1F_celltype_concordance.pdf"), p_concordance,
       width = 9, height = 7, device = "pdf")
ggsave(paste0(out_dir, "Fig1F_celltype_concordance.png"), p_concordance,
       width = 9, height = 7, dpi = 300)


# 8. Marker 表达一致性（POSTN/TREM2/ACTA2 在两数据集中都正确定位）
# ============================================================================
cat("\n========== Step 8: Marker concordance ==========\n")

key_markers <- c("ACTA2","MYH11",       # SMC contractile
                 "POSTN","COL1A1",       # Fibroblast
                 "TREM2","CD68",         # Macrophage
                 "CDH5","PECAM1",        # EC
                 "CD3D","CD3E")          # T cell

key_markers <- intersect(key_markers, rownames(combined))

# 按数据集分别画 FeaturePlot
gene_panels <- list()
for (g in key_markers) {
  
  p_our <- FeaturePlot(subset(combined, dataset == "Our_study"),
                       features = g, reduction = "umap",
                       cols = c("lightgrey","red"), 
                       pt.size = 0.2, raster = FALSE) +
    ggtitle(paste(g, "- Our")) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "none",
          axis.text = element_blank(), axis.ticks = element_blank())
  
  p_als <- FeaturePlot(subset(combined, dataset == "Alsaigh_2022"),
                       features = g, reduction = "umap",
                       cols = c("lightgrey","red"),
                       pt.size = 0.2, raster = FALSE) +
    ggtitle(paste(g, "- Alsaigh")) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "none",
          axis.text = element_blank(), axis.ticks = element_blank())
  
  gene_panels[[paste0(g,"_our")]] <- p_our
  gene_panels[[paste0(g,"_als")]] <- p_als
}

p_marker_concord <- wrap_plots(gene_panels, ncol = 4)

ggsave(paste0(out_dir, "Fig1F_marker_concordance.pdf"), p_marker_concord,
       width = 16, height = 4 * ceiling(length(key_markers)/2),
       device = "pdf", limitsize = FALSE)
ggsave(paste0(out_dir, "Fig1F_marker_concordance.png"), p_marker_concord,
       width = 16, height = 4 * ceiling(length(key_markers)/2),
       dpi = 300, limitsize = FALSE)


# ============================================================================
# 9. thology vs Symptom 分组对比
cat("\n========== Step 9: Pathology vs Symptom comparison ==========\n")


# 计算两个 healing markers (POSTN+TREM2) 的 module score
healing_markers <- list(
  fibroblast_healing = c("POSTN","TNC","COL1A1","COL3A1","LOXL2"),
  TREM2_LAM = c("TREM2","CD9","FOLR2","CD163","CTSD"))

combined <- AddModuleScore(combined, features = healing_markers,
                           name = c("FibroHealing","TREM2_LAM"),
                           search = TRUE)

# Plot 1: 按 tissue_type(类似 Alsaigh 设计)
plot_data_tissue <- combined@meta.data %>%
  filter(!is.na(tissue_type)) %>%
  dplyr::select(tissue_type, dataset, FibroHealing1, TREM2_LAM2)

p_tissue_healing <- ggplot(plot_data_tissue,
                           aes(x = tissue_type, y = FibroHealing1, 
                               fill = tissue_type)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white") +
  facet_wrap(~ dataset) +
  scale_fill_manual(values = c("control" = "#FFE705FF",
                               "plaque" = "#C71000FF")) +
  labs(title = "By tissue type (control vs plaque)\n[Alsaigh-like grouping]",
       y = "Fibroblast healing module score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

# Plot 2: 按 pathology group (仅您的数据)
plot_data_path <- combined@meta.data %>%
  filter(dataset == "Our_study", !is.na(group))

# 排序
plot_data_path$group <- factor(plot_data_path$group, 
                               levels = c("ASYM","SYM"))

p_path_healing <- ggplot(plot_data_path,
                         aes(x = group, y = FibroHealing1, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c(
                               "ASYM"   = "#C71000FF",
                               "SYM"    = "#8BEAFFFF")) +
  labs(title = "By pathology (Stable/Unstable)\n[Our pathology-anchored grouping]",
       y = "Fibroblast healing module score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

# TREM2 module score 同样
p_tissue_lam <- ggplot(plot_data_tissue,
                       aes(x = tissue_type, y = TREM2_LAM2, fill = tissue_type)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white") +
  facet_wrap(~ dataset) +
  scale_fill_manual(values = c("control" = "#FFE705FF",
                               "plaque" = "#C71000FF")) +
  labs(title = "TREM2+ LAM signature by tissue type",
       y = "TREM2+ LAM module score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

p_path_lam <- ggplot(plot_data_path,
                     aes(x = group, y = TREM2_LAM2, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c(
                               "ASYM"   = "#C71000FF",
                               "SYM"    = "#8BEAFFFF")) +
  labs(title = "TREM2+ LAM signature by pathology",
       y = "TREM2+ LAM module score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

# 拼接
fig1F_path_vs_sym <- (p_tissue_healing | p_path_healing) /
  (p_tissue_lam | p_path_lam) +
  plot_annotation(
    title = "Fig 1F. Pathology-anchored grouping reveals healing program (vs tissue-type grouping)",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(paste0(out_dir, "Fig1F_pathology_vs_symptom.pdf"), fig1F_path_vs_sym,
       width = 14, height = 11, device = "pdf")
ggsave(paste0(out_dir, "Fig1F_pathology_vs_symptom.png"), fig1F_path_vs_sym,
       width = 14, height = 11, dpi = 300)

cat("??? Saved Fig1F_pathology_vs_symptom.pdf (核心方法学证据图)\n")


# 10. 统计检验 + 报告

# Wilcoxon 检验
wilcox_tissue <- wilcox.test(
  FibroHealing1 ~ tissue_type,
  data = plot_data_tissue[plot_data_tissue$dataset == "Our_study", ])

# Kruskal-Wallis（多组）
kw_path <- kruskal.test(FibroHealing1 ~ group, data = plot_data_path)

# 报告
sink(paste0(out_dir, "Integration_REPORT.txt"))

cat("=================================================================\n")
cat("ALSAIGH 2022 INTEGRATION REPORT\n")
cat("=================================================================\n\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

cat("--- DATASET SUMMARY ---\n")
cat("Our study:", sum(combined$dataset == "Our_study"), "cells, 12 patients\n")
cat("Alsaigh 2022:", sum(combined$dataset == "Alsaigh_2022"), "cells, 3 patients\n")
cat("Total integrated:", ncol(combined), "cells, 15 patients\n\n")

cat("--- INTEGRATION QUALITY ---\n")
cat("Method: Harmony (theta=2)\n")
cat("Resolution:", 0.5, "\n")
cat("Number of integrated clusters:", length(levels(combined$seurat_clusters)), "\n\n")

cat("--- KEY STATISTICAL RESULTS ---\n\n")
cat("Fibroblast healing module score:\n")
cat("  By tissue_type (Alsaigh-like): Wilcoxon P =", 
    format.pval(wilcox_tissue$p.value), "\n")
cat("  By pathology (our design): Kruskal-Wallis P =", 
    format.pval(kw_path$p.value), "\n\n")

cat("--- INTERPRETATION ---\n")
cat("If pathology grouping yields LOWER P-value than tissue grouping:\n")
cat("  → Pathology-anchored design has BETTER discrimination power\n")
cat("  → This is the methodological advantage of our study\n\n")

cat("--- FIGURES ---\n")
cat("????????? Fig1F_pathology_vs_symptom.pdf (核心方法学证据)\n")
cat("?????? Fig1F_integrated_UMAP.pdf (主图)\n")
cat("?????? Fig1F_celltype_concordance.pdf (一致性证据)\n")
cat("??? Fig1F_marker_concordance.pdf (Marker 一致性)\n")

sink()

# 保存整合对象
saveRDS(combined, paste0(out_dir, "integrated_obj_with_Alsaigh.rds"))

#Pathology vs Symptom 对比图 — 限定在相关细胞群中比较

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
out_dir <- "G:/13sample/Alsaigh_integration/"

# 关键：只在 fibroblast/SMC 群体中评估 fibroblast healing
#       只在 macrophage 群体中评估 TREM2+ LAM
# ============================================================================

# 把您的 cell type 信息也加给 Alsaigh 细胞（基于整合 cluster）
# Cluster 3 = Fibroblasts，Cluster 0/1 = T cells，等等
# 我们用整合 cluster 来分类细胞类型

# 根据 concordance 热图，定义大类
fibroblast_smc_clusters <- c("3","5","7","9","11","12","16","20")  # SMCs+Fibroblasts
macrophage_clusters     <- c("4","6","13","14")                       # Macro+Mono+DC

# 提取相关细胞做分析
combined$broad_type <- "Other"
combined$broad_type[combined$seurat_clusters %in% fibroblast_smc_clusters] <- "Fibroblast_SMC"
combined$broad_type[combined$seurat_clusters %in% macrophage_clusters]     <- "Macrophage"

# 检查数量
cat("Cell counts by broad type:\n")
print(table(combined$broad_type))

# ============================================================================
# Plot 1: Fibroblast healing — 仅在 Fibroblast/SMC 中比较
# ============================================================================

# 数据子集：仅 fibroblast/SMC 细胞
fibro_data <- combined@meta.data %>%
  filter(broad_type == "Fibroblast_SMC")

# Alsaigh 的：tissue_type (control vs plaque)
fibro_alsaigh <- fibro_data %>%
  filter(dataset == "Alsaigh_2022")

p_alsaigh_fibro <- ggplot(fibro_alsaigh,
                          aes(x = tissue_type, y = FibroHealing1, 
                              fill = tissue_type)) +
  geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = c("control" = "#FFE705FF",
                               "plaque" = "#C71000FF")) +
  labs(title = "Alsaigh 2022: Tissue-type grouping",
       subtitle = "Fibroblast healing in Fibroblast/SMC cells",
       y = "Fibroblast healing module score", x = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

# 添加统计 P 值
p_val_alsaigh <- wilcox.test(FibroHealing1 ~ tissue_type, data = fibro_alsaigh)$p.value
p_alsaigh_fibro <- p_alsaigh_fibro +
  labs(caption = paste0("Wilcoxon P = ", format.pval(p_val_alsaigh, digits = 3)))


# 您的：pathology (ASYM/SYM)
fibro_ours <- fibro_data %>%
  filter(dataset == "Our_study", !is.na(group))
fibro_ours$group <- factor(fibro_ours$group, levels = c("ASYM","SYM"))

p_ours_fibro <- ggplot(fibro_ours,
                       aes(x = group, y = FibroHealing1, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = c(
                               "ASYM"   = "#C71000FF",
                               "SYM"    = "#8BEAFFFF")) +
  labs(title = "Our study: Pathology-anchored grouping",
       subtitle = "Fibroblast healing in Fibroblast/SMC cells",
       y = "Fibroblast healing module score", x = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

# 添加统计 P 值
p_val_ours <- kruskal.test(FibroHealing1 ~ group, data = fibro_ours)$p.value
p_ours_fibro <- p_ours_fibro +
  labs(caption = paste0("Kruskal-Wallis P = ", format.pval(p_val_ours, digits = 3)))


# ============================================================================
# Plot 2: TREM2+ LAM — 仅在 Macrophage 中比较
# ============================================================================

mac_data <- combined@meta.data %>%
  filter(broad_type == "Macrophage")

# Alsaigh 的
mac_alsaigh <- mac_data %>%
  filter(dataset == "Alsaigh_2022")

p_alsaigh_lam <- ggplot(mac_alsaigh,
                        aes(x = tissue_type, y = TREM2_LAM2, fill = tissue_type)) +
  geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = c("control" = "#FFE705FF",
                               "plaque" = "#C71000FF")) +
  labs(title = "Alsaigh 2022: Tissue-type grouping",
       subtitle = "TREM2+ LAM signature in Macrophages",
       y = "TREM2+ LAM module score", x = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

p_val_alsaigh_lam <- wilcox.test(TREM2_LAM2 ~ tissue_type, data = mac_alsaigh)$p.value
p_alsaigh_lam <- p_alsaigh_lam +
  labs(caption = paste0("Wilcoxon P = ", format.pval(p_val_alsaigh_lam, digits = 3)))


# 您的
mac_ours <- mac_data %>%
  filter(dataset == "Our_study", !is.na(group))
mac_ours$group <- factor(mac_ours$group, levels = c("ASYM","SYM"))

p_ours_lam <- ggplot(mac_ours,
                     aes(x = group, y = TREM2_LAM2, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = c(
                               "ASYM"   = "#C71000FF",
                               "SYM"    = "#8BEAFFFF")) +
  labs(title = "Our study: Pathology-anchored grouping",
       subtitle = "TREM2+ LAM signature in Macrophages",
       y = "TREM2+ LAM module score", x = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

p_val_ours_lam <- kruskal.test(TREM2_LAM2 ~ group, data = mac_ours)$p.value
p_ours_lam <- p_ours_lam +
  labs(caption = paste0("Kruskal-Wallis P = ", format.pval(p_val_ours_lam, digits = 3)))


# ============================================================================
# 拼接 + 保存
# ============================================================================
fig1F_v2 <- (p_alsaigh_fibro | p_ours_fibro) /
  (p_alsaigh_lam   | p_ours_lam) +
  plot_annotation(
    title = "Fig 1F. Pathology-anchored grouping reveals healing program in cell-type-specific manner",
    subtitle = "Module scores evaluated in target cell populations (Fibroblast/SMC for healing; Macrophages for LAM)",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, color = "darkgrey")))

ggsave(paste0(out_dir, "Fig1F_pathology_vs_symptom_FIXED.pdf"), fig1F_v2,
       width = 13, height = 10, device = "pdf")
ggsave(paste0(out_dir, "Fig1F_pathology_vs_symptom_FIXED.png"), fig1F_v2,
       width = 13, height = 10, dpi = 300)

cat("\n??? Saved Fig1F_pathology_vs_symptom_FIXED.pdf\n")

# 批次验证: (1) 按13样本上色的 tSNE  (2) cluster × sample 细胞数矩阵
#加载包
rm(list = ls())
gc()
setwd("G:/13sample/")
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidydr)
combined131<- readRDS(file="combined131.rds")

# ============================================================================
library(Seurat); library(tidyverse); library(ggplot2)

## 先确认样本ID列名 (常见 orig.ident; 你的对象可能不同, 先看一眼)
cat("meta 列名:\n"); print(colnames(combined131@meta.data))
cat("\n可能的样本列预览:\n")
for (col in c("orig.ident","sample","Sample","samples")) {
  if (col %in% colnames(combined131@meta.data)) {
    cat(col, ":", paste(unique(combined131@meta.data[[col]]), collapse=", "), "\n")
  }
}


# ============================================================================
# 批次 tSNE —
# ============================================================================
library(Seurat); library(tidyverse); library(ggplot2); library(ggrastr)

sample_col <- "sample"
samples_all <- sort(unique(combined131@meta.data[[sample_col]]))
n_sample <- length(samples_all)
sample_cols <- setNames(
  colorRampPalette(c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
                     "#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85",
                     "#FDAF91","#1B9E77","#666666"))(n_sample),
  samples_all)

df <- data.frame(combined131@reductions$tsne@cell.embeddings,
                 sample = combined131@meta.data[[sample_col]])
colnames(df)[1:2] <- c("tSNE_1","tSNE_2")
df <- df[sample(nrow(df)), ]   # 打乱绘制顺序

# ── (1) 主图: 点大小还原 V1 (DimPlot 默认 pt.size≈0.1 视觉, 这里用 0.5 接近 V1 观感) ──
p_main <- ggplot(df, aes(tSNE_1, tSNE_2, color = sample)) +
  rasterise(geom_point(size = 0.5, alpha = 0.8, stroke = 0), dpi = 600) +  # ★ 栅格化, 点同V1大小
  scale_color_manual(values = sample_cols, name = "Sample") +
  labs(title = "Batch mixing across 13 samples", x = "tSNE_1", y = "tSNE_2") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(color = "black"),
        legend.text = element_text(size = 9))

ggsave("G:/13sample/SuppFig_batch_tSNE_bySample_v3.pdf", p_main,
       width = 8, height = 6, device = cairo_pdf)
ggsave("G:/13sample/SuppFig_batch_tSNE_bySample_v3.png", p_main,
       width = 8, height = 6, dpi = 600)

# ── (2) split 版: 每样本一格, 当前样本高亮其余灰底 (看得清每个样本分布) ──
# 关键改进: 每格画"全部细胞灰色背景 + 当前样本彩色高亮", 这样每格都能看清该样本铺在哪
df_bg <- df %>% dplyr::select(tSNE_1, tSNE_2)   # 背景: 所有细胞(无 sample 分面)

p_split <- ggplot() +
  rasterise(geom_point(data = df_bg, aes(tSNE_1, tSNE_2),
                       color = "grey88", size = 0.3, stroke = 0), dpi = 600) +  # 灰背景
  rasterise(geom_point(data = df, aes(tSNE_1, tSNE_2, color = sample),
                       size = 0.4, stroke = 0), dpi = 600) +                    # 高亮当前样本
  scale_color_manual(values = sample_cols, guide = "none") +
  facet_wrap(~ sample, ncol = 5) +                # ★ 5列布局
  labs(title = "Per-sample distribution (each panel: one sample highlighted)",
       x = "tSNE_1", y = "tSNE_2") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey95", color = NA),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

# ★ 画布够大才能看清每格 (13样本5列=3行, 给足尺寸)
ggsave("G:/13sample/SuppFig_batch_tSNE_split_v3.pdf", p_split,
       width = 15, height = 9, device = cairo_pdf)
ggsave("G:/13sample/SuppFig_batch_tSNE_split_v3.png", p_split,
       width = 15, height = 9, dpi = 600)

--------------------------------------------------------------------------------
# FigG/H: donor-level analysis of module scores
# Purpose
#   Correct the unit-of-analysis error in the original cell-level violin plots.
#   Existing cell-level module scores are retained, but they are first aggregated
#   within each donor and target cell population. All statistical tests then use
#   donors, not cells, as independent observations.
#
# Design
#   Fig. 1G: FibroHealing1 in Fibroblast/SMC cells
#   Fig. 1H: TREM2_LAM2 in Macrophages
#   Alsaigh: 3 matched donors, control (PA1-3) versus plaque (AC1-3)
#             -> paired Wilcoxon signed-rank test
#   Our study: 6 ASYM versus 6 SYM donors
#              -> independent exact Wilcoxon rank-sum test
#
# Important terminology
#   This is a donor-level aggregation of cell-level module scores. It should not
#   be described as raw-count pseudobulk differential expression.
################################################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# 1. Paths and analysis settings
# ==============================================================================

input_rds <- "G:/13sample/Alsaigh_integration/integrated_obj_with_Alsaigh.rds"

# A separate output folder is used so the original figures and scripts are not
# overwritten.
out_dir <- file.path(
  "G:/13sample/A 颈动脉单细胞文章/TREM2??? macrophage–reparative fibroblast healing axis",
  "Communications_Biology_submission/Communications_Biology_submission",
  "Fig1GH_donor_level_results"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_iterations <- 10000L
set.seed(20260831)

group_colors <- c(
  "ASYM" = "#FFE705FF",
  "SYM"  = "#C71000FF"
)

tissue_colors <- c(
  "control" = "#A8A8A8",
  "plaque"  = "#C71000FF"
)

# ==============================================================================
# 2. Load the existing integrated object and validate metadata
# ==============================================================================

if (!exists("combined")) {
  combined <- readRDS(input_rds)
}

required_metadata <- c(
  "dataset", "orig.ident", "sample", "group", "tissue_type",
  "seurat_clusters", "FibroHealing1", "TREM2_LAM2"
)

missing_metadata <- setdiff(required_metadata, colnames(combined@meta.data))
if (length(missing_metadata) > 0L) {
  stop(
    "Missing required metadata column(s): ",
    paste(missing_metadata, collapse = ", ")
  )
}

# Keep the same integrated-cluster definitions used by the original Fig. 1 code.
fibroblast_smc_clusters <- c("3", "5", "7", "9", "11", "12", "16", "20")
macrophage_clusters     <- c("4", "6", "13", "14")

md <- combined@meta.data %>%
  mutate(
    cell_id = rownames(combined@meta.data),
    cluster_chr = as.character(seurat_clusters),
    broad_type = case_when(
      cluster_chr %in% fibroblast_smc_clusters ~ "Fibroblast_SMC",
      cluster_chr %in% macrophage_clusters     ~ "Macrophage",
      TRUE                                     ~ "Other"
    )
  )

# In the Alsaigh data, PA1/AC1, PA2/AC2 and PA3/AC3 are matched tissue
# specimens from the same three donors. Do not use orig.ident as six independent
# donor IDs.
alsaigh_specimens <- unique(md$orig.ident[md$dataset == "Alsaigh_2022"])
if (any(!grepl("^(PA|AC)[0-9]+$", alsaigh_specimens))) {
  stop(
    "Unexpected Alsaigh orig.ident format. Expected PA1-3 and AC1-3; found: ",
    paste(sort(alsaigh_specimens), collapse = ", ")
  )
}

md <- md %>%
  mutate(
    donor_id = case_when(
      dataset == "Alsaigh_2022" ~ paste0(
        "Alsaigh_", sub("^[A-Za-z]+", "", as.character(orig.ident))
      ),
      dataset == "Our_study" ~ as.character(sample),
      TRUE ~ NA_character_
    )
  )

if (any(is.na(md$donor_id[md$dataset == "Our_study" &
                          md$group %in% c("ASYM", "SYM")]))) {
  stop("At least one ASYM/SYM cell has a missing donor ID in metadata column 'sample'.")
}

# ==============================================================================
# 3. Aggregate cell-level module scores to one value per donor
# ==============================================================================

aggregate_donor_scores <- function(data, target_type, cohort, group_var,
                                   score_var, signature_label) {
  data %>%
    filter(
      broad_type == target_type,
      dataset == cohort,
      !is.na(.data[[group_var]]),
      is.finite(.data[[score_var]])
    ) %>%
    transmute(
      donor_id,
      specimen_id = as.character(orig.ident),
      group = as.character(.data[[group_var]]),
      cell_score = as.numeric(.data[[score_var]])
    ) %>%
    group_by(donor_id, group) %>%
    summarise(
      specimen_id = paste(sort(unique(specimen_id)), collapse = ";"),
      n_cells = n(),
      score_mean = mean(cell_score),
      score_median = median(cell_score),
      score_sd_cells = sd(cell_score),
      .groups = "drop"
    ) %>%
    mutate(
      cohort = cohort,
      signature = signature_label,
      target_cell_population = target_type
    )
}

fibro_alsaigh <- aggregate_donor_scores(
  md, "Fibroblast_SMC", "Alsaigh_2022", "tissue_type",
  "FibroHealing1", "Fibroblast healing"
) %>%
  filter(group %in% c("control", "plaque"))

fibro_ours <- aggregate_donor_scores(
  md, "Fibroblast_SMC", "Our_study", "group",
  "FibroHealing1", "Fibroblast healing"
) %>%
  filter(group %in% c("ASYM", "SYM"))

mac_alsaigh <- aggregate_donor_scores(
  md, "Macrophage", "Alsaigh_2022", "tissue_type",
  "TREM2_LAM2", "TREM2+ LAM"
) %>%
  filter(group %in% c("control", "plaque"))

mac_ours <- aggregate_donor_scores(
  md, "Macrophage", "Our_study", "group",
  "TREM2_LAM2", "TREM2+ LAM"
) %>%
  filter(group %in% c("ASYM", "SYM"))

donor_scores <- bind_rows(
  fibro_alsaigh,
  fibro_ours,
  mac_alsaigh,
  mac_ours
)

# Report low-cell donors rather than silently removing them. In the current
# object every target-population donor has far more than 20 cells.
low_cell_donors <- donor_scores %>% filter(n_cells < 20L)
if (nrow(low_cell_donors) > 0L) {
  warning(
    "The following donor/signature rows contain fewer than 20 cells and are retained: ",
    paste(low_cell_donors$donor_id, collapse = ", ")
  )
}

# Validate the matched external design: each Alsaigh donor must contribute one
# control and one plaque aggregate for each signature.
validate_alsaigh_pairs <- function(data, panel_name) {
  pair_check <- data %>%
    count(donor_id, group) %>%
    pivot_wider(names_from = group, values_from = n, values_fill = 0L)
  
  needed <- c("control", "plaque")
  if (!all(needed %in% colnames(pair_check)) ||
      any(pair_check$control != 1L) || any(pair_check$plaque != 1L)) {
    stop(panel_name, ": incomplete or duplicated Alsaigh donor pairs.")
  }
}

validate_independent_donors <- function(data, panel_name) {
  donor_group_n <- data %>% count(donor_id) %>% filter(n != 1L)
  if (nrow(donor_group_n) > 0L) {
    stop(panel_name, ": at least one donor occurs in more than one ASYM/SYM group.")
  }
  if (!all(c("ASYM", "SYM") %in% unique(data$group))) {
    stop(panel_name, ": ASYM or SYM donors are missing.")
  }
}

validate_alsaigh_pairs(fibro_alsaigh, "Fig. 1G / Alsaigh")
validate_alsaigh_pairs(mac_alsaigh, "Fig. 1H / Alsaigh")
validate_independent_donors(fibro_ours, "Fig. 1G / Our study")
validate_independent_donors(mac_ours, "Fig. 1H / Our study")

# ==============================================================================
# 4. Donor-level statistics
# ==============================================================================

bootstrap_mean_difference <- function(x1, x0, paired = FALSE,
                                      iterations = bootstrap_iterations) {
  stopifnot(length(x1) > 1L, length(x0) > 1L)
  
  if (paired) {
    stopifnot(length(x1) == length(x0))
    delta <- x1 - x0
    boot_values <- replicate(
      iterations,
      mean(sample(delta, size = length(delta), replace = TRUE))
    )
    estimate <- mean(delta)
  } else {
    boot_values <- replicate(
      iterations,
      mean(sample(x1, size = length(x1), replace = TRUE)) -
        mean(sample(x0, size = length(x0), replace = TRUE))
    )
    estimate <- mean(x1) - mean(x0)
  }
  
  ci <- unname(quantile(boot_values, probs = c(0.025, 0.975), na.rm = TRUE))
  c(estimate = estimate, ci_low = ci[1], ci_high = ci[2])
}

paired_wilcoxon <- function(data, score_col, level0, level1) {
  wide <- data %>%
    select(donor_id, group, score = all_of(score_col)) %>%
    pivot_wider(names_from = group, values_from = score) %>%
    arrange(donor_id)
  
  if (any(!complete.cases(wide[, c(level0, level1)]))) {
    stop("Missing value in a paired donor comparison.")
  }
  
  x0 <- wide[[level0]]
  x1 <- wide[[level1]]
  delta <- x1 - x0
  
  # Base R cannot calculate an exact signed-rank P value in the presence of
  # zero differences or tied absolute differences.
  exact_ok <- all(delta != 0) && !anyDuplicated(abs(delta))
  fit <- wilcox.test(
    x1, x0,
    paired = TRUE,
    exact = exact_ok,
    correct = !exact_ok,
    alternative = "two.sided"
  )
  
  nonzero <- delta != 0
  if (any(nonzero)) {
    rank_abs <- rank(abs(delta[nonzero]), ties.method = "average")
    rank_biserial <- (
      sum(rank_abs[delta[nonzero] > 0]) -
        sum(rank_abs[delta[nonzero] < 0])
    ) / sum(rank_abs)
  } else {
    rank_biserial <- 0
  }
  
  boot <- bootstrap_mean_difference(x1, x0, paired = TRUE)
  
  tibble(
    comparison = paste0(level1, " - ", level0),
    test = ifelse(
      exact_ok,
      "Exact paired Wilcoxon signed-rank",
      "Paired Wilcoxon signed-rank with continuity correction"
    ),
    n_level0 = length(x0),
    n_level1 = length(x1),
    statistic = unname(fit$statistic),
    p_value = fit$p.value,
    rank_biserial = rank_biserial,
    mean_difference = unname(boot["estimate"]),
    mean_difference_ci_low = unname(boot["ci_low"]),
    mean_difference_ci_high = unname(boot["ci_high"])
  )
}

independent_wilcoxon <- function(data, score_col, level0, level1) {
  x0 <- data[[score_col]][data$group == level0]
  x1 <- data[[score_col]][data$group == level1]
  
  # No ties in the pooled donor-level values -> exact Mann-Whitney P value.
  exact_ok <- !anyDuplicated(c(x0, x1))
  fit <- wilcox.test(
    x1, x0,
    paired = FALSE,
    exact = exact_ok,
    correct = !exact_ok,
    alternative = "two.sided"
  )
  
  # Positive values mean that level1 tends to have higher scores than level0.
  mann_whitney_u <- unname(fit$statistic)
  rank_biserial <- 2 * mann_whitney_u / (length(x1) * length(x0)) - 1
  boot <- bootstrap_mean_difference(x1, x0, paired = FALSE)
  
  tibble(
    comparison = paste0(level1, " - ", level0),
    test = ifelse(
      exact_ok,
      "Exact Wilcoxon rank-sum (Mann-Whitney)",
      "Wilcoxon rank-sum with continuity correction"
    ),
    n_level0 = length(x0),
    n_level1 = length(x1),
    statistic = unname(fit$statistic),
    p_value = fit$p.value,
    rank_biserial = rank_biserial,
    mean_difference = unname(boot["estimate"]),
    mean_difference_ci_low = unname(boot["ci_low"]),
    mean_difference_ci_high = unname(boot["ci_high"])
  )
}

run_four_tests <- function(score_col) {
  bind_rows(
    paired_wilcoxon(fibro_alsaigh, score_col, "control", "plaque") %>%
      mutate(panel = "Fig1G", cohort = "Alsaigh_2022",
             signature = "Fibroblast healing"),
    independent_wilcoxon(fibro_ours, score_col, "ASYM", "SYM") %>%
      mutate(panel = "Fig1G", cohort = "Our_study",
             signature = "Fibroblast healing"),
    paired_wilcoxon(mac_alsaigh, score_col, "control", "plaque") %>%
      mutate(panel = "Fig1H", cohort = "Alsaigh_2022",
             signature = "TREM2+ LAM"),
    independent_wilcoxon(mac_ours, score_col, "ASYM", "SYM") %>%
      mutate(panel = "Fig1H", cohort = "Our_study",
             signature = "TREM2+ LAM")
  ) %>%
    relocate(panel, cohort, signature) %>%
    mutate(
      aggregation = score_col,
      p_holm_four_tests = p.adjust(p_value, method = "holm")
    )
}

# Primary analysis: mean cell-level module score within each donor.
primary_tests <- run_four_tests("score_mean")

# Sensitivity analysis only: donor median rather than donor mean. This result is
# exported but is not used to replace the prespecified primary figure.
sensitivity_tests <- run_four_tests("score_median")

descriptive_stats <- donor_scores %>%
  group_by(signature, cohort, group) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    total_cells = sum(n_cells),
    donor_score_mean = mean(score_mean),
    donor_score_sd = sd(score_mean),
    donor_score_median = median(score_mean),
    donor_score_q1 = unname(quantile(score_mean, 0.25)),
    donor_score_q3 = unname(quantile(score_mean, 0.75)),
    .groups = "drop"
  )

# ==============================================================================
# 5. Plot one point per donor
# ==============================================================================

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
}

make_donor_panel <- function(data, test_row, group_levels, colors,
                             panel_title, y_label, paired = FALSE) {
  data <- data %>%
    mutate(group = factor(group, levels = group_levels))
  
  n_by_group <- table(data$group)
  x_labels <- setNames(
    paste0(group_levels, "\n(n=", as.integer(n_by_group[group_levels]), " donors)"),
    group_levels
  )
  
  score_range <- range(data$score_mean, na.rm = TRUE)
  score_span <- diff(score_range)
  if (!is.finite(score_span) || score_span == 0) {
    score_span <- max(abs(score_range), 0.1)
  }
  bracket_y <- score_range[2] + 0.14 * score_span
  tick_y <- bracket_y - 0.035 * score_span
  text_y <- bracket_y + 0.055 * score_span
  
  test_label <- paste0(
    "P = ", format_p(test_row$p_value),
    "; Holm q = ", format_p(test_row$p_holm_four_tests)
  )
  
  p <- ggplot(data, aes(x = group, y = score_mean))
  
  if (paired) {
    p <- p +
      geom_line(
        aes(group = donor_id),
        color = "#777777", linewidth = 0.50, alpha = 0.75
      ) +
      geom_point(
        aes(fill = group),
        shape = 21, size = 3.0, stroke = 0.55, color = "black"
      ) +
      stat_summary(
        aes(group = group),
        fun = median, geom = "crossbar",
        width = 0.48, linewidth = 0.65, color = "black"
      )
  } else {
    p <- p +
      geom_boxplot(
        aes(fill = group),
        width = 0.52, outlier.shape = NA,
        alpha = 0.48, linewidth = 0.55
      ) +
      geom_point(
        aes(fill = group),
        position = position_jitter(width = 0.08, height = 0, seed = 20260831),
        shape = 21, size = 2.8, stroke = 0.55, color = "black"
      )
  }
  
  p +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
             linewidth = 0.45) +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = tick_y,
             linewidth = 0.45) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = tick_y,
             linewidth = 0.45) +
    annotate("text", x = 1.5, y = text_y, label = test_label, size = 3.15) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.28))) +
    labs(
      title = panel_title,
      x = NULL,
      y = y_label
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 10.5, margin = margin(r = 5)),
      axis.line = element_line(linewidth = 0.55, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      plot.margin = margin(7, 7, 5, 7)
    )
}

get_test_row <- function(panel_name, cohort_name) {
  out <- primary_tests %>%
    filter(panel == panel_name, cohort == cohort_name)
  if (nrow(out) != 1L) stop("Could not identify one test row for plotting.")
  out
}

p_g_alsaigh <- make_donor_panel(
  fibro_alsaigh,
  get_test_row("Fig1G", "Alsaigh_2022"),
  c("control", "plaque"), tissue_colors,
  "Alsaigh 2022 (paired)",
  "Donor-level mean\nfibroblast-healing score",
  paired = TRUE
)

p_g_ours <- make_donor_panel(
  fibro_ours,
  get_test_row("Fig1G", "Our_study"),
  c("ASYM", "SYM"), group_colors,
  "Our study (independent)",
  "Donor-level mean\nfibroblast-healing score",
  paired = FALSE
)

p_h_alsaigh <- make_donor_panel(
  mac_alsaigh,
  get_test_row("Fig1H", "Alsaigh_2022"),
  c("control", "plaque"), tissue_colors,
  "Alsaigh 2022 (paired)",
  "Donor-level mean\nTREM2+ LAM score",
  paired = TRUE
)

p_h_ours <- make_donor_panel(
  mac_ours,
  get_test_row("Fig1H", "Our_study"),
  c("ASYM", "SYM"), group_colors,
  "Our study (independent)",
  "Donor-level mean\nTREM2+ LAM score",
  paired = FALSE
)

fig1G <- p_g_alsaigh | p_g_ours
fig1H <- p_h_alsaigh | p_h_ours

fig1GH <- fig1G / fig1H +
  plot_annotation(
    title = "Donor-level analysis of reparative fibroblast and TREM2+ LAM programs",
    subtitle = paste0(
      "Each point represents one donor; Alsaigh samples are paired by donor. ",
      "P values are two-sided donor-level Wilcoxon tests."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5)
    )
  )

# ==============================================================================
# 6. Export figures, donor-level data and complete statistics
# ==============================================================================

ggsave(
  file.path(out_dir, "Fig1G_fibroblast_healing_donor_level.pdf"),
  fig1G, width = 8.8, height = 3.6, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)

ggsave(
  file.path(out_dir, "Fig1H_TREM2_LAM_donor_level.pdf"),
  fig1H, width = 8.8, height = 3.6, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)

ggsave(
  file.path(out_dir, "Fig1GH_donor_level_combined.pdf"),
  fig1GH, width = 8.8, height = 7.2, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)

ggsave(
  file.path(out_dir, "Fig1GH_donor_level_combined.png"),
  fig1GH, width = 8.8, height = 7.2, units = "in",
  dpi = 600, bg = "white"
)

write.csv(
  donor_scores,
  file.path(out_dir, "Fig1GH_donor_level_scores.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  descriptive_stats,
  file.path(out_dir, "Fig1GH_donor_level_descriptive_statistics.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  primary_tests,
  file.path(out_dir, "Fig1GH_primary_donor_mean_Wilcoxon_tests.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

write.csv(
  sensitivity_tests,
  file.path(out_dir, "Fig1GH_sensitivity_donor_median_Wilcoxon_tests.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

methods_text <- c(
  "Fig. 1G/H donor-level analysis",
  "",
  paste0(
    "Cell-level FibroHealing1 and TREM2_LAM2 module scores were summarized ",
    "as the mean score per donor within the prespecified target cell population."
  ),
  paste0(
    "The fibroblast-healing score was evaluated in Fibroblast/SMC cells and ",
    "the TREM2+ LAM score in macrophages, using the same integrated-cluster ",
    "definitions as the original figure."
  ),
  paste0(
    "Matched control and plaque samples from the three Alsaigh donors were ",
    "compared using a two-sided paired Wilcoxon signed-rank test."
  ),
  paste0(
    "ASYM and SYM samples from our cohort (six independent donors per group) ",
    "were compared using a two-sided Wilcoxon rank-sum test."
  ),
  paste0(
    "Rank-biserial correlations and bootstrap 95% confidence intervals for ",
    "the donor-level mean difference were reported."
  ),
  paste0(
    "Holm-adjusted P values across the four displayed comparisons were ",
    "reported alongside unadjusted exact P values."
  ),
  paste0(
    "All donor-level observations are displayed; cells were not treated as ",
    "independent biological replicates."
  )
)

writeLines(
  methods_text,
  con = file.path(out_dir, "Fig1GH_methods_text.txt"),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "Fig1GH_R_sessionInfo.txt")
)

message("Completed donor-level Fig. 1G/H analysis.")
message("Results saved to: ", out_dir)
print(primary_tests)








































