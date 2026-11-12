# Project root: the folder containing SMC/, macro/, Fig3/, Fig5_CellChat/, etc.
# Every path in this script is built from this one setting, so the code no longer
# depends on a single hard-coded drive letter. Either start R/RStudio in the
# project root, or edit the one line below to point at it directly.
# Stored as an R option so it survives the rm(list = ls()) calls between sections.
options(cea.project_root = normalizePath(".", mustWork = FALSE))

# Supplemental Figure 4 
# ==========================================================================================
# 共用准备代码：原路径定位与 00_common.R 完整原文
# ==========================================================================================

source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE))
  source_path <- rstudioapi::getActiveDocumentContext()$path
if (!nzchar(source_path)) stop("Use source('.../01_donor_panels.R', encoding='UTF-8').")
FIG4_HOME <- dirname(source_path)

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


local({
  source_path <- tryCatch(normalizePath(sys.frame(1)$ofile,mustWork=TRUE),error=function(e) "")
  if (!nzchar(source_path) && requireNamespace("rstudioapi",quietly=TRUE)) source_path <- rstudioapi::getActiveDocumentContext()$path
  if (!nzchar(source_path)) stop("Run with source('.../03_trajectory_validation.R',encoding='UTF-8').")
  FIG4_HOME <- dirname(source_path)
  source(file.path(FIG4_HOME,"00_common.R"),local=TRUE,encoding="UTF-8")
  need(c("monocle3","SingleCellExperiment","SummarizedExperiment","igraph"))
  out <- new_run("03_trajectory"); set.seed(cfg$seed)
  dat <- read_source(); scope_table(dat$md,out)
  use <- dat$md$population%in%cfg$populations & dat$md$group%in%cfg$groups
  md <- dat$md[use,c("cell_id","donor","sample_label","group","population")]
  rownames(md) <- md$cell_id
  counts <- dat$counts[,md$cell_id,drop=FALSE]; rm(dat); gc()
  keep <- Matrix::rowSums(counts>0)>=10L
  csv(data.frame(gene=rownames(counts),used_in_trajectory_PCA=keep),"trajectory_feature_filter",out)
  genes <- rownames(counts)[keep]
  base <- monocle3::new_cell_data_set(counts,cell_metadata=md,
                                      gene_metadata=data.frame(gene_short_name=rownames(counts),row.names=rownames(counts)))
  base <- monocle3::preprocess_cds(base,num_dim=50,use_genes=genes,norm_method="log",method="PCA")
  if(!is.null(cfg$alignment_group)) {
    if(!cfg$alignment_group%in%names(md)) stop("Supply a validated technical-batch field; do not substitute clinical group.")
    base <- monocle3::align_cds(base,alignment_group=cfg$alignment_group)
  }
  # Choose a monocyte-supported root region in each partition; not one FCN1 outlier.
  # No monocyte-supported root -> that partition remains unrooted, not forcibly linked.
  choose_roots <- function(cds) {
    part <- as.character(monocle3::partitions(cds)[colnames(cds)])
    near <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)
    near <- near[colnames(cds),1]
    vertex_names <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name
    nodes <- if(all(as.character(near)%in%vertex_names)) as.character(near) else vertex_names[as.integer(near)]
    d <- data.frame(md,partition=part,node=nodes)
    totals <- d %>% count(partition,node,name="n_all")
    mono <- d %>% filter(population=="Monocytes")
    bydonor <- mono %>% count(partition,donor,node,name="n_node") %>%
      group_by(partition,donor) %>% mutate(weight=n_node/sum(n_node)) %>% ungroup()
    # Missing donor-node pairs contribute zero to sums, so high weights cannot come
    # simply from a donor contributing more cells. This is a root heuristic, not a test.
    support <- bydonor %>% group_by(partition,node) %>%
      summarise(donor_balanced_support=sum(weight),n_mono=sum(n_node),
                n_mono_donors=n_distinct(donor),.groups="drop") %>%
      left_join(totals,by=c("partition","node")) %>% mutate(mono_fraction=n_mono/n_all,
                                                            eligible=n_mono>=20L & n_mono_donors>=3L & mono_fraction>=.5) %>%
      arrange(partition,desc(donor_balanced_support),desc(mono_fraction),node)
    selected <- support %>% filter(eligible) %>% group_by(partition) %>%
      mutate(root_rank=row_number()) %>% filter(root_rank<=2L) %>% ungroup()
    list(selected=selected,support=support,partition=part)
  }
  normalized_pt <- function(pt,part) {
    ans <- rep(NA_real_,length(pt))
    for(p in unique(part)) {
      ix <- which(part==p & is.finite(pt)); if(length(ix)<2L) next
      r <- range(pt[ix]); if(diff(r)>0) ans[ix] <- (pt[ix]-r[1])/diff(r)
    }
    ans
  }
  all_cells <- list(); all_roots <- list(); all_support <- list()
  for(i in seq_along(cfg$trajectory_seeds)) {
    seed <- cfg$trajectory_seeds[i]; set.seed(seed)
    message("Trajectory seed ",seed," (",i,"/",length(cfg$trajectory_seeds),")")
    cds <- monocle3::reduce_dimension(base,reduction_method="UMAP",preprocess_method=if(is.null(cfg$alignment_group))"PCA" else "Aligned",
                                      umap.n_neighbors=30L,umap.min_dist=.1,umap.fast_sgd=FALSE,cores=1)
    cds <- monocle3::cluster_cells(cds,reduction_method="UMAP",random_seed=seed)
    cds <- monocle3::learn_graph(cds,use_partition=TRUE,close_loop=FALSE)
    root <- choose_roots(cds)
    all_roots[[i]] <- mutate(root$selected,seed=seed)
    all_support[[i]] <- mutate(root$support,seed=seed)
    primary <- root$selected %>% filter(root_rank==1L)
    # A lack of eligible roots is a result. Do not lower the threshold for desired ordering.
    if(!nrow(primary)) {
      csv(root$support,paste0("root_support_no_eligible_root_seed_",seed),out)
      warning("No monocyte-supported root for seed ",seed,"; pseudotime left undefined.")
      pt <- rep(Inf,ncol(cds)); alt_pt <- pt
    } else {
      cds <- monocle3::order_cells(cds,root_pr_nodes=primary$node)
      pt <- monocle3::pseudotime(cds)[colnames(cds)]
      alternative <- primary
      for(j in seq_len(nrow(alternative))) {
        other <- root$selected %>% filter(partition==alternative$partition[j],root_rank==2L)
        if(nrow(other)) alternative$node[j] <- other$node[1]
      }
      alt <- monocle3::order_cells(cds,root_pr_nodes=alternative$node)
      alt_pt <- monocle3::pseudotime(alt)[colnames(cds)]
    }
    d <- md; d$seed <- seed; d$partition <- root$partition
    d$pseudotime <- as.numeric(pt); d$pt_norm <- normalized_pt(pt,root$partition)
    d$alternative_root_pt_norm <- normalized_pt(alt_pt,root$partition)
    root_pairs <- root$selected %>% group_by(partition) %>% summarise(alternative_available=n()>=2,.groups="drop")
    d <- left_join(d,root_pairs,by="partition")
    all_cells[[i]] <- d
    # Supplementary algorithm QC only; no t-SNE overview is generated.
    pgraph <- monocle3::plot_cells(cds,color_cells_by="population",label_cell_groups=FALSE,
                                   label_leaves=FALSE,label_branch_points=FALSE,label_roots=FALSE,cell_size=.15) +
      labs(title=paste("Exploratory UMAP graph; seed",seed),subtitle="Partitions not forcibly connected; population labels retained") + theme_fig
    save_fig(pgraph,paste0("QC_UMAP_graph_seed_",seed),out,7,5)
    saveRDS(cds,file.path(out,"objects",paste0("cds_seed_",seed,".rds")))
    rm(cds); gc()
  }
  cells <- bind_rows(all_cells); csv(cells,"all_cell_pseudotime_new_seeds",out)
  csv(bind_rows(all_roots),"root_choices",out); csv(bind_rows(all_support),"root_support_all_candidates",out)
  cover <- cells %>% group_by(seed,partition,population,donor,sample_label,group) %>%
    summarise(n_cells=n(),n_finite=sum(is.finite(pt_norm)),finite_fraction=mean(is.finite(pt_norm)),.groups="drop")
  csv(cover,"partition_donor_population_coverage",out)
  primary <- all_cells[[1]]
  # Median pseudotime is compared descriptively within each partition; no cell-level P.
  med <- primary %>% filter(is.finite(pt_norm)) %>%
    group_by(partition,donor,sample_label,group,population) %>%
    summarise(n_cells=n(),median_pt=median(pt_norm),.groups="drop") %>%
    mutate(retained=n_cells>=cfg$min_cells)
  csv(med,"donor_median_pseudotime",out)
  pmed <- ggplot(filter(med,retained),aes(factor(population,levels=cfg$populations),median_pt)) +
    geom_line(data=med %>% filter(retained) %>% group_by(partition,donor) %>% filter(n()>=2L) %>% ungroup(),
              aes(group=donor),color="grey75",linewidth=.3) +
    geom_boxplot(width=.5,outlier.shape=NA,fill="grey90") +
    geom_point(aes(fill=group),shape=21,size=2) + facet_wrap(~partition) +
    scale_fill_manual(values=group_colors) + labs(x=NULL,y="Median within-partition pseudotime",
                                                  title="Exploratory donor-level pseudotime positions",subtitle="Root = monocyte-supported region; not proof of Monocytes -> Macrophages1 -> Macrophages2") +
    theme_fig + theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="bottom")
  if(nrow(filter(med,retained))) save_fig(pmed,"Fig4S_candidate_donor_pseudotime",out,8,5)
  # Root sensitivity: report NA when no genuinely different alternative was available.
  safe_rho <- function(x,y) if(sum(is.finite(x)&is.finite(y))>=10L && sd(x,na.rm=TRUE)>0 && sd(y,na.rm=TRUE)>0)
    cor(x,y,method="spearman",use="complete.obs") else NA_real_
  root_sens <- primary %>% group_by(partition,donor,sample_label) %>%
    summarise(alternative_available=any(alternative_available%in%TRUE),
              rho=if(any(alternative_available%in%TRUE))safe_rho(pt_norm,alternative_root_pt_norm) else NA_real_,.groups="drop")
  csv(root_sens,"alternative_root_donor_consistency",out)
  # Partition IDs can permute between runs. Match by overlap, show mapping and overlap.
  second <- all_cells[[2]]
  pair <- merge(primary[,c("cell_id","donor","sample_label","partition","pt_norm")],
                second[,c("cell_id","partition","pt_norm")],by="cell_id",suffixes=c("_primary","_second"))
  mapping <- pair %>% count(partition_primary,partition_second,name="overlap") %>%
    group_by(partition_primary) %>% mutate(primary_overlap_fraction=overlap/sum(overlap)) %>%
    arrange(desc(overlap),partition_second,.by_group=TRUE) %>% mutate(best_match=row_number()==1L) %>% ungroup()
  csv(mapping,"seed_partition_overlap_mapping",out)
  best <- filter(mapping,best_match) %>% select(partition_primary,partition_second)
  paired_seeds <- inner_join(pair,best,by=c("partition_primary","partition_second"))
  seed_sens <- paired_seeds %>% group_by(partition_primary,partition_second,donor,sample_label) %>%
    summarise(n_common_cells=n(),rho=safe_rho(pt_norm_primary,pt_norm_second),.groups="drop")
  csv(seed_sens,"alternative_seed_donor_consistency",out)
  # Compare old and new only on overlapping cells, within each NEW partition.
  if(file.exists(cfg$old_trajectory)) {
    old <- readRDS(cfg$old_trajectory)@meta.data
    old <- data.frame(cell_id=rownames(old),old_pseudotime=old$pseudotime)
    old$old_pseudotime[!is.finite(old$old_pseudotime)] <- NA_real_
    legacy <- left_join(primary,old,by="cell_id") %>% group_by(partition,donor,sample_label) %>%
      summarise(rho=safe_rho(pt_norm,old_pseudotime),.groups="drop")
    csv(legacy,"old_tSNE_vs_new_UMAP_donor_consistency",out)
  }
  # Donor-balanced marker heatmap. Fixed-width bins; no smoothing across missing bins.
  markers <- list(Inflammatory=c("IL1B","TNF","NLRP3","S100A8","S100A9"),
                  `Lipid/LAM-associated`=c("TREM2","APOE","SPP1","LPL","LIPA"),
                  `Efferocytosis-associated`=c("MERTK","CD163","C1QA"),
                  `Antigen presentation`=c("CIITA","HLA-DRA","CD74"))
  marker_table <- bind_rows(lapply(names(markers),function(nm)data.frame(program=nm,gene=markers[[nm]])))
  marker_table$present <- marker_table$gene%in%rownames(counts); csv(marker_table,"trajectory_marker_QC",out)
  present <- intersect(marker_table$gene,rownames(counts))
  sizes <- Matrix::colSums(counts)
  lognorm <- as.matrix(counts[present,,drop=FALSE] %*% Matrix::Diagonal(x=1e4/sizes))
  lognorm <- log1p(lognorm)
  colnames(lognorm) <- colnames(counts)
  finite <- primary %>% filter(is.finite(pt_norm)) %>%
    mutate(bin=pmin(cfg$bins,1L+floor(pt_norm*cfg$bins)))
  if(!nrow(finite)) stop("No finite pseudotime. QC files saved; no heatmap fabricated.")
  groups <- split(seq_len(nrow(finite)),paste(finite$partition,finite$donor,finite$bin,sep="___"))
  bin_data <- bind_rows(lapply(groups,function(ix) {
    m <- finite[ix,,drop=FALSE]
    data.frame(partition=m$partition[1],donor=m$donor[1],sample_label=m$sample_label[1],group=m$group[1],
               bin=m$bin[1],n_cells=nrow(m),gene=present,
               mean_log1p_CP10K=rowMeans(lognorm[,m$cell_id,drop=FALSE]))
  }))
  bin_data$retained <- bin_data$n_cells>=cfg$bin_cells
  csv(bin_data,"donor_bin_marker_means",out)
  cov <- bin_data %>% distinct(partition,donor,bin,n_cells,retained) %>%
    group_by(partition,bin) %>% summarise(n_donors=sum(retained),n_cells=sum(n_cells[retained]),.groups="drop")
  cov <- tidyr::complete(cov,partition,bin=seq_len(cfg$bins),fill=list(n_donors=0L,n_cells=0L))
  cov$retained <- cov$n_donors>=cfg$bin_donors
  csv(cov,"donor_bin_coverage_including_missing_bins",out)
  heat <- bin_data %>% filter(retained) %>% group_by(partition,bin,gene) %>%
    summarise(mean_log_expression=mean(mean_log1p_CP10K),.groups="drop") %>%
    left_join(cov %>% select(partition,bin,retained),by=c("partition","bin")) %>%
    mutate(mean_log_expression=ifelse(retained,mean_log_expression,NA_real_)) %>%
    select(-retained) %>% tidyr::complete(partition,bin=seq_len(cfg$bins),gene=present) %>%
    group_by(partition,gene) %>% mutate(gene_z=(mean_log_expression-mean(mean_log_expression,na.rm=TRUE))/sd(mean_log_expression,na.rm=TRUE)) %>%
    ungroup() %>% left_join(marker_table,by="gene")
  heat$gene_z[!is.finite(heat$gene_z)] <- NA_real_
  heat$program <- factor(heat$program,levels=names(markers)); heat$gene <- factor(heat$gene,levels=rev(present))
  csv(heat,"donor_balanced_heatmap_source",out)
  limit <- max(abs(heat$gene_z),na.rm=TRUE)
  if(!is.finite(limit)) stop("Insufficient supported bins for standardized heatmap. Coverage table saved.")
  ph <- ggplot(heat,aes((bin-.5)/cfg$bins,gene,fill=gene_z)) + geom_tile(width=1/cfg$bins) +
    facet_grid(program~partition,scales="free_y",space="free_y",labeller=program_labeller) +
    scale_fill_gradient2(low="#326A9C",mid="white",high="#B62932",limits=c(-limit,limit),na.value="grey85",name="Within-gene z") +
    labs(x="Within-partition pseudotime (separate axes, not aligned clocks)",y=NULL,
         title="Donor-balanced marker expression along the exploratory trajectory",
         subtitle=paste0("At least ",cfg$bin_cells," cells/donor/bin and ",cfg$bin_donors," donors/bin; grey = unsupported; no interpolation")) +
    theme_fig + theme(axis.text.y=element_text(face="italic"),strip.text.y=element_text(size=8,angle=0))
  save_fig(ph,"Fig4S_candidate_pseudotime_marker_heatmap",out,9,6.5)
  csv(data.frame(partition=unique(primary$partition)),"observed_primary_partitions",out)
  writeLines(c("This is a sensitivity analysis, not an automatic replacement for the final main figure.",
               "Review partition composition, finite-cell coverage, branch structure and root/seed stability together.",
               "Pseudotime axes from different partitions are not directly comparable.",
               "Monocyte rooting partly imposes low monocyte pseudotime; it is not independent differentiation evidence.",
               "A branched graph does not imply a sequential Macrophages1-to-Macrophages2 transition.",
               "Shared genes used for topology and display make marker trends descriptive, not independent validation.",
               "No technical batch was assumed or regressed out. Assess donor/batch composition before interpretation.",
               "Do not automatically substitute this pseudotime into Fig5; inspect concordance first."),file.path(out,"REVIEW_BEFORE_MAIN_FIGURE.txt"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo_end.txt"))
  writeLines("Trajectory candidates generated; scientific review remains required.",file.path(out,"COMPLETE.txt"))
  message("Step 3 complete: ",out)
})


--------------------------------------------------------------------------------
# ==========================================================================================
# Figure A : Donor/patient median pseudotime

local({
  path <- normalizePath(sys.frame(1)$ofile,mustWork=TRUE); FIG4_HOME <- dirname(path)
  source(file.path(FIG4_HOME,"00_common.R"),local=TRUE,encoding="UTF-8")
  out <- new_run("04_review"); set.seed(cfg$seed)
  tr <- latest_run("03_trajectory"); er <- latest_run("02_enrichment")
  table_in <- function(dir,nm) read.csv(file.path(dir,"tables",paste0(nm,".csv")),check.names=FALSE)
  writeLines(c(trajectory=tr,enrichment=er),file.path(out,"UPSTREAM_RUNS.txt"))
  cover <- table_in(tr,"partition_donor_population_coverage")
  totals <- cover %>% filter(seed==cfg$trajectory_seeds[1]) %>% group_by(partition) %>%
    summarise(n_cells=sum(n_cells),n_finite=sum(n_finite),.groups="drop") %>%
    mutate(fraction=n_cells/sum(n_cells)) %>% arrange(desc(n_cells))
  major <- totals$partition[1]
  csv(totals,"partition_overview",out)
  # Selection is by cell count only. All partitions remain in source tables and full heatmap.
  med <- table_in(tr,"donor_median_pseudotime") %>% filter(retained,partition==major)
  med$population <- factor(med$population,levels=cfg$populations)
  pmed <- ggplot(med,aes(population,median_pt)) +
    geom_line(data=med %>% group_by(donor) %>% filter(n()>=2) %>% ungroup(),aes(group=donor),color="grey75",linewidth=.3) +
    geom_boxplot(width=.5,fill="grey92",outlier.shape=NA) + geom_point(aes(fill=group),shape=21,size=2) +
    scale_fill_manual(values=group_colors,name="Clinical group") +
    labs(x=NULL,y="Donor median pseudotime",title=paste("Largest partition:",major),
         subtitle=sprintf("%.1f%% of focused plaque cells; monocyte-supported root",100*totals$fraction[1])) +
    theme_fig + theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="bottom")
  save_fig(pmed,"Fig4S_donor_PT",out,6.5,4.5)

  rs <- table_in(tr,"alternative_root_donor_consistency") %>% filter(partition==major) %>%
    transmute(donor,sample_label,rho,check="Alternative root")
  ss <- table_in(tr,"alternative_seed_donor_consistency") %>% filter(partition_primary==major) %>%
    transmute(donor,sample_label,rho,check="Alternative seed")
  os <- table_in(tr,"old_tSNE_vs_new_UMAP_donor_consistency") %>% filter(partition==major) %>%
    transmute(donor,sample_label,rho,check="Old vs new trajectory")
  consistency <- bind_rows(rs,ss,os)
  consistency$check <- factor(consistency$check,levels=c("Alternative root","Alternative seed","Old vs new trajectory"))
  csv(consistency,"largest_partition_consistency",out)
  summary <- consistency %>% group_by(check) %>% summarise(n=sum(is.finite(rho)),
                                                           minimum=min(rho,na.rm=TRUE),median=median(rho,na.rm=TRUE),maximum=max(rho,na.rm=TRUE),.groups="drop")
  csv(summary,"largest_partition_consistency_summary",out)
  pr <- ggplot(consistency,aes(check,rho)) + geom_hline(yintercept=0,linetype=2,color="grey65") +
    geom_boxplot(width=.5,outlier.shape=NA,fill="grey92") +
    geom_point(position=position_jitter(width=.06,height=0,seed=cfg$seed),size=2) +
    scale_y_continuous(limits=c(-1,1)) + scale_x_discrete(labels=c("Alternative root"="Alternative\nroot",
                                                                   "Alternative seed"="Alternative\nseed","Old vs new trajectory"="Old vs new\ntrajectory")) +
    labs(x=NULL,y="Within-donor Spearman rho",title="Ordering sensitivity",subtitle="One value per donor; largest partition only") + theme_fig
  
  # ==========================================================================================
  # Figure C : Functional-marker dynamics along within-partition pseudotime
  # 对应原 heat_plot / ph，保留原 marker 顺序及原代码中的所有计算与样式。
  # ==========================================================================================
  
  heat <- table_in(tr,"donor_balanced_heatmap_source")
  heat$program <- factor(heat$program,levels=primary_programs)
  marker_order <- c("IL1B","TNF","NLRP3","S100A8","S100A9","TREM2","APOE","SPP1","LPL","LIPA","MERTK","CD163","C1QA","CIITA","HLA-DRA","CD74")
  heat$gene <- factor(heat$gene,levels=rev(marker_order))
  limit <- max(abs(heat$gene_z),na.rm=TRUE)
  heat_plot <- function(h,all=FALSE) {
    p <- ggplot(h,aes((bin-.5)/cfg$bins,gene,fill=gene_z)) + geom_tile(width=1/cfg$bins) +
      scale_fill_gradient2(low="#326A9C",mid="white",high="#B62932",limits=c(-limit,limit),na.value="grey85",name="Gene z") +
      labs(x="Within-partition pseudotime",y=NULL,title="Donor-balanced marker expression",
           subtitle="Grey bins: insufficient support; no interpolation") + theme_fig +
      theme(axis.text.y=element_text(face="italic",size=8),strip.text.y=element_text(size=8,angle=0))
    if(all) p+facet_grid(program~partition,scales="free_y",space="free_y",labeller=program_labeller) else
      p+facet_grid(program~.,scales="free_y",space="free_y",labeller=program_labeller)
  }
  ph <- heat_plot(filter(heat,partition==major))
  save_fig(heat_plot(heat,TRUE),"Fig4S_all_partition_heatmap",out,10,7)

  save_fig(ph,"Fig4S_marker_PT",out,6.5,6.5)
  
  # ==========================================================================================
  # Figure D 代码 | Donor/patient coverage per pseudotime bin
  pcovdata <- table_in(tr,"donor_bin_coverage_including_missing_bins") %>% filter(partition==major)
  pcov <- ggplot(pcovdata,aes((bin-.5)/cfg$bins,n_donors,fill=retained)) +
    geom_col(width=.9/cfg$bins) + geom_hline(yintercept=cfg$bin_donors,linetype=2) +
    scale_fill_manual(values=c(`FALSE`="grey80",`TRUE`="#3A6FA0"),guide="none") +
    labs(x="Within-partition pseudotime",y="Donors passing cell threshold",title="Coverage along pseudotime",
         subtitle=paste0("At least ",cfg$bin_cells," cells/donor/bin; dashed line = ",cfg$bin_donors," donors")) + theme_fig
  assembled <- ((pmed+pr)/(ph+pcov)) + plot_layout(heights=c(1,1.4)) + plot_annotation(tag_levels="A")
  save_fig(assembled,"Fig4S_trajectory_review",out,11,9)
})

# ==========================================================================================
# Figure E: 12-gene sensitivity LAM program
#  Fig4S_legacy_LAM_gene_set_sensitivity。

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

  pLegacy <- ggplot(filter(scores,program=="Legacy LAM (sensitivity)"),aes(population,score,group=donor)) +
    geom_line(color="grey70") + geom_point(aes(fill=group),shape=21,size=2.5) +
    scale_fill_manual(values=group_colors) + labs(x=NULL,y="Mean gene z",title="Original 12-gene LAM set: score-method sensitivity") +
    theme_fig + theme(legend.position="bottom")
  save_fig(pLegacy,"Fig4S_legacy_LAM_gene_set_sensitivity",out,6,4)
  
})

# ==========================================================================================
# Figure F :Gene-level membership across Our and three published LAM signatures
#  Supp Figure S5 脚本的 p_B（Gene membership heatmap）。

suppressMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

out_dir <- paste0(getOption("cea.project_root"), "/SuppFigS5/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. 定义 4 个 LAM signature ────────────────────────────────────────────
sig_lists <- list(
  "Our (21 genes)" = c(
    "TREM2","CD9","SPP1","APOE","GPNMB","FABP5",
    "LIPA","LPL","CTSD","CTSB","CTSL","CTSZ",
    "ACP5","NPC2","FBP1","CHIT1","CHI3L1","CD63",
    "CD68","CD163","MERTK"),
  
  "Cochain 2018 (9 genes)" = c(
    "TREM2","CD9","SPP1","GPNMB","FABP5","MMP12","CTSD","CD63","LGALS3"),
  
  "Patterson 2023 (14 genes)" = c(
    "TREM2","APOE","LPL","CD9","SPP1","GPNMB","FABP4","FABP5",
    "CTSD","CTSB","CD63","CD36","PLIN2","LIPA"),
  
  "Willemsen 2024 (10 genes)" = c(
    "TREM2","APOE","SPP1","GPNMB","FABP5","LIPA","MMP9","CD9","ACP5","CSTB")
)

cat("Signature sizes:\n")
for (n in names(sig_lists)) {
  cat("  ", n, ":", length(sig_lists[[n]]), "genes\n")
}

# 简化标签 (供 B/C 用)
sig_short <- c(
  "Our (21 genes)"            = "Our",
  "Cochain 2018 (9 genes)"    = "Cochain 2018",
  "Patterson 2023 (14 genes)" = "Patterson 2023",
  "Willemsen 2024 (10 genes)" = "Willemsen 2024"
)

# ─────────────────────────────────────────────────────────────────────────
# Gene membership heatmap
# ─────────────────────────────────────────────────────────────────────────
all_genes <- unique(unlist(sig_lists))
cat("\n=== Total unique genes across 4 signatures:", length(all_genes), "===\n")

gene_matrix <- sapply(sig_lists, function(s) all_genes %in% s) * 1
rownames(gene_matrix) <- all_genes
colnames(gene_matrix) <- sig_short

n_studies   <- rowSums(gene_matrix)
gene_matrix <- gene_matrix[order(-n_studies), ]

heat_df <- as.data.frame(gene_matrix) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Signature", values_to = "Present") %>%
  mutate(Present = factor(ifelse(Present == 1, "Yes", "No"),
                          levels = c("No", "Yes")))

heat_df$Gene      <- factor(heat_df$Gene, levels = rev(rownames(gene_matrix)))
heat_df$Signature <- factor(heat_df$Signature, levels = sig_short)

gene_count_df <- data.frame(
  Gene  = factor(rownames(gene_matrix), levels = rev(rownames(gene_matrix))),
  count = n_studies[rownames(gene_matrix)]
)

p_B <- ggplot(heat_df, aes(x = Signature, y = Gene, fill = Present)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Yes" = "#8B0000", "No" = "#F0F0F0"),
                    name = NULL) +
  geom_text(data = gene_count_df,
            aes(x = 4.6, y = Gene, label = paste0(count, "/4")),
            inherit.aes = FALSE,
            size = 2.8, fontface = "bold", color = "grey30") +
  coord_cartesian(xlim = c(0.5, 5.3), clip = "off") +
  labs(title = "B. Gene-level membership across 4 signatures",
       subtitle = "Sorted by number of studies (right column)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 10,
                                     face = "bold", color = "black"),
        axis.text.y   = element_text(size = 8.5, face = "italic", color = "black"),
        panel.grid    = element_blank(),
        legend.position = "bottom",
        plot.margin   = margin(10, 60, 10, 10))


# ==========================================================================================
# Figure G : GO biological processes of SPP1+ foamy LAMs
# Fig4G_GO_foamy。

suppressMessages({
  library(Seurat); library(tidyverse); library(ggplot2)
  library(patchwork); library(ggrepel); library(ggrastr)
  library(clusterProfiler); library(org.Hs.eg.db)
})

macro <- readRDS(paste0(getOption("cea.project_root"), "/macro/Fig4E/macro_with_pseudotime_3grp.rds"))  # ★ 已剔除 Macro3
out_dir <- paste0(getOption("cea.project_root"), "/macro/Fig4G/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. 仅 foamy LAM + mature LAM ────────────────────────────────────
macro$lam_compare <- dplyr::case_when(
  macro$celltype == "Macrophages1" ~ "Foamy_LAM",
  macro$celltype == "Macrophages2" ~ "Mature_LAM",
  TRUE ~ NA_character_
)
macro_lam <- subset(macro, subset = !is.na(lam_compare))
Idents(macro_lam) <- macro_lam$lam_compare
cat("LAM cells for comparison:\n"); print(table(Idents(macro_lam)))

# ── 2. DEG: foamy vs mature ─────────────────────────────────────────
deg <- FindMarkers(macro_lam, ident.1 = "Foamy_LAM", ident.2 = "Mature_LAM",
                   logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox") %>%
  rownames_to_column("gene") %>%
  mutate(direction = case_when(
    avg_log2FC >  0.5 & p_val_adj < 0.05 ~ "Foamy_up",
    avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Mature_up",
    TRUE ~ "ns"),
    nlog10_padj = -log10(p_val_adj + 1e-300))
cat("DEG counts:\n"); print(table(deg$direction))
write.csv(deg, paste0(out_dir, "DEG_foamy_vs_mature.csv"), row.names = FALSE)

# GO 富集
# ============================================================================
run_go <- function(gene_list, ont = "BP", top_n = 8) {
  gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db) %>% pull(ENTREZID)
  go <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, ont = ont,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                 readable = TRUE)
  if (is.null(go) || nrow(go@result) == 0) return(NULL)
  go@result %>% arrange(p.adjust) %>% head(top_n) %>%
    mutate(Description = str_wrap(Description, width = 38),
           nlog10_padj = -log10(p.adjust),
           GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                                  function(x) as.numeric(x[1])/as.numeric(x[2]))) %>%
    arrange(desc(nlog10_padj))
}

foamy_genes  <- deg %>% filter(direction == "Foamy_up")  %>% pull(gene)
mature_genes <- deg %>% filter(direction == "Mature_up") %>% pull(gene)
cat("\nFoamy up:", length(foamy_genes), " Mature up:", length(mature_genes), "\n")

go_foamy  <- run_go(foamy_genes,  ont = "BP", top_n = 8)

# Fig4G/Figure 4G — Foamy LAM vs Mature LAM DEG + GO.R 
if (!is.null(go_foamy))  write.csv(go_foamy,  paste0(out_dir,"GO_foamy_BP.csv"),  row.names = FALSE)

# FigG: Foamy LAM vs Mature LAM DEG + GO.R 
# ============================================================================
# G2 — 背靠背条形图 (diverging bar)
# ============================================================================
n_top <- 9

make_bar <- function(go_df, fill_col, title) {
  d <- go_df %>% slice_min(p.adjust, n = n_top) %>%
    mutate(Description = str_replace_all(Description, "\n", " "),
           Description = str_wrap(Description, width = 32),   # ★ 长名折行, 控制每行字数
           nlog10 = -log10(p.adjust)) %>%
    arrange(nlog10) %>%
    mutate(Description = factor(Description, levels = Description))
  
  ggplot(d, aes(x = nlog10, y = Description)) +
    geom_col(fill = fill_col, width = 0.78) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = "-log10(pvalue)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title  = element_text(size = 11, hjust = 0),
          axis.text.y = element_text(color = "black", size = 8, lineheight = 0.85),  # 行距紧凑
          axis.text.x = element_text(color = "black"),
          axis.line   = element_line(linewidth = 0.4),
          axis.ticks  = element_line(linewidth = 0.4),
          plot.margin = margin(8, 12, 8, 8))
}

p_left  <- make_bar(go_foamy,  "#E59085", "SPP1+ foamy LAM")

#
ggsave(paste0(out_dir, "Fig4G_GO_foamy.pdf"),  p_left,  width = 7, height = 4, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig4G_GO_foamy.png"),  p_left,  width = 7, height = 4, dpi = 600)

# ==========================================================================================
# Figure H : GO biological processes of TREM2+ mature LAMs
#  Fig4G_GO_mature；使用 G 中原有的 mature_genes、run_go、make_bar。
# ==========================================================================================

# Fig4G/Figure 4G — Foamy LAM vs Mature LAM DEG + GO.R 
go_mature <- run_go(mature_genes, ont = "BP", top_n = 8)

#  Fig4G/Figure 4G — Foamy LAM vs Mature LAM DEG + GO.R 
if (!is.null(go_mature)) write.csv(go_mature, paste0(out_dir,"GO_mature_BP.csv"), row.names = FALSE)
#  Fig4G/Figure 4G — Foamy LAM vs Mature LAM DEG + GO.R
p_right <- make_bar(go_mature, "#8B0000", "TREM2+ mature LAM")

# ★ mature 通路名长, 画布加宽 (width 6→8)

ggsave(paste0(out_dir, "Fig4G_GO_mature.pdf"), p_right, width = 8, height = 4, device = cairo_pdf)
ggsave(paste0(out_dir, "Fig4G_GO_mature.png"), p_right, width = 8, height = 4, dpi = 600)

print(p_left); print(p_right)
