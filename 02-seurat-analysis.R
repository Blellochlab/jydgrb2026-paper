# Seurat Analysis
# Author: Deniz Goekbuget

# Setup
rm(list = ls())
library(Seurat)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)

# Define relative paths for portability
data_dir <- "./data"
plot_dir <- "./plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir)

# Load data
infiles <- list.files(path = data_dir, pattern = "-nodups-seurat.rds", recursive = TRUE, full.names = TRUE)

# -----------------------------------------------------------------------------
# 3.1 Global analysis (Extended Figure 5, 6)
# -----------------------------------------------------------------------------
dat <- readRDS(infiles[1])

# Define condition
dat$condition <- factor(ifelse(grepl("PD-L1", dat$sample), "dKO", "KO"))
dat$condition <- relevel(dat$condition, ref = "KO")
dat$replicate <- ifelse(grepl("#1", dat$sample), "N1", "N2")

# QC filter
high_mito <- quantile(dat[["percent.mt"]]$percent.mt, probs = 0.98)
low_ft <- quantile(dat[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02)
high_ft <- quantile(dat[["nFeature_RNA"]]$nFeature_RNA, probs = 0.98)

# Remove unwanted cells (~2000 cells remaining)
dat <- subset(dat, subset = percent.mt < high_mito & nFeature_RNA > low_ft & nFeature_RNA < high_ft) 

# Process
all.genes <- rownames(dat)
dat <- NormalizeData(dat, verbose=F) %>% 
  FindVariableFeatures(verbose=F) %>% 
  ScaleData(verbose=F, features=all.genes) %>% 
  RunPCA(verbose=F)

dat <- RunUMAP(dat, dims=1:30, reduction='pca')
dat <- FindNeighbors(dat, reduction = "pca", dims = 1:30)
dat <- FindClusters(dat)

# UMAP plots
set.seed(12)
cols <- DiscretePalette(2, shuffle = FALSE)
p1 <- DimPlot(dat, reduction = "umap", shuffle = TRUE, group.by = c("condition"), cols=cols)
ggsave(file.path(plot_dir, "global_umap_genotypes.pdf"), plot = p1, width=90, height=75, units="mm")

cols <- DiscretePalette(21, shuffle = FALSE)
p2 <- DimPlot(dat, reduction = "umap", shuffle = TRUE, group.by = c("seurat_clusters"), cols=cols)
ggsave(file.path(plot_dir, "global_umap_seuratclusters.pdf"), plot = p2, width=90, height=75, units="mm")

p3 <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "seurat_clusters", pt.size = FALSE)
ggsave(file.path(plot_dir, "global_violin_qc_scores_by_cluster.pdf"), plot = p3, width=300, height=75, units="mm")

# Global Il12b Violin
p_il12 <- VlnPlot(dat, features = "Il12b", pt.size = 0)
ggsave(file.path(plot_dir, "global_violin_Il12b.pdf"), plot = p_il12, width=150, height=100, units="mm")

# Fraction by condition (Global)
meta.data <- dat[[]]
counts <- group_by(meta.data, condition, replicate, seurat_clusters) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  group_by(condition, replicate) %>% 
  mutate(total_cells = sum(count), fraction = count / total_cells) %>%
  ungroup()
counts$condition_replicate <- interaction(counts$replicate, counts$condition, sep = "_")

# Save global fractions
write.csv(counts, file.path(plot_dir, "global_cluster_fractions.csv"), row.names = FALSE)

# Dotplot of global marker genes (Excluding low quality cluster 8)
dat_subset <- dat[, dat$seurat_clusters != 8]

global_markers <- list(
  DCs = c("Flt3", "Zbtb46"),
  cDC1s = c("Clec9a", "Xcr1"),
  cDC2s = c("Clec10a", "Mgl2"),
  `Migratory DCs` = c("Ccr7", "Ly75"),
  pDCs = c("Siglech", "Bcl11a"),
  `NK cells` = c("Ncr1", "Klrb1c"),
  `T cells` = c("Cd3g", "Trac"),
  `CD8 T cells` = c("Cd8a"),
  `CD4 T cells` = c("Cd4"),
  Tregs = c("Foxp3", "Tnfrsf4"),
  `Gamma-delta T cell` = c("Tcrg-C1", "Il23r"),
  `Monocytes/Macrophages` = c("Csf1r"),
  `M1-M2 polarization` = c("Ccl3", "Ccl4", "Ly6c2", "Atf3", "Pf4", "Fcrls", "Timp2", "Arg1", "Fn1", "Plac8", "S100a4"),
  ISGs = c("Isg15", "Rsad2", "Ifit2"),
  Neutrophils = c("S100a8", "S100a9"),
  `B cells` = c("Mzb1", "Pou2af1"),
  `Mast cells` = c("Cpa3", "Tpsab1"),
  Proliferation = c("Mki67", "Top2a", "Cdk1")
)

p4 <- DotPlot(object = dat_subset, features=global_markers, cluster.idents=TRUE, scale=FALSE) + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(file.path(plot_dir, "global_dotplot_cluster_markers.pdf"), plot = p4, width=700, height=150, units="mm", limitsize = FALSE)

# Clean up memory slightly without losing core variables
rm(dat_subset, all.genes, meta.data, counts, p1, p2, p3, p4, p_il12)
gc()

# -----------------------------------------------------------------------------
# Subclustering analysis (Figure 7, Ext. Figure 5)
# -----------------------------------------------------------------------------

subcluster_marker_analysis <- function(SeuratObject, SeuratClusters, MarkerList, Name) {
  dat_sub <- subset(SeuratObject, seurat_clusters %in% SeuratClusters)
  
  # Pre-Process
  dat_sub <- NormalizeData(dat_sub, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE, features = rownames(dat_sub)) %>%
    RunPCA(verbose = FALSE)
  
  dat_sub <- RunUMAP(dat_sub, dims = 1:30, reduction = 'pca')
  dat_sub <- FindNeighbors(dat_sub, reduction = "pca", dims = 1:30)
  dat_sub <- FindClusters(dat_sub)
  
  # UMAP plots
  set.seed(12)
  cols_cond <- DiscretePalette(2, shuffle = FALSE)
  p1 <- DimPlot(dat_sub, reduction = "umap", shuffle = TRUE, group.by = c("condition"), cols = cols_cond)
  ggsave(filename = file.path(plot_dir, paste0("UMAP_", Name, "_ByCondition.pdf")), plot = p1, width = 90, height = 75, units = "mm")
  
  cols_clust <- DiscretePalette(length(unique(dat_sub$seurat_clusters)), shuffle = FALSE)
  p2 <- DimPlot(dat_sub, reduction = "umap", shuffle = TRUE, group.by = c("seurat_clusters"), cols = cols_clust)
  ggsave(filename = file.path(plot_dir, paste0("UMAP_", Name, "_BySeuratClusters.pdf")), plot = p2, width = 90, height = 75, units = "mm")
  
  saveRDS(dat_sub, file.path(data_dir, paste0(Name, "_processed.rds")))
  
  # Fractions of cells by cluster and genotype
  meta.data <- dat_sub[[]]
  counts <- group_by(meta.data, condition, replicate, seurat_clusters) %>% 
    summarise(count = n(), .groups = 'drop') %>% 
    group_by(condition, replicate) %>% 
    mutate(total_cells = sum(count), fraction = count / total_cells) %>%
    ungroup()
  counts$condition_replicate <- interaction(counts$replicate, counts$condition, sep = "_")
  
  p3 <- ggplot(counts, aes(x = condition_replicate, y = fraction, fill = condition)) +
    facet_wrap(~ seurat_clusters, scales = "free") +
    geom_bar(stat = 'identity', position = position_dodge(width = 0.9)) + 
    scale_fill_manual(values = c("gray","dodgerblue")) + 
    theme_minimal() +
    labs(x = "Condition and Replicate", y = "Fraction of Cells", fill = "Condition") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  ggsave(filename = file.path(plot_dir, paste0("Barplot_Fractions_", Name, ".pdf")), plot = p3, width = 100, height = 100, units = "mm")
  
  # Dotplots
  p4 <- DotPlot(object = dat_sub, features = MarkerList, cluster.idents = TRUE, scale = FALSE) + 
    theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = file.path(plot_dir, paste0("Dotplot_", Name, "_Markers.pdf")), plot = p4, width = 300 * sum(lengths(MarkerList)) / 20, height = 100 * length(unique(dat_sub$seurat_clusters)) / 5, units = "mm")
  
  # Interleukin Violins
  p5 <- VlnPlot(object = dat_sub, features = c("Il12b","Il15","Il18"), split.by = "condition", pt.size=0) + 
    theme(legend.position = 'right')
  ggsave(filename = file.path(plot_dir, paste0("Violin_", Name, "_Interleukins.pdf")), plot = p5, width = 300 * sum(lengths(MarkerList)) / 20, height = 100 * length(unique(dat_sub$seurat_clusters)) / 5, units = "mm")
  
  return(dat_sub)
}

# DC subclustering execution
dc_clusters <- c(9, 12, 15, 19)

dc_specific_markers <- list(
  DC = c("Flt3", "Zbtb46", "Kit"),
  cDC1 = c("Clec9a", "Xcr1", "Itgae", "Batf3"),
  cDC2 = c("Ccl17", "Mgl2", "Il1r2", "Clec10a"),
  Migratory_DC = c("Ccr7", "Mreg", "Ccl22", "Fscn1", "Ly75"),
  pDC = c("Siglech", "Ccr9", "Pacsin1", "Dntt"),
  Proliferation = c("Mki67", "Top2a", "Cdk1")
)

dc_subset <- subcluster_marker_analysis(SeuratObject = dat,
                                        SeuratClusters = dc_clusters,
                                        MarkerList = dc_specific_markers,
                                        Name = "DCs")