library(DTUrtle)
library(BiocParallel)
library(ggplot2)
library(Seurat)
library(dplyr)

# set appropriate CPU
biocpar <- BiocParallel::MulticoreParam(16)

# import mapping annotation that used during quantification 
tx2gene <- import_gtf(gtf_file  = '.PATH/TO/FILE/gencode.v47.annotation.gtf')
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)
tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

# check that moved transcript names to front 
head(tx2gene, n=5)

# name samples the names of the directories with the alevin quantification results SAMPLE_NAME_alevin

samples <- c("CP1_alevin", "CP2_alevin", "CP3_alevin", "CP6_alevin", "CP7_alevin")
files <- file.path("./PATH/TO/SAMPLE/DIRECTORY", samples, "alevin", "quants_mat.gz")

names(files) <- samples
names(files)

# import transcript levle counts and gene counts 
cts_list <- import_counts(files = files, type = "alevin")
cts_t <- combine_to_matrix(tx_list = cts_list, cell_extensions = samples, cell_extension_side = "prepend")
cts_g <- import_dge_counts(files = files, type = "alevin", tx2gene=tx2gene)
cts_g <- combine_to_matrix(tx_list = cts_g, cell_extensions = samples, cell_extension_side = "prepend")

# create seurat obj with gene_counts for integration and clustering 
seurat_obj <- CreateSeuratObject(counts = cts_g, min.cells = 5, min.features = 100)
seurat_obj$sample <- sapply(strsplit(rownames(seurat_obj@meta.data), "_"), function(x) paste(x[1:2], collapse = "_"))
table(seurat_object$sample)
seurat_list <- SplitObject(seurat_obj, split.by = "sample")
save.image(file = "single_cell_analysis_R.RData")

# folder for plots 
output_dir <- "./OUTPUT_NAME"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# finding feature/count/mt for each sample 
for (i in names(seurat_list)) {
  seu <- seurat_list[[i]]
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  set.seed(42)  
  cells_to_plot <- sample(colnames(seu), size = min(1000, ncol(seu)))
  seu_sub <- subset(seu, cells = cells_to_plot)
  
  vln_plot <- VlnPlot(
    seu_sub, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3, 
    pt.size = 1
  )
  print(vln_plot)
  ggsave(
    filename = file.path(output_dir, paste0(i, "_QC_violin.pdf")),
    plot = vln_plot,
    width = 10,
    height = 4
  )
  
  fs_plot <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(fs_plot)
  ggsave(
    filename = file.path(output_dir, paste0(i, "_FeatureScatter_nCount_vs_nFeature.pdf")),
    plot = fs_plot,
    width = 6,
    height = 5
  )
  
  pdf(file.path(output_dir, paste0(i, "_nCount_RNA_hist.pdf")), width = 8, height = 6)
  hist(seu$nCount_RNA, breaks = 100, main = paste("nCount_RNA distribution -", i), xlab = "Counts per cell")
  dev.off()
  
  seurat_list[[i]] <- seu
}
# chose thresholds of Feature <5000 nCountRNa <30,000 and percent.MT <15 to get rid of dead and likely doublet cells 
# loop over samples for applying filters and normalisation by SCTransform

for (i in names(seurat_list)) {
  seu <- seurat_list[[i]]
  
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA <30000 & percent.mt < 10)
  seu <- SCTransform(seu, verbose = FALSE)
  seurat_list[[i]] <- seu
}

# integrate for clustering anlysis 
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = features
)

seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

saveRDS(seurat_integrated, "seurat_integrated_final.rds")
saveRDS(seurat_list, "seurat_list_SCT_normalized.rds")
saveRDS(features,"integration_features.rds")


DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
ElbowPlot(seurat_integrated, ndims = 50)

#choose dims based on elbow plot 

dims_to_use <- 1:30
seurat_integrated <- FindNeighbors(seurat_integrated, dims = dims_to_use)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.1)
seurat_integrated <- RunUMAP(seurat_integrated, dims = dims_to_use)

plot_by_cluster_ids <- DimPlot(seurat_integrated, label = TRUE, group.by = "seurat_clusters") + ggtitle("Clusters")
ggsave(filename = file.path(output_dir, "DimPlot_clusters.pdf"), plot = plot_by_cluster_ids, width = 8, height = 6)
plot_by_cluster_samples <- DimPlot(seurat_integrated, reduction = "umap", group.by = "sample") + ggtitle("Samples")
ggsave(filename = file.path(output_dir, "DimPlot_samples.pdf"), plot = plot_by_cluster_samples, width = 8, height = 6)

# Find markers to identify the cell types 
seurat_integrated_cell_clusters <- data.frame(cluster = seurat_integrated$seurat_clusters)
Idents(seurat_integrated) <- "seurat_clusters"
all_markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# get top 10 for identification 
top10_markers <- top_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(top10_markers)

write.csv(top10_markers,file = file.path(output_dir, "top10_markers.csv"), row.names = FALSE)

# read in database for identification 
database_cell_markers <-read.csv("./qc_plots/PanglaoDB_markers_27_Mar_2020.tsv", sep = "\t", stringsAsFactors = FALSE)
database_cell_markers <- database_cell_markers %>%
  select(species,official.gene.symbol,cell.type,nicknames,organ)

matching_gene_cells <-top10_markers %>%
  filter(gene %in% database_cell_markers$official.gene.symbol) %>%
  select(cluster, gene)

#get cell types based on the most expressed genes 
matching_genes_with_celltypes <- matching_gene_cells %>%
  inner_join(database_cell_markers, by = c("gene" = "official.gene.symbol")) %>%
  select(gene, cluster, cell.type, organ, species)

# summarise by cluster with the cell types 
summary_by_cluster <- matching_genes_with_celltypes %>%
  group_by(cluster) %>%
  summarise(
    genes = paste(unique(gene), collapse = ", "),
    cell_types = paste(unique(cell.type), collapse = ", "),
    organs = paste(unique(organ), collapse = ", "),
    species = paste(unique(species), collapse = ", ")
  ) %>%
  arrange(as.numeric(as.character(cluster)))

#save file and inspect to assign cell type 

summary_by_cluster
write.csv(
  summary_by_cluster,
  file = file.path(output_dir, "cluster_markers_cell_type.csv"),
  row.names = FALSE
)

# also looked at markers within the clusters to check if agreed with database 
# tumour cells ("KRT14", "KRT17", "KRT18"), 
# fibroblasts ("COL1A1", "COL1A2", "DCN"), 
# endothelial cells ("VWF", "PECAM1", "AQP1")
# immune cells ("PTPRC", "IL7R", "GZMK", "MS4A1") 
# common in ACPs

# can generate feature plots to show chosen markers expression within cell clusters 
p1 <- FeaturePlot(seurat_integrated, 
                  features = c("KRT14", "KRT17", "KRT18", 
                               "COL1A1", "COL1A2", "DCN",
                               "VWF", "PECAM1", "AQP1",
                               "PTPRC", "IL7R", "GZMK", "MS4A1"), 
                  label = TRUE)



# assign final cell types to the groups 
summary_by_cluster <- summary_by_cluster %>%
  mutate(cell_types = case_when(
    cluster == 0 ~  "T_cells",
    cluster == 1 ~ "M2_macrophages",
    cluster == 2 ~ "Fibroblasts",
    cluster == 3 ~ "Endothelial_cells",
    cluster == 4 ~ "B_cells",
    cluster == 5 ~ "Plasma_cells",
    cluster == 6 ~ "Erythroid_cells",
    cluster == 7 ~ "Dendritic_cells",
    cluster == 8 ~ "Tumour_Epithelial_cells",
    cluster == 9 ~ "Other_immune_cells",a
    cluster == 10 ~ "Mast_cells",
    cluster == 11 ~ "Basal_epithelial_cells"
  ))

# save updated csv with cell types 
write.csv(
  summary_by_cluster,
  file = file.path(output_dir, "cluster_cell_types.csv"),
  row.names = FALSE
)

# set assat to visulise expression
DefaultAssay(seurat_integrated)

DimPlot(seurat_integrated, label = FALSE, group.by = "cell_type") + ggtitle("Annotated Cell Types")
# can ggsave(filename = file.path(output_dir, "cell_clusters_cell_type.png"), plot = plot_by_cluster_ids, width = 10, height = 7)

# update default assay by cell type + normalised data and can visualise chosen gene 
DefaultAssay(seurat_integrated) <- "SCT"
Idents(seurat_integrated) <- "cell_type"

FeaturePlot(seurat_integrated, features = "CPLX1", reduction = "umap", label = TRUE) + ggtitle("CPLX1 Expression by Cell Type")

Idents(seurat_integrated) <- "seurat_clusters"

# Subsetting the alevin transcript based on cell type identification + post quality control 
# using DTUrtle to get transcirpt level matrix of seurat object and do DTU based on the types of cell clusters

tissue <- combine_to_matrix(tx_list = cts_list, seurat_obj = seurat_integrated, tx2gene = tx2gene, cell_extension_side = "prepend")

tissue@active.assay
head(tissue@meta.data)

# make default assay transcript matrix
DefaultAssay(tissue) <- "dtutx"

# assign cells to groups for anlysis between clusters
cell_groups <- c(
  "t_cells" = "Immune",
  "macrophages_M2_type" = "Immune",
  "fibroblasts" = "Tumour/Other",
  "endothelial_cells" = "Tumour/Other",
  "tumour_epithelial_cells" = "Tumour/Other",
  "basal_epithelial_cells" = "Tumour/Other", 
  "b_cells" = "Immune",
  "plasma_cells"= "Immune", 
  "erythoid_cells" = "Immune",
  "Dendritic_cells" = "Immune", 
  "mast_cells" = "Immune",
  "other_immune_cells" = "Immune"
)

cell_types_only <- unique(tissue$cell_type)

# get common cells between the alevin output post QC processing
tissue_list <- lapply(cell_types_only, function(ct) {
  cells <- WhichCells(tissue, expression = cell_type == ct)
  Matrix::rowSums(GetAssayData(tissue, assay = "dtutx", slot = "counts")[, cells])
})
tissue_counts <- do.call(cbind, tissue_list)
colnames(tissue_counts) <- cell_types_only

# set up metadata df
pd <- data.frame(
  id = colnames(tissue_counts),
  stringsAsFactors = FALSE
)

pd$group <- cell_groups[pd$id]

pd$group <- factor(pd$group)

all(rownames(tissue_counts) %in% tx2gene$transcript_name)

#check any missing transcripts between map and object  
missing_tx <- rownames(tissue_counts)[!rownames(tissue_counts) %in% tx2gene$transcript_name]
tx2gene_filtered <- tx2gene[tx2gene$transcript_name %in% rownames(tissue_counts), ]
tissue_counts_filtered <- pseudo_bulk_counts[rownames(pseudo_bulk_counts) %in% tx2gene_filtered$transcript_name, ]
all(rownames(tissue_counts_filtered) %in% tx2gene_filtered$transcript_name)

# run differential transcript usage testing 
dtu <- run_drimseq(
  counts = pseudo_bulk_counts_filtered,
  tx2gene = tx2gene_filtered,
  pd = pd_bulk,
  id_col = "id",
  cond_col = "group",
  filtering_strategy = "sc",
  cond_levels = c("Tumour/Other", "Immune"),  # or any two levels you're comparing, 
  BPPARAM = BiocParallel::MulticoreParam(10)
)
dtu_post <- posthoc_and_stager(dturtle = dtu, ofdr = 1, posthoc = 0)

# visualising results can select interested gene to look at the heatmap/proportion plot of isofrom expression 

cplx1_isoform_plot <- FeaturePlot(bulk_test_tissue, features = c("CPLX1-201", "CPLX1-203"), order = TRUE, label = TRUE)
cplx1_barplot <- plot_proportion_barplot(dturtle = dtu_post,  genes = "CPLX1",  meta_gene_id = "gene_id.1")
cplx1_barplot$CPLX1

cplx1_heatmap <- plot_proportion_pheatmap(dturtle = dtu_post, genes = "CPLX1", sample_meta_table_columns = c("sample_id","condition"),                                 include_expression = TRUE, treeheight_col=20)
temp_2$CPLX1

cplx1_track_plot <- plot_transcripts_view(dturtle = dtu_post,   genes = "CPLX1",   gtf = "PATH_TO_FILE/gencode.v47.annotation.gtf",   genome = 'hg38',  one_to_one = TRUE)

# visualise CPLX1-201 and CPLX2-203 expression plots 

FeaturePlot(tissue, features = c("CPLX1-201", "CPLX1-203"), order = TRUE, label = TRUE)
ggsave(filename = file.path(output_dir, "CPLX1-201_203_expression_tx.pdf"), plot = p9, width = 10, height = 6)

