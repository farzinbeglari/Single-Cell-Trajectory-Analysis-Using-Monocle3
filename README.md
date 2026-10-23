# Single-Cell-Trajectory-Analysis-Using-Monocle3
Exploring cellular progression and gene expression dynamics under normoxic and hypoxic conditions in lung cells (GSM10113694 vs. GSM10113696)

library(Seurat)
library(monocle3)
library(ggplot2)
library(SeuratWrappers)

# Data Import
data_directory <- getwd()
samples <- list.dirs(path = data_directory,
                     full.names = F,
                     recursive = F)
Data.list <- list()
for (sample in samples) {
  sample_directory <- file.path(data_directory, sample)
  
  counts <- ReadMtx(mtx = file.path(sample_directory, 'matrix.mtx.gz'),
                    cells = file.path(sample_directory, 'barcodes.tsv.gz'),
                    features = file.path(sample_directory, 'features.tsv.gz'),
                    feature.column = 2)
  
  Data.list[[sample]] <- CreateSeuratObject(counts = counts,
                                            project = sample)
} 

# Integrating and Merging of Single-Cell RNA-seq Datasets (GSM10113694 vs. GSM10113696)
merged.data <- merge(x = Data.list[[1]],
                     y = Data.list[-1],
                     add.cell.ids = names(Data.list))
merged.data$samples <- merged.data$orig.ident
View(merged.data@meta.data)

# Quality Control
merged.data[['percent.mt']] <- PercentageFeatureSet(merged.data, pattern = '^MT-')
View(merged.data@meta.data)
merged.data <- NormalizeData(merged.data)
merged.data <- FindVariableFeatures(merged.data)
merged.data <- JoinLayers(merged.data)

#  Trajectory Analysis of Single-Cell RNA-seq Using Monocle3
cds <- SeuratWrappers::as.cell_data_set(merged.data, group.by = 'samples')
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds, reduction_method = 'UMAP')
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = T)

# Identifying the Root Cells of a Trajectory to Establish the Direction of Cellular Progression Along Pseudotime
root_cells <- colnames(cds)[colData(cds)$samples == 'Normoxic control']
cds <- order_cells(cds,
                   root_cells = root_cells)
plot_cells(cds,
           color_cells_by = 'pseudotime')

# Identification of Pseudotime-Associated Genes
gene_results <- graph_test(cds, neighbor_graph = 'principal_graph')
significant.genes <- gene_results[gene_results$q_value < 0.05,]

Gene results identify genes with significant expression patterns across neighboring cells. These results help identify candidate genes associated with cellular trajectories and potential changes under hypoxic conditions.

Alas2 <- significant.genes[rownames(significant.genes) == 'Alas2',]
Col1a1 <- significant.genes[rownames(significant.genes) == 'Col1a1',]

Alas2 and Col1a1 showed significant spatial autocorrelation along the inferred cellular trajectory (q < 0.05, Moran’s I > 0.5), suggesting that their expression patterns vary across different regions of the trajectory. These findings may be consistent with previous studies investigating molecular changes associated with hypoxic conditions [1].

[1] Li X, Li M, Huang M, Li J, Huang S, Wang B, Gao Y, Wang Z, Shi L. Hypoxic response patterns in lung tissue: an integrated analysis of comparative physiological and transcriptomic studies from Neodon fuscus and Lasiopodomys brandtii. Science of the Total Environment. 2023 Sep 20;892:164537.


# Dataset Download Link
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350399 
