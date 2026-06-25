#analysis 19 dimenssions 1 resolution

##load packages
library(Seurat)
library(tidyverse)
library(future)
plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)
library(ggplot2)
library(leidenbase)
library(devtools)
library(dittoSeq)
library(projectLSI)

set.seed(1234)

#load the files 
untar("filtered_feature_bc_matrix.tar.gz", exdir = "data")
raw_data_1<- Read10X(data.dir = "data")

data<-CreateSeuratObject(counts = raw_data_1[["Gene Expression"]],
                                min.cells = 3,
                                min.features = 200, 
                                project = "TFM", names.delim = "-", names.field = 2)

pdf(file="features and counts.pdf", width = 15, height = 15)
VlnPlot(object = data, features =  "nFeature_RNA", group.by = "orig.ident", pt.size = 0.1) +  ggtitle("nFeature_RNA filtered") + theme(legend.position = "none") + xlab("")
VlnPlot(object = data, features =  "nCount_RNA", group.by = "orig.ident", pt.size = 0.1) +  ggtitle("nCount_RNA filtered") + theme(legend.position = "none") + xlab("")
dev.off()


##filter cells 

genes <- rownames(data[["RNA"]])
mt_genes <- grep("^MT-", genes, value = TRUE)
ribo_genes<-grep("RPS|RPL|MRPL|MRPS" , genes, value = TRUE)

mt_genes  
ribo_genes

data[["percent.mito"]] <- Seurat::PercentageFeatureSet(
  data, pattern = "^MT-"
)

data[["percent.ribo"]] <- Seurat::PercentageFeatureSet(
  data, pattern = "RPS|RPL|MRPL|MRPS"
)

print(VlnPlot(data,
              features = c("nCount_RNA", "nFeature_RNA", "percent.mito"),
              ncol = 3, pt.size = 0.1))

pdf(file="raw counts.pdf", width = 15, height = 15)
VlnPlot(data,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"),
        ncol = 4, pt.size = 0.1)
dev.off()


pdf(file="filter genes.pdf", width = 15, height = 15)

FeatureScatter(object = data, feature1 = "nCount_RNA", 
                        feature2 = "percent.mito")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeatureScatter(object = data, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 FeatureScatter(object = data, feature1 = "nCount_RNA", 
                        feature2 = "percent.ribo")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



VlnPlot(object = data, features = c( "percent.mito", "nFeature_RNA", "percent.ribo"),
        group.by = "orig.ident", pt.size = 0.1)



#I have changed this!
#data <- subset( x = data,subset = nFeature_RNA > 500 & nFeature_RNA < 12000 &
#                  percent.mito < 5)

VlnPlot(object = data, features = c( "percent.mito", "nFeature_RNA", "percent.ribo"),
        group.by = "orig.ident", pt.size = 0.1)

dev.off()


data <- subset(data,
               subset = nFeature_RNA > 4000 &
                 #nFeature_RNA < 10000 &
                 nCount_RNA < 100000 &
                 percent.mito < 10
)

print(VlnPlot(data,
              features = c("nCount_RNA", "nFeature_RNA", "percent.mito"),
              ncol = 3, pt.size = 0.1))




##Merged with the full_annotation information from the HUmanATLAS

#to remove the cell name identifier
for (i in 1:6) {
  colnames(data) <- make.unique(gsub(paste0("-", i, "$"), "", colnames(data)))
}

df1 <- read.csv("/Users/elenacamachoaguilar/Library/CloudStorage/OneDrive-UNIVERSIDADDESEVILLA/Research Projects/2025-TFM-InmaculadaVargasRomero/Analysis/scRNAseq/HuEm_stable_reference_projection_tool/full_annotation (3).csv",
                header = T, sep=";")
df2 <- data@meta.data

common_cells <- intersect(df1$query_cell, colnames(data))

df1$key <- df1$query_cell
pred <- df1$sub_pred_EML
names(pred) <- df1$query_cell
df2$cluster <- pred[rownames(df2)]

df2 -> data@meta.data

head(df2)
tail(df2)


#name the orig.ident as the condition and create a new category with the time 

ident <- c(
  "1" = "BMP2.5_0-30h",
  "2" = "BMP10_0-16h_mTeSR_0-30h",
  "3" = "BMP10_0-30h",
  "4" = "BMP2.5_0-48h",
  "5" = "BMP10_0-16h_mTeSR_0-48h",
  "6" = "BMP10_0-48h_nogin_30-48h"
)

new_idents <- ident[as.character(data$orig.ident)]
names(new_idents) <- colnames(data)
data$orig.ident <- new_idents

unique(data@meta.data$orig.ident)

time <- c(
  "BMP2.5_0-30h"= "30h",
  "BMP10_0-16h_mTeSR_0-30h"= "30h",
  "BMP10_0-30h" = "30h",
  "BMP2.5_0-48h" = "48h",
  "BMP10_0-16h_mTeSR_0-48h"= "48h",
  "BMP10_0-48h_nogin_30-48h"= "48h"
)

#Adding the new columns
data <- AddMetaData(data,
                    metadata = data.frame(time = time[data$orig.ident],
                                          row.names = colnames(data)
                    ))

head(data@meta.data)

##normalization

normalize <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the 10 most highly variable genes with 3000 features

feature3000 <- FindVariableFeatures(normalize, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(feature3000), 10)

# plot variable features with and without labels

pdf(file="features.pdf", width = 15, height = 15)
plot1 <- VariableFeaturePlot(feature3000)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

#check the precense of the genes 

genes_interes <- c("SOX2", "NANOG", "POU5F1", "USP44", "GATA3", "ISL1", "TFAP2A", 
                   "TFAP2C", "CDX2", "HAND1", "SOX17", "NODAL", "CDH1", "CDH2",
                   "SNAI1", "TBXT","MIXL1")

#Check if the genes of interest exist in the pre-filtered matrix
genes_interes %in% rownames(normalize)

#Check if the genes of interest exist in the filtered matrix
genes_interes %in% VariableFeatures(feature3000)


##adding the IMPORTANT genes that were loosed

VariableFeatures(feature3000) <- unique(c(
  VariableFeatures(feature3000),
  genes_interes
))
#Check if the genes of interest were added

genes_interes %in% VariableFeatures(feature3000)

#cell cycle clasiffication 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all.genes <- rownames(normalize)

# Classification of cells according to their cell cycle stage

scale <- CellCycleScoring(feature3000, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

head(scale[[]])

#applying data standardization to remove the influence of cell cycle variation

scale <-  ScaleData( 
  scale,
  vars.to.regress = c("S.Score", "G2M.Score")
)

# Run the standard workflow for visualization and clustering and representations 
#of the variations

PCA <- RunPCA(scale,features = VariableFeatures(object = scale),  nfeatures.print = 10)
print(PCA[["pca"]], dims = 1:5, nfeatures = 10)

pdf("ElbowPlot_26_5.pdf")

DimPlot(PCA, reduction = "pca", dims = 2:3) + NoLegend()
ElbowPlot(PCA)

dev.off()


#Extract the standard deviation of each PC

stdev <- PCA[["pca"]]@stdev

#Calculate the percentage of variance explained
var_explained <- (stdev^2 / sum(stdev^2))

head(var_explained)

pdf("PCA_var_26_5.pdf")
plot(
  var_explained,
  type = "b",
  xlab = "Principal Component",
  ylab = "Percentage of Variance Explained"
)
dev.off()
#Create a bar plot of explained variance

pdf("Barplot_variance_26_5.pdf")

barplot(
  var_explained,
  names.arg = 1:length(var_explained),
  xlab = "Principal Components",
  ylab = "% Variance Explained",
  main = "Variance Explained by Each PC"
)

plot(
  var_explained,
  type = "b",
  pch = 19,
  xlab = "Principal Component",
  ylab = "% Variance Explained",
  main = "Scree Plot"
)

dev.off()

#Calculate cumulative variance explained

pc <- 1:length(var_explained)

explained <- var_explained
cumvar <- cumsum(var_explained)

df <- data.frame(
  PC = pc,
  Explained = explained,
  Cumulative = cumvar
)

pdf("PCA_variance_plot_26_5.pdf", width = 14, height = 7)

ggplot(df, aes(x = PC)) +
  
  geom_bar(aes(y = Explained),
           stat = "identity",
           fill = "forestgreen",
           alpha = 0.7) +
  
  geom_line(aes(y = Cumulative),
            color = "red",
            linewidth = 1.2) +
  
  geom_point(aes(y = Cumulative),
             color = "red",
             size = 2) +
  
  geom_text(aes(y = Cumulative,
                label = paste0(round(Cumulative * 100), "%")),
            vjust = -0.5,
            size = 3) +
  
  scale_x_continuous(breaks = pc) +
  scale_y_continuous(labels = scales::percent) +
  
  labs(
    title = "Explained Variance by Principal Components",
    x = "Principal Components",
    y = "Explained Variance"
  ) +
  
  theme_minimal()

dev.off()

# find neighbors and Clustering (PCA 19, leiden res 1)

sem_19 <- FindNeighbors(PCA, dims = 1:19,  k.param = 10)
sem_19<- FindClusters(
  sem_19,
  algorithm = 4,
  resolution = 1,
  modularity.fxn = 1,
  group.singletons = F,
  n.start = 20,
  n.iter = 50,
  random.seed = 1234,
  leiden_objective_function = "CPM")
sem_19<- RunUMAP(sem_19, dims = 1:19, min.dist=0.1, n.neighbors = 50)
pdf("DimPlot_19_1.pdf", width = 14, height = 7)
DimPlot(sem_19, reduction = "umap", label = TRUE)+
  ggtitle("19 dims resolution 1")
dev.off()

p1 <- DimPlot(sem_19, reduction = "umap", label = TRUE) +
  ggtitle("Clusters")
p2 <- DimPlot(sem_19, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample")

print(p1 + p2)

# find neighbors and Clustering (PCA 19, leiden res 1.5)

sem_19_15 <- FindNeighbors(PCA, dims = 1:19,  k.param = 10)
sem_19_15<- FindClusters(
  sem_19_15,
  algorithm = 4,
  resolution = 1.5,
  modularity.fxn = 1,
  group.singletons = F,
  n.start = 20,
  n.iter = 50,
  random.seed = 1234,
  leiden_objective_function = "CPM")
sem_19_15<- RunUMAP(sem_19_15, dims = 1:19, min.dist=0.1, n.neighbors = 50)
pdf("DimPlot_19_15.pdf", width = 14, height = 7)
DimPlot(sem_19_15, reduction = "umap", label = TRUE)+
  ggtitle("19 dims resolution 1.5")
dev.off()

p1 <- DimPlot(sem_19_15, reduction = "umap", label = TRUE) +
  ggtitle("Clusters")
p2 <- DimPlot(sem_19_15, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample")

print(p1 + p2)

pdf("DimPlot_AtlasLabels_19_15.pdf", width = 14, height = 7)
DimPlot(sem_19_15, reduction = "umap", group.by = "cluster") +
  ggtitle("Human Atlas Labels")
dev.off()

####################################################################
# RUN HARMONY INTEGRATION
####################################################################

library(harmony)

sem_harmony <- RunHarmony(sem_19,
                          group.by.vars = "time",
                          reduction = "pca",
                          reduction.save = "harmony")

sem_harmony <- RunUMAP(sem_harmony,
                       reduction = "harmony",
                       dims = 1:19,
                       min.dist = 0.1,
                       n.neighbors = 50)

sem_harmony <- FindNeighbors(sem_harmony,
                             reduction = "harmony",
                             dims = 1:19,
                             k.param = 10)

sem_harmony <- FindClusters(sem_harmony,
                            resolution = 1)

p1 <- DimPlot(sem_harmony, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample")
p2 <- DimPlot(sem_harmony, reduction = "umap", label = TRUE) +
  ggtitle("Clusters")
print(p1 + p2)

####################################################################
# SELECT WHAT OBJECT TO USE

sem<- sem_19_15
#sem<- sem_19
#sem<-sem_harmony
####################################################################

####################################################################
# PLOT UMAPS GENES

print(FeaturePlot(sem, features = c("SOX2", "CDX2", "HAPLN1", "TBXT"),
                ncol = 2,
                cols = c("lightgrey", "red")))

print(FeaturePlot(sem, features = c("SOX2", "CDX2", "HAPLN1", "TBXT"),
                  ncol = 2,
                  cols = c("lightgrey", "red")))

pdf("FeaturePlot_Pluripotent_NoHarmony.pdf", width = 14, height = 5)
print(FeaturePlot(sem, features = c("SOX2", "NANOG", "USP44"),
                  ncol = 3,
                  cols = c("lightgrey", "red")))
dev.off()

pdf("FeaturePlot_PS_NoHarmony.pdf", width = 14, height = 5)
print(FeaturePlot(sem, features = c("EOMES", "TBXT", "MIXL1"),
                  ncol = 3,
                  cols = c("lightgrey", "red")))
dev.off()

pdf("FeaturePlot_Amnion_NoHarmony.pdf", width = 14, height = 5)
print(FeaturePlot(sem, features = c("HAND1", "CDX2", "ISL1"),
                  ncol = 3,
                  cols = c("lightgrey", "red")))
dev.off()

pdf("FeaturePlot_ExEM_NoHarmony.pdf", width = 14, height = 5)
print(FeaturePlot(sem, features = c("VIM", "SNAI1", "ANXA1"),
                  ncol = 3,
                  cols = c("lightgrey", "red")))
dev.off()

pdf("FeaturePlot_Signals_NoHarmony.pdf", width = 5, height = 14)
print(FeaturePlot(sem, features = c("HAPLN1", "WNT6", "BMP2"),
                  ncol = 1,
                  cols = c("lightgrey", "red")))
dev.off()

DimPlot(sem, 
        reduction = "umap", 
        group.by = "orig.ident") +
  ggtitle("Sample")


####################################################################
# PLOT UMAPS SAMPLES vs LEIDEN
library(ggrepel)
p1 <- DimPlot(sem, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample")


p <- DimPlot(sem,
             reduction = "umap",
             label = FALSE,
             pt.size = 0.5) +
  theme(legend.position = "none") +
  ggtitle("Clusters")

# Add labels with white background
p2<-LabelClusters(p, 
              id = "ident",
              fontface = "bold",
              size = 5,
              box = TRUE,
              fill = "white",
              alpha = 0.7)
print(p1 + p2)
pdf("UMAP_NoHarmony_LeidenRes1p5.pdf", width = 14, height = 7)
p1+p2
dev.off()

DimPlot(sem, reduction = "umap", group.by = "time") +
  ggtitle("Time")

#sem_rep2 <- RunUMAP(sem,
#               reduction = "harmony",
#               dims = 1:19,
#               min.dist = 0.3,
#               n.neighbors = 30,
#               spread = 1.5)

####################################################################
# Calculate cluster proximity using silhouette analysis
# 

library(cluster)

# Extract PCA coordinates (first 19 principal components) for each cell
pca <- Embeddings(sem, "pca")[, 1:19]

# Cluster assignments generated by Seurat
clusters <- Idents(sem)

# Compute the pairwise distance matrix between cells based on PCA coordinates
# Cells with similar PC values have small distances,
# whereas cells with different PC values have large distances
d <- dist(pca)

# Calculate silhouette scores for each cell
# s(i) = (b(i) - a(i)) / max(a(i), b(i))
#
# where:
# a(i) = average distance to other cells within the same cluster
# b(i) = average distance to cells in the nearest neighboring cluster
#
# Silhouette values range from -1 to 1:
#   close to 1  -> well-classified cells
#   around 0    -> cells located at cluster boundaries
#   below 0     -> potentially misclassified cells

sil <- silhouette(
  as.numeric(clusters),
  d
)

# Store silhouette scores in the Seurat metadata
sem$silhouette <- sil[, "sil_width"]

# Calculate the mean silhouette score for each cluster
resumen_sil <- aggregate(
  sil[, "sil_width"],
  by = list(cluster = clusters),
  FUN = mean
)

pdf("silhouette_analysis_noharmony_leiden1p5.pdf",
    width = 14,
    height = 20)

# Classical silhouette plot (Rousseeuw plot)
# Provides a visual assessment of cluster quality
plot(
  sil,
  col = scales::hue_pal()(length(levels(clusters))),
  main = "Silhouette Plot"
)
dev.off()

pdf("silhouette_analysis_optimized_2_noharmony_Leiden1p5.pdf",
    width = 14,
    height = 7)
# Display silhouette scores on the UMAP embedding
# Red cells indicate potentially problematic classifications
FeaturePlot(
  sem,
  features = "silhouette"
) +
  scale_colour_gradient2(
    low = "firebrick3",    # Negative values (potentially misclassified cells)
    mid = "white",         # Boundary cells (score ≈ 0)
    high = "dodgerblue4",  # Well-classified cells
    midpoint = 0
  ) +
  ggtitle("Silhouette Score on UMAP (Red = Potentially Misclassified)")

# Violin plot showing the distribution of silhouette scores per cluster
VlnPlot(
  sem,
  features = "silhouette",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red"
  ) # Reference line at silhouette = 0

dev.off()

#####################################################################
#FIND MARKERS
#find the markers for every cluster compared to all remaining cells 
sem.markers <- FindAllMarkers(
  object = sem,
  only.pos = TRUE,
  features = VariableFeatures(sem),
  test.use = "roc",
  random.seed = 1234
)

write.csv(sem.markers, "markers_noharmony_roc_leiden1p5.csv", row.names = FALSE)

top5 <- sem.markers %>%
  group_by(cluster) %>%
  dplyr::filter(myAUC > 0.65 & avg_log2FC > 0.25) %>%
  slice_max(myAUC, n = 5) %>%
  ungroup()

pdf("heatmap_noharmony_leiden1p5.pdf", width = 14, height = 10)
print(DoHeatmap(sem, 
                features = top5$gene, 
                group.by = "seurat_clusters") +
        ggtitle("Top 5 ROC markers per cluster"))
dev.off()

genes_ordered <- c(
  # Pluripotente
  "POU5F1", "SOX2", "NANOG", "USP44",
  # PS
  "TBXT", "EOMES", "MIXL1",
  # Amnion
  "TFAP2A", "GATA3", "HAND1", "ISL1", "CDX2",
  # ExEM
  "HAPLN1", "VIM", "ANXA1", "SNAI1",
  # OTHER
  "CR1L", "NODAL", "LEFTY1", "HPGD", "NTS"
)

genes_ordered_expanded <- c(
  # Pluripotent / epiblast-like
  "POU5F1", "SOX2", "NANOG", "USP44", "APELA", "OTX2",
  
  # Primitive streak / mesendoderm
  "TBXT", "EOMES", "MIXL1", "MESP1", "KDR", "PDGFRA",
  
  # Amnion / extraembryonic epithelial-like
  "TFAP2A", "TFAP2C", "GATA3", "GATA2", "CDX2",
  "HAND1", "ISL1", "KRT8", "KRT18", "KRT19",
  "ID1", "ID2", "ID3",
  
  # Extraembryonic mesoderm / mesenchymal ECM-like
  "HAPLN1", "VIM", "ANXA1", "SNAI1", "SNAI2",
  "COL1A1", "COL3A1", "LUM", "DCN", "FN1", "IGFBP7",
  
  # Signalling / transitional
  "BMP2", "BMP4", "WNT6",
  "CR1L", "NODAL", "LEFTY1", "LEFTY2", "HPGD", "NTS", "GAL"
)
pdf("heatmap_noharmony_ALLkeymarkers_leiden1p5.pdf", width = 14, height = 10)
print(DoHeatmap(sem, 
                features = genes_ordered_expanded, 
                group.by = "seurat_clusters") +
        ggtitle("Key cell type markers"))
dev.off()

genes_ordered_clean <- c(
  # Pluripotent / epiblast-like
  "POU5F1", "SOX2", "NANOG", "USP44", "APELA",
  
  # Primitive streak / mesoderm
  "TBXT", "EOMES", "MIXL1", "KDR",
  
  # Amnion / extraembryonic epithelial-like
  "TFAP2A", "TFAP2C", "GATA3", "GATA2", "CDX2",
  "HAND1", "ISL1",
  "KRT8", "KRT18", "KRT19",
  "ID1", "ID2", "ID3", "BMP4",
  
  # Extraembryonic mesoderm / mesenchymal
  "HAPLN1", "VIM", "ANXA1", "SNAI1",
  "TAGLN", "CALD1", "CNN1", "MYL9",
  
  # Signalling / transitional
  "BMP2", "WNT6", "CR1L", "NODAL",
  "LEFTY1", "HPGD", "NTS", "GAL"
)
pdf("heatmap_noharmony_ALLkeymarkersCLEAN_leiden1p5.pdf", width = 14, height = 10)
print(DoHeatmap(sem, 
                features = genes_ordered_clean, 
                group.by = "seurat_clusters") +
        ggtitle("Key cell type markers"))
dev.off()

###### SAVE EXPRESSION AS TABLE

genes_use <- genes_ordered_clean  # or genes_ordered_expanded

# Make sure genes are present in scale.data
genes_use <- intersect(genes_use, rownames(sem[["RNA"]]@scale.data))

heatmap_matrix <- sem[["RNA"]]@scale.data[genes_use, ]

# Convert to table
heatmap_table <- as.data.frame(as.matrix(heatmap_matrix))

# Add gene names as a column
heatmap_table$gene <- rownames(heatmap_table)

# Put gene column first
heatmap_table <- heatmap_table[, c("gene", setdiff(colnames(heatmap_table), "gene"))]

# Save
write.csv(heatmap_table, "heatmap_scaled_values_per_cell.csv", row.names = FALSE)

##########################################################################
# PLOT SAMPLE COMPOSITION PER CLUSTERS 

# Get cluster and sample info
df <- sem@meta.data %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(n = n(), .groups = "drop")

# Plot
p <- ggplot(df, aes(x = factor(seurat_clusters), y = n, fill = orig.ident)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster", y = "Number of cells", fill = "Sample") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p)

pdf("barplot_clustersamplecomposition_leiden1p5.pdf", width = 14, height = 10)
print(p)
dev.off()


##########################################################################

##########################################################################
# PLOT ATLAS COMPOSITION PER CLUSTERS 

# Get cluster and sample info
df <- sem@meta.data %>%
  group_by(seurat_clusters, cluster) %>%
  summarise(n = n(), .groups = "drop")

# Plot
p <- ggplot(df, aes(x = factor(seurat_clusters), y = n, fill = cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster", y = "Number of cells", fill = "Sample") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p)

pdf("barplot_clusteratlascomposition_leiden1p5.pdf", width = 14, height = 10)
print(p)
dev.off()




##########################################################################



sem.markers <- FindAllMarkers(
  object = sem,
  only.pos = T,                              
  features = VariableFeatures(sem),
  test.use = "roc", 
  random.seed = 1234
)

write_xlsx(sem.markers, path = "marcadores_leiden_roc_19:1.xlsx")


top5 <- sem.markers %>%
  group_by(cluster) %>%
  dplyr::filter(myAUC > 0.65 & avg_log2FC > 0.25) %>% 
  slice_max(myAUC, n = 5) %>% 
  ungroup()


top10 <- sem.markers %>%
  group_by(cluster) %>%
  dplyr::filter(power > 0.3) %>% 
  slice_max(power, n = 10) %>%  
  ungroup()

pdf("heatmap_19_1_(2)_1_2_6.pdf",  width = 14, height = 10)
DoHeatmap(sem, features = top5$gene, group.by = "seurat_clusters") +  ggtitle("top 5 genes only possitive")
dev.off()

######################################
pdf("umap_per_experiment_19_1_30_05.pdf",  width = 14, height = 7)
DimPlot(sem, reduction = "umap", group.by = "orig.ident",label = TRUE)
dev.off()

#named the clusters after check in the csv documents the markets 

unique(sem@meta.data$orig.ident)


sem_clustering <- sem 

#These tables are used to examine the proportions of each cluster across the different experiments.

tabla <- sem_clustering@meta.data %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop")

tabla_pct <- tabla %>%
  group_by(orig.ident) %>%
  mutate(porcentaje = 100 * n / sum(n))

tabla_final <- tabla_pct %>%
  select(orig.ident, seurat_clusters, porcentaje) %>%
  pivot_wider(
    names_from = seurat_clusters,
    values_from = porcentaje
  )

tabla_pct_2 <- tabla %>%
  group_by(seurat_clusters) %>%
  mutate(porcentaje = 100 * n / sum(n))


tabla_fina_2 <- tabla_pct_2 %>%
  select(orig.ident, seurat_clusters, porcentaje) %>%
  pivot_wider(
    names_from = seurat_clusters,
    values_from = porcentaje
  )


cluster_id <- c(
  "1"  = "amnion",
  "2"  = "amnion",
  "3"  = "amnion",
  "4"  = "ExEM",
  "5"  = "pluri",
  "6"  = "pluri_PS",
  "7"  = "amnion",
  "8"  = "PS_mesoderm",
  "9"  = "amnion",
  "10" = "PS_mesoderm",
  "11" = "pluri",
  "12" = "ExEM",
  "13" = "PS_mesoderm",
  "14" = "ExEM",
  "15" = "pluri",
  "16" = "pluri",
  "17" = "pluri",
  "18" = "amnion",
  "19" = "amnion",
  "20" = "pluri_PS",
  "21" = "pluri",
  "22" = "pluri",
  "23" = "ambiguous",
  "24" = "amnion_ExEM",
  "25" = "ExEM"
)


sem <- AddMetaData(sem,
                    metadata = data.frame(cluster_id = cluster_id[sem$seurat_clusters],
                                          row.names = colnames(data)
                    ))




library(dplyr)
#to calculate the proportion of each cell cluster within each experimental condition (orig.ident).

#Extract metadata and count cells per group
df <- sem@meta.data %>%
  group_by(orig.ident, cluster_id) %>%
  summarise(n = n(), .groups = "drop")

#Compute proportions per sample
df <- df %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n))

pdf("barras.pdf",  width = 14, height = 7)

colores <- c( "mesoderm"        = "#ffcc00", 
              "amnion"     = "#00FFFF",  
              "pluri"           = "#ff33ff", 
              
              "extra_mesoder"  =  "#3A9B7A",  
              "pluri_extra"     = "#B266FF",  
              "pluri_meso"      = "#FF7A7A"  
              )

ggplot(df, aes(x = orig.ident, y = prop, fill = cluster_id)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            size = 3) +
  scale_fill_manual(values = colores) +
  ylab("cell proportion") +
  xlab("Experimet") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



DimPlot(sem, reduction = "umap", label = TRUE, cols = colores, group.by = "cluster_id")



df <- sem@meta.data %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop")

df <- df %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n))


colores_2 <- c(
  # Amnion - blues
  "1"  = "#D6ECFF",
  "2"  = "#A9D8FF",
  "3"  = "#7CC4FF",
  "7"  = "#4DAFFF",
  "9"  = "#1E95E6",
  "18" = "#006BB6",
  "19" = "#004C8C",
  
  # PS / mesoderm - yellows / oranges
  "8"  = "#FFF3B0",
  "10" = "#FFD166",
  "13" = "#F4A261",
  
  # Pluripotent - pinks
  "5"  = "#FFD6E8",
  "11" = "#FFB3D1",
  "15" = "#FF8FBD",
  "16" = "#F25CA2",
  "17" = "#D93682",
  "21" = "#B51F6A",
  "22" = "#8E1554",
  
  # ExEM - greens
  "4"  = "#D8F3DC",
  "12" = "#95D5B2",
  "14" = "#52B788",
  "25" = "#2D6A4F",
  
  # Mixed / transitional / ambiguous - greys
  "6"  = "#CFCFCF",
  "20" = "#AFAFAF",
  "23" = "#7F7F7F",
  "24" = "#5F5F5F"
)



ggplot(df, aes(x = orig.ident, y = prop, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            size = 3) +
  scale_fill_manual(values = colores_2) +
  ylab("cell proportion") +
  xlab("Experimet") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df, aes(x = seurat_clusters, y = prop, fill = seurat_clusters)) +
  geom_bar(stat = "identity", show.legend = FALSE)
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            vjust = -0.5,                             
            size = 2.5) +                            
  facet_wrap(~ orig.ident, ncol = 2) +                
  scale_fill_manual(values = colores_2) +               
  ylab("Cell proportion") +
  xlab("Seurat Clusters") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9),            
    strip.background = element_blank(),              
    strip.text = element_text(face = "bold", size = 11), 
    panel.spacing = unit(1.5, "lines")               
  )

DimPlot(sem, reduction = "umap", label = TRUE, cols = colores_2)

dev.off()

# Make a dataframe from Seurat metadata
df_counts <- sem@meta.data %>%
  count(orig.ident, seurat_clusters, name = "n_cells")
ggplot(df_counts, aes(x = seurat_clusters, y = n_cells, fill = seurat_clusters)) +
  geom_col(width = 0.8) +
  facet_wrap(~ orig.ident) +
  scale_fill_manual(values = colores_2) +
  labs(
    x = "Cluster",
    y = "Number of cells",
    title = "Number of cells per cluster in each experimental condition"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )


##############

pdf("genesPlot_UMAP_19_1.pdf", width = 14, height = 7)

genes <- 
  c("NANOG", "USP44","SOX2","ISL1", "HAND1", "GATA3", "TFAP2A", "TFAP2C", "CDX2",
    "VIM",  "ANXA1","TBXT", "MIXL1", "NODAL", "EOMES","CDH1", "CDH2")

for (g in genes) {
  print(
    FeaturePlot(sem, features =g, pt.size = 0
    ) +
      ggtitle(g)
  )
  
  print(
    
    VlnPlot(sem, features =g, pt.size = 0,  cols = colores_2
    ) +
      ggtitle(g)
  )
}
dev.off()


pdf("DotPlot_markers.pdf", width = 12, height = 8)

sem$cluster_label <- paste0(
  Idents(sem),
  " (",
  sem$cluster_id,
  ")"
)

DotPlot(
  sem,
  features = genes,
  group.by = "cluster_label",
) +
  RotatedAxis()

dev.off()

##3D plot##

sem_19_1_separado_3D <- RunUMAP(
  sem,
  dims = 1:19,
  n.components = 3,
  min.dist = 0.5
)
library(plotly)
umap <- Embeddings(sem_19_1_separado_3D, "umap")

df <- data.frame(
  UMAP_1 = umap[,1],
  UMAP_2 = umap[,2],
  UMAP_3 = umap[,3],
  experiment= sem$seurat_clusters
)

plot_ly(
  data = df,
  x = ~UMAP_1,
  y = ~UMAP_2,
  z = ~UMAP_3,
  symbol = ~experiment,
  colors = colores_2_claude,
  type = "scatter3d",
  mode = "markers", 
  marker = list(size = 3)
)

df_ident <- data.frame(
  UMAP_1 = umap[,1],
  UMAP_2 = umap[,2],
  UMAP_3 = umap[,3],
  experiment= sem$orig.ident
)

plot_ly(
  data = df_ident,
  x = ~UMAP_1,
  y = ~UMAP_2,
  z = ~UMAP_3,
  symbol = ~experiment,
  colors = sample_cols,
  type = "scatter3d",
  mode = "markers", 
  marker = list(size = 3)
)



DimPlot(sem_19_1_separado_3D,
        reduction = "umap",
        group.by  = "cluster_id_claude",
        cols      = colores_celltype) 


DimPlot(sem_19_1_separado_3D,
        reduction = "umap",
        group.by  = "orig.ident",
        cols      = sample_cols) 



amnion_sig_cluster <- subset(sem, subset = cluster_id_claude %in% 'amnion_signalling')





