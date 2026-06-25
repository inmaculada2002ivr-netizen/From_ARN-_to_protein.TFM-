(library(SummarizedExperiment))
library(tximeta) 
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)

## Análisis de expresión diferencial y enriquecimiento 

### Análisis de la calidad de los mapeadores
## Read gene expression for sample1, Elena1.tsv
WT1 <- read.table(file="data/Elena1.tsv",header=T,sep="\t")
gene.expression.WT1 <- WT1$TPM
names(gene.expression.WT1) <- WT1$Gene.ID

## Read gene expression for sample2, Elena_rep2.tsv
BMP_0_48_1 <- read.table(file="data/Elena2.tsv",header=T,sep="\t")
gene.expression.BMP_0_48_1 <- BMP_0_48_1$TPM
names(gene.expression.BMP_0_48_1) <- BMP_0_48_1$Gene.ID

## Read gene expression for sample3, Elena3.tsv
WT2 <- read.table(file="data/Elena3.tsv",header=T,sep="\t")
gene.expression.WT2 <- WT2$TPM
names(gene.expression.WT2) <- WT2$Gene.ID

## Read gene expression for sample4, Elena4.tsv
BMP_0_48_2 <- read.table(file="data/Elena4.tsv",header=T,sep="\t")
gene.expression.BMP_0_48_2 <- BMP_0_48_2$TPM
names(gene.expression.BMP_0_48_2) <- BMP_0_48_2$Gene.ID

## To ensure that genes are in the same order across all samples, we extract
## the gene IDs from each gene expression data frame and obtain a unique list
## containing all genes present in at least one of the samples
gene.ids <- unique(c(
  WT1$Gene.ID,
  WT2$Gene.ID,
  BMP_0_48_1$Gene.ID,
  BMP_0_48_2$Gene.ID
))

# Crear un data frame vacío con gene_id
expr_matrix <- data.frame(gene.ids = gene.ids)

# Añadir columnas para cada muestra
expr_matrix$WT1 <- gene.expression.WT1[gene.ids]
expr_matrix$WT2 <- gene.expression.WT2[gene.ids]
expr_matrix$BMP_0_48_1 <- gene.expression.BMP_0_48_1[gene.ids]
expr_matrix$BMP_0_48_2 <- gene.expression.BMP_0_48_2[gene.ids]

# Reemplazar NAs por 0 (genes no detectados en alguna muestra)
expr_matrix[is.na(expr_matrix)] <- 0

prefix <- "C:/Users/inmav/Desktop/TFM MADOBIS/dato/bulk/data"
fn = 'quant.sf'
files <- c(file.path(prefix,'transcripts_quant_Elena1',fn),
           file.path(prefix,'transcripts_quant_Elena2',fn),
           file.path(prefix,'transcripts_quant_Elena3',fn),
           file.path(prefix,'transcripts_quant_Elena4',fn))

# Create a data frame with file paths and sample names
coldata <- data.frame(
  files = files,  # Column with paths to quantification files
  names= c("WT1", "BMP_0_48_1", "WT2", "BMP_0_48_2"), # Column with sample names
  stringsAsFactors = FALSE # Keep strings as characters, not factors
) 
# Import transcript-level quantification and attach metadata using tximeta
se <- tximeta(coldata)  # Creates a SummarizedExperiment object
# Summarize transcript-level data to gene-level counts/abundances
gse <- summarizeToGene(se)  # Now each row corresponds to a gene
# Extract TPM (Transcripts Per Million) values from the gene-level object
tpm <- assays(gse)$abundance  # Matrix with genes as rows and samples as columns

# Combine gene IDs with TPM values into a single data frame
gene_tpm <- data.frame(
  gene.ids = rownames(gse),  # Gene identifiers
  tpm  # TPM values for each sample
)
##safe the gene id from the salmon quantification in a variable. 
gene.ids_SALMON <- gene_tpm$gene.ids


length(gene.ids)

length(gene.ids_SALMON)

genes_comunes <-intersect(gene.ids,gene.ids_SALMON)

no_comunes<-setdiff(gene.ids,genes_comunes)

length(no_comunes)
expresion_no_comunes <- expr_matrix[expr_matrix$gene.ids%in% no_comunes, ]

write.csv(expresion_no_comunes, "expresión_no_comunes.csv")



draw.pairwise.venn(area1 = length(gene.ids),
                   area2 = length(gene.ids_SALMON),
                   cross.area =  
                     length(genes_comunes),lwd = 3,
                   category = c("STAR","SALMON"),
                   euler.d = TRUE ,scaled = TRUE, col = c("gold2","palevioletred2"),
                   fill = c("gold2","palevioletred2"),alpha = 0.5,cex = 1.2,cat.cex = 1.3, cat.pos = c(-20,20), cat.dist = c(0.05,0.05))


# Creamos un data frame con todos los genes
all_genes <- unique(c(gene.ids, gene.ids_SALMON))

data <- data.frame(
  gene = all_genes,
  STAR = all_genes %in% gene.ids,
  Salmon = all_genes %in% gene.ids_SALMON
)

# Generamos el UpSet plot
upset(
  data,
  intersect = c("STAR", "Salmon"),
  
  # Barras de intersección con números
  base_annotations = list(
    "Intersection size" = intersection_size(
      text = list(
        vjust = -0.5,
        size = 5,
        color = "black"
      ),
      fill = c( "palevioletred2","goldenrod2", "tomato3")
    )
  ),
  
  # Barras laterales (tamaño de cada set)
  set_sizes = upset_set_size(
    geom = geom_bar(fill = c("goldenrod2","palevioletred2"))
  )
) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold")
  )


genes_comunes <- intersect(expr_matrix$gene.ids, gene_tpm$gene.ids)


star_sub <- expr_matrix[expr_matrix$gene.ids %in% genes_comunes, ]
star_sub <- star_sub[match(genes_comunes, star_sub$gene.ids), ]
rownames(star_sub) <- star_sub$gene.ids
star_sub$gene.ids <- NULL
colnames(star_sub) <- paste0(colnames(star_sub), "_STAR")

salmon_sub <- gene_tpm[gene_tpm$gene.ids %in% genes_comunes, ]
salmon_sub <- salmon_sub[match(genes_comunes, salmon_sub$gene.ids), ]
rownames(salmon_sub) <- salmon_sub$gene.ids
salmon_sub$gene.ids <- NULL
colnames(salmon_sub) <- paste0(colnames(salmon_sub), "_Salmon")


full_matrix <- cbind(star_sub, salmon_sub)


log_matrix <- log2(full_matrix + 1)

genes_con_varianza <- apply(log_matrix, 1, var) > 0


log_matrix_filtrada <- log_matrix[genes_con_varianza, ]


cat("Genes iniciales:", nrow(log_matrix), 
    "\nGenes tras eliminar varianza cero:", nrow(log_matrix_filtrada))



pca_res <- prcomp(t(log_matrix_filtrada), scale. = TRUE)



library(ggplot2)

pca_data <- data.frame(
  Sample = colnames(log_matrix_filtrada),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  stringsAsFactors = FALSE)

pca_data$Metodo <- ifelse(grepl("_STAR", pca_data$Sample), "STAR", "Salmon")

pca_data$Condicion[grepl("WT", pca_data$Sample)] <- "mTeser"
pca_data$Condicion[grepl("BMP", pca_data$Sample)] <- "BMP4(10 ng/mL)"


var_pc1 <- round(100 * summary(pca_res)$importance[2,1], 1)
var_pc2 <- round(100 * summary(pca_res)$importance[2,2], 1)


ggplot(pca_data, aes(x = PC1, y = PC2, color = Metodo, shape = Condicion)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = Sample), vjust = -1.2, size = 3, show.legend = FALSE) +
  theme_bw() +theme_bw() +
  theme(
    panel.background = element_rect(fill = "grey90", color = NA)
  )+
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.2)
  )+geom_text(
    aes(label = Sample),
    inherit.aes = FALSE,
    x = pca_data$PC1,
    y = pca_data$PC2,
    color = "black",
    fontface = "bold",
    vjust = -1.2,
    size = 3
  )
  scale_color_manual(values = c(
    "STAR" = "goldenrod2",
    "Salmon" = "palevioletred1"
  )) +
  scale_shape_manual(values = c(
    "mTeser" = 16,
    "BMP4(10 ng/mL)" = 18
  )) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15))+
  labs(
    title = "PCA: Comparison of STAR vs Salmon Mappers",
    subtitle = paste("Based on", nrow(log_matrix_filtrada), "genes with variance > 0"),
    x = paste0("PC1 (", var_pc1, "%)"),
    y = paste0("PC2 (", var_pc2, "%)")
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )
# Scatter plot comparing STAR vs Salmon TPMs
# Merge STAR and Salmon by gene_id
common_matrix <- merge(expr_matrix,  gene_tpm, by="gene.ids", suffixes=c("_STAR","_SALMON"))

samples <- c("WT1","WT2","BMP_0_48_1","BMP_0_48_2") #merge with the same varibale names 
par(mfrow=c(2,2))

#buble for ,representing the the diffenrets scratter plots 

for(s in samples){
  x <- common_matrix[[paste0(s,"_STAR")]]
  y <- common_matrix[[paste0(s,"_SALMON")]]
  
  
  p<-plot(
    log2(x+1), log2(y+1),
    xlab = paste(s, "STAR log2(TPM+1)"),
    ylab = paste(s, "Salmon log2(TPM+1)"),
    pch = 16, col=rgb(0,0,1,0.3),
    main = paste(s, "STAR vs Salmon")
  )
  abline(0,1,col="red")
  cor_val <- cor(x,y,method="spearman")
  legend("topleft", legend=paste("Spearman:", round(cor_val,2)), bty="n")
  print(p)
}
#load the cound matrix (prepDE.py)
data <- read.csv("data/gene_count_matrix.csv", header = T, row.names = "gene_id")
data <- data[,sort(colnames(data))] #sort columns alphabetically 
head(data)
colSums(data)
condition <- c("mTeser","BMP4(10 ng/mL)","mTeser","BMP4(10 ng/mL)") # Create a data frame with sample condition information
my_colData <- as.data.frame(condition)# Make sure row names match column names of the count matrix
rownames(my_colData) <- colnames(data)
my_colData
# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = data, #raw counts matrix
                              colData = my_colData, #sample information 
                              design = ~condition) #experimental desig 
# Run differential expression analysis
dds <- DESeq(dds) # Fit model and estimate dispersion
res <- results(dds) # Extract results
vsd <- vst(dds, blind = FALSE)
# Get normalized counts
#the normalization in DEseq2 in mediated by the median of rations, counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene. 
#fuente: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

normalized_counts <- counts(dds, normalized = T)
head(normalized_counts)

normalized_counts[is.na(normalized_counts)] <- 0

#to merged the raw data and the normalized data in a csv
matriz_raw_normalizado_STAR <- cbind(expr_matrix, normalized_counts)
matriz_raw_normalizado_STAR <- matriz_raw_normalizado_STAR[, colnames(matriz_raw_normalizado_STAR) != "gene.ids"]
colnames(matriz_raw_normalizado_STAR) <- c("WT1_raw", "WT2_raw", "BMP_0_48_1_raw", "BMP_0_48_2_raw", "WT1", "WT2", "BMP_0_48_1", "BMP_0_48_2" )

write.csv(matriz_raw_normalizado_STAR, "matriz_raw_normalizado_STAR.csv")

# Heatmap of sample distances

# Compute Euclidean distances between samples
sampleDists <- dist(t(assay(vsd)))

# Convert to matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Set row names with sample conditions
rownames(sampleDistMatrix) <- vsd$condition

# Create color palette for heatmap
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

# Draw heatmap of sample distances

#install.packages("pheatmap")



pheatmap::pheatmap(sampleDistMatrix,
                   clustering_distance_rows = sampleDists,
                   clustering_distance_cols = sampleDists,
                   col = colors)

# Function for heatmap of most variable genes
variable_gene_heatmap <- function(vsd.obj, num_genes = 500, title = "") {
  # Red-blue color palette
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))
  mr <- ramp(256)[256:1]
  # Define annotation colors
  ann_colors <- list(
    condition = c(
      "BMP4(10 ng/mL)" = "turquoise2",
      "mTeser" = "deeppink1"
    )
  )
  
  # Extract VST (variance-stabilized) counts
  stabilized_counts <- assay(vsd.obj)
  
  # Compute variance per gene
  row_variances <- matrixStats::rowVars(stabilized_counts)
  
  # Select most variable genes
  top_variable_genes <- stabilized_counts[
    order(row_variances, decreasing = TRUE)[1:num_genes], ]
  
  # Center values by gene (mean = 0)
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm = TRUE)
  
  # Sample metadata
  coldata <- as.data.frame(colData(vsd.obj))
  coldata <- coldata[, "condition", drop = FALSE]  # Keep only 'condition'
  
  # Heatmap
  pheatmap::pheatmap(top_variable_genes,
                     color = mr,
                     annotation_col = coldata,
                     annotation_colors = ann_colors,
                     fontsize_col = 8,
                     border_color = NA,
                     cellheight = 12,
                     cellwidth = 20,
                     main = title)
}
variable_gene_heatmap(vsd, num_genes = 40, title = "Top variable genes")

plot_PCA <- function(vsd.obj, group_variable = "condition") {
  
  # PCA desde DESeq2
  pcaData <- plotPCA(vsd.obj,
                     intgroup = group_variable,
                     returnData = TRUE)
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Añadir nombres de muestras
  pcaData$sample <- rownames(pcaData)
  
  # Definir colores
  colors <- c(
    "BMP4(10 ng/mL)" = "turquoise2",
    "mTeser" = "deeppink1"
  )
  
  # Plot PCA con ggplot2 y elipses
  ggplot(pcaData, aes(x = PC1, y = PC2, color = .data[[group_variable]])) +
    geom_point(size = 4)  +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    geom_text_repel(aes(label = sample), color = "black") +
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      title = paste("PCA colored by", group_variable)
    ) +
    theme_classic()
}

# Llamada a la función
plot_PCA(vsd, "condition")


res <- results(dds, contrast = c("condition", "BMP4(10 ng/mL)", "mTeser"))
res_df <- as.data.frame(res)
res_df$Gene.name <- rownames(res_df)
head(res_df)

plot_volcano2 <- function(res,
                          padj_cutoff = 0.01,
                          nlabel = 10,
                          label.by = "padj",
                          gene_column = NULL) {
  
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  # Filtrar NA
  res <- res %>% filter(!is.na(padj))
  
  # Crear columna de tipo de regulación
  res <- res %>% mutate(regulation = case_when(
    padj < padj_cutoff & log2FoldChange > 0 ~ "Activated",
    padj < padj_cutoff & log2FoldChange < 0 ~ "Repressed",
    TRUE ~ "NS"  # No significativo
  ))
  
  
  # Seleccionar genes top para etiquetar
  top_genes <- res %>% filter(regulation != "NS")
  
  if(label.by == "padj"){
    top_genes <- top_genes %>% arrange(padj) %>% head(nlabel)
  } else if(label.by == "log2FoldChange"){
    top_genes <- top_genes %>% arrange(desc(abs(log2FoldChange))) %>% head(nlabel)
  } else {
    stop("label.by debe ser 'padj' o 'log2FoldChange'")
  }
  
  # Colores
  colors <- c("Activated" = "red", "Repressed" = "blue", "NS" = "grey50")
  
  # Volcano plot
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ggrepel::geom_text_repel(data = top_genes,
                             aes(label = .data[[gene_column]]),
                             size = 3,
                             max.overlaps = Inf) +
    labs(x = "Log2 Fold Change",
         y = "-log10(adjusted p-value)",
         title = paste("Volcano Plot (padj <", padj_cutoff, ")")) +
    xlim(c(-10,10)) +
    ylim(c(0,500)) +
    theme_minimal()
}

plot_volcano2(res_df,
              padj_cutoff = 0.01,
              nlabel = 15,
              label.by = "padj",
              gene_column = "Gene.name")