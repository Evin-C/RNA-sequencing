# Program for Exploratory Data Analysis, Differential Expression Analysis, and Overrepresentation analysis.

# Before you start - Do Required Installations and or Load Libraries:

# Install DESeq2 if needed:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

# Install pheatmap if needed:
if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

library(pheatmap)

# Get Colour Palette for plots:
# Install RColorBrewer if needed:
if (!require("RColorBrewer", quietly = TRUE))
  install.packages("RColorBrewer")

library(RColorBrewer)

# Install packages and load library:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# Install and load VennDiagram package:
if (!require("VennDiagram", quietly = TRUE))
  install.packages("VennDiagram")
library(VennDiagram)

# Install clusterProfiler and org.Hs.eg.db if needed:
if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

# Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Install DOSE if needed
if (!require("DOSE", quietly = TRUE)) {
  BiocManager::install("DOSE")
  BiocManager::install("enrichplot")
}
# Load DOSE library
library(DOSE)
# Load enrichplot library
library(enrichplot)

### ------------------------------------ Part I: Exploratory Data Analysis ------------------------------------ ###
# Need the "reformatted_counts.txt" file from running script "09_reformat_featureCounts.sh".

# --------- Make Counts Table ---------#

# Define the reformatted counts file (can also be the path to the file):
counts_file <- "reformatted_counts.txt"

# Read the file using read.table:
counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1)

#Change column names to condition names:
colnames(counts) <- gsub("^.*\\.([^.]+)_mappedReads_sorted\\.bam$", "\\1", colnames(counts))

# Inspect the data:
head(counts)


# -------- Make Metadata Table --------#

# Create metadata table:
coldata <- data.frame(
  sampleName = colnames(counts),
  condition = c(
    "HER2", "HER2", "HER2", 
    "NonTNBC", "NonTNBC", "NonTNBC", 
    "Normal", "Normal", "Normal", 
    "TNBC", "TNBC", "TNBC"
  )
)

# Set sampleName as row names (must match column names in counts):
rownames(coldata) <- coldata$sampleName # Get HER21,HER22,HER23,NonTNBC1, etc. as row names.

# Inspect metadata:
print(coldata)


# -- Construct a DESeqDataSet -- #

# Create DESeqDataSet object:
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~condition
)

# Inspect DESeqDataSet:
dds

# -- Do Variance Stabilizing Transformation (vst) -- # 

# Remove the dependence of the variance on the mean by running DESeq2::vst():
vsd <- vst(dds, blind = TRUE)
# Set blind = TRUE to to ignore the experimental design when stabilizing variance
# (useful for visualizations like PCA).

# ---- Create Principal Component Analysis (PCA) ---- #

# Define custom colors for the conditions:
custom_colors <- c(
  "Normal" = "#39C481",
  "HER2" = "#ADEC66",
  "NonTNBC" = "#8DE0E7",
  "TNBC" = "#2499C0"
)

# PCA is an important quality control step to assess clustering of samples based on their gene expression profiles.
# Use the transformed data (vsd) to generate the PCA plot:
pcaData <- plotPCA(vsd, intgroup= "condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size=4) +
  scale_color_manual(values = custom_colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme_bw() #Remove background

# Save the PCA plot
ggsave("PCA_plot.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# ------- Create Heatmap for Sample Distances ------- #

# To explore colour palettes(optional):
# display.brewer.all()

# Choose a color palette:
heatmap_colors <- colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(255)
# Again for the samples:
annotation_colors <- list(
  condition = c(
    "Normal" = "#39C481",
    "HER2" = "#ADEC66",
    "NonTNBC" = "#8DE0E7",
    "TNBC" = "#2499C0"
  )
)

# Compute sample-to-sample distances:
sampleDists <- dist(t(assay(vsd)))  # Transpose matrix so samples are rows
# The function dist() computes distances between rows or columns of a matrix.
# The default method is Euclidean distance.

sampleDistMatrix <- as.matrix(sampleDists)  # Convert to distance matrix

# Annotate the samples with metadata:
annotation <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

# Save the heatmap
png("Heatmap_sample_distance.png", width = 3600, height = 2700, res = 300)

# Generate heatmap:
pheatmap(sampleDistMatrix, 
         annotation_col = annotation, 
         main = "Sample Distance Heatmap", 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists,
         color = heatmap_colors,
         annotation_colors = annotation_colors,
         cellwidth = 40,
         cellheight = 40)
dev.off()

### -------------------------------- Part II: Differential Expression Analysis -------------------------------- ###

# Perform DESeq2 differential expression analysis.
# This will analyze the gene expression differences across the different conditions:
dds <- DESeq(dds)

# Perform Pairwise Comparisons:

# Normal vs. HER2
res_her2 <- results(dds, contrast = c("condition", "HER2", "Normal"))

# Normal vs. NonTNBC
res_nontnbc <- results(dds, contrast = c("condition", "NonTNBC", "Normal"))

# Normal vs. TNBC
res_tnbc <- results(dds, contrast = c("condition", "TNBC", "Normal"))

# Analyze the Results: Extract the number of Differentially Expressed Genes for each contrast:
# (Using padj < 0.05 significance level)
# Number of DE genes with padj < 0.05 for each contrast;
num_DE_genes_her2 <- sum(res_her2$padj < 0.05, na.rm = TRUE)      # Normal vs HER2
num_DE_genes_nontnbc <-sum(res_nontnbc$padj < 0.05, na.rm = TRUE)   # Normal vs NonTNBC
num_DE_genes_tnbc <-sum(res_tnbc$padj < 0.05, na.rm = TRUE)      # Normal vs TNBC

# Identify Upregulated and Downregulated Genes using the log2-fold change:
# There is a "log2FoldChange" value in result variables (res_her2, res_nontnbc, res_tnbc).

# Normal vs. HER2
# Upregulated and Downregulated DE genes for Normal vs HER2:
upreg_her2 <- sum(res_her2$padj < 0.05 & res_her2$log2FoldChange > 0, na.rm = TRUE)             # Upregulated
downreg_her2 <- sum(res_her2$padj < 0.05 & res_her2$log2FoldChange < 0, na.rm = TRUE)           # Downregulated

# Normal vs. NonTNBC
# Upregulated and Downregulated DE genes for Normal vs NonTNBC:
upreg_nontnbc <- sum(res_nontnbc$padj < 0.05 & res_nontnbc$log2FoldChange > 0, na.rm = TRUE)    # Upregulated
downreg_nontnbc <- sum(res_nontnbc$padj < 0.05 & res_nontnbc$log2FoldChange < 0, na.rm = TRUE)  # Downregulated

# Normal vs. TNBC
# Upregulated and Downregulated DE genes for Normal vs TNBC
upreg_tnbc <- sum(res_tnbc$padj < 0.05 & res_tnbc$log2FoldChange > 0, na.rm = TRUE)             # Upregulated
downreg_tnbc <- sum(res_tnbc$padj < 0.05 & res_tnbc$log2FoldChange < 0, na.rm = TRUE)           # Downregulated

# -------------------------------- Investigate Expression levels of genes ------------------------------- #

# Investigate selected genes of particular interest (based on the original publication).
# Investigate their expression level using normalised counts.
# Selected genes: - ENSG00000204628, receptor for activated C kinase 1 (RACK1)
#                 - ENSG00000172551, mucin like 1 (MUCL1)
#                 - ENSG00000172354, G protein subunit beta 2 (GNB2)

# Define genes of interest to variable (= genes that were mentioned in the publication):
genes_of_interest <- c("ENSG00000204628", "ENSG00000172551", "ENSG00000172354")

# Get normalized counts for genes of interest:
normalized_counts <- counts(dds, normalized = TRUE)

# Display normalized counts for the selected genes
# This will provide the expression values for these genes across all samples, adjusted for sequencing depth.
normalized_counts[genes_of_interest, ]

# ------------------------------- Create Volcano Plots for each comparison ------------------------------ #

# ---------------- Volcano Plot for Normal vs HER2 ----------------- #

# Define custom colors for the volcano plot
keyvals_her2 <- ifelse(
  res_her2$padj < 0.05 & res_her2$log2FoldChange > 1, 'red',                  # Upregulated
  ifelse(res_her2$padj < 0.05 & res_her2$log2FoldChange < -1, 'blue', 'grey') # Downregulated or Not Significant
)

# Assign names to the keyvals (for the legend)
names(keyvals_her2)[keyvals_her2 == 'red']  <- 'Upregulated'
names(keyvals_her2)[keyvals_her2 == 'blue'] <- 'Downregulated'
names(keyvals_her2)[keyvals_her2 == 'grey'] <- 'Not Significant'

# Create Volcano Plot using EnhancedVolcano
EnhancedVolcano(res_her2,
                lab = rownames(res_her2),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = '',
                title = 'HER2 vs Normal',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.5,
                legendPosition = 'top',
                colCustom = keyvals_her2,
                drawConnectors = TRUE,
                caption = NULL,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  xlim(c(-15,25)) +
  ylim(c(-5,100))

ggsave("Volcano_Normal_vs_HER2.png", width = 8, height = 6, dpi = 300)
dev.off()

# --------------- Volcano Plot for Normal vs NonTNBC --------------- #

# Define custom colors for the volcano plot
keyvals_nontnbc <- ifelse(
  res_nontnbc$padj < 0.05 & res_nontnbc$log2FoldChange > 1, 'red',                  # Upregulated
  ifelse(res_nontnbc$padj < 0.05 & res_nontnbc$log2FoldChange < -1, 'blue', 'grey') # Downregulated or Not Significant
)

# Assign names to the keyvals (for the legend)
names(keyvals_nontnbc)[keyvals_nontnbc == 'red']  <- 'Upregulated'
names(keyvals_nontnbc)[keyvals_nontnbc == 'blue'] <- 'Downregulated'
names(keyvals_nontnbc)[keyvals_nontnbc == 'grey'] <- 'Not Significant'

# Create Volcano Plot using EnhancedVolcano
EnhancedVolcano(res_nontnbc,
                lab = rownames(res_nontnbc),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = '',
                title = 'NonTNBC vs Normal',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.5,
                legendPosition = 'top',
                colCustom = keyvals_nontnbc,
                drawConnectors = TRUE,
                caption = NULL,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  xlim(c(-15,25)) +
  ylim(c(-5,100))

ggsave("Volcano_Normal_vs_NonTNBC.png", width = 8, height = 6, dpi = 300)
dev.off()

# ---------------- Volcano Plot for Normal vs TNBC ---------------- #

# Define custom colors for the volcano plot
keyvals_tnbc <- ifelse(
  res_tnbc$padj < 0.05 & res_tnbc$log2FoldChange > 1, 'red',                  # Upregulated
  ifelse(res_tnbc$padj < 0.05 & res_tnbc$log2FoldChange < -1, 'blue', 'grey') # Downregulated or Not Significant
)

# Assign names to the keyvals (for the legend)
names(keyvals_tnbc)[keyvals_tnbc == 'red']  <- 'Upregulated'
names(keyvals_tnbc)[keyvals_tnbc == 'blue'] <- 'Downregulated'
names(keyvals_tnbc)[keyvals_tnbc == 'grey'] <- 'Not Significant'

# Create Volcano Plot using EnhancedVolcano
EnhancedVolcano(res_tnbc,
                lab = rownames(res_tnbc),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = '',
                title = 'TNBC vs Normal',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.5,
                legendPosition = 'top',
                colCustom = keyvals_tnbc,
                drawConnectors = TRUE,
                caption = NULL,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  xlim(c(-15,25)) +
  ylim(c(-5,100))

ggsave("Volcano_Normal_vs_TNBC.png", width = 8, height = 6, dpi = 300)
dev.off()

# ----- Create a heatmap to visualize the expression patterns of the selected genes across all conditions. ----- #

# Save the heatmap
png("Heatmap_gene_expression.png", width = 2700, height = 1000, res = 150)

# Generate heatmap:
pheatmap(normalized_counts[genes_of_interest, ], 
         annotation_col = annotation, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         color = heatmap_colors,
         scale = "row",
         main = "Gene Expression Pattern Heatmap",
         annotation_colors = annotation_colors,
         cellwidth = 80,
         cellheight = 100)
dev.off()

# -- Create Venn Diagram of all sets of differentially expressed (DE) genes across all conditions -- #

# Extract DE genes from results filtered for padj < 0.5:
de_genes_her2 <- rownames(res_her2[!is.na(res_her2$padj) & res_her2$padj < 0.05, ])
de_genes_nontnbc <- rownames(res_nontnbc[!is.na(res_nontnbc$padj) & res_nontnbc$padj < 0.05, ])
de_genes_tnbc <- rownames(res_tnbc[!is.na(res_tnbc$padj) & res_tnbc$padj < 0.05, ])
# Ensure sets are unique
de_genes_her2 <- unique(de_genes_her2)
de_genes_nontnbc <- unique(de_genes_nontnbc)
de_genes_tnbc <- unique(de_genes_tnbc)

# Create the Venn diagram:

# Save the Venn diagram:
png("VennDiagram_DE_genes.png", width = 1000, height = 1000, res = 150)
venn.plot <- draw.triple.venn(
  area1 = length(de_genes_her2),           # Number of DE genes in HER2
  area2 = length(de_genes_nontnbc),        # Number of DE genes in NonTNBC
  area3 = length(de_genes_tnbc),           # Number of DE genes in TNBC
  n12 = length(intersect(de_genes_her2, de_genes_nontnbc)),        # HER2 & NonTNBC
  n13 = length(intersect(de_genes_her2, de_genes_tnbc)),           # HER2 & TNBC
  n23 = length(intersect(de_genes_nontnbc, de_genes_tnbc)),        # NonTNBC & TNBC
  n123 = length(Reduce(intersect, list(de_genes_her2, de_genes_nontnbc, de_genes_tnbc))), # All 3 conditions
  category = c("HER2", "NonTNBC", "TNBC"), # Category names
  fill = c("#ADEC66","#8DE0E7","#2499C0"), # Colors for the sets
  alpha = 0.75,                            # Transparency level
  cex = 2,                                 # Font size for numbers
  cat.cex = 2,                             # Font size for labels
  cat.pos = c(0, 0, 0),                    # Adjust label positions
  cat.dist = c(0.04, 0.04, -0.44),         # Adjust label distances
  euler.d = TRUE,                          # Enable Euler diagram (proportional sizing)
  scaled = TRUE,                           # Enable proportional scaling
  cat.fontfamily = "sans",                 # Change font of labels to default R font (=TT Arial)
  fontfamily = "sans"                      # Change font of numbers to default R font (=TT Arial)
)
dev.off()

### ---------------------------------- Part III: Overrepresentation Analysis ---------------------------------- ###

# --- Extract Ensemble IDs of DE genes (padj < 0.05) for each comparison --- #

#Normal vs. HER2, exclude NA values:
de_genes_her2 <- rownames(res_her2[!is.na(res_her2$padj) & res_her2$padj < 0.05, ])
  # Ensure the IDs are unique
de_genes_her2 <- unique(de_genes_her2)

#Normal vs. NonTNBC, exclude NA values:
de_genes_nontnbc <- rownames(res_nontnbc[!is.na(res_nontnbc$padj) & res_nontnbc$padj < 0.05, ])
  # Ensure the IDs are unique
de_genes_nontnbc <- unique(de_genes_nontnbc)

#Normal vs. TNBC, exclude NA values:
de_genes_tnbc <- rownames(res_tnbc[!is.na(res_tnbc$padj) & res_tnbc$padj < 0.05, ])
# Ensure the IDs are unique
de_genes_tnbc <- unique(de_genes_tnbc)

# ---------- Define universe genes (all measured genes in analysis, =rownames of counts) ---------- #

universe_genes <- rownames(counts)

# ------- Run Overrepresentation Analysis (GO enrichment analysis) ------- #

# --- Normal vs. HER2 --- #
go_results_her2 <- enrichGO(
  gene         = de_genes_her2,      # DE genes
  universe     = universe_genes,     # All measured genes
  OrgDb        = org.Hs.eg.db,       # Annotation database for humans
  keyType      = "ENSEMBL",          # Format of gene IDs
  ont          = "BP",               # GO subontology: "BP" (Biological Process)
  pAdjustMethod = "BH",              # Multiple testing correction method
  pvalueCutoff  = 0.01,              # p-value cutoff
  qvalueCutoff = 0.05,               # q-value cutoff
  readable     = TRUE                # Convert Ensembl IDs to gene symbols
)

# View GO terms
head(go_results_her2)

# Save GO results to a TSV file (can be opened in e.g. Microsoft Excel):
write.table(as.data.frame(go_results_her2), file = "GO_enrichment_her2.tsv", sep = "\t", row.names = FALSE, quote = TRUE)

# --- Normal vs. NonTNBC --- #
go_results_nontnbc <- enrichGO(
  gene         = de_genes_nontnbc,   # DE genes
  universe     = universe_genes,     # All measured genes
  OrgDb        = org.Hs.eg.db,       # Annotation database for humans
  keyType      = "ENSEMBL",          # Format of gene IDs
  ont          = "BP",               # GO subontology: "BP" (Biological Process)
  pAdjustMethod = "BH",              # Multiple testing correction method
  pvalueCutoff  = 0.01,              # p-value cutoff
  qvalueCutoff = 0.05,               # q-value cutoff
  readable     = TRUE                # Convert Ensembl IDs to gene symbols
)

# View GO terms
head(go_results_nontnbc)

# Save GO results to a TSV file (can be opened in e.g. Microsoft Excel):
write.table(as.data.frame(go_results_nontnbc), file = "GO_enrichment_nontnbc.tsv", sep = "\t", row.names = FALSE, quote = TRUE)

# --- Normal vs. TNBC --- #
go_results_tnbc <- enrichGO(
  gene         = de_genes_tnbc,      # DE genes
  universe     = universe_genes,     # All measured genes
  OrgDb        = org.Hs.eg.db,       # Annotation database for humans
  keyType      = "ENSEMBL",          # Format of gene IDs
  ont          = "BP",               # GO subontology: "BP" (Biological Process)
  pAdjustMethod = "BH",              # Multiple testing correction method: "BH" (Benjamini and Hochberg method)
  pvalueCutoff  = 0.01,              # p-value cutoff
  qvalueCutoff = 0.05,               # q-value cutoff
  readable     = TRUE                # Convert Ensembl IDs to gene symbols
)

# View GO terms
head(go_results_tnbc)

# Save GO results to a TSV file (can be opened in e.g. Microsoft Excel):
write.table(as.data.frame(go_results_tnbc), file = "GO_enrichment_tnbc.tsv", sep = "\t", row.names = FALSE, quote = TRUE)

# ------------- Visualize the results ------------- #

# ----- Barplot for Normal vs HER2 ----- #

# Ensure the GO results are in a format that can be used by DOSE
go_results_her2_df <- as.data.frame(go_results_her2)
#Order by 'Count' (descending order) to get the most enriched pathways
go_results_her2_df <- go_results_her2_df[order(go_results_her2_df$Count, decreasing = TRUE), ]
# Retain only the top 10 pathways
top_10_pathways_her2 <- go_results_her2_df[1:10, ]
# Recreate an 'enrichResult' object containing only the top 10 pathways
go_results_her2_top10 <- go_results_her2
go_results_her2_top10@result <- go_results_her2_top10@result[rownames(top_10_pathways_her2), ]

# Create barplot
# Save the barplot to a PNG file
png("GO_enrichment_her2_barplot.png", width = 2400, height = 1800, res = 300)
barplot(go_results_her2_top10, showCategory = 10, title = "Enriched GO Pathways in HER2 vs Normal") + 
  scale_fill_gradient(low = "#d8120e", high = "#ffd926") +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5),   # Adjust spacing of y-axis text
        plot.margin = margin(1, 1, 1, 1.5, "cm"))
 # showCategory: controls how many GO terms to display
dev.off()

# ----- Barplot for Normal vs NonTNBC ----- #

# Ensure the GO results are in a format that can be used by DOSE
go_results_nontnbc_df <- as.data.frame(go_results_nontnbc)
#Order by 'Count' (descending order) to get the most enriched pathways
go_results_nontnbc_df <- go_results_nontnbc_df[order(go_results_nontnbc_df$Count, decreasing = TRUE), ]
# Retain only the top 10 pathways
top_10_pathways_nontnbc <- go_results_nontnbc_df[1:10, ]
# Recreate an 'enrichResult' object containing only the top 10 pathways
go_results_nontnbc_top10 <- go_results_nontnbc
go_results_nontnbc_top10@result <- go_results_nontnbc_top10@result[rownames(top_10_pathways_nontnbc), ]

# Create barplot
# Save the barplot to a PNG file
png("GO_enrichment_nontnbc_barplot.png", width = 2400, height = 1800, res = 300)
barplot(go_results_nontnbc_top10, showCategory = 10, title = "Enriched GO Pathways in NonTNBC vs Normal") + 
  scale_fill_gradient(low = "#d8120e", high = "#ffd926") +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5),   # Adjust spacing of y-axis text
        plot.margin = margin(1, 1, 1, 1.5, "cm"))
# showCategory: controls how many GO terms to display
dev.off()

# ----- Barplot for Normal vs TNBC ----- #

# Ensure the GO results are in a format that can be used by DOSE
go_results_tnbc_df <- as.data.frame(go_results_tnbc)
# Order by 'Count' (descending order) to get the most enriched pathways
go_results_tnbc_df <- go_results_tnbc_df[order(go_results_tnbc_df$Count, decreasing = TRUE), ]
# Retain only the top 10 pathways
top_10_pathways_tnbc <- go_results_tnbc_df[1:10, ]
# Recreate an 'enrichResult' object containing only the top 10 pathways
go_results_tnbc_top10 <- go_results_tnbc
go_results_tnbc_top10@result <- go_results_tnbc_top10@result[rownames(top_10_pathways_tnbc), ]

# Create barplot
# Save the barplot to a PNG file
png("GO_enrichment_tnbc_barplot.png", width = 2400, height = 1800, res = 300)
barplot(go_results_tnbc_top10, showCategory = 10, title = "Enriched GO Pathways in TNBC vs Normal") + 
  scale_fill_gradient(low = "#d8120e", high = "#ffd926") +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5),   # Adjust spacing of y-axis text
        plot.margin = margin(1, 1, 1, 1.5, "cm"))
# showCategory: controls how many GO terms to display
dev.off()

# ----- Dotplot for Normal vs HER2 ----- #

# Create a dotplot
# Save the dotplot to a PNG file
png("GO_enrichment_her2_dotplot.png", width = 2400, height = 1800, res = 300)
dot_plot <- dotplot(
  go_results_her2_top10, 
  showCategory = 10,        # Show top 10 GO terms
  x = "GeneRatio",          # Use GeneRatio on the x-axis
  color = "p.adjust",       # Color points by adjusted p-value
  title = "Enriched GO Pathways in HER2 vs Normal"
)
dot_plot + 
  scale_fill_gradient(low = "#d8120e", high = "#ffd926") # Custom color gradient
dev.off()

# ----- Dotplot for Normal vs NonTNBC ----- #

# Create a dotplot
# Save the dotplot to a PNG file
png("GO_enrichment_nontnbc_dotplot.png", width = 2400, height = 1800, res = 300)
dot_plot <- dotplot(
  go_results_nontnbc_top10, 
  showCategory = 10,        # Show top 10 GO terms
  x = "GeneRatio",          # Use GeneRatio on the x-axis
  color = "p.adjust",       # Color points by adjusted p-value
  title = "Enriched GO Pathways in NonTNBC vs Normal"
)
dot_plot + 
  scale_fill_gradient(low = "#d8120e", high = "#ffd926") # Custom color gradient
dev.off()

# ----- Dotplot for Normal vs TNBC ----- #

# Create a dotplot
# Save the dotplot to a PNG file
png("GO_enrichment_tnbc_dotplot.png", width = 2400, height = 1800, res = 300)
dot_plot <- dotplot(
  go_results_tnbc_top10, 
  showCategory = 10,        # Show top 10 GO terms
  x = "GeneRatio",          # Use GeneRatio on the x-axis
  color = "p.adjust",       # Color points by adjusted p-value
  title = "Enriched GO Pathways in TNBC vs Normal"
)
dot_plot + 
  scale_fill_gradient(low = "#d8120e", high = "#ffd926") # Custom color gradient
dev.off()

# -- Heatmap of gene expression of the top two genes from the top 3 GO pathways for all comparisons -- #
# -- Assumption: There will be an overlap of genes between the comparisons, so there will be less than 2*3*3 genes. -- #

# Function to extract two Ensembl Gene IDs for each of the top 3 GO terms
extract_top_genes <- function(go_results) {
  go_results_df <- as.data.frame(go_results)
  go_results_df <- go_results_df[order(go_results_df$Count, decreasing = TRUE), ]
  
  # Extract the top 3 GO terms
  top_3_pathways <- go_results_df[1:3, ]
  
  # Extract Ensembl Gene IDs for the top 3 GO terms
  top_3_go_terms <- rownames(top_3_pathways)
  selected_genes <- c()  # Create an empty vector to store genes
  
  for (go_term in top_3_go_terms) {
    # Extract the gene symbols for this GO term
    gene_symbols <- go_results@result[go_term, "geneID"]
    
    # Split gene symbols if there are multiple genes in one string
    gene_symbols <- unlist(strsplit(gene_symbols, "/"))
    
    # Convert gene symbols to Ensembl Gene IDs with org.Hs.eg.db
    ensembl_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
    
    # Select only the first two Ensembl IDs
    if (length(ensembl_ids$ENSEMBL) > 1) {
      selected_genes <- c(selected_genes, ensembl_ids$ENSEMBL[1:2])  # Add the first two genes
    } else {
      selected_genes <- c(selected_genes, ensembl_ids$ENSEMBL)  # If only one gene, add it
    }
  }
  
  return(selected_genes)
}

# Extract Ensembl Gene IDs for HER2
selected_genes_her2 <- extract_top_genes(go_results_her2)
# Extract Ensembl Gene IDs for TNBC
selected_genes_tnbc <- extract_top_genes(go_results_tnbc)
# Extract Ensembl Gene IDs for Non-TNBC
selected_genes_nontnbc <- extract_top_genes(go_results_nontnbc)

# Combine all selected genes into a single list
all_top_genes_ensembl <- c(selected_genes_her2, selected_genes_tnbc, selected_genes_nontnbc)

# Remove duplicate Ensembl Gene IDs
all_top_genes_ensembl <- unique(all_top_genes_ensembl)
all_top_genes_ensembl <- intersect(all_top_genes_ensembl, rownames(normalized_counts))

# Save the heatmap
png("Heatmap_top_genes_expression.png", width = 3500, height = 2700, res = 300)
pheatmap(normalized_counts[all_top_genes_ensembl, ], 
         annotation_col = annotation, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         color = heatmap_colors,
         scale = "row",
         main = "Gene Expression Pattern Heatmap",
         annotation_colors = annotation_colors,
         cellwidth = 40,
         cellheight = 45)
dev.off()
