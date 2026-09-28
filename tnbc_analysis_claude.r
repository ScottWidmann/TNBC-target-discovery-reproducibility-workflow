# TCGA-BRCA TNBC Target Discovery Analysis Pipeline
# Author: RNA-seq Analysis for Drug Discovery
# Description: Comprehensive analysis of TCGA-BRCA data to identify TNBC drug targets

# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(ggrepel)
library(VennDiagram)
library(survival)
library(survminer)
library(biomaRt)

# Set working directory and create output folders
#setwd("~/TNBC_Analysis")
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# =====================================================
# PART 1: DATA ACQUISITION FROM TCGA
# =====================================================

# Query TCGA-BRCA RNA-seq data with smaller chunks
query_exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# Robust download with error handling and retry logic
print("Downloading TCGA-BRCA expression data...")
print("This is a large dataset (~5GB) - download may take 15-30 minutes")

download_success <- FALSE
max_attempts <- 3
attempt <- 1

while(!download_success && attempt <= max_attempts) {
  print(paste("Download attempt", attempt, "of", max_attempts))
  
  tryCatch({
    # Try downloading with smaller chunks
    GDCdownload(query_exp, method = "api", files.per.chunk = 50)
    download_success <- TRUE
    print("Download completed successfully!")
  }, error = function(e) {
    print(paste("Download attempt", attempt, "failed:", e$message))
    
    if(attempt < max_attempts) {
      print("Waiting 30 seconds before retry...")
      Sys.sleep(30)
    } else {
      print("All download attempts failed. Trying alternative approach...")
    }
    attempt <<- attempt + 1
  })
}

# If download still fails, try alternative approach
if(!download_success) {
  print("Trying alternative download method...")
  
  # Alternative 1: Download only a subset first
  tryCatch({
    # Query for fewer samples to test
    query_exp_small <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling", 
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal"),
      barcode = getResults(query_exp)$cases[1:100]  # Only first 100 samples
    )
    
    GDCdownload(query_exp_small, method = "api", files.per.chunk = 20)
    query_exp <- query_exp_small  # Use smaller dataset
    download_success <- TRUE
    print("Downloaded subset of data successfully!")
    
  }, error = function(e) {
    print("Subset download also failed. Consider these alternatives:")
    print("1. Try running the download at a different time (less network traffic)")
    print("2. Use a different internet connection")
    print("3. Download manually from GDC portal: https://portal.gdc.cancer.gov/")
    print("4. Use pre-processed data from other sources")
    stop("Unable to download TCGA data. Please try manual download.")
  })
}

# Prepare the data
print("Preparing expression data...")
tryCatch({
  exp_data <- GDCprepare(query_exp)
  print("Data preparation completed!")
}, error = function(e) {
  print(paste("Data preparation failed:", e$message))
  
  # Alternative: try to load if files exist
  print("Checking for existing downloaded files...")
  if(dir.exists("GDCdata")) {
    print("Found existing GDC data directory. Attempting to load...")
    exp_data <- GDCprepare(query_exp)
  } else {
    stop("Data preparation failed and no existing files found.")
  }
})

# Query clinical data
query_clinical <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

# Download clinical data
GDCdownload(query_clinical)
clinical_data <- GDCprepare(query_clinical)

# =====================================================
# PART 2: DATA PREPROCESSING AND SAMPLE CLASSIFICATION
# =====================================================

# Extract expression matrix
exp_matrix <- assay(exp_data, "unstranded")
gene_info <- rowData(exp_data)
sample_info <- colData(exp_data)

# Get clinical information for receptor status
clinical_patient <- clinical_data$clinical_patient_brca
clinical_followup <- clinical_data$clinical_follow_up_v1.5_nte_brca

# Create comprehensive sample annotation
sample_annotation <- data.frame(
  barcode = sample_info$barcode,
  patient = sample_info$patient,
  sample_type = sample_info$sample_type,
  tissue_type = sample_info$tissue_or_organ_of_origin,
  stringsAsFactors = FALSE
)

# Merge with clinical data for receptor status
# Extract patient IDs from barcodes
sample_annotation$patient_id <- substr(sample_annotation$patient, 1, 12)

# Convert clinical_patient to regular data frame to avoid tibble issues
clinical_patient <- as.data.frame(clinical_patient)

# First, let's examine the clinical data structure
print("Available columns in clinical_patient:")
print(colnames(clinical_patient))

# Check for receptor status columns with different names
receptor_cols <- colnames(clinical_patient)[grepl("er_|pr_|her2|receptor", colnames(clinical_patient), ignore.case = TRUE)]
print("Potential receptor status columns:")
print(receptor_cols)

# Create a more flexible approach to find clinical data
available_cols <- colnames(clinical_patient)

# Map common column name variations
col_mapping <- list(
  patient_id = intersect(available_cols, c("bcr_patient_barcode", "submitter_id", "patient_id"))[1],
  er_status = intersect(available_cols, c("er_status_by_ihc", "er_status", "estrogen_receptor_status"))[1],
  pr_status = intersect(available_cols, c("pr_status_by_ihc", "pr_status", "progesterone_receptor_status"))[1],
  her2_status = intersect(available_cols, c("her2_status_by_ihc", "her2_status", "her2_neu_immunohistochemistry_receptor_status"))[1],
  histological_type = intersect(available_cols, c("histological_type", "histologic_diagnosis", "primary_diagnosis"))[1],
  age = intersect(available_cols, c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis", "age"))[1],
  race = intersect(available_cols, c("race_list", "race", "race_category"))[1],
  ethnicity = intersect(available_cols, c("ethnicity", "hispanic_or_latino_ethnicity"))[1]
)

# Remove NULL mappings
col_mapping <- col_mapping[!sapply(col_mapping, is.null)]

print("Column mappings found:")
print(col_mapping)

# Create clinical subset with available columns using base R instead of dplyr
if(!is.null(col_mapping$patient_id)) {
  # Get columns that exist
  existing_cols <- unlist(col_mapping)[unlist(col_mapping) %in% colnames(clinical_patient)]
  
  # Select columns using base R
  clinical_subset <- clinical_patient[, existing_cols, drop = FALSE]
  
  # Rename columns to standard names using base R
  if("bcr_patient_barcode" %in% names(clinical_subset)) {
    names(clinical_subset)[names(clinical_subset) == "bcr_patient_barcode"] <- "patient_id"
  }
  if(!is.null(col_mapping$er_status) && col_mapping$er_status %in% names(clinical_subset)) {
    names(clinical_subset)[names(clinical_subset) == col_mapping$er_status] <- "er_status_by_ihc"
  }
  if(!is.null(col_mapping$pr_status) && col_mapping$pr_status %in% names(clinical_subset)) {
    names(clinical_subset)[names(clinical_subset) == col_mapping$pr_status] <- "pr_status_by_ihc"
  }
  if(!is.null(col_mapping$her2_status) && col_mapping$her2_status %in% names(clinical_subset)) {
    names(clinical_subset)[names(clinical_subset) == col_mapping$her2_status] <- "her2_status_by_ihc"
  }
  
} else {
  print("Warning: Could not find patient ID column in clinical data")
  # Create empty data frame as fallback
  clinical_subset <- data.frame(
    patient_id = character(0),
    er_status_by_ihc = character(0),
    pr_status_by_ihc = character(0),
    her2_status_by_ihc = character(0),
    stringsAsFactors = FALSE
  )
}

sample_annotation <- merge(sample_annotation, clinical_subset, 
                           by = "patient_id", all.x = TRUE)

# Define TNBC samples
sample_annotation$subtype <- "Unknown"
sample_annotation$subtype[
  sample_annotation$er_status_by_ihc == "Negative" &
    sample_annotation$pr_status_by_ihc == "Negative" &
    sample_annotation$her2_status_by_ihc == "Negative"
] <- "TNBC"

sample_annotation$subtype[
  sample_annotation$sample_type == "Solid Tissue Normal"
] <- "Normal"

# Other subtypes
sample_annotation$subtype[
  sample_annotation$er_status_by_ihc == "Positive" |
    sample_annotation$pr_status_by_ihc == "Positive"
] <- "Hormone_Positive"

sample_annotation$subtype[
  sample_annotation$her2_status_by_ihc == "Positive"
] <- "HER2_Positive"

# Print sample distribution
print("Sample distribution by subtype:")
print(table(sample_annotation$subtype, useNA = "always"))

# Filter for samples with known receptor status
valid_samples <- sample_annotation[
  sample_annotation$subtype %in% c("TNBC", "Normal", "Hormone_Positive", "HER2_Positive"),
]

print(paste("Total samples for analysis:", nrow(valid_samples)))
print(paste("TNBC samples:", sum(valid_samples$subtype == "TNBC")))
print(paste("Normal samples:", sum(valid_samples$subtype == "Normal")))

# =====================================================
# PART 3: GENE FILTERING AND NORMALIZATION
# =====================================================

# Filter expression matrix for valid samples
exp_filtered <- exp_matrix[, valid_samples$barcode]

# Remove genes with low expression (keep genes with >10 counts in >10% of samples)
min_samples <- ceiling(0.1 * ncol(exp_filtered))
keep_genes <- rowSums(exp_filtered >= 10) >= min_samples
exp_filtered <- exp_filtered[keep_genes, ]

print(paste("Genes retained after filtering:", nrow(exp_filtered)))

# Add gene symbols
gene_info_filtered <- gene_info[keep_genes, ]
rownames(exp_filtered) <- make.unique(gene_info_filtered$gene_name)

# =====================================================
# PART 4: DIFFERENTIAL EXPRESSION ANALYSIS
# =====================================================

# Prepare DESeq2 object for TNBC vs Normal comparison
sample_info_deseq <- valid_samples
rownames(sample_info_deseq) <- sample_info_deseq$barcode

# TNBC vs Normal analysis
tnbc_normal_samples <- sample_info_deseq[
  sample_info_deseq$subtype %in% c("TNBC", "Normal"),
]
tnbc_normal_exp <- exp_filtered[, tnbc_normal_samples$barcode]

# Create DESeq2 object
dds_tnbc <- DESeqDataSetFromMatrix(
  countData = tnbc_normal_exp,
  colData = tnbc_normal_samples,
  design = ~ subtype
)

# Set reference level
dds_tnbc$subtype <- relevel(dds_tnbc$subtype, ref = "Normal")

# Run DESeq2
dds_tnbc <- DESeq(dds_tnbc)
res_tnbc <- results(dds_tnbc, contrast = c("subtype", "TNBC", "Normal"))

# Convert to data frame and add gene information
res_tnbc_df <- as.data.frame(res_tnbc)
res_tnbc_df$gene_symbol <- rownames(res_tnbc_df)
res_tnbc_df <- res_tnbc_df[!is.na(res_tnbc_df$padj), ]

# =====================================================
# PART 5: TARGET IDENTIFICATION AND PRIORITIZATION
# =====================================================

# Define criteria for drug targets
target_candidates <- res_tnbc_df %>%
  filter(
    padj < 0.05,           # Significant
    log2FoldChange > 1,    # Upregulated in TNBC
    baseMean > 100         # Reasonable expression level
  ) %>%
  arrange(desc(log2FoldChange))

print(paste("Number of upregulated genes in TNBC:", nrow(target_candidates)))

# Get top 50 candidates for detailed analysis
top_targets <- head(target_candidates, 50)

# Save results
write.csv(res_tnbc_df, "results/TNBC_vs_Normal_DEG.csv", row.names = FALSE)
write.csv(target_candidates, "results/TNBC_target_candidates.csv", row.names = FALSE)
write.csv(top_targets, "results/top_50_TNBC_targets.csv", row.names = FALSE)

# =====================================================
# PART 6: VISUALIZATION
# =====================================================

# Volcano plot
p_volcano <- EnhancedVolcano(
  res_tnbc_df,
  lab = res_tnbc_df$gene_symbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'TNBC vs Normal',
  subtitle = 'Differential Gene Expression',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.5,
  legendLabels = c('Not significant',
                   'Log2FC > 1',
                   'p-adj < 0.05',
                   'p-adj < 0.05 & Log2FC > 1'),
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  colConnectors = 'grey50'
)

# View Volcano Plot
p_volcano

# Save Volcano Plot
ggsave("plots/volcano_plot_TNBC_vs_Normal.png", p_volcano, 
       width = 12, height = 10, dpi = 300)

# MA plot
p_ma <- ggplot(res_tnbc_df, aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_point(data = target_candidates, aes(x = log10(baseMean), y = log2FoldChange), 
             color = "red", alpha = 0.8, size = 1.2) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(title = "MA Plot: TNBC vs Normal",
       x = "Log10(Base Mean)",
       y = "Log2(Fold Change)") +
  theme_minimal()

# View MA Plot
p_ma

# Save MA Plot
ggsave("plots/ma_plot_TNBC_vs_Normal.png", p_ma, 
       width = 10, height = 8, dpi = 300)

# Heatmap of top targets
top_genes_exp <- tnbc_normal_exp[top_targets$gene_symbol[1:20], ]
top_genes_exp_norm <- log2(top_genes_exp + 1)

# Prepare annotation for heatmap
annotation_col <- data.frame(
  Subtype = tnbc_normal_samples$subtype,
  row.names = tnbc_normal_samples$barcode
)

colors_subtype <- list(Subtype = c("TNBC" = "red", "Normal" = "blue"))

# Heatmap of top targets
top_genes_exp <- tnbc_normal_exp[top_targets$gene_symbol[1:20], ]
top_genes_exp_norm <- log2(top_genes_exp + 1)

# Prepare annotation for heatmap
annotation_col <- data.frame(
  Subtype = tnbc_normal_samples$subtype,
  row.names = tnbc_normal_samples$barcode
)

colors_subtype <- list(Subtype = c("TNBC" = "red", "Normal" = "blue"))

# Heatmap of top targets
top_genes_exp <- tnbc_normal_exp[top_targets$gene_symbol[1:20], ]
top_genes_exp_norm <- log2(top_genes_exp + 1)

# Prepare annotation for heatmap
annotation_col <- data.frame(
  Subtype = tnbc_normal_samples$subtype,
  row.names = tnbc_normal_samples$barcode
)

colors_subtype <- list(Subtype = c("TNBC" = "red", "Normal" = "blue"))

# Heatmap of top targets
top_genes_exp <- tnbc_normal_exp[top_targets$gene_symbol[1:20], ]
top_genes_exp_norm <- log2(top_genes_exp + 1)

# Prepare annotation for heatmap
annotation_col <- data.frame(
  Subtype = tnbc_normal_samples$subtype,
  row.names = tnbc_normal_samples$barcode
)

colors_subtype <- list(Subtype = c("TNBC" = "red", "Normal" = "blue"))

# Heatmap of top targets
top_genes_exp <- tnbc_normal_exp[top_targets$gene_symbol[1:20], ]
top_genes_exp_norm <- log2(top_genes_exp + 1)

# Prepare annotation for heatmap
annotation_col <- data.frame(
  Subtype = tnbc_normal_samples$subtype,
  row.names = tnbc_normal_samples$barcode
)

colors_subtype <- list(Subtype = c("TNBC" = "red", "Normal" = "blue"))

# Create and display heatmap using ComplexHeatmap
print("Creating heatmap of top 20 target genes with ComplexHeatmap...")

library(ComplexHeatmap)

# Prepare data and annotations
top_genes_exp <- tnbc_normal_exp[top_targets$gene_symbol[1:20], ]
top_genes_exp_norm <- log2(top_genes_exp + 1)

# Apply row scaling (Z-score normalization)
top_genes_exp_scaled <- t(scale(t(top_genes_exp_norm)))

# Check the scaling worked
print(paste("Data range after scaling:", round(min(top_genes_exp_scaled, na.rm = TRUE), 2), 
            "to", round(max(top_genes_exp_scaled, na.rm = TRUE), 2)))

# Create column annotation for sample types
col_anno <- HeatmapAnnotation(
  Type = tnbc_normal_samples$subtype,
  col = list(Type = c("TNBC" = "red", "Normal" = "blue")),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Type = list(
      title = "Type",
      direction = "vertical"
    )
  )
)

# Create the heatmap with properly scaled data
ht <- Heatmap(
  top_genes_exp_scaled,
  name = "Z-score",  # This sets the legend title
  top_annotation = col_anno,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_title = "Top 20 TNBC Target Genes",
  col = colorRampPalette(c("blue", "white", "red"))(50),
  heatmap_legend_param = list(
    title = "Z-score",
    title_position = "topcenter",
    legend_direction = "vertical"
  )
)

# Display the heatmap with stacked legends
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE)

# Save as PNG using Cairo
Cairo::CairoPNG("plots/heatmap_top20_targets.png", width = 900, height = 600)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

print("Heatmap displayed and saved to plots/heatmap_top20_targets.png")

# =====================================================
# PART 7: PATHWAY ENRICHMENT ANALYSIS
# =====================================================

library(stringr)

# Prepare gene lists for enrichment analysis
upregulated_genes <- res_tnbc_df[
  res_tnbc_df$padj < 0.05 & res_tnbc_df$log2FoldChange > 1,
  "gene_symbol"
]

downregulated_genes <- res_tnbc_df[
  res_tnbc_df$padj < 0.05 & res_tnbc_df$log2FoldChange < -1,
  "gene_symbol"
]

# Convert gene symbols to Entrez IDs
upregulated_entrez <- bitr(upregulated_genes, 
                           fromType = "SYMBOL",
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)$ENTREZID

# GO enrichment analysis
go_bp <- enrichGO(gene = upregulated_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# KEGG pathway analysis
kegg_pathways <- enrichKEGG(gene = upregulated_entrez,
                            organism = 'hsa',
                            pvalueCutoff = 0.05)

# Save enrichment results
write.csv(go_bp@result, "results/GO_enrichment_upregulated_TNBC.csv", row.names = FALSE)
write.csv(kegg_pathways@result, "results/KEGG_enrichment_upregulated_TNBC.csv", row.names = FALSE)

# Plot GO enrichment results
p_go <- dotplot(go_bp, showCategory = 20) + 
  ggtitle("GO Biological Processes - Upregulated in TNBC") +
  theme(axis.text.y = element_text(size = 12)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50))

# View GO Enrichment Plot
p_go

# Save GO Enrichment Plot
ggsave("plots/GO_enrichment_dotplot.png", p_go, width = 12, height = 10, dpi = 300)

# Plot KEGG enrichment results
p_kegg <- dotplot(kegg_pathways, showCategory = 20, ) + 
  ggtitle("KEGG Pathways - Upregulated in TNBC") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50))

# View KEGG Enrichment Plot
p_kegg

# Save KEGG Enrichment Plot
ggsave("plots/KEGG_enrichment_dotplot.png", p_kegg, width = 12, height = 10, dpi = 300)

# =====================================================
# PART 8: DRUGGABILITY ASSESSMENT
# =====================================================

# Known druggable gene families (simplified list)
druggable_families <- c(
  "kinase", "phosphatase", "protease", "receptor", "channel", 
  "transporter", "enzyme", "transcription", "epigenetic"
)

# Simple druggability scoring based on gene name patterns
# In practice, you'd use more sophisticated databases like ChEMBL, DrugBank
assess_druggability <- function(gene_names) {
  druggability_score <- rep(0, length(gene_names))
  
  # Kinases
  kinase_pattern <- "K$|KINASE|CDK|PKC|AKT|MAPK|ERK|JNK|PI3K"
  druggability_score[grep(kinase_pattern, gene_names, ignore.case = TRUE)] <- 3
  
  # Receptors
  receptor_pattern <- "R$|RECEPTOR|EGFR|VEGFR|PDGFR|IGF1R"
  druggability_score[grep(receptor_pattern, gene_names, ignore.case = TRUE)] <- 3
  
  # Enzymes
  enzyme_pattern <- "ASE$|ASE[0-9]|ALDH|IDH|PARP|HDAC"
  druggability_score[grep(enzyme_pattern, gene_names, ignore.case = TRUE)] <- 2
  
  # Channels/Transporters
  channel_pattern <- "CHANNEL|TRANSPORTER|SLC|ABC"
  druggability_score[grep(channel_pattern, gene_names, ignore.case = TRUE)] <- 2
  
  # Transcription factors (lower druggability)
  tf_pattern <- "TRANSCRIPTION|FOX|MYC|TP53|STAT"
  druggability_score[grep(tf_pattern, gene_names, ignore.case = TRUE)] <- 1
  
  return(druggability_score)
}

# Add druggability scores to target candidates
target_candidates$druggability_score <- assess_druggability(target_candidates$gene_symbol)

# Prioritize targets based on multiple criteria
target_candidates$priority_score <- (
  target_candidates$log2FoldChange * 0.3 +           # Expression fold change
    (-log10(target_candidates$padj)) * 0.3 +           # Statistical significance
    target_candidates$druggability_score * 0.4         # Druggability
)

# Reorder by priority score
target_candidates <- target_candidates[order(-target_candidates$priority_score), ]

# Save prioritized targets
write.csv(target_candidates, "results/prioritized_TNBC_targets.csv", row.names = FALSE)

print("Creating druggability assessment visualizations...")

# 1. Druggability Score Distribution
druggability_dist <- data.frame(
  score = factor(target_candidates$druggability_score, 
                 levels = 0:3,
                 labels = c("Non-druggable (0)", "Low (1)", "Medium (2)", "High (3)")),
  count = 1
)

p_drug_dist <- ggplot(druggability_dist, aes(x = score, fill = score)) +
  geom_bar(stat = "count", alpha = 0.8) +
  scale_fill_manual(values = c("gray60", "orange", "gold", "darkgreen")) +
  labs(
    title = "Druggability Score Distribution",
    subtitle = "TNBC Target Candidates",
    x = "Druggability Score",
    y = "Number of Genes",
    fill = "Druggability"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# View Druggability Score Distribution
p_drug_dist

# Save Druggability Score Distribution
ggsave("plots/druggability_distribution.png", p_drug_dist, 
       width = 8, height = 6, dpi = 300)

# 2. Scatter plot: Log2FC vs Druggability Score
p_drug_scatter <- ggplot(target_candidates, aes(x = log2FoldChange, y = druggability_score)) +
  geom_point(aes(size = -log10(padj), color = druggability_score), alpha = 0.7) +
  scale_color_gradient(low = "gray", high = "darkgreen", name = "Druggability\nScore") +
  scale_size_continuous(name = "-log10(p-adj)", range = c(1, 6)) +
  labs(
    title = "Expression vs Druggability",
    subtitle = "Size = Statistical Significance",
    x = "Log2 Fold Change (TNBC vs Normal)",
    y = "Druggability Score"
  ) +
  theme_minimal() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue", alpha = 0.5)

# View Scatter plot: Log2FC vs Druggability Score
p_drug_scatter

# Save Scatter plot: Log2FC vs Druggability Score
ggsave("plots/expression_vs_druggability.png", p_drug_scatter, 
       width = 10, height = 8, dpi = 300)

# 3. Top 15 targets ranked by priority score
top_15_targets <- head(target_candidates, 15)

p_top_targets <- ggplot(top_15_targets, aes(x = reorder(gene_symbol, priority_score))) +
  geom_col(aes(y = priority_score, fill = factor(druggability_score)), alpha = 0.8) +
  scale_fill_manual(
    values = c("0" = "gray60", "1" = "orange", "2" = "gold", "3" = "darkgreen"),
    name = "Druggability\nScore",
    labels = c("0" = "Non-druggable", "1" = "Low", "2" = "Medium", "3" = "High")
  ) +
  coord_flip() +
  labs(
    title = "Top 15 TNBC Drug Targets",
    subtitle = "Ranked by Priority Score (Expression + Significance + Druggability)",
    x = "Gene Symbol",
    y = "Priority Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# View Top 15 targets ranked by priority score
p_top_targets

# Save Top 15 targets ranked by priority score
ggsave("plots/top_15_drug_targets.png", p_top_targets, 
       width = 12, height = 8, dpi = 300)

# 4. Create a summary table of top druggable targets
top_druggable <- target_candidates[target_candidates$druggability_score >= 2, ]
top_druggable <- head(top_druggable, 20)

if(nrow(top_druggable) > 0) {
  # Drug target class annotation
  target_class <- character(nrow(top_druggable))
  
  for(i in 1:nrow(top_druggable)) {
    gene <- top_druggable$gene_symbol[i]
    if(grepl("K$|KINASE|CDK|PKC|AKT|MAPK|ERK|JNK|PI3K", gene, ignore.case = TRUE)) {
      target_class[i] <- "Kinase"
    } else if(grepl("R$|RECEPTOR|EGFR|VEGFR|PDGFR|IGF1R", gene, ignore.case = TRUE)) {
      target_class[i] <- "Receptor"
    } else if(grepl("ASE$|ASE[0-9]|ALDH|IDH|PARP|HDAC", gene, ignore.case = TRUE)) {
      target_class[i] <- "Enzyme"
    } else if(grepl("CHANNEL|TRANSPORTER|SLC|ABC", gene, ignore.case = TRUE)) {
      target_class[i] <- "Channel/Transporter"
    } else {
      target_class[i] <- "Other"
    }
  }
  
  top_druggable$target_class <- target_class
  
  # Pie chart of target classes
  class_summary <- table(target_class)
  class_df <- data.frame(
    class = names(class_summary),
    count = as.numeric(class_summary)
  )
  
  p_target_classes <- ggplot(class_df, aes(x = "", y = count, fill = class)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Set3", name = "Target Class") +
    labs(
      title = "Drug Target Classes",
      subtitle = "Distribution of Druggable TNBC Targets"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave("plots/target_classes_pie.png", p_target_classes, 
         width = 8, height = 8, dpi = 300)
  
  # Save detailed druggable targets table
  write.csv(top_druggable[, c("gene_symbol", "log2FoldChange", "padj", 
                              "druggability_score", "priority_score", "target_class")], 
            "results/top_druggable_targets_detailed.csv", row.names = FALSE)
  
  # Detailed visualization of the 20 highly druggable targets
  print("Creating detailed visualization of highly druggable targets...")
  
  # Create a detailed plot showing gene names and scores
  p_detailed_targets <- ggplot(top_druggable, aes(x = reorder(gene_symbol, priority_score), 
                                                  y = priority_score)) +
    geom_col(aes(fill = target_class), alpha = 0.8, width = 0.7) +
    geom_text(aes(label = druggability_score), 
              hjust = -0.1, size = 3.5, fontface = "bold") +
    scale_fill_brewer(palette = "Set2", name = "Target Class") +
    coord_flip() +
    labs(
      title = "Top 20 Highly Druggable TNBC Targets",
      subtitle = "Numbers show druggability scores (2 = Medium, 3 = High)",
      x = "Gene Symbol", 
      y = "Priority Score",
      caption = "Targets ranked by combined expression, significance, and druggability"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  ggsave("plots/detailed_druggable_targets.png", p_detailed_targets, 
         width = 12, height = 10, dpi = 300)
  
  # Create a summary table visualization
  druggability_summary <- top_druggable %>%
    group_by(target_class, druggability_score) %>%
    summarise(count = n(), 
              genes = paste(gene_symbol, collapse = ", "),
              .groups = "drop") %>%
    arrange(desc(druggability_score), target_class)
  
  # Create a heatmap-style summary
  p_summary_heatmap <- ggplot(druggability_summary, 
                              aes(x = factor(druggability_score), y = target_class)) +
    geom_tile(aes(fill = count), color = "white", linewidth = 1) +
    geom_text(aes(label = count), size = 6, fontface = "bold", color = "white") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of\nTargets") +
    labs(
      title = "Druggable Target Summary",
      subtitle = "Count of targets by class and druggability score",
      x = "Druggability Score",
      y = "Target Class"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  
  ggsave("plots/druggable_targets_summary_heatmap.png", p_summary_heatmap, 
         width = 8, height = 6, dpi = 300)
  
  # Print detailed breakdown
  print("Detailed breakdown of highly druggable targets:")
  print(paste(rep("=", 50), collapse = ""))
  
  for(class in unique(top_druggable$target_class)) {
    class_targets <- top_druggable[top_druggable$target_class == class, ]
    cat("\n", class, " (", nrow(class_targets), " targets):\n", sep="")
    
    for(i in 1:nrow(class_targets)) {
      cat("  • ", class_targets$gene_symbol[i], 
          " (Score: ", class_targets$druggability_score[i], 
          ", FC: ", round(class_targets$log2FoldChange[i], 2), ")\n", sep="")
    }
  }
  
  print("\nVisualization files created:")
  print("- plots/detailed_druggable_targets.png")
  print("- plots/druggable_targets_summary_heatmap.png")
}

print("Druggability visualizations saved to plots/ directory")

# View Pie chart of target classes 
p_target_classes

# View Detailed visualization of the 20 highly druggable targets
p_detailed_targets

# View heatmap-style summary
p_summary_heatmap

# =====================================================
# PART 10: SUMMARY REPORT
# =====================================================

# Create a summary report
summary_stats <- list(
  total_samples = nrow(valid_samples),
  tnbc_samples = sum(valid_samples$subtype == "TNBC"),
  normal_samples = sum(valid_samples$subtype == "Normal"),
  genes_analyzed = nrow(exp_filtered),
  significant_genes = sum(res_tnbc_df$padj < 0.05, na.rm = TRUE),
  upregulated_genes = nrow(target_candidates),
  top_targets_druggable = sum(target_candidates$druggability_score >= 2)
)

# Save summary
capture.output(
  {
    cat("TCGA-BRCA TNBC Analysis Summary\n")
    cat("================================\n\n")
    cat("Sample Information:\n")
    cat(sprintf("- Total samples analyzed: %d\n", summary_stats$total_samples))
    cat(sprintf("- TNBC samples: %d\n", summary_stats$tnbc_samples))
    cat(sprintf("- Normal samples: %d\n", summary_stats$normal_samples))
    cat(sprintf("\nGene Analysis:\n"))
    cat(sprintf("- Genes analyzed: %d\n", summary_stats$genes_analyzed))
    cat(sprintf("- Significantly differentially expressed: %d\n", summary_stats$significant_genes))
    cat(sprintf("- Upregulated in TNBC: %d\n", summary_stats$upregulated_genes))
    cat(sprintf("- Potentially druggable targets: %d\n", summary_stats$top_targets_druggable))
  },
  file = "results/analysis_summary.txt"
)
