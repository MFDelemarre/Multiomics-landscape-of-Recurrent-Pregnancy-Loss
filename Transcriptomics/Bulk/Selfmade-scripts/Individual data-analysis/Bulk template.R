# GSE_RPLStudy -> based on
library("rstudioapi", "tidyverse", DESeq2)

Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files(Folder_of_your_choice)

# "GSE_RPLStudy_expr.tsv" "GSE_RPLStudy_meta.tsv"
GSE_RPLStudy_counts_data <- read.delim("GSE_RPLStudy_expr.tsv")
colnames(GSE_RPLStudy_counts_data)

GSE_RPLStudy_meta <- read.delim("GSE_RPLStudy_meta.tsv")
glimpse(GSE_RPLStudy_meta) # To see case or control
print(GSE_RPLStudy_meta[, c("title", "characteristics_ch1","characteristics_ch1.1")])


# To see if they match
colnames(GSE_RPLStudy_counts_data)
print(GSE_RPLStudy_counts_data[["geo_accession"]])

GSE_RPLStudy_sample_info <- data.frame(
  condition = c(rep("Infertile", 10),
                rep("RPL", 9)), # "Infertile" "RPL"
  row.names = c("GSM1587300",
                "GSM1587301",
                "GSM1587302", 
                "GSM1587303",
                "GSM1587304", 
                "GSM1587305",
                "GSM1587306",
                "GSM1587307",
                "GSM1587308",
                "GSM1587309", 
                "GSM1587310", 
                "GSM1587311", 
                "GSM1587312",
                "GSM1587313",
                "GSM1587314", 
                "GSM1587316",
                "GSM1587317",
                "GSM1587318",
                "GSM1587319")
)

# ==================
# or
GSE_RPLStudy_sample_info <- data.frame(
  condition = c("Infertile", 
                "Infertile",
                "Infertile",
                "Infertile",
                "Infertile",
                "Infertile",
                "Infertile",
                "Infertile",
                "Infertile",
                "Infertile",
                "RPL",
                "RPL",
                "RPL", 
                "RPL",
                "RPL", 
                "RPL",
                "RPL",
                "RPL",
                "RPL"), # "Infertile" "RPL"
  row.names = c("GSM1587300",
                "GSM1587301",
                "GSM1587302", 
                "GSM1587303",
                "GSM1587304", 
                "GSM1587305",
                "GSM1587306",
                "GSM1587307",
                "GSM1587308",
                "GSM1587309", 
                "GSM1587310", 
                "GSM1587311", 
                "GSM1587312",
                "GSM1587313",
                "GSM1587314", 
                "GSM1587316",
                "GSM1587317",
                "GSM1587318",
                "GSM1587319")
)
# ==================


GSE_RPLStudy_sample_info$condition <- factor(GSE_RPLStudy_sample_info$condition, levels = c("Infertile", "RPL"))
all(rownames(GSE_RPLStudy_sample_info) == colnames(GSE_RPLStudy_counts_data))

# ============================================= #
# Makes a "DESeqDataSet" from matrix so it can  #
# be used for Statistical                       #
# Differential expression analysis              #
# ============================================= #
GSE_RPLStudy_DESeq_object <- DESeq2::DESeqDataSetFromMatrix(countData = GSE_RPLStudy_counts_data,
                                                            colData = GSE_RPLStudy_sample_info,
                                                            design = ~ condition)
?DESeqDataSetFromMatrix

# ============================================= #
#             Filtering Expressed genes         #
# ============================================= #
GSE_RPLStudy_Expressed_genes <- rowSums(counts(GSE_RPLStudy_DESeq_object)) > 0

GSE_RPLStudy_DESeq_object <- GSE_RPLStudy_DESeq_object[GSE_RPLStudy_Expressed_genes, ]
# ============================================= # 
#           Checking normalization              #
# ============================================= # 
group_colors <- c("Control" = "#4CAF50",              # Soft Green
                  "Implantation_Failure" = "#2196F3", # Soft Blue
                  "RPL" = "#FF5722")                  # Soft Orange
GSE_RPLStudy_colors <- group_colors[GSE_RPLStudy_condition]
boxplot(exprs(GSE165004_eset), 
        outline=FALSE, 
        col= GSE_RPLStudy_colors,
        main = "GSE_RPLStudy: Post-Log2 Normalization Check", 
        las=2, 
        cex.axis = 0.7)
legend("topleft", legend = names(group_colors), fill = group_colors, bty = "n")


# ============================================= # 
#            Statistical Analysis               #
#                 Using the                     #
#      Gamma–Poisson (mixture) distribution     #
#                     or                        #
#         Negative Binomial Distribution        # 
# ============================================= #
GSE_RPLStudy_DESeq_object <- DESeq2::DESeq(GSE_RPLStudy_DESeq_object)
# ==========================================
# STEP 5: Extract the Blueprint (The Results)
# ==========================================
GSE_RPLStudy_results <- DESeq2::results(GSE_RPLStudy_DESeq_object, contrast=c("condition", "RPL", "Infertile"))

# Look at your top significantly changed genes
GSE_RPLStudy_results_ordered <- GSE_RPLStudy_results[order(GSE_RPLStudy_results$padj), ]

# ========================================================== #
# STEP 2: Annotating Dataset
# ========================================================== #
GSE_RPLStudy_results_df <- as.data.frame(GSE_RPLStudy_results_ordered)
# ========================================================== #
# Check if annotation file exists
exists("Annotation")
# ========================================================== #
GSE_RPLStudy_annotated_results <- merge(GSE_RPLStudy_results_df, Annotation, by=0, sort=FALSE)
colnames(GSE_RPLStudy_annotated_results)[1] <- "EntrezID"

# ========================================================== #
GSE_RPLStudy_annotated_results <- GSE_RPLStudy_annotated_results %>%
  mutate(threshold = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))

# ========================================================== #
# Save your results to a .csv
write.csv(GSE_RPLStudy_annotated_results, file="GSE_RPLStudy_annotated_results.csv", row.names=FALSE)
write.table(GSE_RPLStudy_annotated_results, file = "GSE_RPLStudy_annotated_results.csv", sep = ";", dec = ".", row.names=FALSE)
# ========================================================== #







# ============================================= # 
#                                               #
#                 Making Plots                  #
#                                               #
# ============================================= # 


# ========================================================== #
#   Plot Dispersion Estimates / Mean-Difference Plot         #
# ========================================================== #
DESeq2::plotDispEsts(GSE_RPLStudy_DESeq_object, main="GSE_RPLStudy Dispersion Estimates")

# ========================================================== #
# Normalising & PCA
# ========================================================== #
GSE_RPLStudy_vsd <- vst(GSE_RPLStudy_DESeq_object, blind=FALSE)
plotPCA(GSE_RPLStudy_vsd, intgroup="condition")  +  ggtitle("Principal Component Analysis GSE_RPLStudy")

# ========================================================== #
# Histogram
# ========================================================== #
hist(GSE_RPLStudy_results$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "GSE_RPLStudy Frequencies of padj-values")
# ========================================================== #
# Volcano Plot
# ========================================================== #
SE_RPLStudy_annotated_results <- GSE_RPLStudy_annotated_results %>%
  mutate(threshold = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))

ggplot(GSE_RPLStudy_annotated_results,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = threshold)) +
  
  geom_point(size = 1.5, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#00BFFF",   # blue
      "Upregulated"   = "#FF3030",   # red
      "Not significant" = "grey70"
    )
  ) +
  
  labs(
    title = "Volcano Plot of GSE_RPLStudy: Controls vs Results",
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = "Gene status"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  )
