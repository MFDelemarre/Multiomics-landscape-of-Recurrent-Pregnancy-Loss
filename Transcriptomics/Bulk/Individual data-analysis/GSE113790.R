# GSE113790
library("rstudioapi", "tidyverse", DESeq2)

Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files(Folder_of_your_choice)


GSE113790_counts_data <- read.delim("GSE113790_expr.tsv")
colnames(GSE113790_counts_data)

GSE113790_meta <- read.delim("GSE113790_meta.tsv")
glimpse(GSE113790_meta) # To see case or control
print(GSE113790_meta[, c("title", "characteristics_ch1","characteristics_ch1.1")])


# To see if they match
colnames(GSE113790_counts_data)
print(GSE113790_meta[["geo_accession"]])

GSE113790_sample_info <- data.frame(
  condition = c("RPL", "RPL", "RPL", "Control", "Control", "Control"),
  row.names = c("GSM3119483", "GSM3119484", "GSM3119485", "GSM3119486", "GSM3119487", "GSM3119488")
)
view(GSE113790_sample_info)
GSE113790_sample_info$condition <- factor(GSE113790_sample_info$condition, levels = c("Control", "RPL"))
all(rownames(GSE113790_sample_info) == colnames(GSE113790_counts_data))
# ==================#
glimpse(GSE113790_DESeq_object)
GSE113790_DESeq_object <- DESeq2::DESeqDataSetFromMatrix(countData = GSE113790_counts_data,
                                                         colData = GSE113790_sample_info,
                                                         design = ~ condition)

GSE113790_DESeq_object <- DESeq2::DESeq(GSE113790_DESeq_object)


GSE113790_results <- DESeq2::results(GSE113790_DESeq_object, contrast=c("condition", "Control", "RPL"))

# Ordering the genes based on p-value
GSE113790_results_ordered <- GSE113790_results[order(GSE113790_results$padj), ]

# ========================================================== #
# STEP 2: Annotating Dataset
# ========================================================== #
GSE113790_results_df <- as.data.frame(GSE113790_results_ordered)
view(GSE113790_results_df)
# ========================================================== #
# Check if annotation file exists
exists("Annotation")
# ========================================================== #
# Here you can annotate
GSE113790_annotated_results <- merge(GSE113790_results_df, Annotation, by=0, sort=FALSE)
colnames(GSE113790_annotated_results)[1] <- "EntrezID"

# ========================================================== #
GSE113790_annotated_results <- GSE113790_annotated_results %>%
  mutate(threshold = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))

view(GSE113790_annotated_results)
# Save your results to a .csv
write.csv(GSE113790_annotated_results, file="GSE113790_annotated_results.csv", row.names=FALSE)








# ============================================= # 
#                                               #
#                 Making Plots                  #
#                                               #
# ============================================= # 







# ========================================================== #
DESeq2::plotDispEsts(GSE113790_DESeq_object, main="GSE113790 Dispersion Estimates")
# ========================================================== #
# Histogram
hist(GSE113790_results$pvalue, 
     breaks=seq(0, 1, length = 21), 
     col = "grey", border = "white",
     xlab = "Raw P-value", 
     ylab = "Number of Genes", 
     main = "GSE113790 Frequencies of p-values")
# ========================================================== #
GSE65099_subset <- subset(GSE113790_annotated_results, padj < 0.05)
# ========================================================== #
# ggplot
ggplot(GSE113790_annotated_results,
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
    title = "Volcano Plot of GSE113790: Controls vs RPL",
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

str(GSE161969_vsd)
# Normalising & PCA
GSE65099_vsd <- vst(GSE113790_DESeq_object, blind=FALSE)
plotPCA(GSE65099_vsd, intgroup="condition")  +  ggtitle("Principal Component Analysis GSE113790")
