# GSE161969
library("rstudioapi", "tidyverse", DESeq2)

Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files(Folder_of_your_choice)


GSE161969_counts_data <- read.delim("GSE161969_expr.tsv")
colnames(GSE161969_counts_data)

GSE161969_meta <- read.delim("GSE161969_meta.tsv")
glimpse(GSE161969_meta) # To see case or control
print(GSE161969_meta[, c("title", "characteristics_ch1","characteristics_ch1.1")])


# To see if they match
colnames(GSE161969_counts_data)
print(GSE161969_meta[["geo_accession"]])

GSE161969_sample_info <- data.frame(
  condition = c("Control", "Control", "Control", "RPL", "RPL","RPL","RPL"),
  row.names = c("GSM4928875", "GSM4928876", "GSM4928877", "GSM4928878", "GSM4928879", "GSM4928880", "GSM4928881")
)

GSE161969_sample_info$condition <- factor(GSE161969_sample_info$condition, levels = c("Control", "RPL"))
all(rownames(GSE161969_sample_info) == colnames(GSE161969_counts_data))
# ==================#
glimpse(GSE161969_DESeq_object)
GSE161969_DESeq_object <- DESeq2::DESeqDataSetFromMatrix(countData = GSE161969_counts_data,
                                                 colData = GSE161969_sample_info,
                                                 design = ~ condition)

# This single line runs the normalization and statistical testing!
GSE161969_DESeq_object <- DESeq2::DESeq(GSE161969_DESeq_object)


# ==========================================
# STEP 5: Extract the Blueprint (The Results)
# ==========================================
GSE161969_results <- DESeq2::results(GSE161969_DESeq_object, contrast=c("condition", "RPL", "Control"))

# MAYBE DO an log2 fold changes "lfcShrink()"?????

# Look at your top significantly changed genes
GSE161969_results_ordered <- GSE161969_results[order(GSE161969_results$padj), ]


# ========================================================== #
# STEP 2: Annotating Dataset
# ========================================================== #
GSE161969_results_df <- as.data.frame(GSE161969_results_ordered)
# ========================================================== #
# Check if annotation file exists
exists("Annotation")
# ========================================================== #
GSE161969_annotated_results <- merge(GSE161969_results_df, Annotation, by=0, sort=FALSE)
colnames(GSE161969_annotated_results)[1] <- "EntrezID"

# ========================================================== #
# 
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files()  
write.table(GSE161969_annotated_results, file = "GSE161969_annotated_results.tsv", quote = FALSE, sep = "\t", dec = ".", row.names=FALSE)
write.table(GSE161969_annotated_results, file = "GSE161969_annotated_results.csv", quote = FALSE, sep = ";", dec = ".", row.names=FALSE)

# ========================================================== #
DESeq2::plotDispEsts(GSE161969_DESeq_object, main="GSE161969 Dispersion Estimates")
# ========================================================== #
# Histogram
hist(GSE161969_results$pvalue, 
     breaks=seq(0, 1, length = 21), 
     col = "grey", border = "white", 
     xlab = "Raw P-value", 
     ylab = "Number of Genes", 
     main = "GSE161969 Frequencies of p-values")
# ========================================================== #
GSE161969_subset <- subset(GSE161969_annotated_results, padj < 0.05)
# ========================================================== #
GSE161969_annotated_results <- GSE161969_annotated_results %>%
  mutate(threshold = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))
# ========================================================== #

Volcanoplot_GSE161969 <- ggplot(GSE161969_annotated_results,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = threshold)) +
  
  geom_point(size = 1.5, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#00BFFF",   # blue
      "Upregulated"   = "#FF3030",   # red
      "Not significant" = "black"
    )
  ) +
  
  labs(
    title = "Volcano Plot of GSE161969: Controls vs RPL
    (Decidual Tissue)",
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

png("Volcano Plot of GSE161969.png",
    width = 631,
    height = 618)
Volcanoplot_GSE161969
dev.off()
