# GSE65099
library("rstudioapi", "tidyverse", DESeq2)

Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files(Folder_of_your_choice)

# "GSE65099_expr.tsv" "GSE65099_meta.tsv"
GSE65099_counts_data <- read.delim("GSE65099_expr.tsv")
colnames(GSE65099_counts_data)

GSE65099_meta <- read.delim("GSE65099_meta.tsv")
glimpse(GSE65099_meta) # To see case or control
print(GSE65099_meta[, c("title", "characteristics_ch1","characteristics_ch1.1")])


# To see if they match
colnames(GSE65099_counts_data)
print(GSE65099_meta[["geo_accession"]])

GSE65099_sample_info <- data.frame(
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
GSE65099_sample_info <- data.frame(
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


GSE65099_sample_info$condition <- factor(GSE65099_sample_info$condition, levels = c("Infertile", "RPL"))
all(rownames(GSE65099_sample_info) == colnames(GSE65099_counts_data))
# ==================#
GSE65099_DESeq_object <- DESeq2::DESeqDataSetFromMatrix(countData = GSE65099_counts_data,
                                                         colData = GSE65099_sample_info,
                                                         design = ~ condition)

# This single line runs the normalization and statistical testing!
GSE65099_DESeq_object <- DESeq2::DESeq(GSE65099_DESeq_object)

# ==========================================
# STEP 5: Extract the Blueprint (The Results)
# ==========================================
GSE65099_results <- DESeq2::results(GSE65099_DESeq_object, contrast=c("condition", "RPL", "Infertile"))

# Look at your top significantly changed genes
GSE65099_results_ordered <- GSE65099_results[order(GSE65099_results$padj), ]

# ========================================================== #
# STEP 2: Annotating Dataset
# ========================================================== #
GSE65099_results_df <- as.data.frame(GSE65099_results_ordered)
# ========================================================== #
# Check if annotation file exists
exists("Annotation")
# ========================================================== #
GSE65099_annotated_results <- merge(GSE65099_results_df, Annotation, by=0, sort=FALSE)
colnames(GSE65099_annotated_results)[1] <- "EntrezID"
# ========================================================== #
# Save your results to a .csv
write.csv(GSE65099_annotated_results, file="GSE65099_annotated_results.csv", row.names=FALSE)
write.table(GSE65099_annotated_results, file = "GSE65099_annotated_results.csv", sep = ";", dec = ".", row.names=FALSE)




# ============================================= # 
#                                               #
#                 Making Plots                  #
#                                               #
# ============================================= # 





# ========================================================== #
# Dispersion Estimates
DESeq2::plotDispEsts(GSE65099_DESeq_object, main="RNA-seq: GSE65099 Dispersion Estimates")
# ========================================================== #
# Histogram
hist(GSE65099_results$pvalue, 
     breaks=seq(0, 1, length = 21), 
     col = "grey", border = "white",
     xlab = "Nominal P-value", 
     ylab = "Number of Genes", 
     main = "GSE65099 Frequencies of p-values")
# ========================================================== #
GSE65099_subset <- subset(GSE65099_annotated_results, padj < 0.05)
# ========================================================== #
GSE65099_annotated_results <- GSE65099_annotated_results %>%
  mutate(threshold = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))
# ========================================================== #
# ggplot
# Volcano plot
ggplot(GSE65099_annotated_results,
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
  
  
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 5)) +
  
  
  labs(
    title = "Volcano Plot of GSE65099 (Endometrial-tissue)",
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
GSE65099_vsd <- vst(GSE65099_DESeq_object, blind=FALSE)
plotPCA(GSE65099_vsd, intgroup="condition")  +  ggtitle("Principal Component Analysis GSE65099")
