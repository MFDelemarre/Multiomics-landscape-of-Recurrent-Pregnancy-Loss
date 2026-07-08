# GSE183555
library("rstudioapi", "tidyverse", DESeq2)

Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files(Folder_of_your_choice)


GSE183555_counts_data <- read.delim("GSE183555_expr.tsv")
colnames(GSE183555_counts_data)

# =============== #
GSE183555_meta <- read.delim("GSE183555_meta.tsv")
glimpse(GSE183555_meta) # To see case or control
print(GSE183555_meta[, c("title", "characteristics_ch1","characteristics_ch1.1")])
# =============== #

# To see if they match
colnames(GSE183555_counts_data)
print(GSE183555_meta[["geo_accession"]])

GSE183555_sample_info <- data.frame(
  condition = c("Control", 
                "Control", 
                "Control",
                "Control",
                "Control",
                "RPL", 
                "RPL", 
                "RPL", 
                "RPL",
                "RPL"),
  row.names = c("GSM5591871", 
                "GSM5591872",
                "GSM5591873", 
                "GSM5591874",
                "GSM5591875",
                "GSM5591876",
                "GSM5591877",
                "GSM5591878",
                "GSM5591879",
                "GSM5591880")
)
view(GSE183555_sample_info)
GSE183555_sample_info$condition <- factor(GSE183555_sample_info$condition, levels = c("Control", "RPL"))
all(rownames(GSE183555_sample_info) == colnames(GSE183555_counts_data))
# ==================#
glimpse(GSE183555_DESeq_object)
GSE183555_DESeq_object <- DESeq2::DESeqDataSetFromMatrix(countData = GSE183555_counts_data,
                                                         colData = GSE183555_sample_info,
                                                         design = ~ condition)

GSE183555_DESeq_object <- DESeq2::DESeq(GSE183555_DESeq_object)


GSE183555_results <- DESeq2::results(GSE183555_DESeq_object, contrast=c("condition", "RPL", "Control"))

# Ordering the genes based on p-value
GSE183555_results_ordered <- GSE183555_results[order(GSE183555_results$padj), ]

# ========================================================== #
# STEP 2: Annotating Dataset
# ========================================================== #
GSE183555_results_df <- as.data.frame(GSE183555_results_ordered)
# ========================================================== #
# Check if annotation file exists
exists("Annotation")
# ========================================================== #
# Here you can annotate
GSE183555_annotated_results <- merge(GSE183555_results_df, Annotation, by=0, sort=FALSE)
colnames(GSE183555_annotated_results)[1] <- "EntrezID"
# ========================================================== #
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files()  
write.table(GSE183555_annotated_results, file = "GSE183555_annotated_results.tsv", quote = FALSE, sep = "\t", dec = ".", row.names=FALSE)
write.table(GSE183555_annotated_results, file = "GSE183555_annotated_results.csv", quote = FALSE, sep = ";", dec = ".", row.names=FALSE)
# ========================================================== #
#   Plot Dispersion Estimates / Mean-Difference Plot         #
# ========================================================== #

DESeq2::plotDispEsts(GSE183555_DESeq_object, main="GSE183555 Dispersion Estimates")
# ========================================================== #
# Histogram
hist(GSE183555_results$pvalue, 
     breaks=seq(0, 1, length = 21), 
     col = "grey", border = "white",
     xlab = "Raw P-value", 
     ylab = "Number of Genes", 
     main = "GSE183555 Frequencies of p-values")






# Normalising & PCA
GSE183555_vsd <- vst(GSE183555_DESeq_object, blind=FALSE)
plotPCA(GSE65099_vsd, intgroup="condition")  +  ggtitle("Principal Component Analysis GSE183555")



# ========================================================== #
GSE183555_annotated_results <- GSE183555_annotated_results %>%
  mutate(threshold = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))
# ========================================================== #
ggplot(GSE183555_annotated_results,
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
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(0, 5)) +
  labs(
    title = "Volcano Plot of GSE183555: RPL vs Control
    (Endometrial Tissue)",
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

png("Volcano Plot of GSE183555.png",
    width = 631,
    height = 618)
get_last_plot()
dev.off()
