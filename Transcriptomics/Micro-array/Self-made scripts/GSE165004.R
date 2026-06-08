# GSE165004
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(limma)

GSE165004 <- getGEO("GSE165004", GSEMatrix = TRUE, AnnotGPL=TRUE)

GSE165004_eset <- GSE165004[[1]]

GSE165004_exprs_data <- exprs(GSE165004_eset)
GSE165004_meta_data <- pData(GSE165004_eset)
view(GSE165004_meta_data)
# ============================================= # 
# Log2Fold-Change
# ============================================= # 
GSE165004_exprs_data[GSE165004_exprs_data <= 0] <- NaN
exprs(GSE165004_eset) <- log2(GSE165004_exprs_data)
#===# 
annotation(GSE165004_eset)
glimpse(GSE165004_meta_data)
print(GSE165004_meta_data[, c(2, 11)])

print(GSE165004_meta_data$`subject status/group:ch1`)
print(GSE165004_meta_data["subject status/group:ch1"])
count(as.list(GSE165004_meta_data["subject status/group:ch1"]))
# GSE165004_condition <- c("Control", "Control", "Control", "Control", "Control", "Control",
#                         "Control", "Control", "Control", "Control", "Control", "Control",
#                         "Control", "Control", "Control", "Control", "Control", "Control",
#                         "Control", "Control", "Control", "Control", "Control", "Control",
#                         "RPL", "RPL", "RPL", "RPL", "RPL", "RPL",
#                         "RPL", "RPL", "RPL", "RPL", "RPL", "RPL",
#                         "RPL", "RPL", "RPL", "RPL", "RPL", "RPL",
#                         "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure",
#                         "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure",
#                         "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure",
#                         "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure")

GSE165004_condition <- c(rep("Control", 24), rep("RPL", 24), rep("Implantation_Failure", 24))

GSE165004_condition <- factor(GSE165004_condition, levels = c("Control", "Implantation_Failure", "RPL"))

# ============================================= # 
# Checking normalisation
group_colors <- c("Control" = "#4CAF50",              # Soft Green
                  "Implantation_Failure" = "#2196F3", # Soft Blue
                  "RPL" = "#FF5722")                  # Soft Orange
sample_colors <- group_colors[GSE165004_condition]
boxplot(exprs(GSE165004_eset), 
        outline=FALSE, 
        col= sample_colors,
        main = "GSE165004: Post-Log2 Normalization Check", 
        las=2, 
        cex.axis = 0.7)

legend("topleft", legend = names(group_colors), fill = group_colors, bty = "n")


GSE165004_design <- model.matrix(~0 + GSE165004_condition)

colnames(GSE165004_design) <- levels(GSE165004_condition)
view(Transcriptomics)
# ============================================= # 
GSE165004_fit <- lmFit(GSE165004_eset, GSE165004_design)
GSE165004_cont_matrix <- makeContrasts(RPL_vs_Control = RPL - Control,
                                       levels = GSE165004_design)
GSE165004_fit2 <- contrasts.fit(GSE165004_fit, GSE165004_cont_matrix)
GSE165004_fit2 <- eBayes(GSE165004_fit2)
GSE165004_limma_results <- topTable(GSE165004_fit2, coef= "RPL_vs_Control", number=Inf, sort.by="P")
view(GSE165004_limma_results)
# ===============# 
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
# ===============# 
write.table(GSE165004_limma_results, file = "GSE165004_RPL_vs_Control_Annotated.csv", sep = ";", dec = ".", row.names=FALSE)







# ============================================= # 
#                                               #
#                 Making Plots                  #
#                                               #
# ============================================= # 


# ========================================================== #
#   Plot Dispersion Estimates / Mean-Difference Plot         #
# ========================================================== #
plotSA(GSE165004_fit2, main="GSE165004: Mean-Variance Trend")

# ========================================================== #
# Histogram
# ========================================================== #
hist(GSE165004_limma_results$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution for GSE165004")

#GSE165004
GSE165004_decideTests <- decideTests(GSE165004_fit2, adjust.method="fdr", p.value=0.05, lfc=0)
# ================ #
GSE165004_T_good <- which(!is.na(GSE165004_fit2$F)) # filter out bad probes
qqt(GSE165004_fit2$t[GSE165004_T_good], GSE165004_fit2$df.total[GSE165004_T_good], main="Moderated t-statistic for GSE165004")

GSE165004_limma_results <- GSE165004_limma_results %>%
  mutate(threshold = case_when(
    P.Value < 0.05 & logFC > 0  ~ "Upregulated",
    P.Value < 0.05 & logFC < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))
str(GSE26787_limma_results)
ggplot(GSE165004_limma_results,
       aes(x = logFC,
           y = -log10(adj.P.Val),
           color = threshold)) +
  
  geom_point(size = 1.5, alpha = 0.7) +
  # Horizontal significance line
  # geom_hline(yintercept = -log10(0.05),
  #           linetype = "dotted",
  #           color = "black") +
  
  # Vertical fold-change cutoffs
  # geom_vline(xintercept = c(-1, 1),
  #           linetype = "dotted",
  #           color = "black") +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#00BFFF",   # blue
      "Upregulated"   = "#FF3030",   # red
      "Not significant" = "black"
    )
  ) +
  scale_x_continuous(limits = c(-1.5, 1.5)) +
  scale_y_continuous(limits = c(0, 35)) +

  labs(
    title = "Volcano Plot of GSE165004: RPL vs ControL
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
# ==========
