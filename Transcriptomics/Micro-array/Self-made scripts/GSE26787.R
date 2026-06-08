if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(limma)

GSE26787 <- getGEO("GSE26787", GSEMatrix = TRUE, AnnotGPL=TRUE)

GSE26787_eset <- GSE26787[[1]]

GSE26787_exprs_data <- exprs(GSE26787_eset)
GSE26787_meta_data <- pData(GSE26787_eset)

?Biobase::fData
# ============== # 
GSE26787_exprs_data_test2 <- GSE26787_eset
# Log2

GSE26787_exprs_data[GSE26787_exprs_data <= 0] <- NaN
exprs(GSE26787_eset) <- log2(GSE26787_exprs_data)
#===# 
annotation(GSE26787_eset)
print(GSE26787_meta_data[, c(2, 11)])

print(GSE26787_meta_data["source_name_ch1"])


GSE26787_condition <- c("Control", "Control", "Control", "Control", "Control",
                        "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", 
                        "RPL", "RPL", "RPL", "RPL", "RPL")

GSE26787_condition <- factor(GSE26787_condition, levels = c("Control", "Implantation_Failure", "RPL"))

GSE26787_design <- model.matrix(~0 + GSE26787_condition)

colnames(GSE26787_design) <- levels(GSE26787_condition)
#===# 
GSE26787_fit <- lmFit(GSE26787_eset, GSE26787_design)
#===# 
GSE26787_cont_matrix <- makeContrasts(RPL_vs_Control = RPL - Control,
                                      levels = GSE26787_design)

GSE26787_fit2 <- contrasts.fit(GSE26787_fit, GSE26787_cont_matrix)
# ========= # 
GSE26787_fit2 <- eBayes(GSE26787_fit2)

GSE26787_limma_results <- topTable(GSE26787_fit2, coef= "RPL_vs_Control", number=Inf, sort.by="P")

# ===============# 
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
# ===============# 
write.table(GSE26787_limma_results, file = "GSE26787_RPL_vs_Control_Annotated.csv", sep = ";", dec = ".", row.names=FALSE)







# ============================================= # 
#                                               #
#                 Making Plots                  #
#                                               #
# ============================================= # 





# ========================================================== #
#   Plot Dispersion Estimates / Mean-Difference Plot         #
# ========================================================== #
plotSA(GSE26787_fit2, main="GSE26787: Mean-Variance Trend")
# ===============# 


GSE26787_t_good <- which(!is.na(GSE26787_fit2$F)) # filter out bad probes
qqt(GSE26787_fit2$t[GSE26787_t_good], GSE26787_fit2$df.total[GSE26787_t_good], main="GSE26787: Moderated t-statistic")

# ========================================================== #
# Volcano Plot
# ========================================================== #
GSE26787_limma_results_test <- GSE26787_limma_results %>%
  mutate(threshold = case_when(
    adj.P.Val < 0.05 & logFC > 0  ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))
str(GSE26787_limma_results)
ggplot(GSE26787_limma_results_test,
       aes(x = logFC,
           y = -log10(adj.P.Val),
           color = threshold)) +
  
  geom_point(size = 1.5, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#00BFFF",   # blue
      "Upregulated"   = "#FF3030",   # red
      "Not significant" = "black"
    )
  ) +
  
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(0, 20)) +
  labs(
    title = "Volcano Plot of GSE26787: RPL vs Controls 
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
