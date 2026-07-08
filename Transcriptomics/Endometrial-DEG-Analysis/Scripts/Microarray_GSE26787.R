if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(limma)

# ============================================= # 
#                 Loading dataset               #
# ============================================= # 

GSE26787 <- getGEO("GSE26787", GSEMatrix = TRUE, AnnotGPL=TRUE)

GSE26787_eset <- GSE26787[[1]]

GSE26787_exprs_data <- exprs(GSE26787_eset)
GSE26787_meta_data <- pData(GSE26787_eset)

rm(GSE26787)
# ============================================= # 
#        Log2Fold-Change of RPL and other       #
# ============================================= # 

GSE26787_exprs_data[GSE26787_exprs_data <= 0] <- NaN
exprs(GSE26787_eset) <- log2(GSE26787_exprs_data)
# ============================================= # 
#             Annotation/Grouping               #
#               of RPL and other                #
# ============================================= # 
annotation(GSE26787_eset)
print(GSE26787_meta_data[, c(2, 11)])

print(GSE26787_meta_data["source_name_ch1"])


GSE26787_condition <- c("Control", "Control", "Control", "Control", "Control",
                        "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", "Implantation_Failure", 
                        "RPL", "RPL", "RPL", "RPL", "RPL")

GSE26787_condition <- factor(GSE26787_condition, levels = c("Control", "Implantation_Failure", "RPL"))
# ============================================= # 
#             Statistical Analysis              #
#                   Using the                   #
#           Linear Regression Model             #
# ============================================= #
# Making a matrix based on the design           #
# That is condition RPL/Condition/Infertile     #
#     'y' will be outcome to see if             #
#     'condition' has an effect                 #
# ============================================= #
# Intercept = 0 (to see if there is a effect)   #
# ============================================= #
GSE26787_design <- model.matrix(~0 + GSE26787_condition)
colnames(GSE26787_design) <- levels(GSE26787_condition)
# ============================================= #
# Linear model                                  #
# ============================================= #
GSE26787_fit <- lmFit(GSE26787_eset, GSE26787_design)
# ============================================= #
# Making a Contrast matrix                      #
# That is based on the design condition         #
#     'y' will be outcome to see if             #
#     'condition' has an effect                 #
# ============================================= #
GSE26787_cont_matrix <- makeContrasts(RPL_vs_Control = RPL - Control,
                                      levels = GSE26787_design)
# ============================================= #
#     Computing estimated coefficients for      #
#               Stated Contrasts                #   
# ============================================= #
GSE26787_fit2 <- contrasts.fit(GSE26787_fit, GSE26787_cont_matrix)
GSE26787_fit2 <- eBayes(GSE26787_fit2)
# ============================================= #
#     Making a table of Top Genes from          #
#               Linear Model Fit                #   
# ============================================= #
GSE26787_limma_results <- topTable(GSE26787_fit2, coef= "RPL_vs_Control", number=Inf, sort.by="P")
# ============================================= #
#               Exporting                       #
# ============================================= #
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files()  
write.table(GSE26787_limma_results, file = "GSE26787_annotated_results.tsv", quote = FALSE, sep = "\t", dec = ".", row.names=FALSE)
write.table(GSE26787_limma_results, file = "GSE26787_annotated_results.csv", quote = FALSE, sep = ";", dec = ".", row.names=FALSE)





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
png("Volcano Plot of GSE26787.png",
   width = 631,
   height = 618)
get_last_plot()
dev.off()
