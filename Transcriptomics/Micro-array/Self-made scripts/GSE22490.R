# GSE22490
library(GEOquery)
library(limma)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(limma)
library(ggplot2)

GSE22490 <- getGEO("GSE22490", GSEMatrix = TRUE, AnnotGPL=TRUE)

GSE22490_eset <- GSE22490[[1]]

GSE22490_exprs_data <- exprs(GSE22490_eset)
GSE22490_meta_data <- pData(GSE22490_eset)

summary(GSE22490_exprs_data)
quantile(
  GSE22490_exprs_data,
  probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
  na.rm = TRUE
)

sum(GSE22490_exprs_data <= 0, na.rm = TRUE)


sum(GSE22490_fit2$sigma == 0)


hist(as.vector(GSE22490_exprs_data),
     breaks = 100,
     main = "GSE22490: Expression values",
     xlab = "Expression")

GSE22490_exprs_data[GSE22490_exprs_data <= 0] <- NaN
exprs(GSE22490_eset) <- log2(GSE22490_exprs_data)
# ============================================= # 
#             Annotation/Grouping               #
#               of RPL and other                #
# ============================================= # 
annotation(GSE22490_eset)
str(GSE22490_meta_data)
print(GSE22490_meta_data["source_name_ch1"])
GSE22490_condition <- c(rep("Control", 6), rep("RPL", 4))
GSE22490_condition <- factor(GSE22490_condition, levels = c("Control", "RPL"))
# ============================================= # 
#           Checking normalization              #
# ============================================= # 
group_colors <- c("Control" = "#4CAF50",              
                  "Implantation_Failure" = "#2196F3", 
                  "RPL" = "#FF5722")                  
GSE22490_colors <- group_colors[GSE22490_condition]
boxplot(exprs(GSE22490_eset), 
        outline=FALSE, 
        col = GSE22490_colors,
        main = "GSE22490: Post-Log2 Normalization Check", 
        las=2, 
        cex.axis = 0.7)
legend("topleft", legend = names(group_colors), fill = group_colors, bty = "n")
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

GSE22490_design <- model.matrix(~0 + GSE22490_condition)
colnames(GSE22490_design) <- levels(GSE22490_condition)
# ============================================= #
# Linear model                                  #
# ============================================= #
GSE22490_fit <- lmFit(GSE22490_eset, GSE22490_design)
# ============================================= #
# Making a Contrast matrix                      #
# That is based on the design condition         #
#     'y' will be outcome to see if             #
#     'condition' has an effect                 #
# ============================================= #
GSE22490_cont_matrix <- makeContrasts(RPL_vs_Control = RPL - Control,
                                       levels = GSE22490_design)
# ============================================= #
#     Computing estimated coefficients for      #
#               Stated Contrasts                #   
# ============================================= #
GSE22490_fit2 <- contrasts.fit(GSE22490_fit, GSE22490_cont_matrix)
GSE22490_fit2 <- eBayes(GSE22490_fit2)

# ============================================= #
#     Making a table of Top Genes from          #
#               Linear Model Fit                #   
# ============================================= #
GSE22490_limma_results <- topTable(GSE22490_fit2, coef= "RPL_vs_Control", number=Inf, sort.by="P")
# ============================================= #
#               Exporting                       #
# ============================================= #
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
# ============================================= #
write.table(GSE22490_limma_results, file = "GSE22490_Control_vs_RPL_Annotated.csv", sep = ";", dec = ".", row.names=FALSE)

# ============================================= # 
#                                               #
#                 Making Plots                  #
#                                               #
# ============================================= #





# ========================================================== #
# Mean variance trend or the Sigma-vs-A plot                 #
# ========================================================== #
plotSA(GSE22490_fit2, main="GSE22490: Mean-Variance Trend")

# ========================================================== #
# Histogram
# ========================================================== #
hist(GSE_RPLStudy_results$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "GSE_RPLStudy Frequencies of padj-values")
# ========================================================== #
# Volcano Plot
# ========================================================== #
GSE22490_limma_results <- GSE22490_limma_results %>%
  mutate(threshold = case_when(
    P.Value < 0.05 & logFC > 0  ~ "Upregulated",
    P.Value < 0.05 & logFC < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))

ggplot(GSE22490_limma_results,
       aes(x = logFC,
           y = -log10(P.Value),
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
    title = "Volcano Plot of GSE22490: Controls vs RPL",
    x = "log2 fold change",
    y = expression(-log[10]("p-value")),
    color = "Gene status"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  )
