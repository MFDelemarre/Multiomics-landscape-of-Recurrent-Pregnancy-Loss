# GSE_RPLStudy
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(limma)

# ============================================= # 
#                 Loading dataset               #
# ============================================= # 


GSE_RPLStudy <- getGEO("GSE_RPLStudy", GSEMatrix = TRUE, AnnotGPL=TRUE)

GSE_RPLStudy_eset <- GSE_RPLStudy[[1]]

GSE_RPLStudy_exprs_data <- exprs(GSE_RPLStudy_eset)
GSE_RPLStudy_meta_data <- pData(GSE_RPLStudy_eset)

# ============================================= # 
#             Annotation/Grouping               #
#               of RPL and other                #
# ============================================= # 

annotation(GSE_RPLStudy_eset)
glimpse(GSE_RPLStudy_meta_data)
print(GSE_RPLStudy_meta_data[, c(2, 11)])

print(GSE_RPLStudy_meta_data$`subject status/group:ch1`)
print(GSE_RPLStudy_meta_data["subject status/group:ch1"])

GSE_RPLStudy_condition <- c(rep("Control", 24), rep("RPL", 24), rep("Implantation_Failure", 24))
GSE_RPLStudy_condition <- factor(GSE_RPLStudy_condition, levels = c("Control", "Implantation_Failure", "RPL"))





# Check whether the GEO dataset was already on a log scale before applying log2().
hist(as.vector(GSE_RPLStudy_exprs_data),
     breaks = 100,
     main = "GSE_RPLStudy_exprs_data: Expression values",
     xlab = "Expression")

# ============================================= # 
#        Log2Fold-Change of RPL and other       #
# ============================================= # 
GSE_RPLStudy_exprs_data[GSE_RPLStudy_exprs_data <= 0] <- NaN
exprs(GSE_RPLStudy_eset) <- log2(GSE_RPLStudy_exprs_data)

# ============================================= # 
#          Boxplot to see normalization         #
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
#        Log2Fold-Change of RPL and other       #
# ============================================= # 
GSE_RPLStudy_exprs_data[GSE_RPLStudy_exprs_data <= 0] <- NaN
exprs(GSE_RPLStudy_eset) <- log2(GSE_RPLStudy_exprs_data)
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
GSE_RPLStudy_design <- model.matrix(~0 + GSE_RPLStudy_condition)
colnames(GSE_RPLStudy_design) <- levels(GSE_RPLStudy_condition)
# ============================================= #
# Linear model                                  #
# ============================================= #

GSE_RPLStudy_fit <- lmFit(GSE_RPLStudy_eset, 
                          GSE_RPLStudy_design)
# ============================================= #
# Making a Contrast matrix                      #
# That is based on the design condition         #
#     'y' will be outcome to see if             #
#     'condition' has an effect                 #
# ============================================= #
GSE_RPLStudy_cont_matrix <- makeContrasts(RPL_vs_ControL = RPL - Control,
                                          levels = GSE_RPLStudy_design)
# ============================================= #
#     Computing estimated coefficients for      #
#               Stated Contrasts                #   
# ============================================= #
GSE_RPLStudy_fit2 <- contrasts.fit(GSE_RPLStudy_fit, GSE_RPLStudy_cont_matrix)
GSE_RPLStudy_fit2 <- eBayes(GSE_RPLStudy_fit2)
# ============================================= #
#     Making a table of Top Genes from          #
#               Linear Model Fit                #   
# ============================================= #
GSE_RPLStudy_limma_results <- topTable(GSE_RPLStudy_fit2, coef= "RPL_vs_ControL", number=Inf, sort.by="P")
# ============================================= #
#               Exporting                       #
# ============================================= #
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
# ============================================= #
write.table(GSE_RPLStudy_limma_results, file = "GSE_RPLStudy_RPL_vs_ControL_Annotated.csv", sep = ";", dec = ".", row.names=FALSE)







# ============================================= # 
#                                               #
#                 Making Plots                  #
#                                               #
# ============================================= # 






plotSA(GSE_RPLStudy_fit2, main="GSE_RPLStudy: Mean-Variance Trend")

# ===============# 
hist(GSE_RPLStudy_limma_results$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution for GSE_RPLStudy")


GSE_RPLStudy_decideTests <- decideTests(GSE_RPLStudy_fit2, adjust.method="fdr", p.value=0.05, lfc=0)
# ================ #
GSE_RPLStudy_T_good <- which(!is.na(GSE_RPLStudy_fit2$F)) # filter out bad probes
qqt(GSE_RPLStudy_fit2$t[GSE_RPLStudy_T_good], GSE_RPLStudy_fit2$df.total[GSE_RPLStudy_T_good], main="Moderated t-statistic for GSE_RPLStudy")


# ggplot
# Making the ggplot
GSE_RPLStudy_limma_results <- GSE_RPLStudy_limma_results %>%
  mutate(threshold = case_when(
    P.Value < 0.05 & logFC > 0  ~ "Upregulated",
    P.Value < 0.05 & logFC < 0  ~ "Downregulated",
    TRUE                              ~ "Not significant"
  ))

ggplot(GSE_RPLStudy_limma_results,
       aes(x = logFC,
           y = -log10(P.Value),
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
  
  labs(
    title = "Volcano Plot of GSE_RPLStudy: RPL vs Controls",
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
# ========== #
