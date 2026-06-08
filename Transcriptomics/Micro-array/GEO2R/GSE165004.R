# GSE165004
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(limma)

GSE165004 <- getGEO("GSE165004", GSEMatrix = TRUE, AnnotGPL=TRUE)

GSE165004_eset <- GSE165004[[1]]

GSE165004_exprs_data <- exprs(GSE165004_eset)
GSE165004_meta_data <- pData(GSE165004_eset)

?Biobase::fData
# ============== # 
# Log2

GSE165004_exprs_data[GSE165004_exprs_data <= 0] <- NaN
exprs(GSE165004_eset) <- log2(GSE165004_exprs_data)
#===# 
annotation(GSE165004_eset)
glimpse(GSE165004_meta_data)
print(GSE165004_meta_data[, c(2, 11)])

print(GSE165004_meta_data$`subject status/group:ch1`)
print(GSE165004_meta_data["subject status/group:ch1"])

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

GSE165004_design <- model.matrix(~0 + GSE165004_condition)

colnames(GSE165004_design) <- levels(GSE165004_condition)
#===# 
GSE165004_fit <- lmFit(GSE165004_eset, GSE165004_design)
#===# 
GSE165004_cont_matrix <- makeContrasts(ControL_vs_RPL = Control - RPL,
                                       levels = GSE165004_design)

GSE165004_fit2 <- contrasts.fit(GSE165004_fit, GSE165004_cont_matrix)
# ========= # 
GSE165004_fit2 <- eBayes(GSE165004_fit2)


GSE165004_limma_results <- topTable(GSE165004_fit2, coef= "ControL_vs_RPL", number=Inf, sort.by="P")
view(GSE165004_limma_results)
# ===============# 
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
# ===============# 
write.table(GSE165004_limma_results, file = "GSE165004_ControL_vs_RPL_Annotated.csv", sep = ";", dec = ".", row.names=FALSE)

# ===============# 
hist(GSE165004_limma_results$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution for GSE165004")

#GSE165004
GSE165004_decideTests <- decideTests(GSE165004_fit2, adjust.method="fdr", p.value=0.05, lfc=0)
# ================ #
GSE165004_T_good <- which(!is.na(GSE165004_fit2$F)) # filter out bad probes
qqt(GSE165004_fit2$t[GSE165004_T_good], GSE165004_fit2$df.total[GSE165004_T_good], main="Moderated t-statistic for GSE165004")

# volcano plot (log P-value vs log fold change)
colnames(GSE165004_fit2) # list contrast names
ct <- 1        # choose contrast of interest

GSE165004_decideTests <- decideTests(GSE165004_fit2, adjust.method="fdr", p.value=0.05, lfc=0)

volcanoplot(GSE165004_fit2, 
            coef="RPL_vs_Control", 
            main=colnames(GSE165004_fit2)[1],
            pch=20,
            highlight=length(which(GSE165004_decideTests[,"RPL_vs_Control"]!=0)), 
            names=rep('+', nrow(GSE165004_fit2)))

print(GSE165004_fit2$genes$ID)
?volcanoplot
