# GEO_ID_BULK
# ========================= #
library(dplyr)
library(purrr)
# ========================= #
exists("GSE65099_annotated_results")
exists("GSE161969_annotated_results")
exists("GSE113790_annotated_results")
exists("GSE183555_annotated_results")
# ========================= #
view(GSE65099_annotated_results)
str(GSE183555_annotated_results)
# ========================= #
view(Transcriptomics_Bulk)
# ========================= #
# Bulk-Pregnant
Transcriptomics_Bulk %>% 
  filter(Pregnant_or_NonPregnant == "Pregnant") %>%
  distinct(Dataset_Name_ID, .keep_all = TRUE) %>% 
  select(Dataset_Name_ID, `Phase (when the sample was taken)`, Sample_Source )


# Bulk-non-Pregnant
Transcriptomics_Bulk %>% 
  filter(Pregnant_or_NonPregnant != "Pregnant") %>%
  distinct(Dataset_Name_ID, .keep_all = TRUE) %>% 
  select(Dataset_Name_ID, `Phase (when the sample was taken)`, Sample_Source )

# ========================= #
# Micro_Array-Pregnant
Transcriptomics_Micro_Array %>% 
  filter(Pregnant_or_NonPregnant == "Pregnant") %>%
  distinct(Dataset_Name_ID, .keep_all = TRUE) %>% 
  select(Dataset_Name_ID, `Phase (when the sample was taken)`, Sample_Source )

# Micro_Array-non-Pregnant
Transcriptomics_Micro_Array %>% 
  filter(Pregnant_or_NonPregnant != "Pregnant") %>%
  distinct(Dataset_Name_ID, .keep_all = TRUE) %>% 
  select(Dataset_Name_ID, `Phase (when the sample was taken)`, Sample_Source )
# ========================= #
# Endometrial
# ========================= #
Endometrial_Bulk_GSE65099 <- GSE65099_annotated_results %>% 
  select(GeneID = EntrezID, 
         GSE65099_logFC =log2FoldChange,
         GSE65099_pvalue  = pvalue)

Endometrial_Bulk_GSE183555 <- GSE113790_annotated_results %>% 
  select(GeneID = EntrezID,
         GSE183555_logFC = log2FoldChange, 
         GSE183555_pvalue  = pvalue)
# ========================= #
Endometrial_microarray_GSE26787 <- GSE26787_limma_results %>% 
  select(GeneID = Gene.ID,
         GSE26787_logFC = logFC, 
         GSE26787_pvalue = P.Value)

# option 2
Endometrial_microarray_GSE26787 <- GSE26787_limma_results %>% 
  group_by(Gene.ID) %>%
  arrange(P.Value) %>% 
  slice_head() %>%         
  ungroup() %>%
  select(GeneID = Gene.ID, 
         GSE26787_logFC = logFC,
         GSE26787_pvalue = P.Value) 
# ========================= #
Endometrial_microarray_GSE165004 <- GSE165004_limma_results %>%
  mutate(GeneID = as.character(LOCUSLINK_ID)) %>% 
  select(GeneID,
         GSE165004_logFC = logFC, 
         GSE165004_pvalue = P.Value)

Endometrial_microarray_GSE165004 <- GSE165004_limma_results %>% 
  mutate(GeneID = as.character(LOCUSLINK_ID)) %>% 
  group_by(GeneID) %>%
  arrange(P.Value) %>% # Sort by p-value first!
  slice_head()         %>%       
  ungroup() %>%
  select(GeneID, 
         GSE165004_logFC = logFC, 
         GSE165004_pvalue = P.Value)

# ========================= #
endometrial_metaanalysis_table <- list(
  Endometrial_Bulk_GSE65099,
  Endometrial_Bulk_GSE183555, 
  Endometrial_microarray_GSE26787, 
  Endometrial_microarray_GSE165004) %>%
  purrr::reduce(full_join, 
         by = "GeneID") #%>% 
  #distinct(GeneID, .keep_all = TRUE)
str(endometrial_metaanalysis_table)

endometrial_metaanalysis_table <- list(
  Endometrial_Bulk_GSE65099,
  Endometrial_Bulk_GSE183555, 
  Endometrial_microarray_GSE26787, 
  Endometrial_microarray_GSE165004) %>%
  purrr::reduce(full_join, 
                by = "GeneID") # %>% 
  
# ========================= #
Decidual_Bulk_GSE161969 <- GSE161969_annotated_results %>% 
  select(Symbol, EntrezID,
         GSE161969_logFC = log2FoldChange, 
         GSE161969_padj  = padj)

Decidual_Bulk_GSE113790 <- GSE113790_annotated_results %>% 
  select(Symbol, EntrezID,
         GSE113790_logFC = log2FoldChange, 
         GSE113790_padj  = padj)
# ========================= #
# Placenta/trophoblastic tissue
# ========================= #
Placenta_microarray_GSE22490 <- GSE22490_limma_results %>% 
  select(Gene.ID,
         GSE22490_logFC = logFC, 
         GSE22490_padj  = adj.P.Val)

# ========================= #
endometrial_fisher_meta_analysis <- endometrial_metaanalysis_table %>% 
  # To see the direction
  mutate( 
    DEGs_direction_GSE65099  = sign(GSE65099_logFC),
    DEGs_direction_GSE183555 = sign(GSE183555_logFC),
    DEGs_direction_GSE26787  = sign(GSE26787_logFC),
    DEGs_direction_GSE165004 = sign(GSE165004_logFC)
  ) %>% 
  rowwise() %>%
  # Remove NA
  mutate(
    valid_DEGs_directions = list(
      na.omit(
        c(
          DEGs_direction_GSE65099,
          DEGs_direction_GSE183555,
          DEGs_direction_GSE26787,
          DEGs_direction_GSE165004
          )
        )
      ),
    is_directionally_consistent = length(unique(valid_DEGs_directions)) <= 1,
    nominal_p_values = list(
      na.omit(
        c(GSE65099_pvalue, 
          GSE183555_pvalue, 
          GSE26787_pvalue, 
          GSE165004_pvalue))),
    k_studies = length(nominal_p_values)) %>% 
  ungroup() %>% 
  # filtering for genes found in at least 2 studies
  filter(k_studies >= 2 & is_directionally_consistent) %>%  # %>% 
  rowwise() %>% 
  mutate(
    fisher_stat = -2 * sum(log(nominal_p_values)),
    # Chi-squared test with 2Xk degrees of freedom (k being the number of studies)
    meta_p_value = pchisq(fisher_stat, df = 2 * k_studies, lower.tail = FALSE)
  ) %>%
  ungroup() %>% 
  mutate(meta_FDR = p.adjust(meta_p_value, method = "BH")) %>%
  arrange(meta_FDR)

