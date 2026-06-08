library("dplyr")

weighted_stouffer_liptak <- function(meta_analysis_dataframe, sample_sizes, gene_id_column = "GeneID") {
  
  studies <- intersect(
    names(sample_sizes),
    sub("_logFC$", "", grep("_logFC$", names(meta_analysis_dataframe), value = TRUE)) 
  )
  
  weights <- sqrt(sample_sizes[studies])
  
  meta_Z <- apply(meta_analysis_dataframe, 1, function(row) {
    
    z_scores <- sapply(studies, function(study) {
      
      p_value <- as.numeric(row[[paste0(study, "_pvalue")]])
      log2_foldchange <- as.numeric(row[[paste0(study, "_logFC")]])
      
      if (is.na(p_value) || is.na(log2_foldchange))
        return(NA_real_)
      safe_p_value <- pmax(p_value, 1e-16)
      sign(log2_foldchange) * (qnorm(1 - safe_p_value / 2))
    })
    
    valid <- !is.na(z_scores)
    
    if (!any(valid))
      return(NA_real_)
    
    sum(weights[valid] * z_scores[valid]) /
      sqrt(sum(weights[valid]^2))
  })
  
  meta_analysis_dataframe %>%
    dplyr::select(all_of(gene_id_column)) %>%
    mutate(meta_Z = meta_Z,
           meta_p = 2 * pnorm(abs(meta_Z), lower.tail = FALSE),
           meta_FDR = p.adjust(meta_p, method = "BH"))
}