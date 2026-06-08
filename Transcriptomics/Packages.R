Packages_needed <- c("tidyverse",
                     "GEOquery",
                     "R.utils",
                     "DESeq2",
                     "rstudioapi",
                     "limma",
                     "STRINGdb",
                     "RCy3",
                     "clusterProfiler",
                     "org.Hs.eg.db")

install_and_load <- function(package_names) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  Available_on_BiocManager <- BiocManager::available()
  
  for (package_name in package_names) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      if (package_name %in% Available_on_BiocManager) {
        BiocManager::install(package_name, ask = FALSE, update = FALSE)
      } else {
        install.packages(package_name)
      }
    }
  }
  invisible(lapply(package_names, library, character.only = TRUE))
  message("Packages are installed & Loaded!")
}

install_and_load(Packages_needed)
select <- dplyr::select
filter <- dplyr::filter
