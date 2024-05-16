#BiocManager::install(version = "3.18")  # Adjust the version as necessary
## Function to install and load packages using BiocManager
install_and_load <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (package in packages) {
    # Inform the user that the package check is starting
    message(paste("Checking and installing/loading package:", package))
    
    # Check if the package is already installed and loaded
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      # If not installed, try to install the package
      message(paste("Attempting to install the package:", package))
      tryCatch({
        suppressMessages(suppressWarnings(
          BiocManager::install(package, dependencies = TRUE, ask = FALSE)
        ))
      }, error = function(e) {
        # Provide a warning if the installation fails
        warning(paste("Failed to install", package, ":", e$message))
      })
      # After attempting installation, check if the package can now be loaded
      if (!require(package, character.only = TRUE, quietly = TRUE)) {
        # Warn if the package still cannot be loaded
        warning(paste("Failed to load", package))
      } else {
        # Confirm successful installation and loading
        message(paste(package, "installed and loaded successfully"))
      }
    } else {
      # Inform that the package was already installed and loaded
      message(paste(package, "is already installed and loaded"))
    }
  }
}

# List of required packages including both CRAN and Bioconductor
all_packages <- c(
  "BiocManager", "gplots", "RColorBrewer", "httr", "tidyverse", #Should inlcude "tidyverse" but this is a very large package and not always necasry
  "igraph", "fgsea", "stringr", "poolr", "jsonlite", "data.table",
  "DOSE", "htmltools", "rvest", "plotly", "WebGestaltR", "openxlsx",
  "msigdbr", "GOSemSim", "enrichplot", "tidytree", "rrvgo", "pRoloc",
  "treemapify", "DT", "org.Hs.eg.db", "org.Mm.eg.db", "clusterProfiler", 
  "ReactomePA", "biomaRt", "glue", "rtracklayer", "IRanges","parallel","XVector",
  "GenomicRanges","rtracklayer","flux"
)


