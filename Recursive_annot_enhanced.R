
# Helper function for QC to filter splice junctions based on read size and number of unique reads.
filter_splicejunct_QC <- function(sj_df, mx_sj_overhang = 20, nr_of_uniqreads = 5) {
  sj_df %>%
    filter(V9 >= mx_sj_overhang, V7 >= nr_of_uniqreads)
}
sj_annot <- function(sj_tab_folder        = cwd,
                     defined_introns_file = str_glue("{cwd}/defined_introns.bed"),
                     sj_annot_out_folder  = str_glue("{cwd}/sj_annot_out_folder")){
  
  dir.create(sj_annot_out_folder, recursive = TRUE, showWarnings = TRUE)
  
  # extract all splice junction files
splice_files <- list.files(path = sj_tab_folder, pattern = "*.out.tab", recursive = TRUE, full.names = TRUE)
# process all splice junction files
splice_data <- lapply(splice_files, function(file) {
  fread(file, header = FALSE) %>%
    filter_splicejunct_QC()
}) %>%
  bind_rows() %>%
  mutate(n = 1) %>%
  group_by(V1, V2, V3, V4, V5, V6) %>%
  summarise_all(sum) %>%
  filter(n >= length(splice_files))

# Write the filtered splice junctions to an output file
filterd_sj_outfile <- str_glue("{sj_annot_out_folder}/SJ.complete.out.tab")
write_tsv(splice_data, filterd_sj_outfile, col_names = FALSE)

# Intersect the splice junction file with the intron file using system command
intersect_command <- str_glue("bedtools intersect -wa -wb -a {filterd_sj_outfile} -b {defined_introns_file} > {sj_annot_out_folder}/sj_annot_all.bed")
system(intersect_command)

# Read the intersection output and perform further analysis
junction_data <- read_tsv(str_glue("{sj_annot_out_folder}/sj_annot_all.bed"), col_names = FALSE)
colnames(junction_data)  <- column_names <- c(
  "chromosome_measured",
  "start_position_measured",
  "end_position_measured",
  "strand_measured",
  "intron_motif",
  "annotated",
  "unique_mapping_count",
  "multi_mapping_count",
  "maximum_overhang",
  "junction_id",
  "gtf_refrnce_chromosome",
  "start_position_gtf_refrnce",
  "end_position_gtf_refrnce",
  "feature_type",
  "gene_id",
  "strand_gtf_refrnce"
)

junction_data <- mutate(junction_data, length = end_position_gtf_refrnce - start_position_gtf_refrnce)

junction_data <- junction_data %>%
  mutate(strand_measured = case_when(
    strand_measured == 1 ~ "+",
    strand_measured == 2 ~ "-",
    TRUE ~ as.character(strand_measured)
  ))

# Process and save final output
results_no_rec <- junction_data %>%
  filter(strand_measured == strand_gtf_refrnce, start_position_measured == start_position_gtf_refrnce, end_position_measured == end_position_gtf_refrnce) %>%
  write_tsv(str_glue("{sj_annot_out_folder}/NoRecIntronJunction.bed"))

results_rec <- junction_data %>%
  filter(strand_measured == strand_gtf_refrnce, start_position_measured != start_position_gtf_refrnce, end_position_measured != end_position_gtf_refrnce) %>%
  write_tsv(str_glue("{sj_annot_out_folder}/RecIntronJunction.bed"))

# Output the first few lines of the processed data for verification
message("Output of the first few lines of the processed data for verification")
print(head(results_no_rec))
print(head(results_rec))
}

#This part of the script is to run the functions if this scirpt is called in the commandline
# Check if the call stack has more than 1 frame, indicating sourcing
if (sys.nframe() == 0) {
  
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
    "ReactomePA", "biomaRt", "glue"
  )
  
  install_and_load(all_packages)
  
# Load necessary libraries
cwd <- getwd()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("No input file specified.", call. = FALSE)}
  
sj_annot(sj_tab_folder        = args[1],
         defined_introns_file = args[2],
         sj_annot_out_folder  = args[3])  # Set default arguments and read them
}
sj_annot("/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/sj_tab_files/", str_glue("/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/Intron_selection/Drosophila_melanogaster.BDGP6.90.intron.bed"), str_glue("{cwd}/sj_annot_out_folder"))
