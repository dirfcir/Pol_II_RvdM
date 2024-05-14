# Helper function for QC to filter splice junctions based on read size and number of unique reads.
filter_splicejunct_QC <- function(sj_df, mx_sj_overhang = 20, nr_of_uniqreads = 5) {
  sj_df %>%
    filter(V9 >= mx_sj_overhang, V7 >= nr_of_uniqreads)
}

sj_annot <- function(sj_tab_folder        = cwd,
                     defined_introns_file = str_glue("{cwd}/defined_introns.bed"),
                     sj_annot_out_folder  = str_glue("{cwd}/sj_annot_out_folder")) {
  
  dir.create(sj_annot_out_folder, recursive = TRUE, showWarnings = TRUE)
  
  message("Extracting splice junction files...")
  # Extract all splice junction files
  splice_files <- list.files(
    path = sj_tab_folder,
    pattern = "*.out.tab",
    recursive = TRUE,
    full.names = TRUE
  )
  
  message("Processing splice junction files...")
  # Process all splice junction files
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
  
  message("Running bedtools intersect...")
  # Intersect the splice junction file with the intron file using system command
  intersect_command <- str_glue("bedtools intersect -wa -wb -a {filterd_sj_outfile} -b {defined_introns_file} > {sj_annot_out_folder}/sj_annot_all.bed")
  system(intersect_command)
  
  message("Reading intersection output...")
  # Read the intersection output and perform further analysis
  junction_data <- read_tsv(str_glue("{sj_annot_out_folder}/sj_annot_all.bed"), col_names = FALSE)
  colnames(junction_data) <- c(
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
  
  message("Processing and saving final output...")
  # Process and save final output
  common_start_end <- junction_data %>%
    filter(
      strand_measured == strand_gtf_refrnce,
      abs(start_position_measured - start_position_gtf_refrnce) <= 1,
      abs(end_position_measured - end_position_gtf_refrnce) <= 1
    )
  
  common_start <- junction_data %>%
    filter(
      strand_measured == strand_gtf_refrnce,
      abs(start_position_measured - start_position_gtf_refrnce) <= 1,
      end_position_measured - end_position_gtf_refrnce <= -2 
    )
  
  common_end <- junction_data %>%
    filter(
      strand_measured == strand_gtf_refrnce,
      start_position_measured - start_position_gtf_refrnce >= 2, 
      abs(end_position_measured - end_position_gtf_refrnce) <= 1
    )
  
  full_diff <- junction_data %>%
    filter(
      strand_measured == strand_gtf_refrnce,
      start_position_measured - start_position_gtf_refrnce >= 2, 
      end_position_measured - end_position_gtf_refrnce <= -2  
    )
  
  results_no_recursive_sj <- common_start_end 
  results_recursive_sj <- bind_rows(common_start, common_end, full_diff)
  
  # Output the first few lines of the processed data for verification
  message("Output of the first few lines of the processed data for verification")
  print(head(results_no_recursive_sj[, c("chromosome_measured", "start_position_measured", "end_position_measured", "strand_measured", "unique_mapping_count", "junction_id", "intron_motif", "annotated")]))
  print(head(results_recursive_sj[, c("chromosome_measured", "start_position_measured", "end_position_measured", "strand_measured", "unique_mapping_count", "junction_id", "intron_motif", "annotated")])) 
  
  # Final output
  write_tsv(
    results_no_recursive_sj[, c("chromosome_measured", "start_position_measured", "end_position_measured", "strand_measured", "junction_id", "intron_motif")], 
    str_glue("{sj_annot_out_folder}/no_recursiv_intron_sj.bed"), col_names = FALSE
  )
  
  write_tsv(
    results_recursive_sj[, c("chromosome_measured", "start_position_measured", "end_position_measured", "strand_measured", "junction_id", "intron_motif")], 
    str_glue("{sj_annot_out_folder}/recursiv_intron_sj.bed"), col_names = FALSE
  )
}

# This part of the script is to run the functions if this script is called in the command line
if (sys.nframe() == 0) {
  args <- commandArgs()
  cwd <- getwd()
  
  find_folder_or_file <- function(files_to_find, max_steps = 4, folder = FALSE) {
    current_dir <- normalizePath(cwd, winslash = "/")
    for (i in seq_len(max_steps + 1)) {
      # List all files and directories recursively in the current directory
      all_entries <- list.files(path = current_dir, recursive = TRUE, full.names = TRUE)
      
      # Filter to include only those entries that match the pattern
      matches <- all_entries[grepl(files_to_find, all_entries)]
      
      if (identical(matches, character(0))) {
        # Move up one directory level
        current_dir <- dirname(current_dir)
      } else {
        if (folder) {
          path_to_folder <- dirname(matches[1])
          return(path_to_folder)
        } else {
          paths_to_files <- matches
          return(paths_to_files)
        }
      }
    }
    stop(paste("Directory", folder_to_find, "not found"))
  }
  
  cat("Usage: Rscript Recursive_annot_enhanced.R [options]\n",
      "Options:\n",
      "  sj_tab_folder                  Path to folder containing SJ tab files (default:", find_folder_or_file("/sj_tab_files", folder = TRUE), "\n", 
      "  defined_introns_file           Path to file with defined introns (default:", find_folder_or_file("*intron.bed"), "\n",
      "  sj_annot_out_folder            Path to output folder for annotations of splice junctions (default:", file.path(cwd, "sj_annot_out_folder/"), "\n",
      "\nExample:\n",
      "  Rscript Recursive_annot_enhanced.R /path/to/sj_tab_folder/ /path/to/defined_introns.bed /path/to/intron.bed /path/to/output_folder/\n\n")
  
  # Checking number of arguments
  if (length(args) != 3) {
    cat("CAREFUL: NOT ALL ARGUMENTS ARE PROVIDED.\n Proceeding with default values:\n",
        "sj_tab_folder: ", find_folder_or_file("/sj_tab_files", folder = TRUE), "\n", 
        "defined_introns_file:", find_folder_or_file("*intron.bed"), "\n",
        "sj_annot_out_folder:", file.path(cwd, "sj_annot_out_folder/"), "\n")
    
    args[1] <- find_folder_or_file("/sj_tab_files", folder = TRUE)
    args[2] <- find_folder_or_file("*intron.bed")
    args[3] <- file.path(cwd, "sj_annot_out_folder/")
  } else {
    message("All arguments provided. Proceeding with execution...")
  }
  
  required_scripts <- find_folder_or_file("R_scripts/required")
  
  message("Sourcing required scripts...")
  for (file in required_scripts) {
    message("Sourcing file: ", file)
    source(file)
  }
  
  install_and_load(all_packages)
  message("All packages installed and loaded")
  
  sj_annot(sj_tab_folder        = args[1],
           defined_introns_file = args[2],
           sj_annot_out_folder  = args[3])  # Set default arguments and read them
}
