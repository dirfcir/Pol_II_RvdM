# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Function to select introns for slope calculation based on various criteria          #
# Args:                                                                               #
#   intron_info_file: Data frame containing information on introns                    #
# Returns:                                                                            #
#   Data frame filtered based on specific criteria, for downstream slope calculation. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
select_introns_for_slopecalc <- function(intron_info_file) 
  {
  intron_info_file <- intron_info_file %>%
    filter(overlap_outside_reference != TRUE) %>%      # Exclude introns overlapping outside reference
    filter(supported_by_samples == 1) %>%              # Include only introns supported by 100% of samples
    filter(strand_measured != 0) %>%                   # Include only introns where the strand annotation is conclusive +/-
    filter((start_match_reference == TRUE & end_match_reference == TRUE) |  # Fully match reference intro
             (start_match_reference == TRUE & end_match_reference == FALSE) | # Match start OR end position of reference intro
             (start_match_reference == FALSE & end_match_reference == TRUE))
  return(intron_info_file)  # Return the filtered intron info file
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function to fit a linear model to the coverage data of the RNAseq to an intron.                                              #
# Args:                                                                                                        #
#   chromosome_intron_coverage: List of coverage data for introns.                                                        #
#   component: Component of the linear model to return (default is "all").                                     #
# Returns:                                                                                                     #
#   Data frame with lm info on y-intercept, x-intercept, slope, p-value, and sample filename for each intron.  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
compute_lm <- function(coverage_per_intron, 
                       component = "all") 
  {
  if (is.null(coverage_per_intron)) { return(NULL) }  # Return NULL if the input is NULL
  
  result <- t(sapply(1:length(coverage_per_intron), function(sample_index) 
    {
    fit <- lm(as.numeric(coverage_per_intron[[sample_index]]) ~ seq(1:length(coverage_per_intron[[sample_index]])))  # Fit linear model seq() to give positions to coverage values
    summary_fit <- summary(fit)  # Get summary of the fit
    y_intercept <- fit$coefficients[1]  # Extract y-intercept
    x_intercept <- fit$fitted.values[[length(coverage_per_intron[[sample_index]])]]  # Extract x-intercept
    slope <- summary_fit$coefficients[2, 1]  # Extract slope
    p_value <- summary_fit$coefficients[2, 4]  # Extract p-value for the slope
    sample_filename <- as.character(metadata(coverage_per_intron[[sample_index]])$sample_filename)  # Extract sample filename
    
    if (component == "y_intercept") {
      return(c(y_intercept, p_value, sample_filename))
    } else if (component == "x_intercept") {
      return(c(x_intercept, p_value, sample_filename))
    } else if (component == "slope") {
      return(c(slope, p_value, sample_filename))
    } else if (component == "all") {
      return(c(y_intercept, x_intercept, slope, p_value, sample_filename))
    }
  }))
  
  if (component == "all") {
    colnames(result) <- c("y_intercept", "x_intercept", "slope", "p_value", "sample_filename")
  } else {
    colnames(result) <- c("value", "p_value", "sample_filename")
  }
  
  return(as.data.frame(result))  # Return the result as a dataframe
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #                  
# Function to get coverage data for a specific chromosome                               #
# Args:                                                                                 #
#   sample_index: Index of the coverage data list.                                      #
#   chromosome: Chromosome name to filter coverage data.                                #
#   bed_graph: List of coverage data.                                                   #
# Returns:                                                                              #
#   Coverage data for the specified chromosome with metadata (filename and chromosome). #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #                  
get_coverage_per_chromosome <- function(sample_index, 
                                        chromosome, 
                                        bed_graph) 
  {
  chromosome_data_for_coverage_calc <- bed_graph[[sample_index]][seqnames(bed_graph[[sample_index]]) == chromosome]  # Subset data for the chromosome
  chromosome_coverage_data <- coverage(chromosome_data_for_coverage_calc, weight = chromosome_data_for_coverage_calc$score)[[chromosome]]  # Compute coverage
  
  # Extract the sample_filename from the metadata
  sample_filename <- unique(mcols(chromosome_data_for_coverage_calc)$sample_filename)
  metadata(chromosome_coverage_data)$sample_filename <- sample_filename
  metadata(chromosome_coverage_data)$chromosome <- chromosome
  
  return(chromosome_coverage_data)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function to filter introns by chromosome and strand                       #
# Args:                                                                     #
#   intron_info_file: Data frame containing intron information.             #
#   chromosome: Chromosome name to filter introns.                          #
#   strand_sense: Strand sense to filter introns ("all", "plus", "minus").  #
# Returns:                                                                  #
#   Data frame filtered by chromosome and strand.                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
filter_introns_by_chromosome_and_strand <- function(intron_info_file, chromosome, strand_sense = '') 
  {
  if (strand_sense == 'all') {
    introns_per_strand_and_chromosome <- intron_info_file[(intron_info_file$chromosome == chromosome),]
  }
  
  if (strand_sense == 'minus') {
    introns_per_strand_and_chromosome <- intron_info_file[(intron_info_file$chromosome_measured == chromosome) & (intron_info_file$strand_measured == "-"),]
  }
  if (strand_sense == 'plus') {
    introns_per_strand_and_chromosome <- intron_info_file[(intron_info_file$chromosome_measured == chromosome) & (intron_info_file$strand_measured == "+"),]
  }
  return(introns_per_strand_and_chromosome)  # Return filtered introns
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function to import bedGraph files in parallel                   #
# Args:                                                           #
#   bedgraph_filepaths: Vector of file paths to bedGraph files.   #
#   nr_of_cores: Number of cores to use for parallel processing.  #
# Returns:                                                        #
#   List of imported bedGraph files with filenames in metadata.   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
paralel_import_bed_graphs <- function(bedgraph_filepaths, 
                                      nr_of_cores) 
  {
  imported_bed_graphs <- mclapply(bedgraph_filepaths, function(bedgraph_filepath) 
    {
    imported_bed_graph <- import.bedGraph(bedgraph_filepath)  # Import bedGraph file
    mcols(imported_bed_graph)$sample_filename <- basename(bedgraph_filepath)  # Add sample filename to metadata
    return(imported_bed_graph)
  }, mc.cores = nr_of_cores)  # Use multiple cores for parallel processing
  return(imported_bed_graphs)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function to get intron coverage data for a specific chromosome  #
# Args:                                                           #
#   chromosome: Chromosome name to filter coverage data.          #
#   bed_graphs_data_fwd: List of forward strand bedGraph data.    #
#   bed_graphs_data_rev: List of reverse strand bedGraph data.    #
#   intron_info: Data frame containing information on introns     #
#   nr_of_cores: Number of cores to use for parallel processing.  #
# Returns:                                                        #
#   List of intron coverage data for the specified chromosome     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
get_intron_coverage_per_chromosome <- function(chromosome,
                                               bed_graphs_data_fwd,
                                               bed_graphs_data_rev, 
                                               intron_info, nr_of_cores = 6) 
  { 
  # Subselect introns for chromosome and strand
  introns_per_strand_and_chromosome_fwd <- filter_introns_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "plus")
  introns_per_strand_and_chromosome_rev <- filter_introns_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "minus")
  
  # Get coverage data for forward and reverse strands
  coverage_data_fwd <- lapply(1:length(bed_graphs_data_fwd), get_coverage_per_chromosome, chromosome, bed_graphs_data_fwd)        
  coverage_data_rev <- lapply(1:length(bed_graphs_data_rev), get_coverage_per_chromosome, chromosome, bed_graphs_data_rev)
  
  size_of_annot_fwd <- nrow(introns_per_strand_and_chromosome_fwd)
  size_of_annot_rev <- nrow(introns_per_strand_and_chromosome_rev)
  
  # Forward strand coverage data with metadata assignment
  intron_coverage_data_fwd <- lapply(1:length(coverage_data_fwd), function(sample_index) 
    {
    mcmapply(function(start_position, end_position) 
      {
      cov_window <- window(coverage_data_fwd[[sample_index]], start_position, end_position)
      metadata(cov_window)$sample_filename <- metadata(coverage_data_fwd[[sample_index]])$sample_filename
      metadata(cov_window)$chromosome <- metadata(coverage_data_fwd[[sample_index]])$chromosome
      return(cov_window)
    }, introns_per_strand_and_chromosome_fwd$start_position_measured,
      introns_per_strand_and_chromosome_fwd$end_position_measured,
      mc.cores = nr_of_cores)
  })
  
  # Reverse strand coverage data with metadata assignment
  intron_coverage_data_rev <- lapply(1:length(coverage_data_rev), function(sample_index) 
    {
    mcmapply(function(start_position, end_position) 
      {
      cov_window <- rev(window(coverage_data_rev[[sample_index]], start_position, end_position))
      metadata(cov_window)$sample_filename <- metadata(coverage_data_rev[[sample_index]])$sample_filename
      metadata(cov_window)$chromosome <- metadata(coverage_data_fwd[[sample_index]])$chromosome
      return(cov_window)
    }, introns_per_strand_and_chromosome_rev$start_position_measured,
      introns_per_strand_and_chromosome_rev$end_position_measured,
      mc.cores = nr_of_cores)
  })
  
  # Combine forward and reverse strand data
  if (length(intron_coverage_data_fwd[[1]]) == 0) intron_coverage_data_fwd <- NULL else intron_coverage_data_fwd <- lapply(1:size_of_annot_fwd, function(y) lapply(intron_coverage_data_fwd, function(x) x[[y]]))
  if (length(intron_coverage_data_rev[[1]]) == 0) intron_coverage_data_rev <- NULL else intron_coverage_data_rev <- lapply(1:size_of_annot_rev, function(y) lapply(intron_coverage_data_rev, function(x) x[[y]]))
  
  chromosome_intron_coverage_data <- c(intron_coverage_data_fwd, intron_coverage_data_rev)  # Combine forward and reverse strand data
  return(chromosome_intron_coverage_data)  # Return combined data
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Function to fit linear models to intron coverage data for a specific chromosome                                      #
# Args:                                                                                                                #
#   chromosome: Chromosome name to filter coverage data.                                                               #
#   intron_coverage_data: List of intron coverage data for different chromosomes.                                      #
#   intron_info: Data frame containing intron information.                                                             #
#   nr_of_cores: Number of cores to use for parallel processing (default is 6).                                        #
#   component_from_lm: Component of the linear model to return (default is "all").                                     #
# Returns:                                                                                                             #
#   Data frame with intron information combined with linear model results (y-intercept, x-intercept, slope, p-value).  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # # # # # # #
lm_fit_intron_coverage_per_chromosome <- function(chromosome, 
                                                  intron_coverage_data, 
                                                  intron_info,
                                                  nr_of_cores = 6, 
                                                  component_from_lm = "all") 
  {
  # Subselect introns for chromosome and strand
  introns_per_strand_and_chromosome_fwd <- filter_introns_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "plus")
  introns_per_strand_and_chromosome_rev <- filter_introns_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "minus")
  
  chromosome_intron_coverage_data <- intron_coverage_data[[chromosome]]
  chromosome_intron_coverage_lm_info <- mclapply(chromosome_intron_coverage_data, compute_lm, component = component_from_lm, mc.cores = nr_of_cores)  # Fit linear models in parallel
  
  cat(paste(format(object.size(chromosome_intron_coverage_lm_info), units = "MB"), "\n"))  # Print the size of the object
  
  chromosome_intron_coverage_lm_info[chromosome_intron_coverage_lm_info == "NULL"] <- NA  # Replace NULLs with NA
  
  chromosome_intron_coverage_lm_info <- do.call(rbind, chromosome_intron_coverage_lm_info)  # Combine results into a dataframe
  
  introns_per_chromosome <- rbind(introns_per_strand_and_chromosome_fwd, introns_per_strand_and_chromosome_rev)  # Combine forward and reverse strand introns
  
  chromosome_intron_coverage_lm_info <- cbind(introns_per_chromosome, chromosome_intron_coverage_lm_info)  # Combine intron info with linear model results
  
  return(chromosome_intron_coverage_lm_info)  # Return the combined data
}

# Function to execute all other functions and calculate the slope of the intron coverage data
calculate_slopes <- function(bedgraph_folder, 
                             annoted_intron_sj_folder, 
                             outfolder_slope_calculations, 
                             chromosome_names, 
                             nr_of_cores_for_parallel_computing) 
{
  cat("Starting slope calculation process...\n")
  
  # Check if the output_directory exists, and create it if it doesn't
  if (!dir.exists(outfolder_slope_calculations)) {
    dir.create(outfolder_slope_calculations, recursive = TRUE)
    cat("Created output directory.\n")
  } else {
    cat("Output directory already exists.\n")
  }
  
  cat("Reading bedGraph file paths...\n")
  bedgraph_fwd_file_pattern <- str_glue("{bedgraph_folder}*Unique.str2.out.bg")
  bedgraph_fwd_filepaths <- Sys.glob(bedgraph_fwd_file_pattern)
  
  bedgraph_rev_file_pattern <- str_glue("{bedgraph_folder}*Unique.str1.out.bg")
  bedgraph_rev_filepaths <- Sys.glob(bedgraph_rev_file_pattern)
  cat("BedGraph file paths read.\n")
  
  cat("Reading intron information...\n")
  intron_info_df <- read_tsv(str_glue("{annoted_intron_sj_folder}intron_sj_annotated_with_reference.tsv"))
  selected_intron_info_df <- select_introns_for_slopecalc(intron_info_df)
  cat("Intron information read and filtered.\n")
  
  cat("Importing bedGraph data in parallel...\n")
  bed_graphs_data_fwd <- paralel_import_bed_graphs(bedgraph_filepaths = bedgraph_fwd_filepaths, nr_of_cores = nr_of_cores_for_parallel_computing)
  bed_graphs_data_rev <- paralel_import_bed_graphs(bedgraph_filepaths = bedgraph_rev_filepaths, nr_of_cores = nr_of_cores_for_parallel_computing)
  cat("BedGraph data imported.\n")
  
  cat("Getting intron coverage data for each chromosome...\n")
  intron_coverage_data <- lapply(chromosome_names, get_intron_coverage_per_chromosome, bed_graphs_data_fwd = bed_graphs_data_fwd, bed_graphs_data_rev = bed_graphs_data_rev, intron_info = selected_intron_info_df, nr_of_cores = nr_of_cores_for_parallel_computing)
  names(intron_coverage_data) <- chromosome_names
  cat("Intron coverage data collected.\n")
  intron_coverage_data_all <- do.call(rbind, intron_coverage_data)
  save(intron_coverage_data_all, file = str_glue("{outfolder_slope_calculations}intron_coverage_data_all.RData"))
  cat("Intron coverage data saved.\n")
  
  cat("Fitting linear models to intron coverage data...\n")
  intron_coverage_slope_lms <- lapply(chromosome_names, lm_fit_intron_coverage_per_chromosome, intron_coverage_data = intron_coverage_data, intron_info = selected_intron_info_df, nr_of_cores = nr_of_cores_for_parallel_computing, component_from_lm = "all")
  names(intron_coverage_slope_lms) <- chromosome_names
  intron_coverage_slope_lms_row_bound <- do.call(rbind, intron_coverage_slope_lms)
  cat("Linear models fitted.\n")
  
  intron_coverage_slope_lms_row_bound <- intron_coverage_slope_lms_row_bound %>%
    mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
    relocate(adjusted_p_value, .after = p_value)  
  cat("Adjusted p_values BH method.\n")
  
  
  write.table(intron_coverage_slope_lms_row_bound, file = str_glue("{outfolder_slope_calculations}intron_coverage_slope_lms.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
  cat("Slope calculation process completed and results saved in",outfolder_slope_calculations,".\n")
}