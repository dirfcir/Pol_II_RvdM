#########################################################
### FOR KEEPS FOR KEEPS FOR KEEPS KOR KEEPS FOR KEEPS ###
#########################################################
# Consolidated function to compute intercepts, slopes, and their p-values
compute_lm_component <- function(x1, component = "slope") {
  if (is.null(x1)) { return(NULL) }
  
  result <- t(sapply(1:length(x1), function(z) {
    fit <- lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]])))
    summary_fit <- summary(fit)
    y_intercept <- fit$coefficients[1]
    x_intercept <- fit$fitted.values[[length(x1[[z]])]]
    slope <- summary_fit$coefficients[2, 1]
    p_value <- summary_fit$coefficients[2, 4]  # p-value for the slope (same for the model)
    
    if (component == "y_intercept") {
      return(c(y_intercept, p_value))
    } else if (component == "x_intercept") {
      return(c(x_intercept, p_value))
    } else if (component == "slope") {
      return(c(slope, p_value))
    } else if (component == "all") {
      return(c(y_intercept, x_intercept, slope, p_value))
    }
  }))
  
  if (component == "all") {
    colnames(result) <- c("y_intercept", "x_intercept", "slope", "p_value")
  } else {
    colnames(result) <- c("value", "p_value")
  }
  
  return(as.data.frame(result))
}

get_coverage_per_chromosome <- function(x, chromosome, cov_data) {
  chromosome_for_calc_coverage <- cov_data[[x]][seqnames(cov_data[[x]]) == chromosome]
  chromosome_coverage_data <- coverage(chromosome_for_calc_coverage, weight = chromosome_for_calc_coverage$score)[[chromosome]]
  return(chromosome_coverage_data)
}

compute_slope <- function(x1, filt = TRUE) {            
  if (is.null(x1)) { return(NULL) }
  slope <- sapply(1:length(x1), function(z) { sum <- summary(lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]])))); return(sum[4][[1]][[2]]) })
  return(slope)                  
}

filter_by_chromosome_and_strand <- function(intron_info_file, chromosome, strand_sense = '') {
  
  if (strand_sense == 'all'){
    introns_per_strand_and_chromosome <- intron_info_file[(intron_info_file$chromosome==chromosome),]                
  }
  
  if (strand_sense == 'minus') {
    introns_per_strand_and_chromosome <- intron_info_file[(intron_info_file$chromosome_measured == chromosome) & (intron_info_file$strand_measured == "-"),]
  }
  if (strand_sense == 'plus') {
    introns_per_strand_and_chromosome <- intron_info_file[(intron_info_file$chromosome_measured == chromosome) & (intron_info_file$strand_measured == "+"),]
  }
  return(introns_per_strand_and_chromosome)
}

paralel_import_bed_graphs <- function(bedgraph_filepaths, cluster) {
  imported_bed_graphs <- parLapply(cluster, bedgraph_filepaths, function(bedgraph_filepath) {
    imported_bed_graph <- import.bedGraph(bedgraph_filepath)
    mcols(imported_bed_graph)$filename <- basename(bedgraph_filepath)

    #seqlengths(imported_bed_graph) <- chr_size[levels(seqnames(imported_bed_graph))]                       
    return(imported_bed_graph)
    
  }) 
  return(imported_bed_graphs)
}

get_intron_coverages_per_chromosome <- function(chromosome, bed_graphs_data_fwd, bed_graphs_data_rev, intron_info, cl) { 
  
  # Sub select introns for chromosome and strand
  introns_per_strand_and_chromosome_fwd <- filter_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "plus")
  introns_per_strand_and_chromosome_rev <- filter_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "minus")
  
  # Split by chromosome
  coverage_data_fwd <- lapply(1:length(bed_graphs_data_fwd), get_coverage_per_chromosome, chromosome, bed_graphs_data_fwd)        
  coverage_data_rev <- lapply(1:length(bed_graphs_data_rev), get_coverage_per_chromosome, chromosome, bed_graphs_data_rev)
  
  # Extract genomic region of interest
  size_of_annot_fwd <- nrow(introns_per_strand_and_chromosome_fwd)
  size_of_annot_rev <- nrow(introns_per_strand_and_chromosome_rev)
  
  intron_coverage_data_fwd <- lapply(1:length(coverage_data_fwd), function(x) {                                                       
    mcmapply(function(start_position, end_position) {window(coverage_data_fwd[[x]], start_position, end_position)}, 
             introns_per_strand_and_chromosome_fwd$start_position_measured, 
             introns_per_strand_and_chromosome_fwd$end_position_measured,
             mc.cores = 6)                                      
  })
  
  intron_coverage_data_rev <- lapply(1:length(coverage_data_rev), function(x) {                                                       
    mcmapply(function(start_position, end_position) {rev(window(coverage_data_rev[[x]], start_position, end_position))},
             introns_per_strand_and_chromosome_rev$start_position_measured, 
             introns_per_strand_and_chromosome_rev$end_position_measured,
             mc.cores = 6)       
  })
  
  if (length(intron_coverage_data_fwd[[1]]) == 0) intron_coverage_data_fwd <- NULL else intron_coverage_data_fwd <- lapply(1:size_of_annot_fwd, function(y) lapply(intron_coverage_data_fwd, function(x) x[[y]]))
  
  if (length(intron_coverage_data_rev[[1]]) == 0) intron_coverage_data_rev <- NULL else intron_coverage_data_rev <- lapply(1:size_of_annot_rev, function(y) lapply(intron_coverage_data_rev, function(x) x[[y]]))
  
  intron_coverage_data <- c(intron_coverage_data_fwd, intron_coverage_data_rev)
  return(intron_coverage_data)
}

# Calculate the specified component of intronic regions
lm_fit_intron_coverage_per_chromosome <- function(chromosome, intron_coverage_data,intron_info, nr_of_cores = 6, cl = cl, component_from_lm = "all") {
  # Sub select introns for chromosome and strand
  introns_per_strand_and_chromosome_fwd <- filter_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "plus")
  introns_per_strand_and_chromosome_rev <- filter_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "minus")
  intron_coverage_lm_data <- mclapply(intron_coverage_data, compute_lm_component, component = component_from_lm, mc.cores = nr_of_cores)
  
  cat(paste(format(object.size(intron_coverage_lm_data), units = "MB"), "\n"))
  
  intron_coverage_lm_data[intron_coverage_lm_data == "NULL"] <- NA
  
  intron_coverage_lm_data <- do.call(rbind, intron_coverage_lm_data)
  
  introns_per_strand_and_chromosome <- rbind(introns_per_strand_and_chromosome_fwd, introns_per_strand_and_chromosome_rev)
  introns_per_strand_and_chromosome_test <<- introns_per_strand_and_chromosome
  intron_coverage_lm_data_test <<- intron_coverage_lm_data
  intron_coverage_lm_data <- cbind(introns_per_strand_and_chromosome, intron_coverage_lm_data)
  
  return(intron_coverage_lm_data)
}
