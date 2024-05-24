select_introns_for_slopecalc <- function(intron_info_file){
  intron_info_file <- intron_info_file %>%
    filter(overlap_outside_reference != TRUE) %>%
    filter(supported_by_samples == 1)         %>%
    filter(strand_measured != 0)              %>%
    filter((start_match_reference == TRUE & end_match_reference == TRUE) 
           | 
             (start_match_reference == TRUE & end_match_reference == FALSE) 
           | 
             (start_match_reference == FALSE & end_match_reference == TRUE)
    )
  return(intron_info_file)
}

#########################################################
### FOR KEEPS FOR KEEPS FOR KEEPS KOR KEEPS FOR KEEPS ###
#########################################################
# Consolidated function to compute intercepts, slopes, and their p-values
compute_lm <- function(x1, component = "slope") {
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

filter_introns_by_chromosome_and_strand <- function(intron_info_file, chromosome, strand_sense = '') {
  
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
  imported_bed_graphs <- parLapply(cluster, bedgraph_filepaths, function(d) {
    imported_bed_graph <- import.bedGraph(d)
    #seqlengths(imported_bed_graph) <- chr_size[levels(seqnames(imported_bed_graph))]                       
    return(imported_bed_graph)
    
  }) 
  return(imported_bed_graphs)
}

get_intron_coverage_per_chromosome <- function(chromosome, bed_graphs_data_fwd, bed_graphs_data_rev, intron_info, nr_of_cores = 6) { 
  
  # Subselect introns for chromosome and strand
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
             mc.cores = nr_of_cores)                                      
  })
  
  intron_coverage_data_rev <- lapply(1:length(coverage_data_rev), function(x) {                                                       
    mcmapply(function(start_position, end_position) {rev(window(coverage_data_rev[[x]], start_position, end_position))},
             introns_per_strand_and_chromosome_rev$start_position_measured, 
             introns_per_strand_and_chromosome_rev$end_position_measured,
             mc.cores = nr_of_cores)       
  })
  
  if (length(intron_coverage_data_fwd[[1]]) == 0) intron_coverage_data_fwd <- NULL else intron_coverage_data_fwd <- lapply(1:size_of_annot_fwd, function(y) lapply(intron_coverage_data_fwd, function(x) x[[y]]))
  
  if (length(intron_coverage_data_rev[[1]]) == 0) intron_coverage_data_rev <- NULL else intron_coverage_data_rev <- lapply(1:size_of_annot_rev, function(y) lapply(intron_coverage_data_rev, function(x) x[[y]]))
  
  chromosome_intron_coverage_data <- c(intron_coverage_data_fwd, intron_coverage_data_rev)
  return(chromosome_intron_coverage_data)
}

# Calculate the specified component of intronic regions
lm_fit_intron_coverage_per_chromosome <- function(chromosome, intron_coverage_data,intron_info, nr_of_cores = 6, component_from_lm = "all") {
  # Sub select introns for chromosome and strand
  introns_per_strand_and_chromosome_fwd <- filter_introns_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "plus")
  introns_per_strand_and_chromosome_rev <- filter_introns_by_chromosome_and_strand(intron_info_file = intron_info, chromosome = chromosome, strand_sense = "minus")
  
  chromosome_intron_coverage_lm_info <- mclapply(intron_coverage_data, compute_lm, component = component_from_lm, mc.cores = nr_of_cores)
  
  cat(paste(format(object.size(chromosome_intron_coverage_lm_info), units = "MB"), "\n"))
  
  chromosome_intron_coverage_lm_info[chromosome_intron_coverage_lm_info == "NULL"] <- NA
  
  chromosome_intron_coverage_lm_info <- do.call(rbind, chromosome_intron_coverage_lm_info)
  
  introns_per_chromosome <- rbind(introns_per_strand_and_chromosome_fwd, introns_per_strand_and_chromosome_rev)
  
  chromosome_intron_coverage_lm_info <- cbind(introns_per_chromosome, chromosome_intron_coverage_lm_info)
  
  return(chromosome_intron_coverage_lm_info)
}


get_coverage_per_chromosome_and_strand <- function(out_bg_folder       = "", 
                                                   defined_introns_file = "",
                                                   sj_annot_folder      = "",
                                                   outfolder_slope_calc = "",
                                                   chromosome_names     = ""){}
out_bg_folder <- "/cellfile/datapublic/cdebes/cdebes/workspace/scripts/dmelanogaster_mut/"
sj_annot_folder      <- "/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/sj_annot_out_folder/"
outfolder_slope_calc <- "/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/outfolder_slope_calc"
source("/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/misc/chromosome_lengths_d_melanogaster")
chrsize <- chrsize_drosophila
nr_of_cluster_for_paralel_computing <- 6
cl <- makeCluster(6, outfile=str_glue("{outfolder_slope_calc}/cluster.log"))
clusterExport(cl, ls())
clusterEvalQ(cl, install_and_load(all_packages))

bedgraph_fwd_file_pattern <- str_glue("{out_bg_folder}*Unique.str2.out.bg")
bedgraph_fwd_filepaths <- Sys.glob(bedgraph_fwd_file_pattern)
bedgraph_rev_file_pattern <- str_glue("{out_bg_folder}*Unique.str1.out.bg")
bedgraph_rev_filepaths <- Sys.glob(bedgraph_rev_file_pattern)

intron_info_df <- read_tsv(str_glue("{sj_annot_folder}intron_sj_annotated.tsv"))
selected_intron_info_df <- select_introns_for_slopecalc(intron_info_df)

#extracted_sj_annotation_bed_files  <- extracted_sj_annotation_bed_files[,c(1:4,6)]
#colnames(extracted_sj_annotation_bed_files) <- c('chr','start','end','strand', 'id')

bed_graphs_data_fwd <- paralel_import_bed_graphs(bedgraph_filepaths = bedgraph_fwd_filepaths, cluster = cl)
bed_graphs_data_rev <- paralel_import_bed_graphs(bedgraph_filepaths = bedgraph_rev_filepaths, cluster = cl)

intron_coverage_data <- lapply(chromosome_names, get_intron_coverage_per_chromosome, bed_graphs_data_fwd = bed_graphs_data_fwd, bed_graphs_data_rev = bed_graphs_data_rev, intron_info = selected_intron_info_df, chrsize = chrsize, cl = cl)
intron_coverage_data_all <-do.call(rbind,intron_coverage_data)
save(gen_data,file=str_glue("/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/intron_coverage_data_all.RData"))
return(intron_coverage_data)
# Load the saved data
selected_intron_info_df <- select_introns_for_slopecalc(intron_info_df)

check_intron_coverage_2L <- get_intron_coverage_per_chromosome(chromosome = "2L", bed_graphs_data_fwd = bed_graphs_data_fwd, bed_graphs_data_rev = bed_graphs_data_rev, intron_info = selected_intron_info_df)

lm_fit <- lm_fit_intron_coverage_per_chromosome(chromosome = "2L", intron_coverage_data = check_intron_coverage_2L, intron_info = selected_intron_info_df, nr_of_cores = 6, component_from_lm = "all")



slope_lms <- lapply(names(chrsize), lm_fit_intron_coverage, intron_coverage_data = gen_data, intron_info = selected_intron_info_df, nr_of_cores = 6, component_from_lm = "all")
intercept_start<-do.call(rbind,gen_auc)
save(intercept_start,file=paste0(path,"intercept_start.RData"))

gen_auc <- lapply(names(chrsize),run.chr.x.intercept.strand.spe, cov_data_fwd, cov_data_rev, gff_b, chrsize, "all", cl)
intercept_end<-do.call(rbind,gen_auc)
save(intercept_end,file=paste0(path,"intercept_end.RData"))

gen_auc <- lapply(names(chrsize_drosophila),run_chr_slope_strand_spe, cov_data_fwd[1:2], cov_data_rev[1:2], gff_b, chrsize_drosophila, "all", cl)
slope<-do.call(rbind,gen_auc)
save(slope,file=paste0(path,"slope.RData"))

stopCluster(cl)

}
