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

calculate_slopes <- function(out_bg_folder        = "", 
                             defined_introns_file = "",
                             sj_annot_folder      = "",
                             outfolder_slope_calc = "",
                             chrsize              = ""){}
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

intron_coverage_data <- lapply(names(chrsize), get_intron_coverages_per_chromosome, bed_graphs_data_fwd = bed_graphs_data_fwd, bed_graphs_data_rev = bed_graphs_data_rev, intron_info = selected_intron_info_df, chrsize = chrsize, cl = cl)
intron_coverage_data_all <-do.call(rbind,intron_coverage_data)
save(gen_data,file=str_glue("/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/intron_coverage_data_all.RData"))
# Load the saved data

check_intron_coverage_2L <- get_intron_coverages_per_chromosome(chromosome = "2L", bed_graphs_data_fwd = bed_graphs_data_fwd, bed_graphs_data_rev = bed_graphs_data_rev, intron_info = selected_intron_info_df, cl = cl)

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
