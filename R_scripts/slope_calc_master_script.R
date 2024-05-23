calculate_slopes <- function(out_bg_folder        = "", 
                             defined_introns_file = "",
                             sj_annot_folder      = "",
                             outfolder_slope_calc = "",
                             chrsize              = ""){}
out_bg_folder <- "/cellfile/datapublic/cdebes/cdebes/workspace/scripts/dmelanogaster_mut/"
sj_annot_folder      <- "/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/sj_annot_out/"
outfolder_slope_calc <- "/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/outfolder_slope_calc"
source("/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/misc/chromosome_lengths_d_melanogaster")
chrsize <- chrsize_drosophila
nr_of_cluster_for_paralel_computing <- 6
cl <- makeCluster(6, outfile=str_glue("{outfolder_slope_calc}/cluster.log"))
clusterExport(cl, ls())
clusterEvalQ(cl, install_and_load(all_packages))
rm(aa)
bedgraph_fwd_file_pattern <- str_glue("{out_bg_folder}*Unique.str2.out.bg")
bedgraph_fwd_filepaths <- Sys.glob(coverage_fwd_file_pattern)
bedgraph_rev_file_pattern <- str_glue("{out_bg_folder}*Unique.str1.out.bg")
bedgraph_rev_filepaths <- Sys.glob(coverage_fwd_file_pattern)

norec<-read.table(str_glue("{sj_annot_folder}no_recursiv_intron_sj.bed"))
rec  <-read.table(str_glue("{sj_annot_folder}recursiv_intron_sj.bed"))
extracted_sj_annotation_bed_files  <- rbind.data.frame(rec,norec)
extracted_sj_annotation_bed_files  <- extracted_sj_annotation_bed_files[,c(1:4,6)]
colnames(extracted_sj_annotation_bed_files) <- c('chr','start','end','strand', 'id')

bed_graphs_data_fwd <- paralel_import_bed_graphs(bedgraph_filepaths = bedgraph_fwd_filepaths, cluster = cl)
bed_graphs_data_rev <- paralel_import_bed_graphs(bedgraph_filepaths = coverage_rev_filepaths, cluster = cl)

gen_data <- lapply(names(chrsize), get_cov_strand_spe, bed_graphs_data_fwd, bed_graphs_data_rev, gff_b, chrsize, "all", cl)
gen_data<-do.call(rbind,gen_data)
save(gen_data,file=paste0(path,"gen_data.RData"))
  
  gen_auc <- lapply(names(chrsize),run.chr.y.intercept.strand.spe, cov_data_fwd, cov_data_rev, gff_b, chrsize, "all", cl)
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
