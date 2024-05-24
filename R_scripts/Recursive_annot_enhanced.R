# Helper function for QC to filter splice junctions based on read size and number of unique reads.
filter_splicejunct_QC <- function(sj_df, mx_sj_overhang = 20, nr_of_uniqreads = 5) {
  sj_df %>%
    filter(V9 >= mx_sj_overhang, V7 >= nr_of_uniqreads)
}

annotate_splice_junctions <- function(sj_tab_folder        = cwd,
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
  print(splice_files)
  message("Processing splice junction files...")
  # Process all splice junction files
  splice_data <- lapply(splice_files, function(file) {
    fread(file, header = FALSE) %>%
      filter_splicejunct_QC()
  }) %>%
    bind_rows() %>%
    mutate(supported_by_samples = 1) %>%
    group_by(V1, V2, V3, V4, V5, V6) %>%
    summarise(
      supported_by_samples = sum(supported_by_samples),
      max_unique_mapping_count = max(V7),
      min_multi_mapping_count = min(V8),
      max_overhang = max(V9),
    ) %>%
    mutate(supported_by_samples = supported_by_samples/length(splice_files))

  # Write the filtered splice junctions to an output file
  filterd_sj_outfile_name <- str_glue("{sj_annot_out_folder}/all_samples_sj.tsv")
  write_tsv(splice_data, filterd_sj_outfile_name, col_names = FALSE)
  
  message("Running bedtools intersect...")
  # Intersect the splice junction file with the intron file using system command
  intersect_command <- str_glue("bedtools intersect -wa -wb -a {filterd_sj_outfile_name} -b {defined_introns_file} > {sj_annot_out_folder}/intersect_sj_measured_vs_reference.tsv")
  system(intersect_command)
  
  message("Reading intersection output...")
  # Read the intersection output and perform further analysis
  junction_data <- read_tsv(str_glue("{sj_annot_out_folder}/intersect_sj_measured_vs_reference.tsv"), col_names = FALSE)
  colnames(junction_data) <- c(
    "chromosome_measured",
    "start_position_measured",
    "end_position_measured",
    "strand_measured",
    "intron_motif",
    "annotated",
    "supported_by_samples",
    "max_unique_mapping_count",
    "min_multi_mapping_count",
    "max_overhang",
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
  
  #Add Flags
  
  # Create a function to check certain condition
  check_intron_within_intron <- function(start, end, df) {
    any(df$start_position_measured > start & df$end_position_measured < end)
  }
  
  junction_data <- junction_data %>%
    mutate(inside_reference_intron = start_position_measured - start_position_gtf_refrnce >= 2 
                                     & 
                                      end_position_measured - end_position_gtf_refrnce <= -2
                                     &
                                      strand_measured == strand_gtf_refrnce
           )
  
  junction_data <- junction_data %>%
    mutate(end_match_reference = strand_measured == strand_gtf_refrnce
                                 &
                                  abs(end_position_measured - end_position_gtf_refrnce) <= 1
    )
  
  junction_data <- junction_data %>%
    mutate(start_match_reference = strand_measured == strand_gtf_refrnce
                                   &
                                    abs(start_position_measured - start_position_gtf_refrnce) <= 1
           )
  
  junction_data <- junction_data %>%
    mutate(overlap_outside_reference = strand_measured == strand_gtf_refrnce
                                       &
                                        (start_position_measured - start_position_gtf_refrnce <= -2 | end_position_measured - end_position_gtf_refrnce >= 2)
                                )
  
  junction_data <- junction_data %>%
   rowwise() %>%
    mutate(measured_intron_within_this_intron = check_intron_within_intron(start_position_measured, end_position_measured, junction_data)) %>%
    ungroup()

  # junction_data <- junction_data %>%
    #mutate(full_match_reference = strand_measured == strand_gtf_refrnce
                                #  &
                                 #  abs(start_position_measured - start_position_gtf_refrnce) <= 1
                                 # &
                                   #abs(end_position_measured - end_position_gtf_refrnce) <= 1
    #)
  message("Processing and saving final output...")
  # Process and save final output
  # Output the first few lines of the processed data for verification
  message("Output of the first few lines of the processed data for verification")
  print(head(junction_data))
  # Final output
  write_tsv(
    junction_data, 
    str_glue("{sj_annot_out_folder}/intron_sj_annotated_with_reference.tsv"), col_names = TRUE
  )
  write_tsv(
    junction_data[,c("chromosome_measured",	"start_position_measured",	"end_position_measured",	"strand_measured",	"intron_motif",	"supported_by_samples",	"max_unique_mapping_count",	"min_multi_mapping_count",	"max_overhang", "inside_reference_intron",	"start_match_reference",	"end_match_reference",	"overlap_outside_reference", "measured_intron_within_this_intron")], 
    str_glue("{sj_annot_out_folder}/intron_sj_annotated.tsv"), col_names = TRUE
  )

}
 