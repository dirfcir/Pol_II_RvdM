args <- commandArgs(trailingOnly = TRUE)
cwd <- getwd()


# Function to display help message
display_help <- function() {
  cat("Usage: Rscript execute_all.R [options]\n",
      "Options:\n",
      "  sj_tab_folder                      Path to folder containing SJ tab files (default:",cwd,"/data/sj_tab_files/Drosophila)\n", 
      "  gtf_reference_file                 Path to the GTF input file. (default: ",cwd,"/data/gtf_reference_files/Drosophila_melanogaster.BDGP6.90.gtf)\n",
      "  bedgraph_folder                    Path to folder containing the bedGraph files (default: ",cwd,"/bedgraph_folder)\n",
      "  defined_introns_folder             Path to file with files created from defning introns in a gtf file (default: ",cwd,"/processed_gtf/)\n",
      "  annoted_intron_sj_folder           Path to folder containing SJ tab files (default: ",cwd, "/output/sj_annot_out_folder)\n",
      "  chromosome_names                   Path to file containing names of the chromosomes for the species (default: ",cwd,"/data/misc/chromosome_lengths_d_melanogaster)\n",
      "  nr_of_cores_for_paralel_computing  Number of cores to use for parallel computing (default: 3)\n",
      "  output_folder                      Path to the designated output folder (default: ",cwd,"/output_slope_calculations)",
      "  r_scripts_folder                   Path to folder containing R scripts (default: ",cwd,"/R_scripts)\n",
      "\nExample:\n",
      "  Rscript execute_all.R /path/to/sj_tab_folder/ /path/to/gtf_file.gtf /path/to/bedgraph_folder/ /path/to/defined_introns.bed /path/to/sj_annot_out_folder/ /path/to/chromosome_names /path/to/R_scripts/\n\n")
  quit(save = "no", status = 0, runLast = FALSE)
}

# Check if help is requested
if (length(args) == 1 && args[1] %in% c("help", "--help", "-h")) {
  display_help()
}

# Checking number of arguments
if (length(args) != 8) {
  cat("Usage: Rscript execute_all.R [options]\n",
      "Options:\n",
      "  sj_tab_folder                      Path to folder containing SJ tab files (default: ",cwd,"/data/sj_tab_files/Drosophila)\n", 
      "  gtf_reference_file                 Path to the GTF input file. (default: ",cwd,"/data/gtf_reference_files/Drosophila_melanogaster.BDGP6.90.gtf)\n",
      "  bedgraph_folder                    Path to folder containing the bedGraph files (default: ",cwd,"/bedgraph_folder)\n",
      "  defined_introns_folder             Path to file with defined introns (default: ",cwd,"/output_slope_calculations/processed_gtf/)\n",
      "  annoted_intron_sj_folder           Path to folder containing SJ tab files (default: ",cwd, "/output_slope_calculations/sj_annot_out_folder)\n",
      "  chromosome_names                   Path to file containing names of the chromosomes for the species (default: ",cwd,"/data/misc/chromosome_lengths_d_melanogaster)\n",
      "  nr_of_cores_for_paralel_computing  Number of cores to use for parallel computing (default: 3)\n",
      "  output_folder                      Path to the designated output folder (default: ",cwd,"/output_slope_calculations)\n",
      "  r_scripts_folder                   Path to folder containing R scripts (default: ",cwd,"/R_scripts)\n",
      "\nExample:\n",
      "  Rscript execute_all.R /path/to/sj_tab_folder/ /path/to/gtf_file.gtf /path/to/bedgraph_folder/ /path/to/defined_introns.bed /path/to/sj_annot_out_folder/ /path/to/chromosome_names /path/to/R_scripts/\n\n")
  
  cat("CAREFUL: Not all arguments are provided. Proceeding with default values:\n",
      "  sj_tab_folder: ",cwd,"/data/sj_tab_files/Drosophila/\n", 
      "  gtf_reference_file:",cwd,"/data/gtf_reference_files/Drosophila_melanogaster.BDGP6.90.gtf\n",
      "  bedgraph_folder: ",cwd,"/data/bedgraph_files/Drosophila/\n",
      "  defined_introns_folder: ",cwd,"/output_slope_calculations/processed_gtf/\n",
      "  annoted_intron_sj_folder: ",cwd, "/output_slope_calculations/sj_annot_out_folder/\n",
      "  chromosome_names: ",cwd,"/data/misc/chromosome_lengths_d_melanogaster\n",
      "  nr_of_cores_for_paralel_computing: 3\n",
      "  output_folder: ",cwd,"/output_slope_calculations/\n",
      "  r_scripts_folder: ",cwd,"/R_scripts/\n")
      
  
  args[1] <- paste0(cwd, "/data/sj_tab_files/Drosophila/")
  args[2] <- paste0(cwd, "/data/gtf_reference_files/Drosophila_melanogaster.BDGP6.90.gtf")
  args[3] <- paste0(cwd, "/data/bedgraph_files/Drosophila/") 
  args[4] <- paste0(cwd, "/output_slope_calculations/processed_gtf/")
  args[5] <- paste0(cwd, "/output_slope_calculations/sj_annot_out_folder/")
  args[6] <- paste0(cwd, "/data/misc/chromosome_lengths_d_melanogaster")
  args[7] <- 3
  args[8] <- paste0(cwd, "/output_folder_slope_calculations/")
  args[9] <- paste0(cwd, "/R_scripts/")
} else {
  message("All arguments provided. Proceeding with execution...")
  # Normal execution with provided arguments
}

source(paste0(args[9],"/required/load_packages.R"))
source(paste0(args[9],"/Recursive_annot_enhanced.R"))
source(paste0(args[9],"/slope_calc_master_script.R"))

message("Loading and installing required packages...")
suppressMessages(suppressWarnings(install_and_load(all_packages)))
message("All packages installed and loaded")

source(args[6]) #Define chromosome names

message("Defining intron regions from a gtf file...")
system2("./execute_scripts/intron_selection_code.sh", args = c(args[2], args[4]))
message("Defined intron regions from a gtf file")

message("Annotating splice junctions...")
annotate_splice_junctions(sj_tab_folder          = args[1],
                          defined_introns_folder = args[4],
                          sj_annot_out_folder    = args[5])

calculate_slopes(bedgraph_folder                    = args[3], 
                 annoted_intron_sj_folder           = args[5],
                 outfolder_slope_calculations       = args[3],
                 chromosome_names                   = chromosome_names,
                 nr_of_cores_for_parallel_computin  = as.numeric(args[7]))


