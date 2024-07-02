args <- commandArgs(trailingOnly = TRUE)
cwd <- getwd()

# Function to display help message
display_help <- function() {
  cat("Usage: Rscript slope_calc_master_script.R [options]\n",
      "Options:\n",
      "  bedgraph_folder                    Path to folder containing the bedGraph files (default: attempt to find file from cwd ./bedgraph_folder/)\n",
      "  annoted_intron_sj_folder           Path to folder containing SJ tab files (default: attempt to find folder from cwd ./output/sj_annot_out_folder)\n",
      "  outfolder_slope_calculations       Path to output file for slope calculations (default:attempt to create file from cwd ./data/outfolder_slope_calc)\n",
      "  chromosome_names                   Path to file containing names of the chromosomes for the species (default: attempt to find file from cwd ./data/misc/chromosome_lengths_d_melanogaster)\n",
      "  nr_of_cores_for_paralel_computing  Number of cores to use for paralel computing (default: 3)",
      "  r_scripts_folder                   Path to folder containing R scripts (default:attempt to find folder from cwd R_scripts)\n",
      "\nExample:\n",
      "  Rscript execute_slope_calc /path/to/outfolder_slope_calc/ /path/to/bedgraph_folder/ nr_of_cores_for_paralel_computing /path/to/chromosome_names /path/to/sj_annot_out_folder/ /path/to/R_scripts/\n\n")
  quit(save = "no", status = 0, runLast = FALSE)
}

# Check if help is requested
if (length(args) == 1 && args[1] %in% c("help", "--help", "-h")) {
  display_help()
}

# Checking number of arguments
if (length(args) != 6) {
  cat("Usage: Rscript slope_calc_master_script.R [options]\n",
      "Options:\n",
      "  bedgraph_folder                    Path to folder containing the bedGraph files (default: ",cwd,"/bedgraph_folder/)\n",
      "  annoted_intron_sj_folder           Path to folder containing SJ tab files (default: attempt to find folder from cwd ",cwd,"/output/sj_annot_out_folder)\n",
      "  outfile_slope_calc                 Path to output file for slope calculations (default: ",cwd,"/data/outfolder_slope_calc\n",
      "  chromosome_names                   Path to file containing names of the chromosomes for the species (default: attempt to find file from cwd ",cwd,"/data/misc/chromosome_lengths_d_melanogaster)\n",
      "  nr_of_cores_for_paralel_computing  Number of cores to use for paralel computing (default: 3)",
      "  r_scripts_folder                   Path to folder containing R scripts (default: ",cwd,"/Pol_II_RvdM/R_scripts/\n",
      "\nExample:\n",
      "  Rscript execute_slope_calc /path/to/outfolder_slope_calc/ /path/to/bedgraph_folder/ nr_of_cores_for_paralel_computing /path/to/chromosome_names /path/to/sj_annot_out_folder/ /path/to/R_scripts/\n\n")
  
  cat("CAREFULL: Not all arguments provided. Proceeding with default values:\n",
      "  bedgraph_folder: /cellfile/datapublic/cdebes/cdebes/workspace/scripts/dmelanogaster_mut/", # should be this in the future ,cwd,"/data/bedgraph_folder/)\n",
      "  outfile_slope_calc: ",cwd,"/data/outfolder_slope_calc\n",
      "  annoted_intron_sj_folder:",cwd,"/data/output/sj_annot_out_folder)\n",
      "  r_scripts_folder: ",cwd,"/R_scripts/\n"
      )
        # Filter out files, keep only directories
  args[1] <- "/cellfile/datapublic/cdebes/cdebes/workspace/scripts/dmelanogaster_mut/" # should be this in the future paste0(cwd,"/data/bedgraph_folder/")
  args[2] <- paste0(cwd,"/data/sj_annot_out_folder/")
  args[3] <- paste0(cwd,"/data/output/sj_annot_out_folder")
  args[4] <- paste0(cwd,"/data/misc/chromosome_lengths_d_melanogaster")
  args[5] <- 3
  args[6] <- paste0(cwd,"/R_scripts/")
} else {
  message("All arguments provided. Proceeding with execution...")
  # Normal execution with provided arguments 
}
source(paste0(args[6],"/required/load_packages.R"))
source(paste0(args[6],"/slope_calc_master_script.R"))

cat("Loading and installing required packages...")
suppressMessages(suppressWarnings(install_and_load(all_packages)))

source(args[4])

calculate_slopes(bedgraph_folder                    = args[1], 
                 annoted_intron_sj_folder           = args[2],
                 outfolder_slope_calculations       = args[3],
                 chromosome_names                   = chromosome_names,
                 nr_of_cores_for_parallel_computin  = args[5])
                             

