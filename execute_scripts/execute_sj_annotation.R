args <- commandArgs(trailingOnly = TRUE)
cwd <- cwd

# Function to display help message
display_help <- function() {
  cat("Usage: Rscript Recursive_annot_enhanced.R [options]\n\n",
      "Arguments:\n",
      "  sj_tab_folder                  Path to folder containing SJ tab files (default: attempt to find folder from cwd /sj_tab_files)\n", 
      "  defined_introns_folder         Path to file with defined introns (default: attempt to find file from cwd /processed_gtf/)\n",
      "  sj_annot_out_folder            Path to output folder for annotations of splice junctions from cwd (default: attempt to create folder /sj_annot_out_folder/)\n",
      "  r_scripts_folder               Path to folder containing R scripts (default:attempt to find folder from cwd R_scripts)\n",
      
      "\nDescription:\n",
      "  This script processes SJ tab files created by STAR and a defined introns file to annotate splice junctions.\n",
      "  It performs several steps including extracting, processing, and annotating the data, and generates output files with the results.\n\n",
      "Steps:\n",
      "  1. Extract splice junction files.\n",
      "  2. Process splice junction files and filter based on quality control criteria.\n",
      "  3. group the processed data.\n",
      "  4. Run bedtools intersect to compare with the defined introns file.\n",
      "  5. Annotate and save the final output.\n\n",
      "Output:\n",
      "  The script generates the following files in the specified output folder:\n",
      "  - [PATH_TO_OUTPUT_FOLDER]/all_samples_sj.tsv\n",
      "  - [PATH_TO_OUTPUT_FOLDER]/intersect_sj_measured_vs_reference.tsv\n",
      "  - [PATH_TO_OUTPUT_FOLDER]/intron_sj_annotated_with_reference.tsv\n",
      "  - [PATH_TO_OUTPUT_FOLDER]/intron_sj_annotated.tsv\n\n",
      "Examples:\n",
      "  Rscript Recursive_annot_enhanced.R /path/to/sj_tab_folder/ /path/to/defined_introns.bed /path/to/output_folder/\n\n")
  quit(save = "no", status = 0, runLast = FALSE)
}

# Check if help is requested
if (length(args) == 1 && args[1] %in% c("help", "--help", "-h")) {
  display_help()
}

cat("Usage: Rscript Recursive_annot_enhanced.R [options]\n",
    "Options:\n",
    "  sj_tab_folder                  Path to folder containing SJ tab files (default: cwd/sj_tab_files\n", 
    "  defined_introns_folder         Path to file with defined introns (default: cwd/processed_gtf/\n",
    "  sj_annot_out_folder            Path to output folder for annotations of splice junctions (default: cwd/sj_annot_out_folder/\n",
    "  r_scripts_folder               Path to folder containing R scripts (default:  cwd/R_scripts)\n",
    "\nExample:\n",
    "  Rscript Recursive_annot_enhanced.R /path/to/sj_tab_folder/ /path/to/processed_gtf_folder/ /path/to/output /path/to/R_scripts/\n\n")

# Checking number of arguments
if (length(args) != 4) {
  cat("CAREFUL: NOT ALL ARGUMENTS ARE PROVIDED.\n Proceeding with default values:\n",
      "sj_tab_folder: ",cwd,"/sj_tab_files\n", 
      "defined_introns_file: ",cwd,"/processed_gtf/\n",
      "sj_annot_out_folder: ",cwd,"sj_annot_out_folder/\n",
      "r_scripts_folder: ",cwd,"R_scripts\n")

  args[1] <- paste0(cwd,"/sj_tab_files/Drosophila/")
  args[2] <- paste0(cwd,"/processed_gtf/")
  args[3] <- paste0(cwd,"sj_annot_out_folder/")
  args[4] <- paste0(cwd,"/R_scripts/")
} else {
  message("All arguments provided. Proceeding with execution...")
}

source(paste0(args[4],"/required/load_packages.R"))
source(paste0(args[4],"/Recursive_annot_enhanced.R"))
install_and_load(all_packages)
message("All packages installed and loaded")

annotate_splice_junctions(sj_tab_folder        = args[1],
                          defined_introns_file = args[2],
                          sj_annot_out_folder  = args[3])  # Set default arguments and read them

#annotate_splice_junctions(sj_tab_folder = "/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/sj_tab_files/", defined_introns_file = str_glue("/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/Intron_selection/Drosophila_melanogaster.BDGP6.90.intron.bed"))


