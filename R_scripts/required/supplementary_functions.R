find_folder_or_file <- function(files_to_find, max_steps = 4, folder = F) {
  
  current_dir <- normalizePath(getwd(), winslash = "/")
  for (i in 1:(max_steps + 1)) {
    # Check if the target folder exists in the current directory or its subdirectories
    # List all files and directories recursively in the current directory
    all_entries <- list.files(path = current_dir, recursive = TRUE, full.names = TRUE)
    
    # Filter to include only those entries that match the pattern
    matches <- all_entries[grepl(files_to_find, all_entries)]
    
    if (identical(matches, character(0))){
      
      # Move up one directory level
      current_dir <- dirname(current_dir)
    } else {
      if (folder){
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