find_folder <- function(folder_to_find, max_steps = 2) {
  current_dir <- normalizePath(getwd(), winslash = "/")
  
  for (i in 1:(max_steps + 1)) {
    # Check if the target folder exists in the current directory or its subdirectories
    matches <- list.files(paste0(current_dir,"/",folder_to_find, sep =""), recursive = TRUE, full.names = T)
    if (identical(matches, character(0))){
      # Move up one directory level
      parent_dir <- dirname(current_dir)
      stop(paste("Directory", folder_to_find, "not found within", max_steps, "steps"))
    } else {
      return(matches)
    }
  }
  stop(paste("Directory", folder_to_find, "not found within", max_steps, "steps"))
}
