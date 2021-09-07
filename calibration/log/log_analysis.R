library(tidyverse)

process_log <- function(log_path, num_file, num_iter) {
  con <- file(log_path, "r")
  log_list <- vector(mode = "list", length = num_file)
  
  for (i in 1:num_file) {
    log_lines <- vector(mode = "character", length = num_iter)
    for (j in 1:num_iter) {
      log_lines[j] <- readLines(con, n = 1)
    }
    log_list[[i]] <- log_lines
  }
  close(con)
  return(log_list)
}

min_rmse <- function(log_list) {
  mins <- vector(mode = "numeric", length(log_list))
  for (i in seq_along(log_list)) {
    log_lines <- log_list[[i]]
    rmse <- str_split(log_lines, pattern = "RMSE: ") %>% 
      map_chr(2) %>% 
      str_trim() %>% 
      as.numeric()
    mins[i] <- rmse %>%
      na.omit() %>% 
      min()
  }
  return(mins)
}

log_dir <- "./calibration/log/"

log_wt <- process_log(str_c(log_dir, "exp9_wt.txt"), 34, 30)
log_ko <- process_log(str_c(log_dir, "exp9_ko.txt"), 33, 30)

rmse_val_wt <- min_rmse(log_wt)
rmse_val_ko <- min_rmse(log_ko)

