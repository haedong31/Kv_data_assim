library(tidyverse)

read_log <- function(log_path, num_file, num_iter) {
  con <- file(log_path, "r")
  log_list <- vector(mode = "list", length = num_file)
  
  for (i in 1:num_file) {
    log_lines <- vector(mode = "character", length = num_iter)
    for (j in 1:num_iter) {
      log_lines[j] <- readLines(con, n = 1)
    }
    log_list[[i]] <- log_lines
  }
  on.exit(close(con))
  return(log_list)
}

find_min_rmse <- function(log_list) {
  mins <- vector(mode = "numeric", length(log_list))
  file_names <- vector(mode = "character", length(log_list))
  for (i in seq_along(log_list)) {
    log_lines <- log_list[[i]]
    log_split <- str_split(log_lines, pattern = ": ")
    
    # file name
    file_name <- log_split %>% 
      map_chr(1) %>% 
      str_split(pattern = " ") %>% 
      map_chr(3) %>% 
      str_remove(".xlsx")
    file_names[i] <- file_name[1]
    
    # min RMSE
    rmse <- log_split %>% 
      map_chr(2) %>% 
      str_trim() %>% 
      as.numeric()
    mins[i] <- rmse %>%
      na.omit() %>% 
      min()
  }
  rmse_df <- tibble(File = file_names, RMSE = mins)
  return(rmse_df)
}

log_dir <- "./calibration/log/"

## create data frame for bar plot -----
num_files <- 34
num_iters <- 30
group <- "wt";
log1_name <- "exp45"
log2_name <- "exp46"
log3_name <- "exp47"

log1 <- read_log(str_c(log_dir,log1_name,"_",group,".txt"), num_files, num_iters)
log2 <- read_log(str_c(log_dir,log1_name,"_",group,".txt"), num_files, num_iters)
log3 <- read_log(str_c(log_dir,log3_name,"_",group,".txt"), num_files, num_iters)

rmse_df1 <- find_min_rmse(log1)
rmse_df2 <- find_min_rmse(log2)
rmse_df3 <- find_min_rmse(log3)

write_csv(rmse_df1, str_c(log_dir,log1_name,"_",group,".csv"))
write_csv(rmse_df2, str_c(log_dir,log2_name,"_",group,".csv"))
write_csv(rmse_df3, str_c(log_dir,log3_name,"_",group,".csv"))
