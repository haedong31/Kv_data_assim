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
num_files2 <- 33
group_name <- "wt"
group_name2 <- "ko"
log1_name <- "exp27"
log2_name <- "exp27"
# log3_name <- "exp26"
# log4_name <- "exp23"

log1 <- read_log(str_c(log_dir, log1_name, "_", group_name, ".txt"), num_files, 30)
log2 <- read_log(str_c(log_dir, log2_name, "_", group_name2, ".txt"), num_files2, 30)
# log3 <- read_log(str_c(log_dir, log3_name, "_", group_name, ".txt"), num_files, 1)
# log4 <- read_log(str_c(log_dir, log4_name, "_", group_name, ".txt"), num_files, 1)

rmse_df1 <- find_min_rmse(log1)
rmse_df2 <- find_min_rmse(log2)
# rmse_val3 <- find_min_rmse(log3)
# rmse_val4 <- find_min_rmse(log4)

write_csv(rmse_df1, str_c(log_dir, log1_name, "_", group_name, ".csv"))
write_csv(rmse_df2, str_c(log_dir, log2_name, "_", group_name2, ".csv"))

# file_names <- seq(1, num_files) %>% as.character()
# extra_idx <- c(18, 19, 29, 31, 32)
# rmse_val1[extra_idx]

## create data frame for multiple bar plot -----
bar_df <- tibble(
  name = file_names,
  mtd1 = rmse_val1,
  mtd2 = rmse_val2,
  mtd3 = rmse_val3,
  mtd4 = rmse_val4)

# export for MATLAB
write_csv(bar_df, './calibration/bar_graph_norm.csv')

bar_df <- bar_df %>% 
  pivot_longer(c('mtd1','mtd2'), names_to = 'mtd', values_to = 'rmse')
  
ggplot(data = bar_df, mapping = aes(x = name, y = rmse, fill = mtd)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_hline(yintercept = mean(rmse_val1), col = 'red') +
  geom_hline(yintercept = mean(rmse_val2), col = 'green') +
  scale_fill_discrete(name = "Model",
                      breaks = c("mtd1", "mtd2"),
                      labels = c("Method 1", "Method 2"))
