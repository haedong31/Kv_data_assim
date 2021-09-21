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
  close(con)
  return(log_list)
}

find_min_rmse <- function(log_list) {
  mins <- vector(mode = "numeric", length(log_list))
  for (i in seq_along(log_list)) {
    log_lines <- log_list[[i]]
    rmse <- str_split(log_lines, pattern = ": ") %>% 
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

## WT -----
num_files <- 34
log1_name <- "exp16"
log2_name <- "exp17"
log3_name <- "exp18"
# wt_file_names <- c('15o26002','15o26005','15o26008','15o26014','15o26017',
#                    '15o26020','15o26023','15o26031','15o29002','15o29009',
#                    '15o29012','15o29015','15o29021','15o29024','15o29027',
#                    '15o20002','15o20005','15o20010','15o20015','15o21003',
#                    '15o21008','15o22002','15o22008','15n10002','15n10006',
#                    '15n10009','15n10012','15n23002','15n23005','15n23008',
#                    '15n23011','15n23014','15n23019','15n23022')
file_names <- seq(1, num_files) %>% as.character()

log1 <- read_log(str_c(log_dir, log1_name, "_wt.txt"), num_files, 30)
log2 <- read_log(str_c(log_dir, log2_name, "_wt.txt"), num_files, 30)
log3 <- read_log(str_c(log_dir, log3_name, "_wt.txt"), num_files, 1)

rmse_val1 <- find_min_rmse(log1)
rmse_val2 <- find_min_rmse(log2)
rmse_val3 <- find_min_rmse(log3)

bar_df <- tibble(
  name = file_names,
  mdl1 = rmse_val1,
  mdl2 = rmse_val2,
  mdl3 = rmse_val3)

# export for MATLAB
write_csv(bar_df, './calibration/bar_graph_exp16-17-18.csv')

bar_df <- bar_df %>% 
  pivot_longer(c('mdl1','mdl2'), names_to = 'mdl', values_to = 'rmse')
  
ggplot(data = bar_df, mapping = aes(x = name, y = rmse, fill = mdl)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_hline(yintercept = mean(rmse_val1), col = 'red') +
  geom_hline(yintercept = mean(rmse_val2), col = 'green') +
  scale_fill_discrete(name = "Model",
                      breaks = c("mdl1", "mdl2"),
                      labels = c("Model 1", "Model 2"))

## Mgat1KO -----
num_files <- 33
log1_name <- "exp16"
log2_name <- "exp17"
log3_name <- "exp18"
file_names <- seq(1, num_files) %>% as.character()

log1 <- read_log(str_c(log_dir, log1_name, "_ko.txt"), num_files, 30)
log2 <- read_log(str_c(log_dir, log2_name, "_ko.txt"), num_files, 30)
log3 <- read_log(str_c(log_dir, log3_name, "_ko.txt"), num_files, 1)

rmse_val1 <- find_min_rmse(log1)
rmse_val2 <- find_min_rmse(log2)
rmse_val3 <- find_min_rmse(log3)

bar_df <- tibble(
  name = file_names,
  mdl1 = rmse_val1,
  mdl2 = rmse_val2,
  mdl3 = rmse_val3)

# export for MATLAB
write_csv(bar_df, './calibration/bar_graph_exp16-17-18_ko.csv')
