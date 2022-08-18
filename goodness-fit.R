library(tidyverse)
library(readxl)
library(PAmeasures)

file_group <- "wt"
exp_num <- "exp51"

calib_dir <- file.path("calibration", str_c("calib_",exp_num))
data_dir <- file.path("calibration", "mgat1ko_data", str_c(file_group,"-preprocessed"))

matching_tbl <- read_excel(
  file.path("calibration", "mgat1ko_data", 
                          str_c("matching-table-", file_group, ".xlsx")))
file_names <- matching_tbl$trace_file_name_4half %>% na.omit()

# prepare for iteration
num_volts <- 9
num_files <- length(file_names)

col_names1 <- c("t", str_c("y",1:num_volts))
col_names2 <- str_c("yhat",1:num_volts)
col_names3 <- "file_name"

a1 <- 2
a2 <- 3
a1_idx <- vector("numeric", length = num_volts)
a2_idx <- vector("numeric", length = num_volts)
for (i in 1:num_volts) {
  col_names3 <- c(col_names3, str_c("R2_v",i), str_c("L2_v",i))
  a1_idx[i] <- a1 + (i-1)*2
  a2_idx[i] <- a2 + (i-1)*2
}
col_names3 <- c(col_names3, "mean_R2", "mean_L2")

for (i in seq_along(file_names)) {
  exp_data <- read_excel(
    file.path(data_dir,file_names[i]),
    col_names = FALSE)
  colnames(exp_data) <- col_names1
  
  yhat <- read_excel(
    file.path(calib_dir, str_c(file_group, "_yhat"), file_names[i]),
    col_names = FALSE)
  colnames(yhat) <- col_names2
  
  r <- vector("list", length = length(col_names3))
  names(r) <- col_names3
  
  r[[1]] <- file_names[i]
  for (j in 1:num_volts) {
    m <- pam.nlm(exp_data[[j+1]], yhat[[j]]) # list of two
    r[[a1_idx[j]]] <- as.numeric(m[[1]])
    r[[a2_idx[j]]] <- as.numeric(m[[2]])
  }
  r[[1+2*num_volts+1]] <- mean(unlist(r[a1_idx]))
  r[[1+2*num_volts+2]] <- mean(unlist(r[a2_idx]))
  
  if (i == 1) {
    fit_df <- tibble(!!!r)
  } else {
    fit_df <- add_row(fit_df,!!!r)
  }
}

write_excel_csv(fit_df, 
                file.path(calib_dir,str_c("r2l2_measures_",file_group,".csv")))
