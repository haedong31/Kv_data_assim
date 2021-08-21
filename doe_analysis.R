library(tidyverse)
library(FrF2)

# design matrix
dgn <- read_csv("./")
# response
res <- read_csv("./sensitivity_analysis/ik1_res.csv", col_names = FALSE)
res1 <- res$X1
res2 <- res$X2

