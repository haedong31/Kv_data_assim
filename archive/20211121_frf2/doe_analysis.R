library(tidyverse)
library(unrepx)
library(FrF2)

# design matrix
set.seed(821)

# design 1: individually for each current
target_resol <- 6

num_param_kto <- 14
num_param_kslow1 <- 10
num_param_k1 <- 10

param_kto <- str_c("p", seq(1, num_param_kto))
param_kslow1 <- str_c("p", seq(1, num_param_kslow1))
param_k1 <- str_c("p", seq(1, num_param_k1))

dgn1 <- FrF2(resolution = target_resol, nfactors = num_param_kto, factor.names = param_kto)
dgn2 <- FrF2(resolution = target_resol, nfactors = num_param_kslow1, factor.names = param_kslow1)
dgn3 <- FrF2(resolution = target_resol, nfactors = num_param_k1, factor.names = param_k1)

# response of ik1 (average of currents at holding and clamp voltages)
res <- read_csv("./sensitivity_analysis/ik1_res.csv", col_names = FALSE)
res1 <- res$X1
res2 <- res$X2

dgn3_res1 <- add.response(dgn3, res1)
dgn3_res2 <- add.response(dgn3, res2)

# estimate main effects
dgn3_anv1 <- lm(dgn3_res1)
dgn3_anv2 <- lm(dgn3_res2)

half_mes1 <- dgn3_anv1$coefficients[1:num_param_k1+1]
half_mes2 <- dgn3_anv2$coefficients[1:num_param_k1+1]

mes1 <- 2*half_mes1
mes2 <- 2*half_mes2

MEPlot(dgn3_res1)
MEPlot(dgn3_res2)

hnplot(mes1, method = "Lenth", half = TRUE, ID = TRUE)
