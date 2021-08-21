library(tidyverse)
library(FrF2)

set.seed(821)

# design 1: individually for each current
target_resol <- 6

num_param1 <- 17
num_param2 <- 13
num_param3 <- 3
num_param4 <- 3
num_param5 <- 4
num_param6 <- 10

param1 <- str_c("p", seq(1, num_param1))
param2 <- str_c("p", seq(1, num_param2))
param3 <- str_c("p", seq(1, num_param3))
param4 <- str_c("p", seq(1, num_param4))
param5 <- str_c("p", seq(1, num_param5))
param6 <- str_c("p", seq(1, num_param6))

dgn1 <- FrF2(resolution = target_resol, nfactors = num_param1, factor.names = param1)
dgn2 <- FrF2(resolution = target_resol, nfactors = num_param2, factor.names = param2)
dgn3 <- FrF2(resolution = target_resol, nfactors = num_param3, factor.names = param3)
dgn4 <- FrF2(resolution = target_resol, nfactors = num_param4, factor.names = param4)
dgn5 <- FrF2(resolution = target_resol, nfactors = num_param5, factor.names = param5)
dgn6 <- FrF2(resolution = target_resol, nfactors = num_param6, factor.names = param6)

write_csv(dgn1, "./sensitivity_analysis/ikto_dgn_matrix.csv")
write_csv(dgn2, "./sensitivity_analysis/ikslow1_dgn_matrix.csv")
write_csv(dgn3, "./sensitivity_analysis/ikslow2_dgn_matrix.csv")
write_csv(dgn4, "./sensitivity_analysis/ikur_dgn_matrix.csv")
write_csv(dgn5, "./sensitivity_analysis/ikss_dgn_matrix.csv")
write_csv(dgn6, "./sensitivity_analysis/iks1_dgn_matrix.csv")

# design 2: whole current model
num_param = num_param1 + num_param2 + num_param3 + num_param4 + num_param5 + num_param6
param <- str_c("p", seq(1, num_param))
dgn <- FrF2(resolution = 5, nfactors = num_param, factor.names = param)
write_csv(dgn, "./sensitivity_analysis/kcurrent_dgn_matrix.csv")
