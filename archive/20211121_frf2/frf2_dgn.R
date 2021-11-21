library(tidyverse)
library(FrF2)

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

write_csv(dgn1, "./sensitivity_analysis/ikto_dgn_matrix.csv")
write_csv(dgn2, "./sensitivity_analysis/ikslow1_dgn_matrix.csv")
write_csv(dgn3, "./sensitivity_analysis/ik1_dgn_matrix.csv")

# design 2: whole current model
num_param <- num_param_kto + num_param_kslow1 + num_param_k1
param <- str_c("p", seq(1, num_param))
dgn <- FrF2(resolution = target_resol, nfactors = num_param, factor.names = param)
write_csv(dgn, "./sensitivity_analysis/kcurrent_dgn_matrix.csv")
