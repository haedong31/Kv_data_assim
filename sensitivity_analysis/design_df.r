library(tidyverse)
library(FrF2)

set.seed(228)

num_param_kto <- 14
num_param_kslow1 <- 10
num_param_kslow2 <- 2
num_param_kss <- 3

param_kto <- str_c("p", seq(1, num_param_kto))
param_kslow1 <- str_c("p", seq(1, num_param_kslow1))
param_kslow2 <- str_c("p", seq(1, num_param_kslow2))
param_kss <- str_c("p", seq(1, num_param_kss))

dgn1 <- FrF2(nruns = 1024, nfactors = num_param_kto, factor.names = param_kto)
dgn2 <- FrF2(nruns = 1024, nfactors = num_param_kslow1, factor.names = param_kslow1)
dgn3 <- FrF2(nruns = 2^num_param_kslow2, nfactors = num_param_kslow2, factor.names = param_kslow2)
dgn4 <- FrF2(nruns = 2^num_param_kss, nfactors = num_param_kss, factor.names = param_kss)

write_csv(dgn1, "./sensitivity_analysis/ikto_dgn.csv")
write_csv(dgn2, "./sensitivity_analysis/ikslow1_dgn.csv")
write_csv(dgn3, "./sensitivity_analysis/ikslow2_dgn.csv")
write_csv(dgn4, "./sensitivity_analysis/ikss_dgn.csv")
