library(tidyverse)
library(FrF2)


## for Ktrace2.m -----
# Ito
fnames <- seq(1, 12) %>% as.character()
fnames <- str_c('p', fnames)
dgn1_to <- FrF2(nfactors = 12, nruns = 1024, factor.names = fnames, seed = 1)
dgn1_to_info <- design.info(dgn1_to)
# run_order <- run.order(dgn1_to)
write_csv(dgn1_to, './dgn-model2-Ito.csv')

# IKslow
fnames <- seq(1, 11) %>% as.character()
fnames <- str_c('p', fnames)
dgn1_slow <- FrF2(nfactors = 11, nruns = 1024, factor.names = fnames, seed = 2)
dgn1_slow_info <- design.info(dgn1_slow)
write_csv(dgn1_slow, './dgn-model2-IKslow.csv')


## for Ktrace3.m -----
# Ito
fnames <- seq(1, 20) %>% as.character()
fnames <- str_c('p', fnames)
dgn2_to <- FrF2(nfactors = 20, nruns =  1024, factor.names = fnames, seed = 3)
dgn2_to_info <- design.info(dgn2_to)
write_csv(dgn2_to, './dgn-model3-Ito.csv')

# IKslow
fnames <- seq(1, 12) %>% as.character()
fnames <- str_c('p', fnames)
dgn2_slow <- FrF2(nfactors = 12, nruns = 1024, factor.names = fnames, seed = 4)
write_csv(dgn2_slow, './dgn-model3-IKslow.csv')
