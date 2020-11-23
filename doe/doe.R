library(tidyverse)
library(unrepx)
library(FrF2)


## experimental design for Ktrace2.m -----
# Ito
fnames <- seq(1, 12) %>% as.character()
fnames <- str_c('p', fnames)
dgn1_to <- FrF2(nfactors = 12, nruns = 1024, factor.names = fnames, seed = 1)
dgn1_to_info <- design.info(dgn1_to)
# run_order <- run.order(dgn1_to)
write_csv(dgn1_to, './doe/dgn-model2-Ito.csv')

# IKslow
fnames <- seq(1, 11) %>% as.character()
fnames <- str_c('p', fnames)
dgn1_slow <- FrF2(nfactors = 11, nruns = 1024, factor.names = fnames, seed = 2)
dgn1_slow_info <- design.info(dgn1_slow)
write_csv(dgn1_slow, './doe/dgn-model2-IKslow.csv')


## experimental design for Ktrace3.m -----
# Ito
fnames <- seq(1, 20) %>% as.character()
fnames <- str_c('p', fnames)
dgn2_to <- FrF2(nfactors = 20, nruns =  1024, factor.names = fnames, seed = 3)
dgn2_to_info <- design.info(dgn2_to)
write_csv(dgn2_to, './doe/dgn-model3-Ito.csv')

# IKslow
fnames <- seq(1, 12) %>% as.character()
fnames <- str_c('p', fnames)
dgn2_slow <- FrF2(nfactors = 12, nruns = 1024, factor.names = fnames, seed = 4)
write_csv(dgn2_slow, './doe/dgn-model3-IKslow.csv')


## model 2 Ito analysis -----
fnames <- seq(1, 12) %>% as.character()
fnames <- str_c('p', fnames)
dgn1_to <- FrF2(nfactors = 12, nruns = 1024, factor.names = fnames, seed = 1)
res_to <- read_csv('./doe/res-model2-Ito.csv')

# main effects - peak
dgn_res_to <- add.response(dgn1_to, res_to$peak)
anv_peak <- lm(dgn_res_to)
summary(anv_peak)

half_eff_peak <- anv_peak$coefficients[2:13]
eff_peak <- 2*half_eff_peak
names(eff_peak) <- fnames

hnplot(eff_peak, method = 'Lenth', half = TRUE, ID = TRUE, main = 'Peak Ito')
MEPlot(dgn_res_to, select = seq(1, 12), main = 'Effct plot of peak Ito')

# main effects - tau
dgn_res_to <- add.response(dgn1_to, res_to$time_const)
anv_tc <- lm(dgn_res_to)
summary(anv_tc)

half_eff_tc <- anv_tc$coefficients[2:13]
eff_tc <- 2*half_eff_tc
names(eff_tc) <- fnames

hnplot(eff_tc, method = 'Lenth', half = TRUE, ID = TRUE, main = 'Tau Ito')
MEPlot(dgn_res_to, select = seq(1, 12), main = 'Effct plot of Tau Ito')

# main effects - SSA
dgn_res_to <- add.response(dgn1_to, res_to$ssa)
anv_ssa <- lm(dgn_res_to)
summary(anv_ssa)

half_eff_ssa <- anv_ssa$coefficients[2:13]
eff_ssa <- 2*half_eff_ssa
names(eff_ssa) <- fnames

hnplot(eff_ssa, method = 'Lenth', half = TRUE, ID = TRUE, main = 'SSA Ito')
MEPlot(dgn_res_to, select = seq(1, 12), main = 'Effct plot of SSA Ito')

# main effects - SSI
dgn_res_to <- add.response(dgn1_to, res_to$ssi)
anv_ssi <- lm(dgn_res_to)
summary(anv_ssi)

half_eff_ssi <- anv_ssi$coefficients[2:13]
eff_ssi <- 2*half_eff_ssi
names(eff_ssi) <- fnames

hnplot(eff_ssi, method = 'Lenth', half = TRUE, ID = TRUE, main = 'SSI Ito')
MEPlot(dgn_res_to, select = seq(1, 12), main = 'Effct plot of SSI Ito')


## model 2 IKslow analysis -----
fnames <- seq(1, 11) %>% as.character()
fnames <- str_c('p', fnames)
dgn1_kslow <- FrF2(nfactors = 11, nruns = 1024, factor.names = fnames, seed = 2)
res_kslow <- read_csv('./doe/res-model2-IKslow.csv')

# main effects - peak
dgn_res_kslow <- add.response(dgn1_kslow, res_kslow$peak)
anv_peak <- lm(dgn_res_kslow)
summary(anv_peak)

half_eff_peak <- anv_peak$coefficients[2:12]
eff_peak <- 2*half_eff_peak
names(eff_peak) <- fnames

hnplot(eff_peak, method = 'Lenth', half = TRUE, ID = TRUE, main = 'Peak IKslow')
MEPlot(dgn_res_kslow, select = seq(1, 11), main = 'Effct plot of peak IKslow')

# main effects - tau
dgn_res_kslow <- add.response(dgn1_kslow, res_kslow$time_const)
anv_tc <- lm(dgn_res_kslow)
summary(anv_tc)

half_eff_tc <- anv_tc$coefficients[2:12]
eff_tc <- 2*half_eff_tc
names(eff_tc) <- fnames

hnplot(eff_tc, method = 'Lenth', half = TRUE, ID = TRUE, main = 'Tau IKslow')
MEPlot(dgn_res_kslow, select = seq(1, 11), main = 'Effct plot of Tau IKslow')

# main effects - SSA
dgn_res_kslow <- add.response(dgn1_kslow, res_kslow$ssa)
anv_ssa <- lm(dgn_res_kslow)
summary(anv_ssa)

half_eff_ssa <- anv_ssa$coefficients[2:12]
eff_ssa <- 2*half_eff_ssa
names(eff_ssa) <- fnames

hnplot(eff_ssa, method = 'Lenth', half = TRUE, ID = TRUE, main = 'SSA IKslow')
MEPlot(dgn_res_kslow, select = seq(1, 11), main = 'Effct plot of SSA IKslow')

# main effects - SSI
dgn_res_kslow <- add.response(dgn1_kslow, res_kslow$ssi)
anv_ssi <- lm(dgn_res_kslow)
summary(anv_ssi)

half_eff_ssi <- anv_ssi$coefficients[2:12]
eff_ssi <- 2*half_eff_ssi
names(eff_ssi) <- fnames

hnplot(eff_ssi, method = 'Lenth', half = TRUE, ID = TRUE, main = 'SSI IKslow')
MEPlot(dgn_res_kslow, select = seq(1, 11), main = 'Effct plot of SSI IKslow')
