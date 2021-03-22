library(tidyverse)
eps <- .Machine$double.eps

obj <- function(param, time_space, y, holdV, P1, Ek=-91.1) {
  yktof <- iktof(param[1:5], time_space, holdV, P1)
  yktos <- iktos(param[6:11], time_space, holdV, P1)
  ykur <- ikur(param[12:17], time_space, holdV, P1)
  ykss <- c(rep(0, length(time_space[[2]])), rep(param[18], length(time_space[[3]])))
  
  yksum <- yktof + yktos + ykur + ykss
  rmse <- sqrt(mean((y - yksum)^2))

  return(rmse)
}

iktof <- function(param, time_space, holdV, P1, Ek=-91.1) {
  # Bondarenko IKtof (A67)
  # 5 parameters
  
  # constants
  act0 <- 0.265563e-2 # ato_f Gating variable for transient outward K+ current
  inact0 <- 0.999977 # ito_f Gating variable for transient outward K+ current
  gmax <- param[5]  # 0.4067
  
  # time space information
  t <- time_space[[1]]
  tH <- time_space[[2]]
  tP1_adj <- time_space[[3]]
  hold_idx <- length(tH)
  
  end <- length(t)
  current_trc <- vector(mode="numeric", end)
  
  # at holding potential
  cv_hold <- cv_ktof(param[1:12], holdV)
  act_hold <- cv_hold[1] - (cv_hold[1] - act0)*exp(-(tH/cv_hold[2]))
  inact_hold <- cv_hold[3] - (cv_hold[3] - inact0)*exp(-(tH/cv_hold[4]))
  current_trc[1:hold_idx] <- gmax*(act_hold^3)*(inact_hold)*(holdV - Ek)
  
  # at P1
  cv_P1 <- cv_ktof(param[1:12], P1)
  act_P1 <- cv_P1[1] - (cv_P1[1] - act0)*exp(-(tP1_adj/cv_P1[2]))
  inact_P1 <- cv_P1[3] - (cv_P1[3] - inact0)*exp(-(tP1_adj/cv_P1[4]))
  current_trc[(hold_idx + 1):end] <- gmax*(act_P1^3)*(inact_P1)*(P1 - Ek)
  
  return(current_trc)
}

cv_ktof <- function(p, V) {
  # characteristic variables of IKtof
  # cv[1]: steady-state activation
  # cv[2]: time constant of activation
  # cv[3]: steady-state inactivation
  # cv[4]: time constant of inactivation
  
  # p0 <- c(30.0, 13.5, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335)
  cv <- vector(mode="numeric", 4) 
  alphaA <- 0.18064*exp(0.03577*(V+p[1]))
  betaA <- 0.3956*exp(-0.06237*(V+p[1]))
  alphaI <- (0.000152*exp(-(V+p[2])/p[4])) / (0.067083*exp(-(V+p[3])/p[4]) + 1)
  betaI <- (0.00095*exp((V+p[3])/p[4])) / (0.051335*exp((V+p[3])/p[4]) + 1)
  
  cv[1] <- alphaA/(alphaA+betaA)
  cv[2] <- 1/(alphaA+betaA)
  cv[3] <- alphaI/(alphaI+betaI)
  cv[4] <- 1/(alphaI+betaI)
  
  return(cv)
}

iktos <- function(param, time_space, holdV, P1, Ek=-91.1) {
    # Bondarenko IKtos (A75)
    # 6 parameters

    act0 <- 0.417069e-3 # ato_s Gating variable for transient outward K+ current
    inact0 <- 0.998543 # ito_s Gating variable for transient outward K+ current
    gmax <- param[6] # 0.0629 in spetal model 0 in apical model

    # time space information
    t <- time_space[[1]]
    tH <- time_space[[2]]
    tP1_adj <- time_space[[3]]
    hold_idx <- length(tH)

    end <- length(t)
    current_trc <- vector(mode="numeric", end)

    # at holding potential
    cv_hold <- cv_ktos(param[1:5], holdV)
    act_hold <- cv_hold[1] - (cv_hold[1] - act0)*exp(-(tH/cv_hold[2]))
    inact_hold <- cv_hold[3] - (cv_hold[3] - inact0)*exp(-(tH/cv_hold[4]))
    current_trc[1:hold_idx] <- gmax*(act_hold)*(inact_hold)*(holdV - Ek)

    # at P1
    cv_P1 <- cv_ktos(param[1:5], P1)
    act_P1 <- cv_P1[1] - (cv_P1[1] - act0)*exp(-(tP1_adj/cv_P1[2]))
    inact_P1 <- cv_P1[3] - (cv_P1[3] - inact0)*exp(-(tP1_adj/cv_P1[4]))
    current_trc[(hold_idx + 1):end] <- gmax*(act_P1)*(inact_P1)*(P1 - Ek)

    return(current_trc)  
}

cv_ktos <- function(p, V) {
    # p0 <- c(22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 270, 1050, 45.2, 5.7)
    cv <- vector(mode="numeric", 4)

    # note similarity of IKtos and IKur
    # same activation gate characteristics
    cv[1] <- 1 / (1 + exp(-(V+p[1])/p[2]))
    cv[2] <- 0.493*exp(-0.0629*V) + 2.058

    # for inactivation gate, cv[3] is same but not for cv[4] 
    cv[3] <- 1 / (1 + exp((V+p[3])/p[4]))
    cv[4] <- p[5] + 1050/(1 + exp((V + p[3])/p[4]))

    return(cv)
}

ikur <- function(param, time_space, holdV, P1, Ek=-91.1) {
    # Bondarenko IKur (A87)
    # 6 parameters

    # constants
    act0 <- 0.417069e-3  # aur Gating variable for 
    inact0 <- 0.998543  # iur Gating variable for ultrarapidly 
    gmax <- param[6] # 0.16

    # time space information
    t <- time_space[[1]]
    tH <- time_space[[2]]
    tP1_adj <- time_space[[3]]
    hold_idx <- length(tH)

    end <- length(t)
    current_trc <- vector(mode="numeric", end)

    # at holding potential
    cv_hold <- cv_kur(param[1:5], holdV)
    act_hold <- cv_hold[1] - (cv_hold[1] - act0)*exp(-(tH/cv_hold[2]))
    inact_hold <- cv_hold[3] - (cv_hold[3] - inact0)*exp(-(tH/cv_hold[4]))
    current_trc[1:hold_idx] <- gmax*(act_hold)*(inact_hold)*(holdV - Ek)

    # at P1
    cv_P1 <- cv_kur(param[1:5], P1)
    act_P1 <- cv_P1[1] - (cv_P1[1] - act0)*exp(-(tP1_adj/cv_P1[2]))
    inact_P1 <- cv_P1[3] - (cv_P1[3] - inact0)*exp(-(tP1_adj/cv_P1[4]))
    current_trc[(hold_idx + 1):end] <- gmax*(act_P1)*(inact_P1)*(P1 - Ek)

    return(current_trc)
}

cv_kur <- function(p, V) {
    # characteristic variables of IKur
    # cv[1]: steady-state activation
    # cv[2]: time constant of activation
    # cv[3]: steady-state inactivation
    # cv[4]: time constant of inactivation
    
    # p0 <- c(22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 1200.0, 170.0, 45.2, 5.7)
    cv <- vector(mode="numeric", 4)
    cv[1] <- 1 / (1 + exp(-(V+p[1])/p[2]))
    cv[2] <- 0.493*exp(-0.0629*V) + 2.058
    cv[3] <- 1 / (1 + exp((V+p[3])/p[4]))
    cv[4] <- p[5] - (170)/(1 + exp((V+p[3])/p[4]))

    return(cv)
}

# L-BFGS-B optimization
df <- read_csv("./benchmarks/ko1.csv")

t <- df$time
th <- t[1:99]
tp1_adj <- t[100:length(t)] - t[100]
y <- df$current/254.3

time_space <- vector(mode="list", 3)
time_space[[1]] <- t
time_space[[2]] <- th
time_space[[3]] <- tp1_adj

p0 <- c(30.0, 13.5, 33.5, 7.0, 0.4067, 
  22.5, 7.7, 45.2, 5.7, 270, 0.0629,
  22.5, 7.7, 45.2, 5.7, 1200, 0.16, 0.05)
low_bd <- c(-50, -60, -80, 1, eps,
  -60, 1, -60, 1, 1, eps,
  -60, 1, -60, 1, 170, eps)
up_bd <- c(60, 60, 60, 15, 2,
  60, 20, 100, 20, 5000, 2,
  60, 20, 100, 20, 5000, 2)

sol <- optim(p0, obj, method="L-BFGS-B", lower=low_bd, upper=up_bd,
  time_space=time_space, y=y, holdV=-70, P1=50)

yktof <- iktof(sol$par[1:5], time_space, -70, 50)
yktos <- iktos(sol$par[6:11], time_space, -70, 50)
ykur <- ikur(sol$par[12:17], time_space, -70, 50)
ykss <- c(rep(0, length(th)), rep(sol$par[18], length(tp1_adj)))
yksum <- yktof + yktos + ykur + ykss

plot(t, y, type="l", col="black", xlab="time (ms)", ylab="current (pA/pF)", lwd=2)
lines(t, yksum, lty="dashed", col="red")
legend("topright", c("Real","Simulated"), col=c("black","red"), lty=1:2)
