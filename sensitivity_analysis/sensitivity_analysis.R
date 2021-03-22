library(lhs)
library(laGP)
eps <- .Machine$double.eps

iktof <- function(param, holdV, P1, time_space, Ek=-91.1) {
  # Bondarenko IKtof (A67)
  # 13 parameters
  
  # constants
  act0 <- 0.265563e-2 # ato_f Gating variable for transient outward K+ current
  inact0 <- 0.999977 # ito_f Gating variable for transient outward K+ current
  gmax <- param[13]  # 0.4067
  
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
  alphaA <- p[5]*exp(p[6]*(V+p[1]))
  betaA <- p[7]*exp(-p[8]*(V+p[1]))
  alphaI <- (p[9]*exp(-(V+p[2])/p[4])) / (p[10]*exp(-(V+p[3])/p[4]) + 1)
  betaI <- (p[11]*exp((V+p[3])/p[4])) / (p[12]*exp((V+p[3])/p[4]) + 1)
  
  cv[1] <- alphaA/(alphaA+betaA)
  cv[2] <- 1/(alphaA+betaA)
  cv[3] <- alphaI/(alphaI+betaI)
  cv[4] <- 1/(alphaI+betaI)
  
  return(cv)
}

output1 <- function(time_space, current_trc) {
  # outputs for IKtof, IKtos, and IKur
  # output 1: peak time
  # output 2: peak current
  # output 3: tau, peak reduced by 1 - exp(-1) (~63%)
  # output 4: current at 1/3 of interval [peak, tau]
  # output 5: current at 2/3 of interval [peak, tau]
  # output 6: current at last
  # output 7: first time hit output 6 value
  out <- vector(mode="numeric", 7)
  
  t <- time_space[[1]]
  th <- time_space[[2]]
  hold_idx <- length(th)
  holdt <- t[hold_idx]
  
  # valid current trace shape?
  peak_time_idx <- which.max(current_trc)
  peak_time <- t[peak_time_idx]
  peak <- current_trc[peak_time_idx]
  
  check_pt1 <- any(current_trc < 0) # negative current
  check_pt2 <- peak_time < holdt # can't generate current properly
  check_pt3 <- var(current_trc[1:holdt]) > 1e-5 # not stable at holding potential
  
  if(check_pt1 || check_pt2 || check_pt3) {return(NA)}	
  
  # output 1
  out[1] <- peak_time - holdt
  
  # output 2
  out[2] <- peak
  
  # output 3
  tau_current <- peak*exp(-1)
  tau_idx <- which.min(abs(current_trc - tau_current))
  tau <- t[tau_idx]
  out[3] <- tau - holdt
  
  # output 4
  out45_jump_size <- floor((peak_time_idx + tau_idx)/3)
  out[4] <- current_trc[peak_time_idx + out45_jump_size]
  
  # output 5
  out[5] <- current_trc[peak_time_idx + out45_jump_size*2]
  
  # output 6
  last_idx <- length(t)
  last_current <- current_trc[last_idx]
  out[6] <- last_current
  
  # output 7
  o7_idx <- which.min(abs(current_trc[tau_idx:last_idx] - last_current))
  o7_idx <- o7_idx + (tau_idx - 1)
  out[7] <- t(o7_idx) - holdt
  
  return(out)
}

# parameters ranges
p1r <- c(-50, 60)
p2r <- c(-60, 60)
p3r <- c(-80, 60)
p4r <- c(eps, 15)
p5r <- c(eps, 1.5)
p6r <- c(eps, 0.5)
p7r <- c(eps, 1.5)
p8r <- c(eps, 1)
p9r <- c(eps, 0.001)
p10r <- c(eps, 0.15)
p11r <- c(eps, 0.003)
p12r <- c(eps, 0.5)
gr <- c(eps, 2)

# input arguments
hold_volt <- -70
volt <- seq(-50, 50, by=10)
num_sample_per_volt <- 100

holdt <- 450
endt <- 5000

th <- seq(0, holdt)
tp <- seq(holdt + 1, endt)
tp_adj <- tp - tp[1]
t <- c(th, tp)

time_space <- vector(mode="list", 3)
time_space[[1]] <- t
time_space[[2]] <- th
time_space[[3]] <- tp_adj

# initial data with LHS
nvar <- 13
X <- randomLHS(length(volt)*num_sample_per_volt, nvar)
Y <- matrix(NA, ncol=7, nrow=nrow(X)) # 7: number of outputs

volt_idx <- 1
for(i in 1:nrow(X)) {
  # assign parameters
  param <- vector(mode="numeric", length=nvar)
  param[1] <- X[i,1]*diff(p1r) + p1r[1]
  param[2] <- X[i,2]*diff(p2r) + p2r[1]
  param[3] <- X[i,3]*diff(p3r) + p3r[1]
  param[4] <- X[i,4]*diff(p4r) + p4r[1]
  param[5] <- X[i,5]*diff(p5r) + p5r[1]
  param[6] <- X[i,6]*diff(p6r) + p6r[1]
  param[7] <- X[i,7]*diff(p7r) + p7r[1]
  param[8] <- X[i,8]*diff(p8r) + p8r[1]
  param[9] <- X[i,9]*diff(p9r) + p9r[1]
  param[10] <- X[i,10]*diff(p10r) + p10r[1]
  param[11] <- X[i,11]*diff(p11r) + p11r[1]
  param[12] <- X[i,12]*diff(p12r) + p12r[1]
  param[13] <- X[i,13]*diff(gr) + gr[1]

  # generate IKtof trace
  if(i %% num_sample_per_volt == 0) {volt_idx <- volt_idx + 1}
  yktof <- iktof(param, hold_volt, volt[volt_idx], time_space)
  
  # calculate outputs
  Y[i,] <- output1(time_space, yktof)
}
