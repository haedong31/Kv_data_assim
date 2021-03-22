library(lhs)
library(laGP)

source("./sensitivity_analysis/kcurrents.R")
source("./sensitivity_analysis/outputs.R")
eps <- .Machine$double.eps

scale_param <- function(unit_param) {
  # parameters ranges
  p1r <- c(-50, 60)
  p2r <- c(-60, 60)
  p3r <- c(-80, 60)
  p4r <- c(1, 15)
  p5r <- c(eps, 1.5)
  p6r <- c(eps, 0.5)
  p7r <- c(eps, 1.5)
  p8r <- c(eps, 1)
  p9r <- c(eps, 0.001)
  p10r <- c(eps, 0.15)
  p11r <- c(eps, 0.003)
  p12r <- c(eps, 0.5)
  gr <- c(eps, 2)
  
  param <- vector(mode="numeric", length=nvar)
  param[1] <- unit_param[1]*diff(p1r) + p1r[1]
  param[2] <- unit_param[2]*diff(p2r) + p2r[1]
  param[3] <- unit_param[3]*diff(p3r) + p3r[1]
  param[4] <- unit_param[4]*diff(p4r) + p4r[1]
  param[5] <- unit_param[5]*diff(p5r) + p5r[1]
  param[6] <- unit_param[6]*diff(p6r) + p6r[1]
  param[7] <- unit_param[7]*diff(p7r) + p7r[1]
  param[8] <- unit_param[8]*diff(p8r) + p8r[1]
  param[9] <- unit_param[9]*diff(p9r) + p9r[1]
  param[10] <- unit_param[10]*diff(p10r) + p10r[1]
  param[11] <- unit_param[11]*diff(p11r) + p11r[1]
  param[12] <- unit_param[12]*diff(p12r) + p12r[1]
  param[13] <- unit_param[13]*diff(gr) + gr[1]
  
  return(param)
}

# input arguments
hold_volt <- -70
volt <- 50
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
X <- randomLHS(1300, nvar)
Y <- matrix(NA, ncol=7, nrow=nrow(X)) # 7: number of outputs]
na_row_idx <- NULL

for(i in 1:nrow(X)) {
  # assign parameters
  param <- scale_param(X[i,])

  # generate IKtof trace
  yktof <- iktof(param, hold_volt, volt, time_space)
  
  # calculate outputs
  o <- output1(time_space, yktof)
  if(any(is.na(o))) {na_row_idx <- append(na_row_idx, i)}
    
  Y[i,] <- o
}

# drop NA
X <- X[-na_row_idx, ]
Y <- Y[-na_row_idx, ]

# fit
out <- Y[,2]
gpi <- newGPsep(X, out, d=0.1, g=var(out)/10, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=c(eps,2), tmax=c(10, var(out)))
