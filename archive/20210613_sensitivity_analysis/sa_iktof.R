library(lhs)
library(laGP)

source("./sensitivity_analysis/kcurrents.R")
source("./sensitivity_analysis/outputs.R")
eps <- .Machine$double.eps

# sensitivity analysis arguments
nvar <- 13
N <- 10000
G <- 30
grid <- seq(0, 1, length=G)
XX <- matrix(NA, ncol=nvar, nrow=N)

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

main_effect <- function(gpi) {
  m <- q1 <- q2 <- matrix(NA, ncol=nvar, nrow=G)
  for(j in 1:nvar) {
    for(i in 1:G) {
      XX[,j] <- grid[i]
      XX[,-j] <- randomLHS(N, nvar-1)
      
      p <- predGPsep(gpi, XX, lite=TRUE, nonug=TRUE)
      
      m[i,j] <- mean(p$mean)
      q1[i,j] <- mean(qnorm(0.05, p$mean, sqrt(p$s2)))
      q2[i,j] <- mean(qnorm(0.95, p$mean, sqrt(p$s2)))
    }
  }

  return(list(m=m, q1=q1, q2=q2))
}

first_sensitivity <- function(gpi) {
  M <- randomLHS(N, nvar)
  pM <- predGPsep(gpi, M, lite=TRUE, nonug=TRUE)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2

  Mprime <- randomLHS(N, nvar)
  S <- EE2j <- rep(NA, nvar)
  for(j in 1:nvar) {
    Mjprime <- Mprime
    Mjprime[,j] <- M[,j]
    
    pMprime <- predGPsep(gpi, Mjprime, lite=TRUE, nonug=TRUE)
    EE2j[j] <- (t(pM$mean) %*% pMprime$mean)/(N - 1)
    S[j] <- (EE2j[j] - Ey^2)/Vary
  }

  return(S)
}

## data preparation -----
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

## sensitivity analysis -----
# output 1
out1 <- Y[,1]
gpi <- newGPsep(X, out1, d=0.1, g=var(out1)/10, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=c(eps,2), tmax=c(10,var(out1)))

# main effects
m1 <- q11 <- q21 <- matrix(NA, ncol=nvar, nrow=G)
for(j in 1:nvar) {
  for(i in 1:G) {
    XX[,j] <- grid[i]
    XX[,-j] <- randomLHS(N, nvar-1)
    
    p <- predGPsep(gpi, XX, lite=TRUE, nonug=TRUE)
    
    m1[i,j] <- mean(p$mean)
    q11[i,j] <- mean(qnorm(0.05, p$mean, sqrt(p$s2)))
    q12[i,j] <- mean(qnorm(0.95, p$mean, sqrt(p$s2)))
  }
}

# first-order sensitivity
M <- randomLHS(N, nvar)
pM <- predGPsep(gpi, M, lite=TRUE, nonug=TRUE)
Ey <- mean(pM$mean)
Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2

Mprime <- randomLHS(N, nvar)
S <- EE2j <- rep(NA, nvar)
for(j in 1:nvar) {
  Mjprime <- Mprime
  Mjprime[,j] <- M[,j]
  
  pMprime <- predGPsep(gpi, Mjprime, lite=TRUE, nonug=TRUE)
  EE2j[j] <- (t(pM$mean) %*% pMprime$mean)/(N - 1)
  S[j] <- (EE2j[j] - Ey^2)/Vary
}

deleteGPsep(gpi)

# output 2 
out2 <- Y[,2]
gpi <- newGPsep(X, out2, d=0.1, g=var(out2)/100, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=c(eps,2), tmax=c(10,var(out2)/10))

me2 <- main_effect(gpi)
S2 <- first_sensitivity(gpi)

deleteGPsep(gpi)

## visualize
m <- me2$m
q1 <- me2$q1
q2 <- me2$q2
plot(0, xlab="grid", ylab="main effect", xlim=c(0,1), ylim=range(c(q1,q2)), type="n")
for(j in 1:nvar) {
  lines(grid, m[,j], col=j, lwd=2)
}
legend("bottomright", paste0("x",1:nvar), fill=1:nvar, horiz=TRUE, cex=0.75)

