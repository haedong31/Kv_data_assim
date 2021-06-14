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

# protocol
hold_volt <- -70
volt <- 50
hold_len <- 0.45 * 1e+3
end_len <- 5 * 1e+3

# time space
hold_t <- seq(0, hold_len)
pulse_t <- seq(hold_len + 1, end_len)
pulse_t_adj <- pulse_t - pulse_t[1]
t <- c(hold_t, pulse_t)

time_space <- vector(mode="list", 3)
time_space[[1]] <- t
time_space[[2]] <- hold_t
time_space[[3]] <- pulse_t_adj

scale_param <- function(unit_param) {
  # parameters ranges
  p1r <- c(-70, 50)
  p2r <- c(-70, 50)
  p3r <- c(-70, 50)
  p4r <- c(eps, 10)
  p5r <- c(eps, 50)
  p6r <- c(eps, 50)
  p7r <- c(eps, 1)
  p8r <- c(eps, 100)
  p9r <- c(eps, 1000)
  p10r <- c(eps, 50)
  p11r <- c(0, 1)
  p12r <- c(0, 0.5)
  p13r <- c(0, 0.5)
  
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
  param[13] <- unit_param[13]*diff(p13r) + p13r[1]
  
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

# p0 <- c(33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956,
#         0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387)
# yktof <- ikto(p0, hold_volt, volt, time_space)

##### main effects -----
# initial data with LHS
X <- randomLHS(2000, nvar)
Y <- matrix(NA, ncol=7, nrow=nrow(X)) # 7: number of outputs]
na_row_idx <- NULL

for(i in 1:nrow(X)) {
  # assign parameters
  param <- scale_param(X[i,])
  
  # generate current trace
  ykslow1 <- ikslow1(param, hold_volt, volt, time_space)
  
  # calculate outputs
  o <- current_output(time_space, ykslow1)
  if(any(is.na(o))) {na_row_idx <- append(na_row_idx, i)}
  
  Y[i,] <- o
}

# drop NA
X <- X[-na_row_idx, ]
Y <- Y[-na_row_idx, ]

# output 1
out1 <- Y[,1]
gpi <- newGPsep(X, out1, d=0.01, g=var(out1)/100, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=c(eps,2), tmax=c(10,var(out1)/10))
me1 <- main_effect(gpi)
deleteGPsep(gpi)

# output 2 
out2 <- Y[,2]
gpi <- newGPsep(X, out2, d=0.01, g=var(out2)/100, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=c(eps,2), tmax=c(10,var(out2)/10))
me2 <- main_effect(gpi)

# visualization
m <- me1$m
q1 <- me1$q1
q2 <- me1$q2
plot(0, xlab="grid", ylab="main effect", xlim=c(0,1), ylim=range(c(q1,q2)), type="n")
for(j in 1:nvar) {
  lines(grid, m[,j], col=j, lwd=2)
}
legend("bottomright", paste0("x",1:nvar), fill=1:nvar, horiz=TRUE, cex=0.75)

m <- me2$m
q1 <- me2$q1
q2 <- me2$q2
plot(0, xlab="grid", ylab="main effect", xlim=c(0,1), ylim=range(c(q1,q2)), type="n")
for(j in 1:nvar) {
  lines(grid, m[,j], col=j, lwd=2)
}
legend("bottomright", paste0("x",1:nvar), fill=1:nvar, horiz=TRUE, cex=0.75)

# save results
save(Y, me1, me2, file = './sensitivity_analysis/ikslow1_me12.RData')
