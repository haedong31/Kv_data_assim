## validate K+ currents functions -----
source("./sensitivity_analysis/kcurrents.R")

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

p0 <- c(30.0, 13.5, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335, 0.4067, -91.1)
yktof <- iktof(p0, -70, 50, time_space)
plot(t, yktof, type="l")

p0 <- c(22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 270, 1050, 45.2, 5.7, 0.0629, -91.1)
yktos <- iktos(p0, -70, 50, time_space)  
plot(t, yktos, type="l")

p0 <- c(22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 1200.0, 170.0, 45.2, 5.7, 0.16, -91.1)
ykur <- ikur(p0, -70, 50, time_space)
plot(t, ykur, type="l")

p0 <- c(22.5, 7.7, 39.3, 0.0862, 13.17, 0.05, -91.1)
ykss <- ikss(p0, -70, 50, time_space)
plot(t, ykss, type="l")

# outputs
peak_time_idx <- which.max(yktof)
peak_time <- t[peak_time_idx]
peak <- yktof[peak_time_idx]

any(yktof < 0) || (peak_time < th)

# calculate tau
peak*exp(-1)
tau_idx <- which.min(abs(yktof - peak*exp(-1)))
t[tau_idx] - tail(th, 1)

# convergence value of yktof
p0 <- c(30.0, 13.5, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335, 0.4067, -91.1)
cv_ktof <- cv_ktof(p0, 50)
p0[13]*cv_ktof[1]^3*cv_ktof[3]
