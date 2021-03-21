## validate K+ currents functions -----
source("./kcurrent_functions/iktof.R")

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
