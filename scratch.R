library(tidyverse)
library(lhs)
library(laGP)

eps <- sqrt(.Machine$double.eps)
df <- read_csv("ko1.csv")

gr <- c(eps, 2)
pr <- c(170, 5000)

XU <- maximinLHS(100, 2)

zyM <- IK