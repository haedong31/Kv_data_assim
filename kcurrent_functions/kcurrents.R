iktof <- function(param, holdv, pulse, time_space) {
    # Bondarenko IKtof (A67)
    # 14 parameters

    # constants
    act0 <- 0.265563e-2 # ato_f Gating variable for transient outward K+ current
    inact0 <- 0.999977 # ito_f Gating variable for transient outward K+ current
    gmax <- param[13]  # 0.4067
    Ek <- param[14] # resting potential

    # time space information
    t <- time_space[[1]]
    tH <- time_space[[2]]
    tP1_adj <- time_space[[3]]
    hold_idx <- length(tH)

    end <- length(t)
    current_trc <- vector(mode<-"numeric", end)

    # at holding potential
    cv_hold <- cv_ktof(param[1:12], holdv)
    act_hold <- cv_hold[1] - (cv_hold[1] - act0)*exp(-(tH/cv_hold[2]))
    inact_hold <- cv_hold[3] - (cv_hold[3] - inact0)*exp(-(tH/cv_hold[4]))
    current_trc[1:hold_idx] <- gmax*(act_hold^3)*(inact_hold)*(holdv - Ek)

    # at P1
    cv_P1 <- cv_ktof(param[1:12], pulse)
    act_P1 <- cv_P1[1] - (cv_P1[1] - act0)*exp(-(tP1_adj/cv_P1[2]))
    inact_P1 <- cv_P1[3] - (cv_P1[3] - inact0)*exp(-(tP1_adj/cv_P1[4]))
    current_trc[(hold_idx + 1):end] <- gmax*(act_P1^3)*(inact_P1)*(pulse - Ek)

    return(current_trc)
}

iktos <- function(param. holdv, pulse, time_space) {
    # Bondarenko IKtos (A75)
    # 11 parameters + 1 conductance + 1 Ek

    act0 <- 0.417069e-3 # ato_s Gating variable for transient outward K+ current
    inact0 <- 0.998543 # ito_s Gating variable for transient outward K+ current
    gmax <- param(12) # 0.0629 in spetal model 0 in apical model
    Ek <- param(13) # resting potential

    # time space information
    t <- time_space[[1]]
    tH <- time_space[[2]]
    tP1_adj <- time_space[[3]]
    hold_idx <- length(tH)

    end <- length(t)
    current_trc <- vector(mode<-"numeric", end)

    # at holding potential
    cv_hold <- cv(param(1:11), holdV)
    act_hold <- cv_hold[1] - (cv_hold[1] - act0)*exp(-(tH/cv_hold[2]))
    inact_hold <- cv_hold[3] - (cv_hold[3] - inact0)*exp(-(tH/cv_hold[4]))
    current_trc(1:hold_idx) <- gmax*(act_hold)*(inact_hold)*(holdV - Ek)

    # at P1
    cv_P1 <- cv(param(1:11), P1)
    act_P1 <- cv_P1[1] - (cv_P1[1] - act0)*exp(-(tP1_adj/cv_P1[2]))
    inact_P1 <- cv_P1[3] - (cv_P1[3] - inact0)*exp(-(tP1_adj/cv_P1[4]))
    current_trc((hold_idx + 1):end) <- gmax*(act_P1)*(inact_P1)*(P1 - Ek)

    return(current_trc)
}

cv_ktof <- function(p, V) {
    # characteristic variables of IKtof
    # cv[1]: steady-state activation
    # cv[2]: time constant of activation
    # cv[3]: steady-state inactivation
    # cv[4]: time constant of inactivation
    
    # p0 <- c(30.0, 13.5, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335)
    cv <- vector(mode<-"numeric", 4) 
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

cv_ktos <- function(p, V) {
    # p0 <- c(22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 270, 1050, 45.2, 5.7)
    cv <- vector(mode="numeric", 4)

    # note similarity of IKtos and IKur
    # same activation gate characteristics
    cv[1] <- 1 / (1 + exp(-(V+p[1])/p[2]))
    cv[2] <- p[3]*exp(-p[4]*V) + p[5]

    # for inactivation gate, cv[3] is same but not for cv[4] 
    cv[3] <- 1 / (1 + exp((V+p[6])/p[7]))
    cv[4] <- p[8] + p[9]/(1 + exp((V + p(10))/p(11)))

    return(cv)
}