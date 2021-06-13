hh_model <- function(t, ss0, ss, tau) {
    y <- ss - (ss - ss0)*exp(-t/tau)
    return(y)
}

ikto <- function(p, hold_volt, volt, time_space, Ek=-91.1) {
    # see 2014 Bondarenko and 2020 Bondarenko
    # 17 parameters 

    # conductance and phosphorylation fraction
    f_eacv <- p[15] # 0.2087704319 Ikto fraction phophorylated
    gmax <- p[16] # 0.14067
    gmaxp <- p[16]*(1.0-p[17]) # 0.14067*(1.0-0.387)

    # initial state values
    act0 <- 0.4139033547E-02
    inact0 <- 0.9999623535E+00
    actp0 <- 0.8637307739E-03
    inactp0 <- 0.9999881755E+00
    
    # time space information
    t <- time_space[[1]]
    hold_t <- time_space[[2]]
    pulse_t <- time_space[[3]]
    hold_idx <- length(hold_t)

    end <- length(t)
    current_trace <- vector(mode="numeric", end)

    # at holding potential
    gv_hold <- ikto_gating(p, hold_volt)

    act_hold <- hh_model(hold_t, act0, gv_hold[1], gv_hold[3])
    inact_hold <- hh_model(hold_t, inact0, gv_hold[2], gv_hold[4])      

    actp_hold <- hh_model(hold_t, actp0, gv_hold[5], gv_hold[7])
    inactp_hold <- hh_model(hold_t, inactp0, gv_hold[6], gv_hold[8])

    current_trace[1:hold_idx] <- (1-f_eacv)*gmax*(act_hold^3)*(inact_hold)*(hold_volt-Ek) + 
        f_eacv*gmaxp*(actp_hold^3)*(inactp_hold)*(hold_volt-Ek)

    # at pulse voltage
    gv_pulse <- ikto_gating(p, volt)

    act_pulse <- hh_model(pulse_t, act0, gv_pulse[1], gv_pulse[3])
    inact_pulse <- hh_model(pulse_t, inact0, gv_pulse[2], gv_pulse[4])

    actp_pulse <- hh_model(pulse_t, actp0, gv_pulse[5], gv_pulse[7])
    inactp_pulse <- hh_model(pulse_t, inactp0, gv_pulse[6], gv_pulse[8])

    current_trace[(hold_idx+1):end] <- (1-f_eacv)*gmax*(act_pulse^3)*(inact_pulse)*(volt-Ek) + 
        f_eacv*gmaxp*(actp_pulse^3)*(inactp_pulse)*(volt-Ek)

    return(current_trace)
}

ikto_gating <- function(p, volt) {
    # p0 <- c(33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956,
    #         0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387)

    gv <- vector(mode="numeric", 8)
    
    # for non-phosphorylated
    alpha1 <- p[9]*exp(p[7]*(volt+p[1]))
    beta1 <- p[10]*exp(-p[8]*(volt+p[1]))
    
    alpha2_temp1 <- p[11]*exp((volt+p[2])/(-1.0*p[6]))
    alpha2_temp2 <- p[12]*exp((volt+p[2]+p[3])/(-1.0*p[6]))
    alpha2 <- alpha2_temp1/(1.0+alpha2_temp2)

    beta2_temp1 <- p[13]*exp((volt+p[2]+p[3])/p[6])
    beta2_temp2 <- p[14]*exp((volt+p[2]+p[3])/p[6])
    beta2 <- beta2_temp1/(1.0+beta2_temp2)

    gv[1] <- alpha1/(alpha1+beta1)
    gv[2] <- alpha2/(alpha2+beta2)
    gv[3] <- 1/(alpha1+beta1)
    gv[4] <- 1/(alpha2+beta2)

    # for phosphorylated
    alpha1p <- p[9]*exp(p[7]*(volt+p[1]-p[4]))
    beta1p <- p[10]*exp(-p[8]*(volt+p[1]-p[4]))
    
    alpha2_temp1p <- p[11]*exp((volt+p[2]-p[5])/(-1.0*p[6]))
    alpha2_temp2p <- p[12]*exp((volt+p[2]+p[3]-p[5])/(-1.0*p[6]))
    alpha2p <- alpha2_temp1p/(1.0+alpha2_temp2p)
    
    beta2_temp1p <- p[13]*exp((volt+p[2]+p[3]-p[5])/p[6])
    beta2_temp2p <- p[14]*exp((volt+p[2]+p[3]-p[5])/p[6])
    beta2p <- beta2_temp1p/(1.0+beta2_temp2p)

    gv[5] <- alpha1p/(alpha1p+beta1p)
    gv[6] <- alpha2p/(alpha2p+beta2p)
    gv[7] <- 1/(alpha1p+beta1p)
    gv[8] <- 1/(alpha2p+beta2p)

    return(gv)
}
