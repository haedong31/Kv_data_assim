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

output2 <- function(time_space, current_trc) {
  # outputs for IKss
	# output 1: time when reaching SS current
	# output 2: SS current
	out <- vector(mode="numeric", 2)

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

	return(out)
}
