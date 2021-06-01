function [amp, time_stats] = trace_stat4exp(t, current_trc, hold_idx)
    % check validty of trace shape
    [peak, peak_idx] = max(current_trc);
    tau_magnitude = peak*exp(-1);
    [~, tau_idx] = min(abs(tau_magnitude - current_trc));

    check_pt1 = any(isnan(current_trc));
    check_pt2 = any(current_trc < 0);
    check_pt3 = var(current_trc(1:hold_idx)) > 1; % not stable at hold_volt
    check_pt4 = peak_idx < hold_idx; % not stable at hold_volt or too flat at pulse
    check_pt5 = tau_idx < hold_idx; 

    if(check_pt1 || check_pt2 || check_pt3 || check_pt4 || check_pt5)
        % -1: code for indicating not valid current trace shape
        amp = -1;
        tau = -1;
    end

    % amplitude
    amp = peak;

    % stats relavant to time
    time_stats = zeros(3,1);

    % time at 1/3 of interval [peak, tau]
    tau2_idx_step = floor((tau_idx - peak_idx)/3);
    tau2_idx = peak_idx + tau2_idx_step;
    time_stats(1) = t(tau2_idx);
    
    % tau
    time_stats(2) = t(tau_idx);

    % time at 2/3 of interval [tau, end_point]
    tau3_idx_step = floor((length(t) - tau_idx)/3);
    tau3_idx = tau_idx + 2*tau3_idx_step;
    time_stats(3) = t(tau3_idx);
end
