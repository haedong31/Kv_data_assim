function [amp, tau] = trace_stat(t, current_trc, hold_idx)
    % check validty of trace shape
    [peak, peak_idx] = max(current_trc);
    tau_magnitude = peak*exp(-1);
    [~, tau_idx] = min(abs(tau_magnitude - current_trc));

    check_pt1 = any(isnan(current_trc));
    check_pt2 = any(current_trc < 0);
    check_pt3 = var(current_trc(1:hold_idx)) > 1e-3; % not stable at hold_volt
    check_pt4 = peak_idx < hold_idx; % not stable at hold_volt or too flat at pulse
    check_pt5 = tau_idx < hold_idx; 

    if(check_pt1 || check_pt2 || check_pt3 || check_pt4 || check_pt5)
        % -1: code for indicating not valid current trace shape
        amp = -1;
        tau = -1;
    else
        % amplitude
        amp = peak;

        % tau
        tau = t(tau_idx);
    end
end
