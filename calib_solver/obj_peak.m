function [z] = obj_peak(p, hold_volt, volt, time_space, Ek, peak_vals, param_select)
    if param_select == true
        [ykto, ykslow1, ykslow2, ykss, ~] = reduced_model(p, hold_volt, volt, time_space, Ek);
    elseif param_select == false
        [ykto, ykslow1, ykslow2, ykss, ~] = full_model(p, hold_volt, volt, time_space, Ek);
    end

    % get peaks
    hold_idx = length(time_space{2});
    pkto = get_peak(ykto, hold_idx);
    pkslow1 = get_peak(ykslow1, hold_idx);
    pkslow2 = get_peak(ykslow2, hold_idx);
    pkss = get_peak(ykss, hold_idx);

    % objective
    sim_peaks = [pkto, pkslow1, pkslow2, pkss];
    if any(sim_peaks == -1)
        z = 1e+3; % arbtrary big number
    else
        z = sum(abs(sim_peaks - peak_vals));
    end
end

function [peak] = get_peak(current_trc, hold_idx)
    % check validty of trace shape
    [peak, peak_idx] = max(current_trc);
    check_pt1 = any(isnan(current_trc));
    check_pt2 = any(current_trc < 0);
    check_pt3 = var(current_trc(1:hold_idx)) > 1e-3; % not stable at hold_volt
    check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse

    if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
        peak = -1; % code for indicating not valid current trace and peak
    end
end
