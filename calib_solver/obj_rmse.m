function [z] = obj_rmse(p, hold_volt, volt, time_space, Ek, exp_ksum)
    [~, ~, ~, ~, yksum] = full_model(p, hold_volt, volt, time_space, Ek);

    % check validty of trace shape
    hold_idx = length(time_space{2});
    [~, peak_idx] = max(yksum);
    check_pt1 = any(isnan(yksum));
    check_pt2 = any(yksum < 0);
    check_pt3 = var(yksum(1:hold_idx)) > 1e-3; % not stable at hold_volt
    check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse
    
    if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
        z = 1e+3; % arbitrary big number
    else
        z = sqrt(mean((yksum-exp_ksum).^2));
    end
end
