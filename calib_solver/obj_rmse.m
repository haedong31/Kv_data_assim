function [z] = obj_rmse(p, hold_volt, hold_idx, volts, t, yksum, Ek, param_select)
    num_volts = length(volts);
    rmse_list = zeros(num_volts, 1);

    for i = 1:num_volts
        volt_idx = i;
        
        time_space = cell(1, 3);
        time_space{1} = t;
        time_space{2} = t(1:hold_idx(volt_idx));
        pulse_t = t(hold_idx(volt_idx)+1:end);
        pulse_t_adj = pulse_t - pulse_t(1);
        time_space{3} = pulse_t_adj;
        
        rmse_list(volt_idx) = rmse(p, hold_volt, volts(volt_idx), time_space, Ek, yksum(:, volt_idx), param_select); 
    end
    z = sum(rmse_list);
end

function [z] = rmse(p, hold_volt, volt, time_space, Ek, yksum, param_select)
    if param_select == true
        [~, ~, ~, ~, yksum_hat] = reduced_model(p, hold_volt, volt, time_space, Ek);
    elseif param_select == false        
        [~, ~, ~, ~, yksum_hat] = full_model(p, hold_volt, volt, time_space, Ek);
    end
    
    % check validty of trace shape
    hold_idx = length(time_space{2});
    [~, peak_idx] = max(yksum_hat);
    check_pt1 = any(isnan(yksum_hat));
    check_pt2 = any(yksum_hat < 0);
    check_pt3 = var(yksum_hat(1:hold_idx)) > 0.01; % not stable at hold_volt
    check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse
    
    if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
        z = 1e+3; % arbitrary big number
    else
        z = sqrt(mean((yksum_hat-yksum).^2));
    end
end
