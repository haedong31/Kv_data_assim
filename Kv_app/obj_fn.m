function [z] = obj_fn(p, hold_volt, volts, time_space, yksum, Ek, norm_select)
    num_volts = length(volts);
    rmse_list = zeros(num_volts, 1);

    for i = 1:num_volts
        volt = volts(i);
        yksum_i = yksum(:, i);
        
        % generate current
        [~, ~, ~, ~, yksum_hat] = reduced_model(p, hold_volt, volt, time_space, Ek);

        % check validty of trace shape
        hold_idx = length(time_space{2});
        [~, peak_idx] = max(yksum_hat);
        check_pt1 = any(isnan(yksum_hat));
        check_pt2 = any(yksum_hat < 0);
        check_pt3 = var(yksum_hat(1:hold_idx)) > 0.001; % not stable at hold_volt
        check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse

        if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
            rmse_list(i) = 1e+3; % arbitrary big number
        else
            running_rmse = sqrt(mean((yksum_i - yksum_hat).^2));
            
            if norm_select == true
                miny = min(yksum_i);
                maxy = max(yksum_i);
                rmse_list(i) = (running_rmse - miny) / (maxy - miny);
            else
                rmse_list(i) = running_rmse;
            end
        end
    end
    z = sum(rmse_list);
end
