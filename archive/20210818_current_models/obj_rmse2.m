function [z] = obj_rmse2(p, hold_volt, volts, time_space, yksum, Ek, norm_select)
    hold_idx = length(time_space{2});
    num_volts = length(volts);

    rmse_list = zeros(num_volts, 1);
    for i = 1:num_volts
        volt = volts(i);
        yksum_i = yksum(:, i);

        % negative noise of yksum_i: convert to 0
        yksum_i(yksum_i < 0) = 0;

        % generate current
        [~, ~, ~, yksum_hat] = kcurrent_model2(p, hold_volt, volt, time_space, Ek);

        % check validty of trace shape
        [~, peak_idx] = max(yksum_hat);
        check_pt1 = any(isnan(yksum_hat));
        check_pt2 = any(yksum_hat < 0);
        check_pt3 = var(yksum_hat(1:hold_idx)) > 0.4812e-4; % not stable at hold_volt
        check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse

        if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
            rmse_list(i) = 1e+3; % arbitrary big number
        else
            running_rmse = sqrt(mean((yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2));
            
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
