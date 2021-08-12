function [z] = obj_multi2(p, hold_volt, volts, time_space, yksum, Ek, norm_select)
    early_portion = 0.1;
    w = 0.5;

    hold_idx = length(time_space{2});
    pulse_t = time_space{3};
    num_volts = length(volts);

    rmse1 = zeros(num_volts, 1);
    rmse2 = zeros(num_volts, 1);
    for i = 1:num_volts
        yksum_i = yksum(:, i);
        
        % generate current
        [~, ~, ~, ~, yksum_hat] = kcurrent_model(p, hold_volt, volts(i), time_space, Ek);

        % check validty of trace shape
        [~, peak_idx] = max(yksum_hat);
        check_pt1 = any(isnan(yksum_hat));
        check_pt2 = any(yksum_hat < 0);
        check_pt3 = var(yksum_hat(1:hold_idx)) > 0.4812e-6; % not stable at hold_volt
        check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse
        
        if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
            % arbitrary big number
            rmse1(i) = 1e+3;
            rmse2(i) = 1e+3;
        else
            early_phase_idx = floor(length(pulse_t)*early_portion);
            
            % model prediction
            yksum_hat1 = yksum_hat((hold_idx + 1):early_phase_idx);
            yksum_hat2 = yksum_hat((hold_idx + 1):end);
            
            % experimental data
            yksum_i1 = yksum_i((hold_idx + 1):early_phase_idx);
            yksum_i2 = yksum_i((hold_idx + 1):end);

            % calculate RMSEs
            rmse1(i) = sqrt(mean((yksum_i1 - yksum_hat1).^2));
            rmse2(i) = sqrt(mean((yksum_i2 - yksum_hat2).^2));

            % normalize?
            if norm_select == true
                miny = min(yksum_i);
                maxy = max(yksum_i);
                
                rmse1(i) = (rmse1(i) - miny) / (maxy - miny);
                rmse2(i) = (rmse2(i) - miny) / (maxy - miny);
            end
        end
    end    
    z = (1-w)*sum(rmse1) + w*sum(rmse2);
end
