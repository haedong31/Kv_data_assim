function [z] = obj_multi(p, hold_volt, hold_idx, volts, t, yksum, Ek, param_select, norm_select)
    num_volts = length(volts);
    rmse_list = zeros(num_volts, 1);
    pts_diff = zeros(num_volts, 1);

    for i = 1:num_volts
        volt = volts(i);
        yksum_i = yksum(:, i);
        
        time_space = cell(1, 3);
        time_space{1} = t;
        time_space{2} = t(1:hold_idx(i));
        pulse_t = t(hold_idx(i)+1:end);
        pulse_t_adj = pulse_t - pulse_t(1);
        time_space{3} = pulse_t_adj;
        
        % generate current
        if param_select == true
            [~, ~, ~, ~, yksum_hat] = reduced_model(p, hold_volt, volt, time_space, Ek);
        elseif param_select == false        
            [~, ~, ~, ~, yksum_hat] = full_model(p, hold_volt, volt, time_space, Ek);
        end

        % check validty of trace shape
        [~, peak_idx] = max(yksum_hat);
        check_pt1 = any(isnan(yksum_hat));
        check_pt2 = any(yksum_hat < 0);
        check_pt3 = var(yksum_hat(1:hold_idx(i))) > 0.01; % not stable at hold_volt
        check_pt4 = peak_idx < hold_idx(i); % not stable at hold_volt of too flat at pulse

        if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
            rmse_list(i) = 1e+3; % arbitrary big number
        else
            % calculate RMSE
            running_rmse = sqrt(mean((yksum_i - yksum_hat).^2));
            yksum_key_pts = key_points(time_space, yksum_i);
            yksum_hat_key_pts = key_points(time_space, yksum_hat);
            
            % calculate differences of key points
            pts_diff(i) = sum(abs(yksum_key_pts - yksum_hat_key_pts));
            
            % normalize?
            if norm_select == true
                miny = min(yksum_i);
                maxy = max(yksum_i);
                rmse_list(i) = (running_rmse - miny) / (maxy - miny);
            else
                rmse_list(i) = running_rmse;                
            end
        end
    end
    z = sum(rmse_list) + sum(pts_diff);
end

function [pts] = key_points(time_space, current_trace)
    % peak and end point
    pts = zeros(2, 1);
    
    hold_idx = length(time_space{2});
    pulse_t_adj = time_space{3};
    
    % early 20%
    early_phase_idx = floor(length(pulse_t_adj)*0.2);
    
    % 1st time stamp from yksum
    current_trace_trunc = current_trace((hold_idx+1):end);
    current_trace_trunc = current_trace_trunc(1:early_phase_idx);

    pts(1) = max(current_trace_trunc);
    pts(2) = current_trace(end);
end
