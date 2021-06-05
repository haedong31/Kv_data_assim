function [z] = obj_multi(p, hold_volt, hold_idx, volts, t, yksum, Ek, param_select)
    num_volts = length(volts);
    eval_list = zeros(num_volts, 1);

    for i = 1:num_volts
        volt = volts(i);
        yksum_i = yksum(:, i);

        % generate time space for given voltage
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
            eval_list(i) = 1e+3; % arbitrary big number
        else
            % RMSE betwenn two curves
            rmse = sqrt(mean((yksum_i - yksum_hat).^2));

            % MSE of trace statistics
            key_stats = key_points_diff(t, yksum_i, yksum_hat);
            rms_key_stats = sqrt(mean(key_stats(:, 1) - key_stats(:, 2).^2));

            eval_list(i) = rmse + rms_key_stats;
        end   
    end
    z = sum(eval_list);
end

function [stats] = key_points_diff(t, yksum, yksum_hat)
    pts = zeros(4, 1);
    pts_hat = zeros(4, 1);
    
    % 1st time stamp from yksum
    [peak, peak_idx] = max(yksum);
    
    % truncate time space and current traces
    t_trunc = t(peak_idx:end);
    t_trunc = t_trunc - t_trunc(1);
    yksum_trunc  = yksum(peak_idx:end);
    yksum_hat_trunc = yksum_hat(peak_idx:end);
    
    % other three time stamps from yksum
    [~, pt2_idx] = min(abs(peak*0.75) - yksum_trunc);
    [~, pt3_idx] = min(abs(peak*exp(-1) - yksum_trunc));
    pt4_idx = length(t_trunc);

    % important current magnitudes in yksum
    pts(1) = peak;
    pts(2) = yksum_trunc(pt2_idx);
    pts(3) = yksum_trunc(pt3_idx);
    pts(4) = yksum_trunc(pt4_idx);

    % important current magnitudes in yksum_hat
    pts_hat(1) = yksum_hat(peak_idx);
    pts_hat(2) = yksum_hat_trunc(pt2_idx);
    pts_hat(3) = yksum_hat_trunc(pt3_idx);
    pts_hat(4) = yksum_hat_trunc(pt4_idx);

    stats = [pts; pts_hat];
end

% function [normalized] = trace_normalize(current_trace)
% 
% end
