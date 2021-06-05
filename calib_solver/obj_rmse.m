function [z] = obj_rmse(p, hold_volt, hold_idx, volts, t, yksum, Ek, param_select, norm_select)
    num_volts = length(volts);
    rmse_list = zeros(num_volts, 1);

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

        % normalize?
        if norm_select == true
            yksum_i = trace_normalize(yksum_i);
            yksum_hat = trace_normalize(yksum_hat);
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
            rmse_list(i) = sqrt(mean((yksum_i - yksum_hat).^2));
        end
    end
    z = sum(rmse_list);
end

function [normalized] = trace_normalize(current_trace)
    % min-max normalization
    miny = min(current_trace);
    maxy = max(current_trace);

    normalized = (current_trace - miny) ./ (maxy - miny);
end
