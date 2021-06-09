function [z] = obj_multi2(p, hold_volt, hold_idx, volts, t, yksum, Ek, param_select, norm_select)
    w = 0.6;
    num_volts = length(volts);
    trace_diff = zeros(num_volts, 1);
    early_trace_diff = zeros(num_volts, 1);

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
            % arbitrary big number
            trace_diff(i) = 1e+3;
            early_trace_diff(i) = 1e+3;
        else
            % calculate similarity for entire domain
            rmse1 = sqrt(mean((yksum_i - yksum_hat).^2));
            
            % calculate similarity for early 20% 
            yksum_i_trunc = early_phase(time_space, yksum_i);
            yksum_hat_trunc = early_phase(time_space, yksum_hat);
            rmse2 = sqrt(mean((yksum_i_trunc - yksum_hat_trunc).^2));

            % normalize?
            if norm_select == true
                miny = min(yksum_i);
                maxy = max(yksum_i);
                
                trace_diff(i) = (rmse1 - miny) / (maxy - miny);
                early_trace_diff(i) = (rmse2 - miny) / (maxy - miny);
            else
                trace_diff(i) = rmse1;
                early_trace_diff(i) = rmse2;
            end
        end
    end    
    z = (1-w)*sum(trace_diff) + w*sum(early_trace_diff);
end

function [current_trace_trunc] = early_phase(time_space, current_trace)    
    hold_idx = length(time_space{2});
    pulse_t_adj = time_space{3};
    
    % early 20%
    early_phase_idx = floor(length(pulse_t_adj)*0.1);
    
    % 1st time stamp from yksum
    current_trace_trunc = current_trace((hold_idx+1):end);
    current_trace_trunc = current_trace_trunc(1:early_phase_idx);
end
