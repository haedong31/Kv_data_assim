function [z] = obj_multi(p, hold_volt, hold_idx, volts, t, yksum, Ek, param_select)
    num_volts = length(volts);
    eval_list = length(num_volts, 1);

    for i = 1:num_volts
        volt = volts(i);

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
        hold_idx = length(time_space{2});
        [~, peak_idx] = max(yksum_hat);
        check_pt1 = any(isnan(yksum_hat));
        check_pt2 = any(yksum_hat < 0);
        check_pt3 = var(yksum_hat(1:hold_idx)) > 0.01; % not stable at hold_volt
        check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse

        if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
            eval_list(i) = 1e+3; % arbitrary big number
        else
            rmse = sqrt(mean((yksum - yksum_hat).^2));
            
        end   
    end
    z = sum(eval_list);
end

function [e] = trace_stat(t, current_trace)
        % statistics of current trace
    % stat 1: peak time
	% stat 2: peak current value
	% stat 3: tau, peak reduced by exp(-1) (~63%)
	% stat 4: current at 1/3 of interval [peak, tau]
	% stat 5: current at 2/3 of interval [peak, tau]
	% stat 6: current at last
    
    stats = zeros(6, 1);

    % truncate trace
    [peak, peak_idx] = max(current_trace);
    current_trace_trunc = current_trace(peak_idx:end);
    t_trunc = t(peak_idx:end);
    t_trunc = t_trunc - t_trunc(1);
    
    % stat 1
    stats(1) = t(peak_idx);
    
    % stat 2
    stats(2) = peak;
    
    % stat 3
    [~, tau_idx] = min(abs(peak*exp(-1) - current_trace_trunc));
    stats(3) = t_trunc(tau_idx);

    % stat 4
    [~, stat45_jump_size] = min(abs(t_trunc(tau_idx)/3 - t_trunc));
    stats(4) = current_trace_trunc(stat45_jump_size);

    % stat 5
    stats(5) = current_trace_trunc(stat45_jump_size*2);

    % stat 6
    stats(6) = current_trace(end);
end
