function [z] = obj_current_stats(p, hold_volt, hold_idx, volts, t, yksum, Ek)
    num_volts = length(volts);

    mse_list = zeros(num_volts, 1);
    for i = 1:num_volts
        volt_idx = i;
        volt = volts(volt_idx);

        time_space = cell(1,3);
        time_space{1} = t;
        time_space{2} = t(1:hold_idx(volt_idx));
        pulse_t = t(hold_idx(volt_idx)+1:end);
        pulse_t_adj = pulse_t - pulse_t(1);
        time_space{3} = pulse_t_adj;

        mse_list(i) = stat_mse(p, hold_volt, volt, time_space, yksum(:, volt_idx), Ek);
    end
    z = mean(mse_list);
end

function [mse] = stat_mse(p, hold_volt, volt, time_space, yksum, Ek)
    [~, ~, ~, ~, yksum_hat] = reduced_model(p, hold_volt, volt, time_space, Ek);

    % time information
    t = time_space{1};
    hold_t = time_space{2};
    hold_idx = length(hold_t);

    % check validty of trace shape
    [peak, peak_idx] = max(yksum_hat);
    yksum_hat_trunc = yksum_hat(peak_idx:end);
    [~, tau_idx] = min(abs(peak*exp(-1) - yksum_hat_trunc));

    check_pt1 = any(isnan(yksum_hat));
    check_pt2 = any(yksum_hat < 0);
    check_pt3 = var(yksum_hat(1:hold_idx)) > 0.001; % not stable at hold_volt
    check_pt4 = peak_idx < hold_idx; % not stable at hold_volt or too flat at pulse
    check_pt5 = tau_idx < hold_idx; 

    if (check_pt1 || check_pt2 || check_pt3 || check_pt4 || check_pt5)
        mse = 1e+4;
    else
        stats = trace_stat(t, yksum);
        stats_hat = trace_stat(t, yksum_hat);
    
        mse = sqrt(mean((stats - stats_hat).^2));
    end
end
