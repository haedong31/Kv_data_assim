function [z] = obj_rmse(p, hold_volt, volts, time_space, yksum, Ek)
    hold_idx = length(time_space{2});
    num_volts = length(volts);

    rmse_list = zeros(num_volts, 1);    
    for i = 1:num_volts
        volt = volts(i);
        yksum_i = yksum(:, 1);

        [yksum_hat, ~] = kcurrent_model(p, hold_volt, volt, time_space, Ek);

        running_rmse = sqrt(mean((yksum_i((hold_idx+1):end) - yksum_hat((hold_idx+1):end)).^2));

        rmse_list(i) = running_rmse;
    end
    z = sum(rmse_list);
end
