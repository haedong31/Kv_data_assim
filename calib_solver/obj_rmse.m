function [z] = obj_rmse(p, model_struct, time_space, volt_space)
    volts = 
    time_space = protocol{3};
    hold_idx = length(time_space{2});
    num_volts = length(volts);

    rmse_list = zeros(num_volts, 1);    
    for i = 1:num_volts
        volt = volts(i);
        yksum_i = yksum(:, 1);

        [yksum, ~] = kcurrent_model(p0, model_struct, protocol);

        running_rmse = sqrt(mean((yksum_i((hold_idx+1):end) - yksum_hat((hold_idx+1):end)).^2));

        rmse_list(i) = running_rmse;
    end
    z = sum(rmse_list);
end
