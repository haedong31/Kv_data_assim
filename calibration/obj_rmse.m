function z = obj_rmse(p, kcurrent_model, model_struct, volt_space, time_space, yksum)
    hold_idx = time_space{4};
    volts = volt_space{2};
    num_volts = length(volts);

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{3} = time_space;
    protocol{4} = volt_space{3};

    rmse_list = zeros(num_volts, 1);
    for i = 1:num_volts
        yksum_i = yksum(:, i);
        protocol{2} = volts(i);

        [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol);
        rmse_list(i) = sqrt(mean((yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2));
    end
    z = sum(rmse_list);
end
