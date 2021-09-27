function z = obj_rmse_norm(p, kcurrent_model, model_struct, volt_space, time_space, yksum)
    hold_idx = time_space{4};
    end_idx = time_space{5};
    
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

        miny = min(yksum_i);
        maxy = max(yksum_i);

        [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol);
        se = (yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2;
        running_rmse = sqrt(mean(se));

        rmse_list(i) = (running_rmse-miny)/(maxy-miny);
    end
    z = sum(rmse_list);
end
