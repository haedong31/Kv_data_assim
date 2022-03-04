function sol = run_calib(y, model_info, protocol_info, init)            
    % optimization options
    options = optimoptions(@fmincon, ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxFunctionEvaluations',1e+6, ...
        'SpecifyObjectiveGradient',false);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];

    if init == 0
        ktof_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, ...
            0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
        ktos_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, ...
            803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
        kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, ...
            803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
        kslow2_default = [5334, 4912, 0.05766];
        kss_default = [0.0862, 1235.5, 13.17, 0.0428];
    else
        init
    end
end

function z = obj_rmse(p, kcurrent_model, model_struct, volt_space, time_space, y)
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
        y_i = y(:, i);
        protocol{2} = volts(i);

        [y_hat, ~] = kcurrent_model(p, model_struct, protocol);
        rmse_list(i) = sqrt(mean((y_i((hold_idx + 1):end) - y_hat((hold_idx + 1):end)).^2));
    end
    z = sum(rmse_list);
end
