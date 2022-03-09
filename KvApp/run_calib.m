function [sol,fval] = run_calib(y, model_info, protocol_info, init_pts, lb, ub)            
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

    opt_fun = @(p) obj_rmse(p, y, model_info, protocol_info);
    [sol,fval] = fmincon(opt_fun, init_pts, A, b, Aeq, beq, lb, ub, nonlcon, options);
end

function z = obj_rmse(p, y, model_info, protocol_info)
    time_space = protocol_info{4};
    hold_idx = length(time_space{2});
    
    volts = protocol_info{3};
    num_volts = length(volts);

    rmse_list = NaN(num_volts, 1);
    for i = 1:num_volts
        y_i = y(:, i);
        protocol_info{3} = volts(i);

        [y_hat, ~] = kcurrent_model(p, model_info, protocol_info);
        rmse_list(i) = sqrt(mean((y_i((hold_idx + 1):end) - y_hat((hold_idx + 1):end)).^2));
    end
    z = sum(rmse_list);
end

function [yksum, comp_currents] = kcurrent_model(p, model_info, protocol_info)
    global param_kslow1
    
    num_currents = length(model_info);
    current_names = cell(num_currents, 1);
    
    for i = 1:num_currents
        current_names{i} = model_info(i).name;
    end
    
    % declare shared parameters of ikslow1 as global variable
    matching_idx = strcmp(current_names, 'ikslow1');
    if any(matching_idx)    
        num_kslow1_param = 11;
        kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, ...
            6.1, 18.0, 2.058, 803.0, 0.16];
        
        tune_idx1 = model_info(matching_idx).idx1;
        tune_idx2 = model_info(matching_idx).idx2;
        fixed_idx = setdiff(1:num_kslow1_param, tune_idx1);

        param_kslow1 = NaN(num_kslow1_param, 1);
        param_kslow1(tune_idx1) = p(tune_idx2);
        param_kslow1(fixed_idx) = kslow1_default(fixed_idx);
    end

    time_space = protocol_info{4};
    yksum = zeros(length(time_space{1}), 1);
    comp_currents = cell(num_currents, 1);
    for i = 1:num_currents
        comp_currents{i} = gen_matching_current(p, model_info(i), protocol_info);
        yksum = yksum + comp_currents{i};
    end
end

function current_trace = gen_matching_current(p, model_info, protocol_info)
    global param_kslow1

    % current model info
    current_name = model_info.name;
    tune_idx1 = model_info.idx1;
    tune_idx2 = model_info.idx2;

    % protocol info
    hold_volt = protocol_info{1};
    ek = protocol_info{2};
    volt = protocol_info{3};
    time_space = protocol_info{4};

    % generate current
    switch current_name
        case 'iktof'
            % generate ikto
            num_param = 13;

            ktof_default = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
                0.000152, 0.067083, 0.00095, 0.051335, 0.3846];
            fixed_idx = setdiff(1:num_param, tune_idx1);
            
            param = NaN(num_param, 1);
            param(tune_idx1) = p(tune_idx2);
            param(fixed_idx) = ktof_default(fixed_idx);
            
            current_trace = iktof(param, hold_volt, volt, time_space, ek);
        case 'ikslow1'
            % generate ikslow1
            current_trace = ikslow1(param_kslow1, hold_volt, volt, time_space, ek);
        case 'ikslow2'
            % generate ikslow2
            num_param = 11;
            kslow2_default = [4912, 5334, 0.16];
            shared_idx = [1:7, 9];
            uniq_idx = setdiff(1:num_param, shared_idx);

            param = NaN(num_param, 1);
            param(shared_idx) = param_kslow1(shared_idx);

            uniq_param = kslow2_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param(uniq_idx) = uniq_param;

            current_trace = ikslow2(param, hold_volt, volt, time_space, ek);
        case 'ikss'
            % generate ikss
            num_param = 7;
            kss_default = [0.0862, 1235.5, 13.17, 0.0611];
            shared_idx1 = 1:3;
            shared_idx2 = [1, 3, 4];
            uniq_idx = setdiff(1:num_param, shared_idx1);

            param = NaN(num_param, 1);
            param(shared_idx1) = param_kslow1(shared_idx2);

            uniq_param = kss_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param(uniq_idx) = uniq_param;

            current_trace = ikss(param, hold_volt, volt, time_space, ek);
    end
end
