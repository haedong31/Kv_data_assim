function [yksum, comp_currents] = kcurrent_model(p, model_struct, protocol_info)
    global param_kslow1
    
    num_currents = length(model_struct);
    current_names = cell(num_currents, 1);
    
    for i = 1:num_currents
        current_names{i} = model_struct(i).name;
    end
    
    % declare shared parameters of ikslow1 as global variable
    matching_idx = strcmp(current_names, 'ikslow1');
    if any(matching_idx)    
        num_kslow1_param = 13;
        kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, ...
            2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
        
        tune_idx1 = model_struct(matching_idx).idx1;
        tune_idx2 = model_struct(matching_idx).idx2;
        fixed_idx = setdiff(1:num_kslow1_param, tune_idx1);

        param_kslow1 = zeros(num_kslow1_param, 1);
        param_kslow1(tune_idx1) = p(tune_idx2);
        param_kslow1(fixed_idx) = kslow1_default(fixed_idx);
    end

    time_space = protocol_info{3};
    yksum = zeros(length(time_space{1}), 1);
    comp_currents = cell(num_currents, 1);
    for i = 1:num_currents
        comp_currents{i} = gen_matching_current(p, model_struct(i), protocol_info);
        yksum = yksum + comp_currents{i};
    end
end

function [current_trace] = gen_matching_current(p, model_info, protocol_info)
    global param_kslow1

    % current model info
    current_name = model_info.name;
    tune_idx1 = model_info.idx1;
    tune_idx2 = model_info.idx2;

    % protocol info
    hold_volt = protocol_info{1};
    volt = protocol_info{2};
    time_space = protocol_info{3};
    ek = protocol_info{4};

    % generate current
    switch current_name
        case "ikto"
            % generate ikto
            num_param = 17;
            kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
                0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
            fixed_idx = setdiff(1:num_param, tune_idx1);
            
            param = zeros(num_param, 1);
            param(tune_idx1) = p(tune_idx2);
            param(fixed_idx) = kto_default(fixed_idx);
            
            current_trace = ikto(param, hold_volt, volt, time_space, ek);
        case "ikslow1"
            % generate ikslow1
            current_trace = ikslow1(param_kslow1, hold_volt, volt, time_space, ek);
        case "ikslow2"
            % generate ikslow2
            num_param = 11;
            kslow2_default = [5334, 4912, 0.05766];
            shared_idx = 1:8;
            uniq_idx = setdiff(1:num_param, shared_idx);

            param = zeros(num_param, 1);
            param(shared_idx) = param_kslow1(shared_idx);
            param(uniq_idx) = kslow2_default;
            param(tune_idx1) = p(tune_idx2);

            current_trace = ikslow2(param, hold_volt, volt, time_space, ek);
        case "ikur"
            % generate ikur
            num_param = 11;
            kur_default = [270, 1050, 0];
            shared_idx = 1:8;
            uniq_idx = setdiff(1:num_param, shared_idx);

            param = zeros(num_param, 1);
            param(shared_idx) = param_kslow1(shared_idx);
            param(uniq_idx) = kur_default;
            param(tune_idx1) = p(tune_idx2);

            current_trace = ikur(param, hold_volt, volt, time_space, ek);
        case "ikss"
            % generate ikss
            num_param = 7;
            kss_default = [0.0862, 1235.5, 13.17, 0.0428];
            shared_idx1 = 1:3;
            shared_idx2 = [1, 3, 4];
            uniq_idx = setdiff(1:num_param, shared_idx1);

            param = zeros(num_param, 1);
            param(shared_idx1) = param_kslow1(shared_idx2);
            param(uniq_idx) = kss_default;
            param(tune_idx1) = p(tune_idx2);

            current_trace = ikss(param, hold_volt, volt, time_space, ek);
        case "ik1"
            % generate ik1
            num_param = 10;
            k1_default = [59.215, 5.476, 594.31, 4.753, ...
                1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];
            fixed_idx = setdiff(1:num_param, tune_idx1);

            param = zeros(num_param, 1);
            param(tune_idx1) = p(tune_idx2);
            param(fixed_idx) = k1_default(fixed_idx);

            current_trace = ik1(param, hold_volt, volt, time_space, ek);
    end
end
