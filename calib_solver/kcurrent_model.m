function [yksum, comp_currents] = kcurrent_model(p, model_struct, protocol_info)
    current_names = model_struct.name;
    num_currents = length(current_names);
    
    % declare shared parameters of ikslow1 as global variable
    matching_idx = strcmp(current_names, 'ikslow1');
    if any(matching_idx)
        global param_kslow1
        kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, ...
            2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
        num_kslow1_param = length(kslow1_default);
        
        tune_idx1 = model_struct(matching_idx).idx1;
        tune_idx2 = model_struct(matching_idx).idx2;
        fixed_idx = setdiff(1:num_kslow1_param, tune_idx1);

        param_kslow1 = zeros(num_kslow1_param, 1);
        param_kslow1(tune_idx1) = p(tune_idx2);
        param_kslow1(fixed_idx) = kslow1_default(fixed_idx);
    end

    comp_currents = cell(num_currents, 1);
    for i = 1:num_currents
        comp_currents{i} = gen_matching_current(p, model_struct(i), protocol_info);
    end

    yksum = sum(comp_currents, 2);
end

function [current_trace] = gen_matching_current(p, model_info, protocol_info)
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
            kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
                0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
            num_param = length(kto_default);
            fixed_idx = setdiff(1:num_param, tune_idx1);
            
            param = zeros(num_param, 1);
            param(tune_idx1) = p(tune_idx2);
            param(fixed_idx) = kto_default(fixed_idx);
            
            current_trace = ikto(param, hold_volt, volt, time_space, ek);
        case "ikslow1"
            % generate ikslow1
            global param_kslow1
            num_kslow1_param = 13;

        case "ikslow2"
            % generate ikslow2
            global param_kslow1
            num_kslow2_param = 11;

        case "ikur"
            % generate ikur
            global param_kslow1
            num_kur_param = 11;

        case "ikss"
            % generate ikss
            global param_kslow1
            num_kss_param = 11;

        case "ik1"
            % generate ik1
            num_k1_param = 4;

    end
end
