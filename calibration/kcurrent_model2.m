function [yksum, comp_currents] = kcurrent_model2(p, model_struct, protocol_info)
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
        case 'ikto'
            % generate ikto
            num_param = 16;
            kto_default = [33.0, 40.0, 33.0, 15.5, 20.0, 7.7, 5.7, ... 
                0.03577, 0.06237, 7.0, 0.18064, 0.3956, 0.0226585, 0.00095, 0.051335, 0.14067];
            fixed_idx = setdiff(1:num_param, tune_idx1);
            
            param = zeros(num_param, 1);
            param(tune_idx1) = p(tune_idx2);
            param(fixed_idx) = kto_default(fixed_idx);
            
            current_trace = ikto(param, hold_volt, volt, time_space, ek);
        case 'ikslow1'
            % generate ikslow1
            current_trace = ikslow1(param_kslow1, hold_volt, volt, time_space, ek);
        case 'ikslow2'
            % generate ikslow2
            num_param = 11;
            kslow2_default = [5334, 4912, 0.05766];
            shared_idx = 1:8;
            uniq_idx = setdiff(1:num_param, shared_idx);

            param = zeros(num_param, 1);
            param(shared_idx) = param_kslow1(shared_idx);

            uniq_param = kslow2_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param(uniq_idx) = uniq_param;

            current_trace = ikslow2(param, hold_volt, volt, time_space, ek);
        case 'ikss'
            % generate ikss
            num_param = 7;
            kss_default = [0.0862, 1235.5, 13.17, 0.0428];
            shared_idx1 = 1:3;
            shared_idx2 = [1, 3, 4];
            uniq_idx = setdiff(1:num_param, shared_idx1);

            param = zeros(num_param, 1);
            param(shared_idx1) = param_kslow1(shared_idx2);

            uniq_param = kss_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param(uniq_idx) = uniq_param;

            current_trace = ikss(param, hold_volt, volt, time_space, ek);
        case 'ikur'
            % generate ikur
            num_param = 11;
            kur_default = [270, 1050, 0];
            shared_idx = 1:8;
            uniq_idx = setdiff(1:num_param, shared_idx);

            param = zeros(num_param, 1);
            param(shared_idx) = param_kslow1(shared_idx);

            uniq_param = kur_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param(uniq_idx) = uniq_param;

            current_trace = ikur(param, hold_volt, volt, time_space, ek);            
        case 'ik1'
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

function [current_trc] = ikto(p, hold_volt, volt, time_space, Ek)
    % constants & initial values
    gmax = p(16); % 0.14067

    act0 = 0.4139033547E-02;
    inact0 = 0.9999623535E+00;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding 
    gv_hold = ikto_gating_variables(p, hold_volt);

    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));      

    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(hold_volt-Ek);

    % current equation at pulse voltage
    gv_pulse = ikto_gating_variables(p, volt);

    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));

    current_trc((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-Ek);
end

function [gv] = ikto_gating_variables(p, V)
    % p0 = [33.0, 40.0, 33.0, 15.5, 20.0, 7.7, 5.7, ...
    %     0.03577, 0.06237, 7.0, 0.18064, 0.3956, 0.0226585, 0.00095, 0.051335, 0.14067]

    gv = zeros(4,1);
    
    % for Ikto
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(6))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(7))); % iss

    alpha1 = p(11).*exp(p(8).*(V+p(3)));
    beta1 = p(12).*exp(-p(9).*(V+p(3)));
    gv(3) = 1./(alpha1+beta1);
    
    alpha2_temp1 = exp((V+p(4))./(-1.0.*p(10)));
    alpha2_temp2 = exp((V+p(4)+p(5))./(-1.0.*p(10)));
    alpha2 = p(13).*(alpha2_temp1./alpha2_temp2);

    beta2_temp1 = p(14).*exp((V+p(4)+p(5))./p(10));
    beta2_temp2 = p(15).*exp((V+p(4)+p(5))./p(10));
    beta2 = beta2_temp1./(beta2_temp2+1.0);

    gv(4) = 1./(alpha2+beta2);
end

function [current_trc] = ikslow1(p, hold_volt, volt, time_space, Ek)
    % constants & initial values
    gmax = p(11); % 0.05766
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;
    
    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding 
    gv_hold = ikslow1_gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));        
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = ikslow1_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - Ek);
end

function [gv] = ikslow1_gating_variables(p, V)
    % gv(1:3) = gv(1:3) in Ikslow2
    % gv(1:3) = gv(1:3) in Ikur
    % gv(1) = gv(1) in Ikss
    % p0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];

    gv = zeros(4, 1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    gv(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    gv(4) = p(9)-p(10)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function [current_trc] = ikslow2(p, hold_volt, volt, time_space, Ek)
    % 11 parameters; {p(11): gmax} 
    % see 2020 Bondarenko

    % constants & initial values
    gmax = p(11); % 0.05766
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding
    gv_hold = ikslow2_gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = ikslow2_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - Ek);
end

function [gv] = ikslow2_gating_variables(p, V)
    % gv(1:3) = gv(1:3) in Ikslow1
    % gv(4) = p(9) - p(1)*[gv(2) in Ikslow1]
    % p0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
    % {p(9): p1, p(10): p2}

    gv = zeros(4,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    gv(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    gv(4) = p(9) - p(10)./(1.0+exp((p(2)+V)./p(5)));
end

function [current_trc] = ikss(p, hold_volt, volt, time_space, Ek)
    % 7 parameters; {p(7): gmax}

    % constants & initial values
    gmax = p(7); % 0.0428
    act0 = 0.5091689794e-03;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);
    
    current_trc = zeros(length(t), 1);

    % current equation at holding
    gv_hold = ikss_gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(2));
    current_trc(1:hold_idx) = gmax.*act_hold.*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = ikss_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(2));
    current_trc((hold_idx + 1):end) = gmax.*act_pulse.*(volt - Ek);
end

function [gv] = ikss_gating_variables(p, V)
    % gv(1) = gv(1) in Ikslow1
    % p0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428]

    gv = zeros(2,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    gv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function [current_trc] = ikur(p, hold_volt, volt, time_space, Ek)
    gmax = p(11); % 0
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);
    
    current_trc = zeros(length(t), 1);

    % current equation at holding
    gv_hold = ikur_gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = ikur_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - Ek);    
end

function [gv] = ikur_gating_variables(p, V)
    % gv(1:3) = gv(1:3) in Ikslow1
    % p0 = [270, 1050];
    % p0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 270, 1050, 0];

    gv = zeros(4,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    gv(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    gv(4) = p(9) + p(10)./(1.0+exp((p(2)+V)./p(5)));
end

function [current_trc] = ik1(p, hold_volt, volt, time_space, Ek)
    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    % pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);
    
    % current equation at holding
    [alpha_k1, beta_k1] = ik1_transition_rates(p, hold_volt, Ek);
    current_trc(hold_idx) = 0.27*sqrt(5400/5400) .* (alpha_k1/(alpha_k1 + beta_k1)) .* (hold_volt - Ek);

    [alpha_k1, beta_k1] = ik1_transition_rates(p, volt, Ek);
    current_trc((hold_idx+1):end) = 0.27*sqrt(5400/5400) .* (alpha_k1/(alpha_k1 + beta_k1)) .* (hold_volt - Ek);
end

function [alpha_k1, beta_k1] = ik1_transition_rates(p, v, Ek)
    % p0 = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143]
    alpha_k1 = p(5) ./ (1.0 + exp(p(6).*(v-Ek-p(1))));
    beta_k1 = (p(7)*exp(p(8)*(v-Ek+p(2))) + exp(p(9)*(v-Ek-p(3)))) ./ ...
        (1.0 + exp(-p(10)*(v-Ek+p(4))));
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
