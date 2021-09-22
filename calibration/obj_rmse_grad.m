function [f, g] = obj_rmse_grad(p, model_struct, volt_space, time_space, yksum)
    global param_kto
    global param_kslow1
    global param_kslow2
    global param_kss
    global param_kur
    global param_k1

    kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
        0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
    kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, ...
        2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
    kslow2_default = [5334, 4912, 0.05766];
    kss_default = [0.0862, 1235.5, 13.17, 0.0428];
    kur_default = [270, 1050, 0];
    k1_default = [59.215, 5.476, 594.31, 4.753, ...
        1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

    % calibration parameters
    for i = 1:length(model_struct)
        % current model info
        current_name = model_struct(i).name;
        tune_idx1 = model_struct(i).idx1;
        tune_idx2 = model_struct(i).idx2;

        switch current_name
        case 'ikto'
            % generate ikto
            num_param = 17;
            param_kto = NaN(num_param, 1);

            fixed_idx = setdiff(1:num_param, tune_idx1);
            param_kto(tune_idx1) = p(tune_idx2);
            param_kto(fixed_idx) = kto_default(fixed_idx);
        case 'ikslow1'
            % generate ikslow1
            num_kslow1_param = 13;
            param_kslow1 = NaN(num_kslow1_param, 1);
            
            fixed_idx = setdiff(1:num_kslow1_param, tune_idx1);
            param_kslow1(tune_idx1) = p(tune_idx2);
            param_kslow1(fixed_idx) = kslow1_default(fixed_idx);
        case 'ikslow2'
            % generate ikslow2
            num_param = 11;
            param_kslow2 = NaN(num_param, 1);
            
            shared_idx = 1:8;
            param_kslow2(shared_idx) = param_kslow1(shared_idx);
            
            uniq_idx = setdiff(1:num_param, shared_idx);
            uniq_param = kslow2_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param_kslow2(uniq_idx) = uniq_param;
        case 'ikss'
            % generate ikss
            num_param = 7;
            param_kss = NaN(num_param, 1);
            
            shared_idx1 = 1:3;
            shared_idx2 = [1, 3, 4];
            param_kss(shared_idx1) = param_kslow1(shared_idx2);
            
            uniq_idx = setdiff(1:num_param, shared_idx1);
            uniq_param = kss_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param_kss(uniq_idx) = uniq_param;
        case 'ikur'
            % generate ikur
            num_param = 11;
            param_kur = NaN(num_param, 1);
            
            shared_idx = 1:8;
            param_kur(shared_idx) = param_kslow1(shared_idx);
            
            uniq_idx = setdiff(1:num_param, shared_idx);
            uniq_param = kur_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param_kur(uniq_idx) = uniq_param;
        case 'ik1'
            % generate ik1
            num_param = 10;
            param_k1 = NaN(num_param, 1);

            fixed_idx = setdiff(1:num_param, tune_idx1);
            param_k1(tune_idx1) = p(tune_idx2);
            param_k1(fixed_idx) = k1_default(fixed_idx);
        end
    end

    hold_idx = time_space{4};    
    volts = volt_space{2};
    num_volts = length(volts);

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{3} = time_space;
    protocol{4} = volt_space{3};

    rmse_list = NaN(num_volts, 1);
    grad_list = NaN(model_struct(end).idx2(end), num_volts);
    for i = 1:num_volts
        yksum_i = yksum(:, i);
        protocol{2} = volts(i);

        [yksum_hat, state_vars_list, trans_rates_list] = kcurrent_model(model_struct, protocol);
        err = yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end);
        mserr = mean(err.^2);

        rmse_list(i) = sqrt(mserr);
        grad_list(:, i) = kcurrent_grad(model_struct, protocol, mserr, err, state_vars_list, trans_rates_list);
    end
    f = sum(rmse_list);
end

function [yksum_hat, state_vars_list, trans_rates_list] = kcurrent_model(model_struct, protocol_info)
    % protocol info
    hold_volt = protocol_info{1};
    volt = protocol_info{2};
    time_space = protocol_info{3};
    ek = protocol_info{4};
    
    num_currents = length(model_struct);
    yksum_hat = zeros(length(time_space{1}), 1);
    state_vars_list = cell(num_currents, 1);
    trans_rates_list = cell(num_currents, 1);
    
    % current model info
    for i = 1:num_currents
        % current model info
        switch model_struct(i).name
        case 'ikto'
            % generate ikto
            [current_trace, state_vars, trans_rates] = ikto(hold_volt, volt, time_space, ek);
        case 'ikslow1'
            % generate ikslow1
            [current_trace, state_vars, trans_rates] = ikslow1(hold_volt, volt, time_space, ek);
        case 'ikslow2'
            % generate ikslow2
            [current_trace, state_vars, trans_rates] = ikslow2(hold_volt, volt, time_space, ek);
        case 'ikss'
            % generate ikss
            [current_trace, state_vars, trans_rates] = ikss(hold_volt, volt, time_space, ek);
        case 'ikur'
            % generate ikur
            [current_trace, state_vars, trans_rates] = ikur(hold_volt, volt, time_space, ek);            
        case 'ik1'
            % generate ik1
            [current_trace, state_vars, trans_rates] = ik1(hold_volt, volt, time_space, ek);
        end
        yksum_hat = yksum_hat + current_trace;
        state_vars_list{i} = state_vars;
        trans_rates_list{i} = trans_rates;
    end
end

function [current_trc, state_vars, trans_rates] = ikto(hold_volt, volt, time_space, Ek)
    global param_kto
    p = param_kto;

    % constants & initial values
    f_eacv = p(15); % 0.2087704319 Ikto fraction phophorylated
    gmax = p(16); % 0.14067
    gmaxp = p(16)*(1.0-p(17)); % 0.14067*(1.0-0.387)

    act0 = 0.4139033547E-02;
    inact0 = 0.9999623535E+00;
    actp0 = 0.8637307739E-03;
    inactp0 = 0.9999881755E+00;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % transition rates
    tr_hold = ikto_trans_rates(hold_volt);
    tr_pulse = ikto_trans_rates(volt);

    % gating variables
    ass_hold = tr_hold(1)./(tr_hold(1)+tr_hold(2));
    atau_hold = 1./(tr_hold(1)+tr_hold(2));
    iss_hold = tr_hold(3)./(tr_hold(3)+tr_hold(4));
    itau_hold = 1./(tr_hold(3)+tr_hold(4));

    assp_hold = tr_hold(5)./(tr_hold(5)+tr_hold(6)); 
    ataup_hold = 1./(tr_hold(6)+tr_hold(6));
    issp_hold = tr_hold(7)./(tr_hold(7)+tr_hold(8)); 
    itaup_hold = 1./(tr_hold(7)+tr_hold(8)); 

    ass_pulse = tr_pulse(1)./(tr_pulse(1)+tr_pulse(2));
    atau_pulse = 1./(tr_pulse(1)+tr_pulse(2));
    iss_pulse = tr_pulse(3)./(tr_pulse(3)+tr_pulse(4));
    itau_pulse = 1./(tr_pulse(3)+tr_pulse(4));

    assp_pulse = tr_pulse(5)./(tr_pulse(5)+tr_pulse(6)); 
    ataup_pulse = 1./(tr_pulse(6)+tr_pulse(6));
    issp_pulse = tr_pulse(7)./(tr_pulse(7)+tr_pulse(8)); 
    itaup_pulse = 1./(tr_pulse(7)+tr_pulse(8)); 

    % state variables
    act_hold = hh_model(hold_t, act0, ass_hold, atau_hold);
    inact_hold = hh_model(hold_t, inact0, iss_hold, itau_hold);      
    
    actp_hold = hh_model(hold_t, actp0, assp_hold, ataup_hold);
    inactp_hold = hh_model(hold_t, inactp0, issp_hold, itaup_hold); 
    
    act_pulse = hh_model(pulse_t, act0, ass_pulse, atau_pulse);
    inact_pulse = hh_model(pulse_t, inact0, iss_pulse, itau_pulse);
    
    actp_pulse = hh_model(pulse_t, actp0, assp_pulse, ataup_pulse);
    inactp_pulse = hh_model(pulse_t, inactp0, issp_pulse, itaup_pulse);
    
    % current equation
    current_trc(1:hold_idx) = (1-f_eacv).*gmax.*(act_hold.^3).*(inact_hold).*(hold_volt-Ek) + ...
        f_eacv.*gmaxp.*(actp_hold.^3).*(inactp_hold).*(hold_volt-Ek);
    current_trc((hold_idx+1):end) = (1-f_eacv).*gmax.*(act_pulse.^3).*(inact_pulse).*(volt-Ek) + ...
        f_eacv.*gmaxp.*(actp_pulse.^3).*(inactp_pulse).*(volt-Ek);

    state_vars = [act_pulse, inact_pulse, actp_pulse, inactp_pulse];
    trans_rates = tr_pulse;
end

function tr = ikto_trans_rates(V)
    global param_kto
    p = param_kto;

    tr = NaN(8,1);
    
    %
    % nonphosphorylated
    %
    % activation
    alpha1 = p(9).*exp(p(7).*(V+p(1)));
    beta1 = p(10).*exp(-p(8).*(V+p(1)));
    
    % inactivation
    alpha2_temp1 = p(11).*exp((V+p(2))./(-1.0.*p(6)));
    alpha2_temp2 = p(12).*exp((V+p(2)+p(3))./(-1.0.*p(6)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(13).*exp((V+p(2)+p(3))./p(6));
    beta2_temp2 = p(14).*exp((V+p(2)+p(3))./p(6));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    tr(1) = alpha1;
    tr(2) = beta1;
    tr(3) = alpha2;
    tr(4) = beta2;

    %
    % phosphorylated
    %
    % activation
    alpha1p = p(9).*exp(p(7).*(V+p(1)-p(4)));
    beta1p = p(10).*exp(-p(8).*(V+p(1)-p(4)));
    
    % inactivation
    alpha2_temp1p = p(11).*exp((V+p(2)-p(5))./(-1.0.*p(6)));
    alpha2_temp2p = p(12).*exp((V+p(2)+p(3)-p(5))./(-1.0.*p(6)));
    alpha2p = alpha2_temp1p./(1.0+alpha2_temp2p);
    
    beta2_temp1p = p(13).*exp((V+p(2)+p(3)-p(5))./p(6));
    beta2_temp2p = p(14).*exp((V+p(2)+p(3)-p(5))./p(6));
    beta2p = beta2_temp1p./(1.0+beta2_temp2p);

    tr(5) = alpha1p;
    tr(6) = beta1p;
    tr(7) = alpha2p;
    tr(8) = beta2p;
end

function [current_trc, state_vars, trans_rates] = ikslow1(hold_volt, volt, time_space, Ek)
    global param_kslow1
    p = param_kslow1;
    
    % constants & initial values
    f_eacv = p(11); % 0.9214774521 Ikslow1 fraction of nonphospholatedl
    gmax = p(12); % 0.05766
    gmaxp = p(13); % 0.07496
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;
    
    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding 
    tr_hold = ikslow1_trans_rates(hold_volt);
    act_hold = hh_model(hold_t, act0, tr_hold(1), tr_hold(3));
    inact_hold = hh_model(hold_t, inact0, tr_hold(2), tr_hold(4));        
    current_trc(1:hold_idx) = (gmax*f_eacv + gmaxp*(1-f_eacv)).*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    tr_pulse = ikslow1_trans_rates(volt);
    act_pulse = hh_model(pulse_t, act0, tr_pulse(1), tr_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, tr_pulse(2), tr_pulse(4));
    current_trc((hold_idx + 1):end) = (gmax*f_eacv + gmaxp*(1-f_eacv)).*(act_pulse).*(inact_pulse).*(volt - Ek);

    state_vars = [act_pulse, inact_pulse];
    trans_rates = tr_pulse;
end

function tr = ikslow1_trans_rates(V)
    global param_kslow1
    p = param_kslow1;

    tr = zeros(4, 1);
    tr(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    tr(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    tr(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    tr(4) = p(9)-p(10)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function [current_trc, state_vars, trans_rates] = ikslow2(hold_volt, volt, time_space, Ek)
    global param_kslow2
    p = param_kslow2;
    
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
    tr_hold = ikslow2_trans_rates(hold_volt);
    act_hold = hh_model(hold_t, act0, tr_hold(1), tr_hold(3));
    inact_hold = hh_model(hold_t, inact0, tr_hold(2), tr_hold(4));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    tr_pulse = ikslow2_trans_rates(volt);
    act_pulse = hh_model(pulse_t, act0, tr_pulse(1), tr_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, tr_pulse(2), tr_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - Ek);
    
    state_vars = [act_pulse, inact_pulse];
    trans_rates = tr_pulse;
end

function tr = ikslow2_trans_rates(V)
    global param_kslow2
    p = param_kslow2;

    tr = zeros(4,1);
    tr(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    tr(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    tr(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    tr(4) = p(9) - p(10)./(1.0+exp((p(2)+V)./p(5)));
end

function [current_trc, state_vars, trans_rates] = ikss(hold_volt, volt, time_space, Ek)
    global param_kss
    p = param_kss;

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
    tr_hold = ikss_trans_rates(hold_volt);
    act_hold = hh_model(hold_t, act0, tr_hold(1), tr_hold(2));
    current_trc(1:hold_idx) = gmax.*act_hold.*(hold_volt - Ek);

    % current equation at pulse voltage
    tr_pulse = ikss_trans_rates(volt);
    act_pulse = hh_model(pulse_t, act0, tr_pulse(1), tr_pulse(2));
    current_trc((hold_idx + 1):end) = gmax.*act_pulse.*(volt - Ek);

    state_vars = act_pulse;
    trans_rates = tr_pulse;
end

function tr = ikss_trans_rates(V)
    global param_kss;
    p = param_kss;

    tr = zeros(2,1);
    tr(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    tr(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function [current_trc, state_vars, trans_rates] = ikur(hold_volt, volt, time_space, Ek)
    global param_kur
    p = param_kur;
    
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
    gv_hold = ikur_trans_rates(hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = ikur_trans_rates(volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - Ek); 
    
    state_vars = [];
    trans_rates = [];
end

function tr = ikur_trans_rates(V)
    global param_kur
    p = param_kur;

    tr = zeros(4,1);
    tr(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    tr(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    tr(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    tr(4) = p(9) + p(10)./(1.0+exp((p(2)+V)./p(5)));
end

function [current_trc, state_vars, trans_rates] = ik1(hold_volt, volt, time_space, Ek)
    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    % pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);
    
    % current equation at holding
    [alpha_k1, beta_k1] = ik1_transition_rates(hold_volt, Ek);
    current_trc(hold_idx) = 0.27*sqrt(5400/5400) .* (alpha_k1/(alpha_k1 + beta_k1)) .* (hold_volt - Ek);

    [alpha_k1, beta_k1] = ik1_transition_rates(volt, Ek);
    current_trc((hold_idx+1):end) = 0.27*sqrt(5400/5400) .* (alpha_k1/(alpha_k1 + beta_k1)) .* (hold_volt - Ek);

    state_vars = [];
    trans_rates = [];
end

function [alpha_k1, beta_k1] = ik1_transition_rates(v, Ek)
    global param_k1
    p = param_k1;
    
    % p0 = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143]
    alpha_k1 = p(5) ./ (1.0 + exp(p(6).*(v-Ek-p(1))));
    beta_k1 = (p(7)*exp(p(8)*(v-Ek+p(2))) + exp(p(9)*(v-Ek-p(3)))) ./ ...
        (1.0 + exp(-p(10)*(v-Ek+p(4))));
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end

function g = kcurrent_grad(model_struct, protocol, mserr, err, state_vars_list, trans_rates_list)
    dfdm = 1/(2*sqrt(mserr));

    g = NaN(model_struct(end).idx2(end), 1);
    for i = 1:length(model_struct)
        current_name = model_struct(i).name;
        tune_idx2 = model_struct(i).idx2;
        state_vars = state_vars_list{i};
        trans_rates = trans_rates_list{i};

        switch current_name
        case 'ikto'
            g(tune_idx2) = grad_ikto(protocol, dfdm, err, state_vars, trans_rates);
        case 'ikslow1'
            g(tune_idx2) = grad_ikslow1(protocol, dfdm, err, state_vars, trans_rates);
        case 'ikslow2'
            g(tune_idx2) = grad_ikslow2(protocol, dfdm, err, state_vars, trans_rates);
        case 'ikss'
            g(tune_idx2) = grad_ikss(protocol, dfdm, err, state_vars, trans_rates);
        otherwise
            g(tune_idx2) = [];
        end
    end
end

function g = grad_ikto(protocol, dfdm, err, state_vars, trans_rates)
    % tune_idx1_kto = [1, 2, 6, 10, 13, 14, 15, 16, 17];
    
    global param_kto
    p = param_kto;
    g = NaN(9, 1);

    act0 = 0.4139033547E-02;
    inact0 = 0.9999623535E+00;
    actp0 = 0.8637307739E-03;
    inactp0 = 0.9999881755E+00;

    f_eacv = param_kto(15);
    gmax = param_kto(16);
    gmaxp = param_kto(16)*(1.0-param_kto(17));
    
    time_space = protocol{3};
    t = time_space{3};
    n = length(t);
    volt = protocol{2};
    ek = protocol{4};

    alpha1 = trans_rates(1);
    beta1 = trans_rates(2);
    alpha2 = trans_rates(3);
    beta2 = trans_rates(4);
    
    alpha1p = trans_rates(5);
    beta1p = trans_rates(6);
    alpha2p = trans_rates(7);
    beta2p = trans_rates(8);

    ass = alpha1/(alpha1+beta1);
    atau = 1/(alpha1+beta1);
    iss = alpha2/(alpha2+beta2);
    itau = 1/(alpha2+beta2);

    assp = alpha1p/(alpha1p+beta1p);
    ataup = 1/(alpha1p+beta1p);
    issp = alpha2p/(alpha2p+beta2p);
    itaup = 1/(alpha2p+beta2p);

    % derivatives
    dmde = 2*err;
    dedy = -1;
    dfdy = dfdm*dmde*dedy;
    
    dyda = 3*gmax*(state_vars(:, 1).^2).*(state_vars(:, 2))*(1-f_eacv)*(volt-ek);
    dydi = gmax*state_vars(:, 1).^3*(1-f_eacv)*(volt-ek);
    dydap = 3*gmaxp*(state_vars(:, 3).^2).*(state_vars(:, 4))*f_eacv*(volt-ek);
    dydip = gmaxp*state_vars(:, 4).^3*f_eacv*(volt-ek);
    
    dadss = 1 - exp(-t./atau);
    dadtau = -(ass - act0)*(t.*exp(-t./atau))/(atau^2);
    didss = 1- exp(-t./itau);
    didtau = -(iss - inact0)*(t.*exp(-t./itau))/(itau^2);

    dapdss = 1 - exp(-t./ataup);
    dapdtau = -(assp - actp0)*(t.*exp(-t./ataup))/(ataup^2);
    dipdss = 1- exp(-t./itaup);
    dipdtau = -(issp - inactp0)*(t.*exp(-t./itaup))/(itaup^2);
    
    dssdalpha1 = beta1/(alpha1+beta1)^2;
    dssdbeta1 = -alpha1/(alpha1+beta1)^2;
    dtaudalpha1 = -1/(alpha1+beta1)^2;
    dtaudbeta1 = -1/(alpha1+beta1)^2;
    dssdalpha2 = beta2/(alpha2+beta2)^2;
    dssdbeta2 = -alpha2/(alpha2+beta2)^2;
    dtaudalpha2 = -1/(alpha2+beta2)^2;
    dtaudbeta2 = -1/(alpha2+beta2)^2;

    dssdalpha1p = beta1p/(alpha1p+beta1p)^2;
    dssdbeta1p = -alpha1p/(alpha1p+beta1p)^2;
    dtaudalpha1p = -1/(alpha1p+beta1p)^2;
    dtaudbeta1p = -1/(alpha1p+beta1p)^2;
    dssdalpha2p = beta2p/(alpha2p+beta2p)^2;
    dssdbeta2p = -alpha2p/(alpha2p+beta2p)^2;
    dtaudalpha2p = -1/(alpha2p+beta2p)^2;
    dtaudbeta2p = -1/(alpha2p+beta2p)^2;
    
    % p1
    dalpha1dp1 = p(7)*p(9)*exp(p(7)*(volt+p(1)));
    dbeta1dp1 = -p(8)*p(10)*exp(-p(8)*(volt+p(1)));
    dalpha1pdp1 = p(7)*p(9)*exp(p(7)*(volt+p(1)-p(4)));
    dbeta1pdp1 = -p(8)*p(10)*exp(-p(8)*(volt+p(1)-p(4)));
    g(1) = (1/n)*sum(dfdy.*(dyda.*dadss*dssdalpha1*dalpha1dp1 + ...
        dyda.*dadss*dssdbeta1*dbeta1dp1 + ...
        dyda.*dadtau*dtaudalpha1*dalpha1dp1 + ...
        dyda.*dadtau*dtaudbeta1*dbeta1dp1 + ...
        dydap.*dapdss*dssdalpha1p*dalpha1pdp1 + ...
        dydap.*dapdss*dssdbeta1p*dbeta1pdp1 + ...
        dydap.*dapdtau*dtaudalpha1p*dalpha1pdp1 + ...
        dydap.*dapdtau*dtaudbeta1p*dbeta1pdp1));

    % p2
    dalpha2dp2_temp1 = p(11)*exp(2*(volt+p(2)+p(3))/p(6)-(volt+p(2))/p(6));
    dalpha2dp2_temp2 = p(6)*(exp((volt+p(2)+p(3))/p(6))+p(12))^2;
    dalpha2dp2 = -dalpha2dp2_temp1/dalpha2dp2_temp2;
    dbeta2dp2_temp1 = p(13)*exp((volt+p(2)+p(3))/p(6));
    dbeta2dp2_temp2 = p(6)*(p(14)*exp((volt+p(2)+p(3))/p(6))+1.0)^2;
    dbeta2dp2 = dbeta2dp2_temp1/dbeta2dp2_temp2;
    dalpha2pdp2_temp1 = p(11)*exp(2*(volt+p(2)+p(3)-p(5))/p(6)-(volt+p(2)-p(5))/p(6));
    dalpha2pdp2_temp2 = p(6)*(exp((volt+p(2)+p(3)-p(5))/p(6))+p(12))^2;
    dalpha2pdp2 = - dalpha2pdp2_temp1/dalpha2pdp2_temp2;
    dbeta2pdp2_temp1 = p(13)*exp((volt+p(2)+p(3)-p(5))/p(6));
    dbeta2pdp2_temp2 = p(6)*(p(14)*exp((volt+p(2)+p(3)-p(5))/p(6))+1.0)^2;
    dbeta2pdp2 = dbeta2pdp2_temp1/dbeta2pdp2_temp2;
    g(2) = (1/n)*sum(dfdy.*(dydi.*didss*dssdalpha2*dalpha2dp2 + ...
        dydi.*didss*dssdbeta2*dbeta2dp2 + ...
        dydi.*didtau*dtaudalpha2*dalpha2dp2 + ...
        dydi.*didtau*dtaudbeta2*dbeta2dp2 + ...
        dydip.*dipdss*dssdalpha2p*dalpha2pdp2 + ...
        dydip.*dipdss*dssdbeta2p*dbeta2pdp2 + ...
        dydip.*dipdtau*dtaudalpha2p*dalpha2pdp2 + ...
        dydip.*dipdtau*dtaudbeta2p*dbeta2pdp2));

    % p6
    dalpha2dp6_temp1 = p(11)*((volt+p(2))*exp((volt+p(2)+p(3))/p(6))-p(12)*p(3))*exp((volt+p(2)+p(3))/p(6)-(volt+p(2))/p(6));
    dalpha2dp6_temp2 = p(6)^2*(exp((volt+p(2)+p(3))/p(6))+p(12))^2;
    dalpha2dp6 = dalpha2dp6_temp1/dalpha2dp6_temp2;
    dbeta2dp6_temp1 = p(13)*(volt+p(2)+p(3))*exp((volt+p(2)+p(3))/p(6));
    dbeta2dp6_temp2 = p(6)^2*(p(14)*exp((volt+p(2)+p(3))/p(6))+1.0)^2;
    dbeta2dp6 = -dbeta2dp6_temp1/dbeta2dp6_temp2;
    dalpha2pdp6_temp1 = p(11)*((volt+p(2)-p(5))*exp((volt+p(2)+p(3)-p(5))/p(6))-p(12)*p(3))*exp((volt+p(2)+p(3)-p(5))/p(6)-(volt+p(2)-p(5))/p(6));
    dalpha2pdp6_temp2 = p(6)^2*(exp((volt+p(2)+p(3)-p(5))/p(6))+p(12))^2;
    dalpha2pdp6 = dalpha2pdp6_temp1/dalpha2pdp6_temp2;
    dbeta2pdp6_temp1 = p(13)*(volt+p(2)+p(3)-p(5))*exp((volt+p(2)+p(3)-p(5))/p(6));
    dbeta2pdp6_temp2 = p(6)^2*(p(14)*exp((volt+p(2)+p(3)-p(5))/p(6))+1.0)^2;
    dbeta2pdp6 = -dbeta2pdp6_temp1/dbeta2pdp6_temp2;
    g(3) = (1/n)*sum(dfdy.*(dydi.*didss*dssdalpha2*dalpha2dp6 + ...
        dydi.*didss*dssdbeta2*dbeta2dp6 + ...
        dydi.*didtau*dtaudalpha2*dalpha2dp6 + ...
        dydi.*didtau*dtaudbeta2*dbeta2dp6 + ...
        dydip.*dipdss*dssdalpha2p*dalpha2pdp6 + ...
        dydip.*dipdss*dssdbeta2p*dbeta2pdp6 + ...
        dydip.*dipdtau*dtaudalpha2p*dalpha2pdp6 + ...
        dydip.*dipdtau*dtaudbeta2p*dbeta2pdp6));
    
    % p10
    dbeta1dp10 = exp(-p(8)*(volt+p(1)));
    dbeta1pdp10 = exp(-p(8)*(volt+p(1)-p(4)));
    g(4) = (1/n)*sum(dfdy.*(dyda.*dadss*dssdbeta1*dbeta1dp10 + ...
        dyda.*dadtau*dtaudbeta1*dbeta1dp10 + ...
        dydap.*dapdss*dssdbeta1p*dbeta1pdp10 + ...
        dydap.*dapdtau*dtaudbeta1p*dbeta1pdp10));

    % p13
    dbeta2dp13 = (exp((volt+p(2)+p(3))/p(6)))/(p(14)*exp((volt+p(2)+p(3))/p(6))+1.0);
    dbeta2pdp13 = (exp((volt+p(2)+p(3)-p(5))/p(6)))/(p(14)*exp((volt+p(2)+p(3)-p(5))/p(6))+1.0);
    g(5) = (1/n)*sum(dfdy.*(dydi.*didss*dssdbeta2*dbeta2dp13 + ...
        dydi.*didtau*dtaudbeta2*dbeta2dp13 + ...
        dydip.*dipdss*dssdbeta2p*dbeta2pdp13 + ...
        dydip.*dipdtau*dtaudbeta2p*dbeta2pdp13));

    % p14
    dbeta2dp14_temp1 = p(13)*exp(2*(volt+p(2)+p(3))/p(6));
    dbeta2dp14_temp2 = (p(14)*exp((volt+p(2)+p(3))/p(6))+1.0)^2;
    dbeta2dp14 = -dbeta2dp14_temp1/dbeta2dp14_temp2;
    dbeta2pdp14_temp1 = p(13)*exp(2*(volt+p(2)+p(3)-p(5))/p(6));
    dbeta2pdp14_temp2 = (p(14)*exp((volt+p(2)+p(3)-p(5))/p(6))+1.0)^2;
    dbeta2pdp14 = -dbeta2pdp14_temp1/dbeta2pdp14_temp2;
    g(6) = (1/n)*sum(dfdy.*(dydi.*didss*dssdbeta2*dbeta2dp14 + ...
        dydi.*didtau*dtaudbeta2*dbeta2dp14 + ...
        dydip.*dipdss*dssdbeta2p*dbeta2pdp14 + ...
        dydip.*dipdtau*dtaudbeta2p*dbeta2pdp14));

    % p15
    k1_feacv = gmax*(state_vars(:, 1).^3).*(state_vars(:, 2));
    k2_feacv = gmaxp*(state_vars(:, 3).^3).*(state_vars(:, 4));
    dydfeacv = (-k1_feacv+k2_feacv)*(volt-ek);
    g(7) = (1/n)*sum(dfdy.*dydfeacv);

    % p16
    k1_gmax = (1-f_eacv)*(state_vars(:, 1).^3).*(state_vars(:, 2));
    k2_gmax = f_eacv*(1-p(17))*(state_vars(:, 3).^3).*(state_vars(:, 4));
    dydgmax = (k1_gmax+k2_gmax)*(volt-ek);
    g(8) = (1/n)*sum(dfdy.*dydgmax);

    % p17
    dydp17 = -f_eacv*gmax*(state_vars(:, 3).^3).*(state_vars(:, 4))*(volt-ek);
    g(9) = (1/n)*sum(dfdy.*dydp17);
end

function g = grad_ikslow1(protocol, dfdm, err, state_vars, trans_rates)
    % tune_idx1_kslow1 = [1, 2, 3, 4, 5, 8, 9, 11, 12, 13];
    
    global param_kslow1
    p = param_kslow1;
    g = NaN(10, 1);

    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;
    
    f_eacv = p(11);
    gmax = p(12);
    gmaxp = p(13);

    time_space = protocol{3};
    t = time_space{3};
    n = length(t);
    volt = protocol{2};
    ek = protocol{4};
    
    ass = trans_rates(1);
    iss = trans_rates(2);
    atau = trans_rates(3);
    itau = trans_rates(4);

    % derivatives
    dmde = 2*err;
    dedy = -1;
    dfdy = dfdm*dmde*dedy;

    dyda = (gmax*f_eacv+gmaxp*(1-f_eacv))*state_vars(:, 2)*(volt-ek);
    dydi = (gmax*f_eacv+gmaxp*(1-f_eacv))*state_vars(:, 1)*(volt-ek);

    dadass = 1.0-exp(-t./atau);
    dadatau = -(ass - act0)*(t.*exp(-t./atau))/(atau^2);
    didiss = 1.0-exp(-t./itau);
    diditau = -(iss - inact0)*(t.*exp(-t./itau))/(itau^2);

    % p1
    dassdp1 = (exp((p(1)+volt)/p(4)))/(p(4)*(exp((p(1)+volt)/p(4))+1.0)^2);
    g(1) = (1/n)*sum(dfdy.*(dyda.*dadass*dassdp1));

    % p2
    dissdp2 = -exp((p(2)+volt)/p(5))/(p(5)*(exp((p(2)+volt)/p(5))+1.0)^2);
    ditaudp2 = (p(10)*exp((p(2)+volt)/p(5)))/(p(5)*(exp((p(2)+volt)/p(5))+1.0)^2);
    g(2) = (1/n)*sum(dfdy.*dydi.*(didiss*dissdp2 + diditau*ditaudp2));

    % p3
    dataudp3_temp1 = p(6)*(p(7)*exp(p(7)*(volt+p(3)))-p(7)*exp(-p(7)*(volt+p(3))));
    dataudp3_temp2 = (exp(p(7)*(volt+p(3))) + exp(-p(7)*(volt+p(3))))^2;
    dataudp3 = -dataudp3_temp1/dataudp3_temp2;
    g(3) = (1/n)*sum(dfdy.*dyda.*dadatau*dataudp3);

    % p4
    dassdp4 = -((p(1)+volt)*exp((p(1)+volt)/p(4)))/(p(4)^2*(exp((p(1)+volt)/p(4))+1.0)^2);
    g(4) = (1/n)*sum(dfdy.*dyda.*dadass*dassdp4);

    % p5
    dissdp5 = ((p(2)+volt)*exp((p(2)+volt)/p(5)))/(p(5)^2*(exp((p(2)+volt)/p(5))+1.0)^2);
    ditaudp5 = -(p(10)*(p(2)+volt)*exp((p(2)+volt)/p(5)))/(p(5)^2*(exp((p(2)+volt)/p(5))+1)^2);
    g(5) = (1/n)*sum(dfdy.*dydi.*(didiss*dissdp5 + diditau*ditaudp5));

    % p8
    dataudp8 = 1;
    g(6) = (1/n)*sum(dfdy.*dyda.*dadatau*dataudp8);

    % p9
    ditaudp9 = 1;
    g(7) = (1/n)*sum(dfdy.*dydi.*diditau*ditaudp9);

    % p11
    dydfeacv = (gmax-gmaxp)*state_vars(:, 1).*state_vars(:, 2)*(volt-ek);
    g(8) = (1/n)*sum(dfdy.*dydfeacv);

    % p12
    dydgmax = f_eacv.*state_vars(:, 1).*state_vars(:, 2)*(volt-ek);
    g(9) = (1/n)*sum(dfdy.*dydgmax);

    % p13
    dydgmaxp = (1-f_eacv)*state_vars(:, 1).*state_vars(:, 2)*(volt-ek);
    g(10) = (1/n)*sum(dfdy.*dydgmaxp);
end

function g = grad_ikslow2(p, protocol, dfdm, err, state_vars, trans_rates)
    global param_kslow2
    g = [];
end
function g = grad_ikss(p, protocol, dfdm, err, state_vars, trans_rates)
    global param_kss
    g = [];
end
