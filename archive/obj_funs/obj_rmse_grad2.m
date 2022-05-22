function [f, g] = obj_rmse_grad2(p, model_struct, volt_space, time_space, yksum)
    global param_kto
    global param_kslow1
    global param_kslow2
    global param_kss
    global param_kur
    global param_k1

    kto_default = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.3846];
    kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, 6.1, 18.0, 2.058, 803.0, 0.16];
    kslow2_default = [4912, 5334, 0.16];
    kss_default = [0.0862, 1235.5, 13.17, 0.0611];
    kur_default = [270, 1050, 0];
    k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

    num_currents = length(model_struct);        
    current_names = cell(num_currents, 1);
    for i = 1:num_currents
        current_names{i} = model_struct(i).name;
    end

    % declare shared parameters of ikslow1 as global variable
    matching_idx = strcmp(current_names, 'ikslow1');
    if any(matching_idx)    
        num_kslow1_param = 11;
        
        tune_idx1 = model_struct(matching_idx).idx1;
        tune_idx2 = model_struct(matching_idx).idx2;
        fixed_idx = setdiff(1:num_kslow1_param, tune_idx1);

        param_kslow1 = zeros(num_kslow1_param, 1);
        param_kslow1(tune_idx1) = p(tune_idx2);
        param_kslow1(fixed_idx) = kslow1_default(fixed_idx);
    end

    % calibration parameters
    for i = 1:num_currents
        % current model info
        current_name = model_struct(i).name;
        tune_idx1 = model_struct(i).idx1;
        tune_idx2 = model_struct(i).idx2;

        switch current_name
        case 'ikto'
            num_param = 13;
            param_kto = NaN(num_param, 1);

            fixed_idx = setdiff(1:num_param, tune_idx1);
            param_kto(tune_idx1) = p(tune_idx2);
            param_kto(fixed_idx) = kto_default(fixed_idx);
        case 'ikslow2'
            num_param = 11;
            param_kslow2 = NaN(num_param, 1);
            
            shared_idx = 1:8;
            param_kslow2(shared_idx) = param_kslow1(shared_idx);
            
            uniq_idx = setdiff(1:num_param, shared_idx);
            uniq_param = kslow2_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param_kslow2(uniq_idx) = uniq_param;
        case 'ikss'
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
            num_param = 11;
            param_kur = NaN(num_param, 1);
            
            shared_idx = 1:8;
            param_kur(shared_idx) = param_kslow1(shared_idx);
            
            uniq_idx = setdiff(1:num_param, shared_idx);
            uniq_param = kur_default;
            uniq_param(tune_idx1) = p(tune_idx2);
            param_kur(uniq_idx) = uniq_param;
        case 'ik1'
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
    rmse_grad_list = NaN(model_struct(end).idx2(end), num_volts);
    for i = 1:num_volts
        yksum_i = yksum(:, i);
        protocol{2} = volts(i);

        [yksum_hat, state_vars_list, trans_rates_list] = kcurrent_model(model_struct, protocol);
        err = yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end);

        rmse_list(i) = sqrt(mean(err.^2));
        rmse_grad_list(:, i) = rmse_grad(model_struct, protocol, err, state_vars_list, trans_rates_list);
    end
    f = sum(rmse_list);
    g = sum(rmse_grad_list, 2);
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

function [current_trc, state_vars, trans_rates] = ikto(hold_volt, volt, time_space, ek)
    global param_kto
    p = param_kto;

    % constants & initial values
    gmax = p(13);
    act0 = 0.4139033547E-02;
    inact0 = 0.9999623535E+00;
    
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

    ass_pulse = tr_pulse(1)./(tr_pulse(1)+tr_pulse(2));
    atau_pulse = 1./(tr_pulse(1)+tr_pulse(2));
    iss_pulse = tr_pulse(3)./(tr_pulse(3)+tr_pulse(4));
    itau_pulse = 1./(tr_pulse(3)+tr_pulse(4));

    % state variables
    act_hold = hh_model(hold_t, act0, ass_hold, atau_hold);
    inact_hold = hh_model(hold_t, inact0, iss_hold, itau_hold);      
    
    act_pulse = hh_model(pulse_t, act0, ass_pulse, atau_pulse);
    inact_pulse = hh_model(pulse_t, inact0, iss_pulse, itau_pulse);
    
    % current equation
    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(hold_volt-ek);
    current_trc((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-ek);

    state_vars = [act_pulse, inact_pulse];
    trans_rates = tr_pulse;
end

function tr = ikto_trans_rates(V)
    global param_kto
    p = param_kto;
    tr = NaN(8,1);
    
    % activation
    alpha1 = p(7).*exp(p(5).*(V+p(1)));
    beta1 = p(8).*exp(-p(6).*(V+p(1)));
    
    % inactivation
    alpha2_temp1 = p(9).*exp((V+p(2))./(-1.0.*p(4)));
    alpha2_temp2 = p(10).*exp((V+p(2)+p(3))./(-1.0.*p(4)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(11).*exp((V+p(2)+p(3))./p(4));
    beta2_temp2 = p(12).*exp((V+p(2)+p(3))./p(4));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    tr(1) = alpha1;
    tr(2) = beta1;
    tr(3) = alpha2;
    tr(4) = beta2;
end

function [current_trc, state_vars, trans_rates] = ikslow1(hold_volt, volt, time_space, ek)
    global param_kslow1
    p = param_kslow1;
    
    % constants & initial values
    gmax = p(11);
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
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt-ek);

    % current equation at pulse voltage
    tr_pulse = ikslow1_trans_rates(volt);
    act_pulse = hh_model(pulse_t, act0, tr_pulse(1), tr_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, tr_pulse(2), tr_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt-ek);

    state_vars = [act_pulse, inact_pulse];
    trans_rates = tr_pulse;
end

function tr = ikslow1_trans_rates(V)
    global param_kslow1
    p = param_kslow1;

    tr = NaN(4, 1);
    tr(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    tr(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    tr(3) = p(7)./(exp(p(6)*(V+p(3))) + exp(-p(6)*(V+p(3))))+p(9); % taua
    tr(4) = p(10) - p(8)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function [current_trc, state_vars, trans_rates] = ikslow2(hold_volt, volt, time_space, Ek)
    global param_kslow2
    p = param_kslow2;
    
    % constants & initial values
    gmax = p(11);
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

    tr = NaN(4,1);
    tr(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    tr(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    tr(3) = p(7)./(exp(p(6)*(V+p(3))) + exp(-p(6)*(V+p(3))))+p(9); % taua
    tr(4) = p(10) - p(8)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function [current_trc, state_vars, trans_rates] = ikss(hold_volt, volt, time_space, Ek)
    global param_kss
    p = param_kss;

    % constants & initial values
    gmax = p(7);
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

    tr = NaN(2,1);
    tr(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    tr(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end

function g = rmse_grad(model_struct, protocol, err, state_vars_list, trans_rates_list)
    dfdm = 1/(2*sqrt(mean(err.^2)));

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
    tune_idx1_kto = [1, 2, 4, 5, 7, 11, 13];
    
    global param_kto
    p = param_kto;
    g = NaN(length(tune_idx1_kto), 1);

    act0 = 0.4139033547E-02;
    inact0 = 0.9999623535E+00;
    gmax = param_kto(13);
    
    time_space = protocol{3};
    t = time_space{3};
    n = length(t);
    volt = protocol{2};
    ek = protocol{4};

    alpha1 = trans_rates(1);
    beta1 = trans_rates(2);
    alpha2 = trans_rates(3);
    beta2 = trans_rates(4);

    ass = alpha1/(alpha1+beta1);
    atau = 1/(alpha1+beta1);
    iss = alpha2/(alpha2+beta2);
    itau = 1/(alpha2+beta2);

    % derivatives
    dmde = 2*err;
    dedy = -1;
    dfdy = dfdm*dmde*dedy;
    
    dyda = 3*gmax*(state_vars(:,1).^2).*(state_vars(:,2))*(volt-ek);
    dydi = gmax*state_vars(:,1).^3*(volt-ek);
    
    dadss = 1 - exp(-t./atau);
    dadtau = -(ass - act0)*(t.*exp(-t./atau))/(atau^2);
    didss = 1- exp(-t./itau);
    didtau = -(iss - inact0)*(t.*exp(-t./itau))/(itau^2);
    
    dssdalpha1 = beta1/(alpha1+beta1)^2;
    dssdbeta1 = -alpha1/(alpha1+beta1)^2;
    dtaudalpha1 = -1/(alpha1+beta1)^2;
    dtaudbeta1 = -1/(alpha1+beta1)^2;
    dssdalpha2 = beta2/(alpha2+beta2)^2;
    dssdbeta2 = -alpha2/(alpha2+beta2)^2;
    dtaudalpha2 = -1/(alpha2+beta2)^2;
    dtaudbeta2 = -1/(alpha2+beta2)^2;
    
    % p1
    dalpha1dp1 = p(5)*p(7)*exp(p(5)*(volt+p(1)));
    dbeta1dp1 = -p(6)*p(8)*exp(-p(6)*(volt+p(1)));
    g(1) = (1/n)*sum(dfdy.*(...
        dyda.*(dadss*dssdalpha1*dalpha1dp1 + ...
            dadss*dssdbeta1*dbeta1dp1 + ...
            dadtau*dtaudalpha1*dalpha1dp1 + ...
            dadtau*dtaudbeta1*dbeta1dp1)));

    % p2
    dalpha2dp2_temp1 = p(9)*exp(2*(volt+p(2)+p(3))/p(4)-(volt+p(2))/p(4));
    dalpha2dp2_temp2 = p(4)*(exp((volt+p(2)+p(3))/p(4))+p(10))^2;
    dalpha2dp2 = -dalpha2dp2_temp1/dalpha2dp2_temp2;
    dbeta2dp2_temp1 = p(11)*exp((volt+p(2)+p(3))/p(4));
    dbeta2dp2_temp2 = p(4)*(p(12)*exp((volt+p(2)+p(3))/p(4))+1.0)^2;
    dbeta2dp2 = dbeta2dp2_temp1/dbeta2dp2_temp2;
    g(2) = (1/n)*sum(dfdy.*(...
        dydi.*(didss*dssdalpha2*dalpha2dp2 + ...
            didss*dssdbeta2*dbeta2dp2 + ...
            didtau*dtaudalpha2*dalpha2dp2 + ...
            didtau*dtaudbeta2*dbeta2dp2)));

    % p4
    dalpha2dp4_temp1 = p(9)*((volt+p(2))*exp((volt+p(2)+p(3))/p(4))-p(10)*p(3))*exp((volt+p(2)+p(3))/p(4)-(volt+p(2))/p(4));
    dalpha2dp4_temp2 = p(4)^2*(exp((volt+p(2)+p(3))/p(4))+p(10))^2;
    dalpha2dp4 = dalpha2dp4_temp1/dalpha2dp4_temp2;
    dbeta2dp4_temp1 = p(11)*(volt+p(2)+p(3))*exp((volt+p(2)+p(3))/p(4));
    dbeta2dp4_temp2 = p(4)^2*(p(12)*exp((volt+p(2)+p(3))/p(4))+1.0)^2;
    dbeta2dp4 = -dbeta2dp4_temp1/dbeta2dp4_temp2;
    g(3) = (1/n)*sum(dfdy.*(...
        dydi.*(didss*dssdalpha2*dalpha2dp4 + ...
            didss*dssdbeta2*dbeta2dp4 + ...
            didtau*dtaudalpha2*dalpha2dp4 + ...
            didtau*dtaudbeta2*dbeta2dp4)));
    
    % p5
    dalpha1dp5 = (volt+p(1))*p(7)*exp(p(5)*(volt+p(1)));
    g(4) = (1/n)*sum(dfdy.*(...
        dyda.*(dadss*dssdalpha1*dalpha1dp5 + ...
            dadtau*dtaudalpha1*dalpha1dp5)));

    % p7
    dalpha1dp7 = exp(p(5)*(volt+p(1)));
    g(5) = (1/n)*sum(dfdy.*(...
        dyda.*(dadss*dssdalpha1*dalpha1dp7 + ...
            dadtau*dtaudalpha1*dalpha1dp7)));

    % p11
    dbeta2dp11 = (exp((volt+p(2)+p(3))/p(4)))/(p(12)*exp((volt+p(2)+p(3))/p(4))+1.0);
    g(6) = (1/n)*sum(dfdy.*(...
        dydi.*(didss*dssdbeta2*dbeta2dp11 + ...
            didtau*dtaudbeta2*dbeta2dp11)));

    % p13
    dydgmax = (state_vars(:,1).^3).*state_vars(:,2)*(volt-ek);
    g(7) = (1/n)*sum(dfdy.*dydgmax);
end

function g = grad_ikslow1(protocol, dfdm, err, state_vars, trans_rates)
    tune_idx1_kslow1 = [1, 2, 4, 5, 9, 10, 11];
    
    global param_kslow1
    p = param_kslow1;
    g = NaN(length(tune_idx1_kslow1), 1);

    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;    
    gmax = p(11);

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

    dyda = gmax*state_vars(:, 2)*(volt-ek);
    dydi = gmax*state_vars(:, 1)*(volt-ek);

    dadass = 1.0-exp(-t./atau);
    dadatau = -(ass - act0)*(t.*exp(-t./atau))/(atau^2);
    didiss = 1.0-exp(-t./itau);
    diditau = -(iss - inact0)*(t.*exp(-t./itau))/(itau^2);

    % p1
    dassdp1 = (exp((p(1)+volt)/p(4)))/(p(4)*(exp((p(1)+volt)/p(4))+1.0)^2);
    g(1) = (1/n)*sum(dfdy.*(dyda.*dadass*dassdp1));

    % p2
    dissdp2 = -exp((p(2)+volt)/p(5))/(p(5)*(exp((p(2)+volt)/p(5))+1.0)^2);
    ditaudp2 = (p(8)*exp((p(2)+volt)/p(5)))/(p(5)*(exp((p(2)+volt)/p(5))+1.0)^2);
    g(2) = (1/n)*sum(dfdy.*dydi.*(didiss*dissdp2 + diditau*ditaudp2));

    % p4
    dassdp4 = -((p(1)+volt)*exp((p(1)+volt)/p(4)))/(p(4)^2*(exp((p(1)+volt)/p(4))+1.0)^2);
    g(3) = (1/n)*sum(dfdy.*dyda.*dadass*dassdp4);

    % p5
    dissdp5 = ((p(2)+volt)*exp((p(2)+volt)/p(5)))/(p(5)^2*(exp((p(2)+volt)/p(5))+1.0)^2);
    ditaudp5 = -(p(8)*(p(2)+volt)*exp((p(2)+volt)/p(5)))/(p(5)^2*(exp((p(2)+volt)/p(5))+1)^2);
    g(4) = (1/n)*sum(dfdy.*dydi.*(didiss*dissdp5 + diditau*ditaudp5));

    % p9
    dataudp9 = 1;
    g(5) = (1/n)*sum(dfdy.*dyda.*dadatau*dataudp9);

    % p10
    ditaudp10 = 1;
    g(6) = (1/n)*sum(dfdy.*dydi.*diditau*ditaudp10);

    % p11
    dydgmax = state_vars(:, 1).*state_vars(:, 2)*(volt-ek);
    g(7) = (1/n)*sum(dfdy.*dydgmax);
end

function g = grad_ikslow2(protocol, dfdm, err, state_vars, trans_rates)
    % tune_idx1_kslow2 = [2, 3]; p([10, 11])
     
    global param_kslow2
    p = param_kslow2;
    g = NaN(2, 1);

    inact0 = 0.9980927689;
    gmax = p(11);

    time_space = protocol{3};
    t = time_space{3};
    n = length(t);
    volt = protocol{2};
    ek = protocol{4};
    
    iss = trans_rates(2);
    itau = trans_rates(4);
    
    % derivatives
    dmde = 2*err;
    dedy = -1;
    dfdy = dfdm*dmde*dedy;

    dydi = gmax*state_vars(:, 1)*(volt-ek);
    diditau = -(iss - inact0)*(t.*exp(-t./itau))/(itau^2);

    % p2 (p(10))
    ditaudp2 = 1;
    g(1) = (1/n)*sum(dfdy.*dydi.*diditau*ditaudp2);

    % p3 (p(11))
    dydgmax = state_vars(:, 1).*state_vars(:, 2)*(volt-ek);
    g(2) = (1/n)*sum(dfdy.*dydgmax);
end

function g = grad_ikss(protocol, dfdm, err, state_vars, trans_rates)
    tune_idx1_kss = [1, 2, 3, 4]; % p([4, 5, 6, 7])
    global param_kss
    global param_kslow1
    p = param_kss;
    x = param_kslow1(3);
    g = NaN(length(tune_idx1_kss), 1);

    gmax = p(7);
    act0 = 0.5091689794e-03;

    time_space = protocol{3};
    t = time_space{3};
    n = length(t);
    volt = protocol{2};
    ek = protocol{4};
    
    ass = trans_rates(1);
    atau = trans_rates(2);

    % derivatives
    dmde = 2*err;
    dedy = -1;
    dfdy = dfdm*dmde*dedy;
    dyda = gmax*(volt-ek);
    dadtau = -(ass - act0)*(t.*exp(-t./atau))/(atau^2);

    % p1
    dtaudp1_temp1 = p(5)*(volt+x)*exp(p(4)*(volt+x))*(exp(2*p(4)*(volt+x))-1.0);
    dtaudp1_temp2 = (exp(2*p(4)*(volt+x))+1.0)^2;
    dtaudp1 = -dtaudp1_temp1/dtaudp1_temp2;
    g(1) = (1/n)*sum(dfdy.*dyda.*dadtau.*dtaudp1);

    % p2
    dtaudp2 = 1/(exp(p(4)*(volt+x))+exp(-p(4)*(volt+x)));
    g(2) = (1/n)*sum(dfdy.*dyda.*dadtau.*dtaudp2);

    % p3
    dtaudp3 = 1;
    g(3) = (1/n)*sum(dfdy.*dyda.*dadtau.*dtaudp3);

    % p4
    dydgmax = state_vars(:, 1)*(volt-ek);
    g(4) = (1/n)*sum(dfdy.*dydgmax);
end
