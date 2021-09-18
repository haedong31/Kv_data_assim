function [f, g] = obj_rmse_grad(p, model_struct, volt_space, time_space, yksum)
    hold_idx = time_space{4};
    end_idx = time_space{5};
    
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

        [yksum_hat, state_vars_list, trans_rates_list] = kcurrent_model(p, model_struct, protocol);
        err = yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end);
        mserr = mean(err.^2);

        rmse_list(i) = sqrt(mserr);
        grad_list(:, i) = kcurrent_grad(p, model_struct, protocol, mserr, err, state_vars_list, trans_rates_list);
    end
    f = sum(rmse_list);
end

function [yksum_hat, state_vars_list, trans_rates_list] = kcurrent_model(p, model_struct, protocol_info)
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
            
            [current_trace, state_vars, trans_rates] = ikto(hold_volt, volt, time_space, ek);
        case 'ikslow1'
            % generate ikslow1
            num_kslow1_param = 13;
            param_kslow1 = NaN(num_kslow1_param, 1);
            
            fixed_idx = setdiff(1:num_kslow1_param, tune_idx1);
            param_kslow1(tune_idx1) = p(tune_idx2);
            param_kslow1(fixed_idx) = kslow1_default(fixed_idx);
    
            [current_trace, state_vars, trans_rates] = ikslow1(hold_volt, volt, time_space, ek);
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

            [current_trace, state_vars, trans_rates] = ikslow2(hold_volt, volt, time_space, ek);
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

            [current_trace, state_vars, trans_rates] = ikss(hold_volt, volt, time_space, ek);
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

            [current_trace, state_vars, trans_rates] = ikur(hold_volt, volt, time_space, ek);            
        case 'ik1'
            % generate ik1
            num_param = 10;
            param_k1 = NaN(num_param, 1);

            fixed_idx = setdiff(1:num_param, tune_idx1);
            param_k1(tune_idx1) = p(tune_idx2);
            param_k1(fixed_idx) = k1_default(fixed_idx);

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
    global param_k1
    p = param_k1;
    
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

function g = kcurrent_grad(p, model_struct, protocol, mserr, err, state_vars_list, trans_rates_list)
    global param_kslow1
    global param_kslow2
    global param_kss
    global param_kur
    global param_k1
    
    df1 = 1/(2*sqrt(mserr));
    df2 = 2*mean(err);
    df3 = -1;

    g = NaN(length(p), 1);
    for i = 1:length(model_struct)
        current_name = model_struct(i).name;
        tune_idx2 = model_struct(i).idx2;
        state_vars = state_vars_list{i};
        trans_rates = trans_rates_list{i};

        switch current_name
        case 'ikto'
            g(tune_idx2) = grad_ikto(p, df1, df2, df3, state_vars, trans_rates);
        case 'ikslow1'
            g(tune_idx2) = grad_ikslow1(p, df1, df2, df3, state_vars, trans_rates);
        case 'ikslow2'
            g(tune_idx2) = grad_ikslow2(p, df1, df2, df3, state_vars, trans_rates);
        case 'ikss'
            g(tune_idx2) = grad_ikss(p, df1, df2, df3, state_vars, trans_rates);
        otherwise
            g(tune_idx2) = [];
        end
    end
end

function g = grad_ikto(p, df1, df2, df3, state_vars, trans_rates)
    % tune_idx1_kto = [1, 2, 6, 10, 13, 14, 15, 16, 17];
    
    global param_kto
    g = NaN(9, 1);

    % a

    % i

end

function g = grad_ikslow1(p, df1, df2, df3, state_vars, trans_rates)
    global param_kslow1
    g = [];
end

function g = grad_ikslow2(p, df1, df2, df3, state_vars, trans_rates)
    global param_kslow2
    g = [];
end
function g = grad_ikss(p, df1, df2, df3, state_vars, trans_rates)
    global param_kss
    g = [];
end
