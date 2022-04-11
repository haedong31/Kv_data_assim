function ymx = kcurrent_basic(p, mdl_struct, pdefault, protocol, yksum)
    current_names = {mdl_struct.name};
    idx_info1 = {mdl_struct.idx1};
    idx_info2 = {mdl_struct.idx2};
    
    % define param_kslow1 first just in case it would be used in advance
    matching_idx = strcmp(current_names, 'ikslow1');
    param_kslow1 = pdefault{3};
    if any(matching_idx)
        param_kslow1(idx_info1{matching_idx}) = p(idx_info2{matching_idx});
    end

    ymx = zeros(length(protocol{4}), length(current_names));
    for i = 1:length(current_names)
        switch current_names{i}
            case 'iktof'
                param = pdefault{1};
                param(idx_info1{i}) = p(idx_info2{i});
                ymx(:,i) = iktof(param, protocol, volt);
            case 'iktos'
                num_param = 11;
                shared_idx = [1:7,9];
                uniq_idx = setdiff(1:num_param, shared_idx);

                param = NaN(num_param,1);
                param(shared_idx) = param_kslow1(shared_idx);
                uniq_param = pdefault{2};
                uniq_param(idx_info1{i}) = p(idx_info2{i});
                param(uniq_idx) = uniq_param;
                ymx(:,i) = iktos(param, protocol, volt);
            case 'ikslow1'
                ymx(:,i) = ikslow1(param_kslow1, protocol, volt);
            case 'ikslow2'
                num_param = 11;
                shared_idx = [1:7,9];
                uniq_idx = setdiff(1:num_param, shared_idx);

                param = NaN(num_param,1);
                param(shared_idx) = param_kslow1(shared_idx);
                uniq_param = pdefault{4};
                uniq_param(idx_info1{i}) = p(idx_info2{i});
                param(uniq_idx) = uniq_param;
                ymx(:,i) = ikslow2(param, protocol, volt);
            case 'ikss'
                num_param = 7;
                shared_idx1 = 1:3;
                shared_idx2 = [1,3,4];
                uniq_idx = setdiff(1:num_param, shared_idx1);

                param = NaN(num_param,1);
                param(shared_idx1) = param_kslow1(shared_idx2);
                uniq_param = pdefault{5};
                uniq_param(idx_info1{i}) = p(idx_info2{i});
                param(uniq_idx) = uniq_param;
                ymx(:,i) = ikss(param, protocol, volt);
        end
    end
end

function y = iktof(p, protocol, volt)
    gmax = p(13);
    act0 = 0.4139033547e-02;
    inact0 = 0.9999623535e+00;
    
    holdv = protocol{1};
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding volt
    kv_hold = ikto_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));
    y(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(holdv-ek);

    % current equation at pulse volt
    kv_pulse = ikto_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-ek);    
end

function y = iktos(p, protocol, volt)
    gmax = p(11);
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689e+00;

    holdv = protocol{1};
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));
    y(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(holdv - ek);
    
    % current equation at pulse voltage
    kv_pulse = rectifier_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - ek);    
end

function y = ikslow1(p, protocol, volt)
    gmax = p(11);
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;
    
    holdv = protocol{1};
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding 
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));        
    y(1:hold_idx) = (gmax).*(act_hold).*(inact_hold).*(holdv-ek);

    % current equation at pulse voltage
    kv_pulse = rectifier_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = (gmax).*(act_pulse).*(inact_pulse).*(volt-ek);
end

function y = ikslow2(p, protocol, volt)
    gmax = p(11); 
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;

    holdv = protocol{1};
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));
    y(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(holdv - ek);

    % current equation at pulse voltage
    kv_pulse = rectifier_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - ek);
end

function y = ikss(p, protocol, volt)
    gmax = p(7);
    act0 = 0.5091689794e-03;

    holdv = protocol{1};
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding
    kv_hold = ikss_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(2));
    y(1:hold_idx) = gmax.*act_hold.*(holdv - ek);

    % current equation at pulse voltage
    kv_pulse = ikss_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(2));
    y((hold_idx+1):end) = gmax.*act_pulse.*(volt - ek);
end

function kv = ikto_kinetic_variables(p, V)
    kv = NaN(8,1);
    
    alpha1 = p(7).*exp(p(5).*(V+p(1)));
    beta1 = p(8).*exp(-p(6).*(V+p(1)));
    
    alpha2_temp1 = p(9).*exp((V+p(2))./(-1.0.*p(4)));
    alpha2_temp2 = p(10).*exp((V+p(2)+p(3))./(-1.0.*p(4)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(11).*exp((V+p(2)+p(3))./p(4));
    beta2_temp2 = p(12).*exp((V+p(2)+p(3))./p(4));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    kv(1) = alpha1./(alpha1+beta1);
    kv(2) = alpha2./(alpha2+beta2);
    kv(3) = 1./(alpha1+beta1);
    kv(4) = 1./(alpha2+beta2);
    kv(5) = alpha1;
    kv(6) = beta1;
    kv(7) = alpha2;
    kv(8) = beta2;
end

function kv = rectifier_kinetic_variables(p, V)
    kv = NaN(4,1);
    kv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    kv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    kv(3) = p(7)./(exp(p(6)*(V+p(3))) + exp(-p(6)*(V+p(3))))+p(9); % taua
    kv(4) = p(10) - p(8)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function kv = ikss_kinetic_variables(p, V)
    kv = NaN(2,1);
    kv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    kv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
