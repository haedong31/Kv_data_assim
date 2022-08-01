function bm = ikto_biomarkers2(p, protocol_info)
    % protocol info
    hold_volt = protocol_info{1};
    volt = protocol_info{2};
    ek = protocol_info{3};
    time_space = protocol_info{4};

    % generate current
    c = ikto(p, hold_volt, volt, time_space, ek);
    t = time_space{1};
    
    hold_t = time_space{2};
    cc = c(length(hold_t):end);
    tt = t(length(hold_t):end);
    pulse_t = length(cc);

    % calculate biomarkers
    bm = NaN(1,6);

    % very early (10 ms)
    bm(1) = cc(10);

    % 25% point
    bm(2) = cc(floor(pulse_t*0.25));
    
    % 50% point
    bm(3) = cc(floor(pulse_t*0.5));
    
    % 75% point
    bm(4) = cc(floor(pulse_t*0.75));

    % peak
    [bm(5), peak_idx] = max(cc(1:floor(pulse_t*0.1)));

    % tau
    [~,tau_idx] = min(abs(cc(peak_idx:end)-bm(5)*exp(-1)));
    tt = tt(peak_idx:end);
    bm(6) = tt(tau_idx)-length(hold_t);
end

function current_trc = ikto(p, hold_volt, volt, time_space, ek)    
    % constants & initial values
    gmax = p(13);
    act0 = 0.265563e-02;
    inact0 = 0.9999623535e+00;

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

    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(hold_volt-ek);

    % current equation at pulse voltage
    gv_pulse = ikto_gating_variables(p, volt);

    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));

    current_trc((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-ek);
end

function gv = ikto_gating_variables(p, V)
    gv = zeros(4,1);
    
    % for Ikto
    alpha1 = p(7).*exp(p(5).*(V+p(1)));
    beta1 = p(8).*exp(-p(6).*(V+p(1)));
    
    alpha2_temp1 = p(9).*exp((V+p(2))./(-1.0.*p(4)));
    alpha2_temp2 = p(10).*exp((V+p(2)+p(3))./(-1.0.*p(4)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(11).*exp((V+p(2)+p(3))./p(4));
    beta2_temp2 = p(12).*exp((V+p(2)+p(3))./p(4));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    gv(1) = alpha1./(alpha1+beta1);
    gv(2) = alpha2./(alpha2+beta2);
    gv(3) = 1./(alpha1+beta1);
    gv(4) = 1./(alpha2+beta2);
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
