function bm = ikss_biomarkers(p, protocol_info)
    % protocol info
    hold_volt = protocol_info{1};
    volt = protocol_info{2};
    ek = protocol_info{3};
    time_space = protocol_info{4};

    % generate current
    c = ikss(p, hold_volt, volt, time_space, ek);
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

function current_trc = ikss(p, hold_volt, volt, time_space, ek)
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
    current_trc(1:hold_idx) = gmax.*act_hold.*(hold_volt - ek);

    % current equation at pulse voltage
    gv_pulse = ikss_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(2));
    current_trc((hold_idx + 1):end) = gmax.*act_pulse.*(volt - ek);
end

function gv = ikss_gating_variables(p, V)
    % gv(1) = gv(1) in Ikslow1
    % p0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428]

    gv = zeros(2,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    gv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
