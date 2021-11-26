function bm = ikslow2_biomarkers(p, protocol_info)
    % protocol info
    hold_volt = protocol_info{1};
    volt = protocol_info{2};
    ek = protocol_info{3};
    time_space = protocol_info{4};

    % generate current
    c = ikslow2(p, hold_volt, volt, time_space, ek);
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

function current_trc = ikslow2(p, hold_volt, volt, time_space, ek)
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
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - ek);

    % current equation at pulse voltage
    gv_pulse = ikslow2_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - ek);
end

function gv = ikslow2_gating_variables(p, V)
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

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end

