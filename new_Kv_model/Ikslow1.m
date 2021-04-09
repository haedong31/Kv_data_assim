function [current_trc] = Ikslow1(p, hold_volt, volt, time_space, Ek)
    % 12 parameters; {p(11): f, p(12): gmax}
    
    % constants & initial values
    f_eacv = p(11); % Ikslow1 fraction of nonphospholatedl; default value: 0.9214774521
    gmax = p(12); % default value: 0.05766
    gmaxp = p(13); % default value: 0.07496
    act0 = 0.417069e-3; 
    inact0 = 0.998543; 
    
    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding 
    gv_hold = gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));        
    current_trc(1:hold_idx) = (gmax*f_eacv + gmaxp*(1-f_eacv)).*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = (gmax*f_eacv + gmaxp*(1-f_eacv)).*(act_pulse).*(inact_pulse).*(volt - Ek);
end

function [gv] = gating_variables(p, V)
    % p0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0]

    gv = zeros(4, 1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0/(1.0+exp((p(2)+V)/p(5))); % iss
    gv(3) = p(6)/(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    gv(4) = p(9)-p(10)/(1.0+exp((p(2)+V)/p(5))); % taui
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
