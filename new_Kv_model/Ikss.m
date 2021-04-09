function [current_trc] = Ikss(p, hold_volt, volt, time_space, Ek)
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
    gv_hold = gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(2));
    current_trc(1:hold_idx) = gmax.*act_hold.*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(2));
    current_trc((hold_idx + 1):end) = gmax.*act_pulse.*(volt - Ek);
end

function [gv] = gating_variables(p, V)
    % gv(1) = gv(1) in Ikslow1 = gv(1) in Ikslow2
    % p0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17]

    gv = zeros(2,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    gv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
