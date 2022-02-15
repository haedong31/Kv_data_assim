function current_trc = ikss(p, hold_volt, volt, time_space, ek)
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
    kv_hold = ikss_kinetic_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, kv_hold(1), kv_hold(2));
    current_trc(1:hold_idx) = gmax.*act_hold.*(hold_volt - ek);

    % current equation at pulse voltage
    kv_pulse = ikss_kinetic_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, kv_pulse(1), kv_pulse(2));
    current_trc((hold_idx + 1):end) = gmax.*act_pulse.*(volt - ek);
end

function kv = ikss_kinetic_variables(p, V)
    kv = zeros(2,1);
    kv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    kv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
