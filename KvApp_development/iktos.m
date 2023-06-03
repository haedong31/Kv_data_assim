function [current_trc,kv_pulse] = iktos(p, hold_volt, volt, time_space, ek)
    % constants & initial values
    gmax = p(11);
    act0 = 0.5091689794E-03;
    inact0 = 0.9980927689E+00;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);
    
    current_trc = NaN(length(t), 1);

    % current equation at holding
    kv_hold = iktos_kinetic_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(hold_t, inact0, kv_hold(2), kv_hold(4));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - ek);
    
    % current equation at pulse voltage
    kv_pulse = iktos_kinetic_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, kv_pulse(2), kv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - ek);    
end

function kv = iktos_kinetic_variables(p, V)
    kv = NaN(4,1);
    kv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    kv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    kv(3) = p(7)./(exp(p(6)*(V+p(3))) + exp(-p(6)*(V+p(3))))+p(9); % taua
    kv(4) = p(10) + p(8)./(1.0+exp((p(2)+V)./p(5))); % taui 
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
