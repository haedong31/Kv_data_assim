function [current_trc,kv_pulse] = iktof(p, hold_volt, volt, time_space, ek)
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

    % current equation at holding 
    kv_hold = ikto_kinetic_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(hold_t, inact0, kv_hold(2), kv_hold(4));      
    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(hold_volt-ek);

    % current equation at pulse voltage
    kv_pulse = ikto_kinetic_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, kv_pulse(2), kv_pulse(4));
    current_trc((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-ek);
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

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
