function [current_trc] = ikto2(p, hold_volt, volt, time_space, Ek)
    % see 2014 Bondarenko and 2020 Bondarenko
    
    % constants & initial values
    gmax = p(16); % 0.14067

    act0 = 0.4139033547E-02;
    inact0 = 0.9999623535E+00;

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

    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(hold_volt-Ek);

    % current equation at pulse voltage
    gv_pulse = gating_variables(p, volt);

    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));

    current_trc((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-Ek);
end

function [gv] = gating_variables(p, V)
    % p0 = [33.0, 40.0, 33.0, 15.5, 20.0, 7.7, 5.7, ...
    %     0.03577, 0.06237, 7.0, 0.18064, 0.3956, 0.0226585, 0.00095, 0.051335, 0.14067]

    gv = zeros(4,1);
    
    % for Ikto
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(6))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(7))); % iss

    alpha1 = p(11).*exp(p(8).*(V+p(3)));
    beta1 = p(12).*exp(-p(9).*(V+p(3)));
    gv(3) = 1./(alpha1+beta1);
    
    alpha2_temp1 = exp((V+p(4))./(-1.0.*p(10)));
    alpha2_temp2 = exp((V+p(4)+p(5))./(-1.0.*p(10)));
    alpha2 = p(13).*(alpha2_temp1./alpha2_temp2);

    beta2_temp1 = p(14).*exp((V+p(4)+p(5))./p(10));
    beta2_temp2 = p(15).*exp((V+p(4)+p(5))./p(10));
    beta2 = beta2_temp1./(beta2_temp2+1.0);

    gv(4) = 1./(alpha2+beta2);
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
