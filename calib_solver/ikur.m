function [current_trc] = ikur(p, hold_volt, volt, time_space, Ek)
    gmax = p(1);

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);
    
    current_trc = zeros(length(t), 1);
end

function [gv] = gating_variables(p, v)
    gv = zeros(4, 1);
    gv(1) = 1.0 ./ (1.0 + exp((-22.5-v)./7.7));
    gv(2) = 1.0 ./ (1.0 + exp((45.2+v)/5.7));
    gv(3) = 6.1 ./ (exp(0.0629*(v+40.0)) + exp(-0.0629*(v+40.0))) + 2.058;
    gv(4) = 270.0 + 1050.0 ./ (1.0 + exp((45.2+v)/5.7));
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
