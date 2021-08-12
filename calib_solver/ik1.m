function [current_trc] = ik1(p, hold_volt, volt, time_space, Ek)
    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    % pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);
    
    % current equation at holding
    [alpha_k1, beta_k1] = transition_rates(p, hold_volt, Ek);
    current_trc(1:hold_idx) = 0.27*sqrt(5400/5400) .* (alpha_k1/(alpha_k1 + beta_k1)) .* (hold_volt - Ek);

    [alpha_k1, beta_k1] = transition_rates(p, volt, Ek);
    current_trc((hold_idx+1):end) = 0.27*sqrt(5400/5400) .* (alpha_k1/(alpha_k1 + beta_k1)) .* (hold_volt - Ek);
end

function [alpha_k1, beta_k1] = transition_rates(p, v, Ek)
    % p0 = [59.215, 5.476, 594.31, 4.753]
    alpha_k1 = 1.02 ./ (1 + exp(0.2385.*(v-Ek-p(1))));
    beta_k1 = (0.8*exp(0.08032*(v-Ek+p(2))) + exp(0.06175*(v-Ek-p(3)))) ./ ...
        (1 + exp(-0.5143*(v-Ek+p(4))));
end
