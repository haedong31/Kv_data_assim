function [current_trc] = IKI(param, holdV, P1, time_space)
    % Bondarenko IKI (A82)
    % 3 parameters + 1 extracellular K+ concentration + 1 Ek
    
    % constants
    extra_Kconcent = param(4);
    Ek = param(5);

    % time space information
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};
    hold_idx = length(tH);

    % p0 = [0.2938+, 210, 0.0896+]
    current_trc = zeros(length(t), 1);

    % at holding potential
    current_trc(1:hold_idx) = param(3).*(extra_Kconcent./(extra_Kconcent + param(2))).*((holdV - Ek)./(1 + exp(param(1).*(holdV - Ek))));

    % at P1
    current_trc((hold_idx + 1):end) = param(3).*(extra_Kconcent./(extra_Kconcent + param(2))).*((P1 - Ek)./(1 + exp(param(1).*(P1 - Ek))));
end