function [current_trc] = IKss(param, holdV, P1, time_space)
    % Bondarenko IKss (A92)
    % 5 parameters + 1 conductance + 1 Ek

    % constants
    act0 = 0.417069e-3;  % aKss; Gating variable for noninactivating steady-state K+ current
    inact0 = 1;  % iKss; Gating variable for noninactivating steady-state K+ current
    gmax = param(6); % 0.05
    Ek = param(7);
 
    % time space information
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};
    hold_idx = length(tH);

    current_trc = zeros(length(t), 1);

    % at holding potential
    cv_hold = cv(param(1:5), holdV);
    act_hold = cv_hold(1) - (cv_hold(1) - act0).*exp(-(tH./cv_hold(2)));
    inact_hold = inact0;
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(holdV - Ek);

    % at P1
    cv_P1 = cv(param(1:5), P1); 
    act_P1 = cv_P1(1) - (cv_P1(1) - act0).*exp(-(tP1_adj./cv_P1(2)));
    inact_P1 = inact0;
    current_trc((hold_idx + 1):end) = gmax.*(act_P1).*(inact_P1).*(P1 - Ek);
end

function [cv] = cv(p, V)
    % characteristic variables of IKss
    % cv(1): (A78) steady-state value of activation gate
    % cv(2): (A95) time constant of activation gate 

    % p0 = [22.5, 7.7, 39.3, 0.0862, 13.17]
    cv = zeros(2, 1);
    % act_Kss
    cv(1) = 1./(1 + exp(-(V + p(1))./p(2)));
    % tasu act_Kss
    cv(2) = p(3).*exp(-p(4).*V) + p(5);
end