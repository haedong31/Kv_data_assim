function [current_trc] = IKs(param, holdV, P1, time_space)
    % Bondarenko IKs (A83)
    % 5 parameters + 1 conductance + 1 Ek
    
    % constants
    n0 = 0.262753e-3;  % nKs; Gating variable for slow delayed-rectifier K+ current
    gmax = param(6); % 0.00575; GKs; Maximum slow delayed-rectifier K+ current conductance:mS/uF
    Ek = param(7);
    
    % time space information
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};
    hold_idx = length(tH);

    current_trc = zeros(length(t), 1);

    % at holding potential
    cv_hold = cv(param(1:5), holdV);
    n = cv_hold(1) - (cv_hold(1) - n0).*exp(-(tH./cv_hold(2)));
    current_trc(1:hold_idx) = gmax.*power(n, 2).*(holdV - Ek);

    % at P1
    cv_P1 = cv(param(1:5), P1);
    n = cv_P1(1) - (cv_P1(1) - n0).*exp(-(tP1_adj./cv_P1(2)));
    current_trc((hold_idx + 1):end) = gmax.*power(n, 2).*(P1 - Ek);
end

function [cv] = cv(p, V)
    % characteristic variables of IKs
    % cv(1): steady-state value of activation gate
    % cv(2): time constant of activation gate

    % p0 = [26.5, 0.128+, 0.038+, 4.81333e-06+, 9.53333e-05+]
    cv = zeros(2, 1);
    a = (p(4).*(V + p(1)))./(1 - exp(-p(2).*(V + p(1))));
    b = p(5).*exp(-p(3).*(V + p(1)));

    cv(1) = a./(a + b);
    cv(2) = 1./(a + b);
end
