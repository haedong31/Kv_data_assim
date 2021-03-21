function [current_trc] = IKtof(param, holdV, P1, time_space)
    % Bondarenko IKtof (A67)
    % 12 parameters;

    % constants
    act0 = 0.265563e-2; % ato_f; Gating variable for transient outward K+ current
    inact0 = 0.999977; % ito_f; Gating variable for transient outward K+ current
    gmax = param(13);  % 0.4067
    Ek = param(14); % resting potential

    % time space information
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};
    hold_idx = length(tH);

    current_trc = zeros(length(t), 1);

    % at holding potential
    cv_hold = cv(param(1:12), holdV);
    act_hold = cv_hold(1) - (cv_hold(1) - act0).*exp(-(tH./cv_hold(2)));
    inact_hold = cv_hold(3) - (cv_hold(3) - inact0).*exp(-(tH./cv_hold(4)));
    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(holdV - Ek);

    % at P1
    cv_P1 = cv(param(1:12), P1);
    act_P1 = cv_P1(1) - (cv_P1(1) - act0).*exp(-(tP1_adj./cv_P1(2)));
    inact_P1 = cv_P1(3) - (cv_P1(3) - inact0).*exp(-(tP1_adj./cv_P1(4)));
    current_trc((hold_idx + 1):end) = gmax.*(act_P1.^3).*(inact_P1).*(P1 - Ek);
end

function [cv] = cv(p, V)
    % characteristic variables of IKtof
    % cv(1): steady-state activation
    % cv(2): time constant of activation
    % cv(3): steady-state inactivation
    % cv(4): time constant of inactivation
    
    % p0 = [30.0, 13.5, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335];
    cv = zeros(4, 1);
    alphaA = p(5).*exp(p(6).*(V+p(1)));
    betaA = p(7).*exp(-p(8).*(V+p(1)));
    alphaI = (p(9).*exp(-(V+p(2))./p(4))) ./ (p(10).*exp(-(V+p(3))./p(4)) + 1);
    betaI = (p(11).*exp((V+p(3))./p(4))) ./ (p(12).*exp((V+p(3))./p(4)) + 1);

    cv(1) = alphaA/(alphaA+betaA);
    cv(2) = 1/(alphaA+betaA);
    cv(3) = alphaI/(alphaI+betaI);
    cv(4) = 1/(alphaI+betaI);
end
