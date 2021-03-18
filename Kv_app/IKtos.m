function [current_trc] = IKtos(param, holdV, P1, time_space)
    % Bondarenko IKtos (A75)
    % 11 parameters + 1 conductance + 1 Ek

    act0 = 0.417069e-3;  % ato_s; Gating variable for transient outward K+ current
    inact0 = 0.998543;  % ito_s; Gating variable for transient outward K+ current
    gmax = param(12); % 0.0629 in spetal model; 0 in apical model
    Ek = param(13); % resting potential

    % time space information
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};
    hold_idx = length(tH);

    current_trc = zeros(length(t), 1);

    % at holding potential
    cv_hold = cv(param(1:11), holdV);
    act_hold = cv_hold(1) - (cv_hold(1) - act0).*exp(-(tH./cv_hold(2)));
    inact_hold = cv_hold(3) - (cv_hold(3) - inact0).*exp(-(tH./cv_hold(4)));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(holdV - Ek);

    % at P1
    cv_P1 = cv(param(1:11), P1);
    act_P1 = cv_P1(1) - (cv_P1(1) - act0).*exp(-(tP1_adj./cv_P1(2)));
    inact_P1 = cv_P1(3) - (cv_P1(3) - inact0).*exp(-(tP1_adj./cv_P1(4)));
    current_trc((hold_idx + 1):end) = gmax.*(act_P1).*(inact_P1).*(P1 - Ek);
end

function [cv] = cv(p, V)
    % p0 = [22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 270, 1050, 45.2, 5.7]
    cv = zeros(4, 1);
    
    % note similarity of IKtos and IKur
    % same activation gate characteristics
    cv(1) = 1 ./ (1 + exp(-(V+p(1))./p(2)));
    cv(2) = p(3).*exp(-p(4).*V) + p(5);

    % for inactivation gate, cv(3) is same but not for cv(4) 
    cv(3) = 1 ./ (1 + exp((V+p(6))./p(7)));
    cv(4) = p(8) + p(9)./(1 + exp((V + p(10))./p(11)));
end
