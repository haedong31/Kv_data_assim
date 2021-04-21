function [current_trc] = IKur(param, holdV, P1, time_space)
    % Bondarenko IKur (87)
    % 11 parameters + 1 conductance + 1 Ek

    % constants
    act0 = 0.417069e-3;  % aur; Gating variable for 
    inact0 = 0.998543;  % iur; Gating variable for ultrarapidly 
    gmax = param(10); % 0.16
    Ek = param(11); % resting potential

    % time space information
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};
    hold_idx = length(tH);

    current_trc = zeros(length(t), 1);

    % at holding potential
    cv_hold = cv(param, holdV);
    act_hold = cv_hold(1) - (cv_hold(1) - act0).*exp(-(tH./cv_hold(3)));
    inact_hold = cv_hold(2) - (cv_hold(2) - inact0).*exp(-(tH./cv_hold(4)));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(holdV - Ek);

    % at P1
    cv_P1 = cv(param, P1);
    act_P1 = cv_P1(1) - (cv_P1(1) - act0).*exp(-(tP1_adj./cv_P1(3)));
    inact_P1 = cv_P1(2) - (cv_P1(2) - inact0).*exp(-(tP1_adj./cv_P1(4)));
    current_trc((hold_idx + 1):end) = gmax.*(act_P1).*(inact_P1).*(P1 - Ek);
end

function [cv] = cv(p, V)
    % characteristic variables of IKur

    % p0 = [22.5, 7.7, 45.2, 7.7, 0.493, 0.0629, 2.058, 1200, 170, 0.16]
    cv = zeros(4, 1);

    % steady-state values
    cv(1) = 1 ./ (1 + exp(-(V+p(1))./p(2)));
    cv(2) = 1 ./ (1 + exp((V+p(3))./p(4)));
    
    % time constants
    cv(3) = p(5).*exp(-p(6).*V) + p(7);
    cv(4) = p(8) - (p(9))./(1 + exp((V+p(3))./p(4)));
end
