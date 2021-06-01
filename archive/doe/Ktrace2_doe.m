function [res_to, res_kslow] = Ktrace2_doe(param, holdV, P1, time_space)
    % using Rasmusson K+ current models

    % constants
    Gto = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    GKslow = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
    ato0 = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    ito0 = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    aKslow0 = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    iKslow0 = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    Ek = -91.1;  % resting potential

    % time space
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};
    hold_idx = length(tH);

    % tuning parameters (decision variables)
    param_Ito = param{1};
    param_IKslow1 = param{2};

    ito_trc = zeros(length(t), 1);
    ikslow_trc = zeros(length(t), 1);

    % Ito; holding potential
    cs_Ito_hold = cs_Ito(param_Ito, holdV);
    act_Ito_hold = cs_Ito_hold(1) - (cs_Ito_hold(1) - ato0).*exp(-(tH./cs_Ito_hold(2)));
    inact_Ito_hold = cs_Ito_hold(3) - (cs_Ito_hold(3) - ito0).*exp(-(tH./cs_Ito_hold(4)));
    ito_trc(1:hold_idx) = Gto.*(act_Ito_hold.^3).*(inact_Ito_hold).*(holdV - Ek);

    % IKslow1; holding potential
    cs_IKslow1_hold = cs_IKslow1(param_IKslow1, holdV);
    act_IKslow1_hold = cs_IKslow1_hold(1) - (cs_IKslow1_hold(1) - aKslow0).*exp(-(tH./cs_IKslow1_hold(2)));
    inact_IKslow1_hold = cs_IKslow1_hold(3) - (cs_IKslow1_hold(3) - iKslow0).*exp(-(tH./cs_IKslow1_hold(4)));
    ikslow_trc(1:hold_idx) = GKslow.*(act_IKslow1_hold).*(inact_IKslow1_hold).*(holdV - Ek);

    % Ito; P1
    cs_Ito_P1 = cs_Ito(param_Ito, P1);
    act_Ito_P1 = cs_Ito_P1(1) - (cs_Ito_P1(1) - ato0).*exp(-(tP1_adj./cs_Ito_P1(2)));
    inact_Ito_P1 = cs_Ito_P1(3) - (cs_Ito_P1(3) - ito0).*exp(-(tP1_adj./cs_Ito_P1(4)));
    ito_trc(hold_idx+1:end) = Gto.*(act_Ito_P1.^3).*(inact_Ito_P1).*(P1 - Ek);

    % IKslow1; P1
    cs_IKslow1_P1 = cs_IKslow1(param_IKslow1, P1);
    act_IKslow1_P1 = cs_IKslow1_P1(1) - (cs_IKslow1_P1(1) - ato0).*exp(-(tP1_adj./cs_IKslow1_P1(2)));
    inact_IKslow1_P1 = cs_IKslow1_P1(3) - (cs_IKslow1_P1(3) - ito0).*exp(-(tP1_adj./cs_IKslow1_P1(4)));
    ikslow_trc(hold_idx+1:end) = GKslow.*(act_IKslow1_P1).*(inact_IKslow1_P1).*(P1 - Ek);

    % responses of interest
    res_to = zeros(1, 4);
    res_kslow = zeros(1, 4);

    % peak
    peak_to = max(ito_trc);
    peak_kslow = max(ikslow_trc);
    res_to(1) = peak_to;
    res_kslow(1) = peak_kslow;

    % time constant
    [~, tau_idx_to] = min(abs(peak_to*exp(-1) - ito_trc));
    [~, tau_idx_kslow] = min(abs(peak_kslow*exp(-1) - ikslow_trc));
    res_to(2) = t(tau_idx_to);
    res_kslow(2) = t(tau_idx_kslow);
    
    % SSA
    res_to(3) = cs_Ito_P1(1);
    res_kslow(3) = cs_IKslow1_P1(1);

    % SSI
    res_to(4) = cs_Ito_P1(3);
    res_kslow(4) = cs_IKslow1_P1(3);
end

function [cs] = cs_Ito(p, V)
    % characteristic statistics of Ito
    % cs(1): steady-state activation
    % cs(2): time constant of activation
    % cs(3): steady-state inactivation
    % cs(4): time constant of inactivation
    
    % p0 = [30.0, 13.5, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335];
    cs = zeros(1, 4);
    alphaA = p(5).*exp(p(6).*(V+p(1)));
    betaA = p(7).*exp(-p(8).*(V+p(1)));
    alphaI = (p(9).*exp(-(V+p(2))./p(4))) ./ (p(10).*exp(-(V+p(3))./p(4)) + 1);
    betaI = (p(11).*exp((V+p(3))./p(4))) ./ (p(12).*exp((V+p(3))./p(4)) + 1);

    cs(1) = alphaA/(alphaA+betaA);
    cs(2) = 1/(alphaA+betaA);
    cs(3) = alphaI/(alphaI+betaI);
    cs(4) = 1/(alphaI+betaI);
end

function [cs] = cs_IKslow1(p, V)
    % characteristic statistics of IKslow
    % cs(1): steady-state activation
    % cs(2): time constant of activation
    % cs(3): steady-state inactivation
    % cs(4): time constant of inactivation
    
    % p0 = [22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 1200.0, 170.0, 45.2, 5.7];
    cs = zeros(1, 4);
    cs(1) = 1 ./ (1 + exp(-(V+p(1))./p(2)));
    cs(2) = p(3).*exp(-p(4).*V) + p(5);
    cs(3) = 1 ./ (1 + exp((V+p(6))./p(7)));
    cs(4) = p(8) - (p(9))./(1 + exp((V+p(10))./p(11)));
end
