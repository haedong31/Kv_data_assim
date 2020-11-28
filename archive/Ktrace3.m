function [Ktrc] = Ktrace3(param, holdV, P1, time_space)
    % Dongping's K+ current models
    % fully parameterized
    % Ito; 20 parameters
    % IKslow; 12 parameters

    % constants
    % Gto = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    % GKslow = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
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
    param_IKslow2 = param{3};
    param_Iss = param{4};
    paramG = param{5};

    Ktrc = zeros(length(t), 4);

    % Ito; holding potential
    cs_Ito_hold = cs_Ito(param_Ito, holdV);
    act_Ito_hold = cs_Ito_hold(1) - (cs_Ito_hold(1) - ato0).*exp(-(tH./cs_Ito_hold(2)));
    inact_Ito_hold = cs_Ito_hold(3) - (cs_Ito_hold(3) - ito0).*exp(-(tH./cs_Ito_hold(4)));
    Ktrc(1:hold_idx, 1) = paramG(1).*(act_Ito_hold.^3).*(inact_Ito_hold).*(holdV - Ek);

    % IKslow1; holding potential
    cs_IKslow1_hold = cs_IKslow1(param_IKslow1, holdV);
    act_IKslow1_hold = cs_IKslow1_hold(1) - (cs_IKslow1_hold(1) - aKslow0).*exp(-(tH./cs_IKslow1_hold(2)));
    inact_IKslow1_hold = cs_IKslow1_hold(3) - (cs_IKslow1_hold(3) - iKslow0).*exp(-(tH./cs_IKslow1_hold(4)));
    Ktrc(1:hold_idx, 2) = paramG(2).*(act_IKslow1_hold).*(inact_IKslow1_hold).*(holdV - Ek);

    % IKslow2; holding potential
    cs_IKslow2_hold = cs_IKslow2(param_IKslow2, holdV);
    act_IKslow2_hold = cs_IKslow2_hold(1) - (cs_IKslow2_hold(1) - aKslow0).*exp(-(tH./cs_IKslow2_hold(2)));
    inact_IKslow2_hold = cs_IKslow2_hold(3) - (cs_IKslow2_hold(3) - iKslow0).*exp(-(tH./cs_IKslow2_hold(4)));
    Ktrc(1:hold_idx, 3) = paramG(3).*(act_IKslow2_hold).*(inact_IKslow2_hold).*(holdV - Ek);
    
    % Ito; P1
    cs_Ito_P1 = cs_Ito(param_Ito, P1);
    act_Ito_P1 = cs_Ito_P1(1) - (cs_Ito_P1(1) - ato0).*exp(-(tP1_adj./cs_Ito_P1(2)));
    inact_Ito_P1 = cs_Ito_P1(3) - (cs_Ito_P1(3) - ito0).*exp(-(tP1_adj./cs_Ito_P1(4)));
    Ktrc(hold_idx+1:end, 1) = paramG(1).*(act_Ito_P1.^3).*(inact_Ito_P1).*(P1 - Ek);

    % IKslow1; P1
    cs_IKslow1_P1 = cs_IKslow1(param_IKslow1, P1);
    act_IKslow1_P1 = cs_IKslow1_P1(1) - (cs_IKslow1_P1(1) - ato0).*exp(-(tP1_adj./cs_IKslow1_P1(2)));
    inact_IKslow1_P1 = cs_IKslow1_P1(3) - (cs_IKslow1_P1(3) - ito0).*exp(-(tP1_adj./cs_IKslow1_P1(4)));
    Ktrc(hold_idx+1:end, 2) = paramG(2).*(act_IKslow1_P1).*(inact_IKslow1_P1).*(P1 - Ek);

    % Ikslow2; P1
    cs_IKslow2_P1 = cs_IKslow2(param_IKslow2, P1);
    act_IKslow2_P1 = cs_IKslow2_P1(1) - (cs_IKslow2_P1(1) - ato0).*exp(-(tP1_adj./cs_IKslow2_P1(2)));
    inact_IKslow2_P1 = cs_IKslow2_P1(3) - (cs_IKslow2_P1(3) - ito0).*exp(-(tP1_adj./cs_IKslow2_P1(4)));
    Ktrc(hold_idx+1:end, 3) = paramG(3).*(act_IKslow2_P1).*(inact_IKslow2_P1).*(P1 - Ek);

    % Iss
    Iss = zeros(length(t), 1);
    Iss(hold_idx+1:end) = param_Iss;
    Ktrc(:, 4) = Iss;
end

function [cs] = cs_Ito(p, V)
    % characteristic statistics of Ito
    % cs(1): steady-state activation
    % cs(2): time constant of activation
    % cs(3): steady-state inactivation
    % cs(4): time constant of inactivation

    cs = zeros(1, 4);
    alphaA = p(5).*exp(p(6).*(V+p(7)));
    betaA = p(8).*exp(p(9).*(V+p(10)));
    alphaI = (p(11).*exp(-(V+p(12))./p(13))) ./ (p(14).*exp(-(V+p(15))./p(16)) + 1);
    betaI = (p(17).*exp((V+p(18))./p(19))) ./ (p(20).*exp((V+p(18))./p(19)) + 1);

    cs(1) = 1 ./ (1+exp(-(V-p(1))./p(2)));
    cs(2) = 1/(alphaA+betaA);
    cs(3) = 1 ./ (1+exp(-(V-p(3))./p(4)));
    cs(4) = 1/(alphaI+betaI);
end

function [cs] = cs_IKslow1(p, V)
    % characteristic statistics of IKslow
    % cs(1): steady-state activation
    % cs(2): time constant of activation
    % cs(3): steady-state inactivation
    % cs(4): time constant of inactivation
    
    cs = zeros(1, 4);
    cs(1) = 1 ./ (1 + exp(-(V+p(1))./p(2)));
    cs(2) = p(7).*exp(p(8).*V) + p(9);
    cs(3) = (0.21./(1+exp(-(V+p(3))./p(4)))) + (0.79./(1+exp(-(V+p(5))./p(6))));
    cs(4) = 500.0 + (p(10))./(1 + exp((V+p(11))./p(12)));
end

function [cs] = cs_IKslow2(p, V)
    % characteristic statistics of IKslow
    % cs(1): steady-state activation
    % cs(2): time constant of activation
    % cs(3): steady-state inactivation
    % cs(4): time constant of inactivation
    
    cs = zeros(1, 4);
    cs(1) = 1 ./ (1 + exp(-(V+p(1))./p(2)));
    cs(2) = p(7).*exp(p(8).*V) + p(9);
    cs(3) = (0.21./(1+exp(-(V+p(3))./p(4)))) + (0.79./(1+exp(-(V+p(5))./p(6))));
    cs(4) = 500.0 + (p(10))./(1 + exp((V+p(11))./p(12)));
end
