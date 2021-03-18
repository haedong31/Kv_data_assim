function [t, I, TR] = analytic_form(holding_p, holding_t, P1, P1_t)
    % constants
    ato0 = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    ito0 = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    aKslow0 = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    iKslow0 = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    Gto = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    GKslow = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
    Ek = -91.1;

    % time space
    t = 0:0.001:P1_t;
    [~, holding_idx] = min(abs(t - holding_t));

    % matrix contains current traces
    I = zeros(length(t), 2);

    TR = zeros(8, 2);
    
    % Ito; holding potential
    TR_hold = tran_rate(holding_p);
    TR(:,1) = TR_hold;
    ass = TR_hold(1)/(TR_hold(1) + TR_hold(2));
    atau = 1/(TR_hold(1) + TR_hold(2));
    iss = TR_hold(3)/(TR_hold(3) + TR_hold(4));
    itau = 1/(TR_hold(3) + TR_hold(4));

    ato = ass - (ass - ato0).*exp(-(t(1:holding_idx)./atau));
    ito = iss - (iss - ito0).*exp(-(t(1:holding_idx)./itau));
    I(1:holding_idx, 1) = Gto.*ato.^3.*ito.*(holding_p-Ek);

    % IKslow; holding potential
    aKslow = TR_hold(5) - (TR_hold(5) - aKslow0).*exp(-(t(1:holding_idx)./TR_hold(6)));
    iKslow = TR_hold(7) - (TR_hold(7) - iKslow0).*exp(-(t(1:holding_idx)./TR_hold(8)));
    I(1:holding_idx, 2) = GKslow.*aKslow.*iKslow.*(holding_p-Ek);
    
    % Ito; P1
    TR_P1 = tran_rate(P1);
    TR(:, 2) = TR_P1;
    ass = TR_P1(1)/(TR_P1(1) + TR_P1(2));
    atau = 1/(TR_P1(1) + TR_P1(2));
    iss = TR_P1(3)/(TR_P1(3) + TR_P1(4));
    itau = 1/(TR_P1(3) + TR_P1(4));

    ato = ass - (ass - ato0).*exp(-(t(holding_idx+1:end)./atau));
    ito = iss - (iss - ito0).*exp(-(t(holding_idx+1:end)./itau));
    I(holding_idx+1:end, 1) = Gto.*ato.^3.*ito.*(P1-Ek);

    % IKslow; P1
    aKslow = TR_P1(5) - (TR_P1(5) - aKslow0).*exp(-(t(holding_idx+1:end)./TR_P1(6)));
    iKslow = TR_P1(7) - (TR_P1(7) - iKslow0).*exp(-(t(holding_idx+1:end)./TR_P1(8)));
    I(holding_idx+1:end, 2) = GKslow.*aKslow.*iKslow.*(P1-Ek);
end

function [TR] = tran_rate(V)
    TR = zeros(1, 8);
    
    % Ito transition rates
    TR(1) = 0.18064.*exp(0.03577.*(V+30.0));
    TR(2) = 0.3956.*exp(-0.06237.*(V+30.0));
    TR(3) = (0.000152.*exp(-(V+13.5)./7.0))./( 0.0067083.*exp(-(V+33.5)./7.0)+1); % typo? 0.067083 in paper
    TR(4) = (0.00095.*exp((V+33.5)./7.0))./( 0.051335.*exp((V+33.5)./7.0)+1);

    % IKslow steady-state activation
    TR(5) = 1./(1+exp(-(V+22.5)./7.7));
    % IKslow activation tau
    TR(6) = 0.493.*exp(-0.0629.*V)+2.058;
    % IKslow steady-state inactivation
    TR(7) = 1./(1+exp((V+45.2)./5.7));
    % IKslow inactivation tau
    TR(8) = 1200.0 - 170.0./(1+exp((V+45.2)./5.7));

    % steady-state activation
    % activation tau
    % steady-state inactivation
    % inactivation tau
end
