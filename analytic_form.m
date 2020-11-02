function [t, A] = analytic_form(holding_p, holding_t, P1, P1_t)
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
    A = zeros(length(t), 2);

    % Ito; holding potential
    CS = ch_stat(holding_p);
    ato = CS(1) - (CS(1) - ato0).*exp(-(t(1:holding_idx)./CS(2)));
    ito = CS(3) - (CS(3) - ito0).*exp(-(t(1:holding_idx)./CS(4)));
    A(1:holding_idx, 1) = Gto.*ato.^3.*ito.*(holding_p-Ek);

    % IKslow; holding potential
    aKslow = CS(5) - (CS(5) - aKslow0).*exp(-(t(1:holding_idx)./CS(6)));
    iKslow = CS(7) - (CS(7) - iKslow0).*exp(-(t(1:holding_idx)./CS(8)));
    A(1:holding_idx, 2) = GKslow.*aKslow.*iKslow.*(holding_p-Ek);
    
    % Ito; P1
    CS = ch_stat(P1);
    ato = CS(1) - (CS(1) - ato0).*exp(-(t(holding_idx+1:end)./CS(2)));
    ito = CS(3) - (CS(3) - ito0).*exp(-(t(holding_idx+1:end)./CS(4)));
    A(holding_idx+1:end, 1) = Gto.*ato.^3.*ito.*(P1-Ek);

    % IKslow; P1
    aKslow = CS(5) - (CS(5) - aKslow0).*exp(-(t(holding_idx+1:end)./CS(6)));
    iKslow = CS(7) - (CS(7) - iKslow0).*exp(-(t(holding_idx+1:end)./CS(8)));
    A(holding_idx+1:end, 2) = GKslow.*aKslow.*iKslow.*(P1-Ek);
end

function [CS] = ch_stat(V)
    CS = zeros(1, 8);
    
    % Ito transition rates
    alpha_a = 0.18064.*exp(0.03577.*(V+30.0));
    beta_a = 0.3956.*exp(-0.06237.*(V+30.0));
    alpha_i = (0.000152.*exp(-(V+13.5)./7.0))./( 0.0067083.*exp(-(V+33.5)./7.0)+1); % typo? 0.067083 in paper
    beta_i = (0.00095.*exp((V+33.5)./7.0))./( 0.051335.*exp((V+33.5)./7.0)+1);
    
    % Ito steady-state activation
    CS(1) = alpha_a/(alpha_a+beta_a);
    % Ito activation tau
    CS(2) = 1/(alpha_a+beta_a);
    % Ito steady-state inactivation
    CS(3) = alpha_i/(alpha_i+beta_i);
    % Ito inactivation tau
    CS(4) = 1/(alpha_i+beta_i);

    % IKslow steady-state activation
    CS(5) = 1./(1+exp(-(V+22.5)./7.7));
    % IKslow activation tau
    CS(6) = 0.493.*exp(-0.0629.*V)+2.058;
    % IKslow steady-state inactivation
    CS(7) = 1./(1+exp((V+45.2)./5.7));
    % IKslow inactivation tau
    CS(8) = 1200.0 - 170.0./(1+exp((V+45.2)./5.7));

    % steady-state activation
    % activation tau
    % steady-state inactivation
    % inactivation tau
end
