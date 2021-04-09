function [t, A] = analytic_form(holdV, holdT, P1, P1T)
    % constants
    ato0 = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    ito0 = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    aKslow0 = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    iKslow0 = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    Gto = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    GKslow = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
    Ek = -91.1;

    % time space
    time_step = 0.001;
    tH = 0:time_step:holdT;
    tP1 = 0:time_step:(P1T-holdT);
    t = [tH, tP1+(holdT+time_step)];

    hold_idx = length(tH);

    % matrix contains current traces
    A = zeros(length(tH)+length(tP1), 2);

    % Ito; holding potential
    CS_hold = ch_stat(holdV);
    ato = CS_hold(1) - (CS_hold(1) - ato0).*exp(-(tH./CS_hold(2)));
    ito = CS_hold(3) - (CS_hold(3) - ito0).*exp(-(tH./CS_hold(4)));
    A(1:hold_idx, 1) = Gto.*ato.^3.*ito.*(holdV-Ek);

    % IKslow; holding potential
    aKslow = CS_hold(5) - (CS_hold(5) - aKslow0).*exp(-(tH./CS_hold(6)));
    iKslow = CS_hold(7) - (CS_hold(7) - iKslow0).*exp(-(tH./CS_hold(8)));
    A(1:hold_idx, 2) = GKslow.*aKslow.*iKslow.*(holdV-Ek);
    
    % Ito; P1
    CS_p1 = ch_stat(P1);
    ato = CS_p1(1) - (CS_p1(1) - ato0).*exp(-(tP1./CS_p1(2)));
    ito = CS_p1(3) - (CS_p1(3) - ito0).*exp(-(tP1./CS_p1(4)));
    A(hold_idx+1:end, 1) = Gto.*ato.^3.*ito.*(P1-Ek);

    % IKslow; P1
    aKslow = CS_p1(5) - (CS_p1(5) - aKslow0).*exp(-(tP1./CS_p1(6)));
    iKslow = CS_p1(7) - (CS_p1(7) - iKslow0).*exp(-(tP1./CS_p1(8)));
    A(hold_idx+1:end, 2) = GKslow.*aKslow.*iKslow.*(P1-Ek);
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
end
