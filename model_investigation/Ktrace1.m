function [A] = Ktrace1(p, holdV, hold_idx, P1, t)
    % constants
    % Gto = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    % GKslow = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
    ato0 = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    ito0 = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    aKslow0 = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    iKslow0 = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    Ek = -91.1;

    % to avoid declaring global variables (MATLAB recommandation)
    fn_input = cell(1, 5);
    fn_input{2} = holdV;
    fn_input{3} = P1;
    fn_input{4} = t;
    fn_input{5} = hold_idx;

    A = zeros(length(t), 4);

    fn_input{1} = p(1:9);
    Ito = generic_Ktrace(fn_input, ato0, ito0, Ek);

    fn_input{1} = p(10:18);
    IKslow1 = generic_Ktrace(fn_input, aKslow0, iKslow0, Ek);

    fn_input{1} = p(19:27);
    IKslow2 = generic_Ktrace(fn_input, aKslow0, iKslow0, Ek);

    Iss = zeros(length(t), 1);
    Iss(hold_idx+1:end) = p(28);

    A(:, 1) = Ito;
    A(:, 2) = IKslow1;
    A(:, 3) = IKslow2;
    A(:, 4) = Iss;
end

function [Itrc] = generic_Ktrace(fn_input, a0, i0, Ek)
    p = fn_input{1};
    holdV = fn_input{2};
    P1 = fn_input{3};
    t = fn_input{4};
    hold_idx = fn_input{5};

    Itrc = zeros(length(t), 1);

    % during holding potential
    CS = ch_stat(p, holdV);
    act = CS(1) - (CS(1) - a0).*exp(-(t(1:hold_idx)./CS(2)));
    inact = CS(3) - (CS(3) - i0).*exp(-(t(1:hold_idx)./CS(4)));
    Itrc(1:hold_idx) = p(9).*act.*inact.*(holdV - Ek);

    % during P1 pulse
    CS = ch_stat(p, P1);
    act = CS(1) - (CS(1) - a0).*exp(-(t(hold_idx+1:end)./CS(2)));
    inact = CS(3) - (CS(3) - i0).*exp(-(t(hold_idx+1:end)./CS(4)));
    Itrc(hold_idx+1:end) = p(9).*act.*inact.*(P1 - Ek);
end

function [CS] = ch_stat(p, V)
    % compute characteristic statistics:
    % 1. steady-state activation
    % 2. time constant of activation
    % 3. steady-state inactivation
    % 4. time constant of inactivation

    CS = zeros(1, 4);

    alpha_a = p(1).*exp(p(2).*V);
    beta_a = p(3).*exp(-p(4).*V);
    alpha_i = p(5).*exp(p(6).*V);
    beta_i = p(7).*exp(-p(8).*V);

    CS(1) = alpha_a/(alpha_a+beta_a);
    CS(2) = 1/(alpha_a+beta_a);
    CS(3) = alpha_i/(alpha_i+beta_i);
    CS(4) = 1/(alpha_i+beta_i);
end
