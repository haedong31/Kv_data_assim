function [current_trc] = IKr(param, holdV, P1, time_space)
    % Bondarenko IKr (A96)
    % 15 parameters + 2 ionic concentrations + 1 conductance

    % constants
    ck10 == 0.992513e-3; % CK1; mERG channel closed state
    ck20 == 0.641229e-3; % CK2; mERG channel closed state
    ok0 == 0.175298e-3; % OK; mERG channel open state
    ik0 == 0.319129e-4; % IK; mERG channel inactivated state
    
    kf = param(14); % 0.023761; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    kb = param(15); % 0.036778; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    Ko = param(16); % 5400; Exracellular K+ concentration:uM
    Nao = param(17); % 140000; Exracellular Na+ concentration:uM
    gmax = param(18); % 0.078; Maximum rapid delayed-rectifier K+ current conductance:mS/uF


end

function [cv] = cv(p, V)
    % p0 = [0.022348, 0.01176, 0.047002, 0.0631, 
    % 0.013733, 0.038198, 6.89e-05, 0.04178, 
    % 5, 0.090821, 0.023391, 0.006497, 0.03268]
    
    cv = zeros(6, 1);
    % A102
    cv(1) = p(1).*exp(p(2).*V);
    % A103
    cv(2) = p(3).*exp(-p(4).*V);
    % A104
    cv(3) = p(5).*exp(p(6).*V);
    % A105
    cv(4) = p(7).*exp(-p(8).*V);
    % A106
    cv(5) = p(10).*exp(p(11).*(V+p(9)));
    % A107
    cv(6) = p(12).*exp(-p(13).*(V+p(9)));
end