
% Steady State Availability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%V= -140:-40;
%[Vhalf K] = simfit_SSA;
%Y_SIM = 1./(1+exp((V - Vhalf)/K));
SSA = load('testSSA_param.txt');
SSA_x = SSA(:,3);
SSA_y = SSA(:,4);

SSA_EXP = load('WT_SSA.txt');
WT_x = SSA_EXP(:, 1);
WT_C = SSA_EXP(:,2);
WT_F = SSA_EXP(:,3);
WT_L = SSA_EXP(:,4);

figure (1);
hold on;
plot (SSA_x, SSA_y, '-k', 'Linewidth', 3 );
plot (WT_x, WT_C, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold off;

% Recovery from UDB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RUDB_Sim = load('testRUDB100_param.txt');
RUDB_x = RUDB_Sim(:,1);
RUDB_y = RUDB_Sim(:,2);

RUDB_Exp = load('WT_rec_exp.txt');
%RUDB_Exp = load('WTLido.txt');
RUDB_Exp_x = RUDB_Exp(:,1);
RUDB_Exp_y = RUDB_Exp(:,2);

figure (2);
axes('XScale','log');
hold on;
plot (RUDB_x, RUDB_y, '-k', 'Linewidth', 3);
plot (RUDB_Exp_x, RUDB_Exp_y, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold off; 

% Tau to Inactivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Sim = load('testTau_param.txt');
Tau_Sim_x = Tau_Sim(:,2);
Tau_Sim_y = Tau_Sim(:,7);

Tau_Exp = load('Tau50.txt');
Tau_Exp_x = Tau_Exp(:,1);
Tau_Exp_y = Tau_Exp(:,2);

figure (3);
hold on;
plot (Tau_Sim_x, Tau_Sim_y, '-k', 'Linewidth', 3);
plot (Tau_Exp_x, Tau_Exp_y, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold off; 

% Recovery from Inactivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RFI_Sim = load('testRFI100_param.txt');
RFI_Sim_x = RFI_Sim(:,1);
RFI_Sim_y = RFI_Sim(:,2);

RFI_Exp = load('WT_RFI100.txt');
RFI_Exp_x = RFI_Exp(:,1);
RFI_Exp_y = RFI_Exp(:,2);

figure (4);
axes('XScale','log');
hold on;
plot (RFI_Sim_x, RFI_Sim_y, '-k', 'Linewidth', 3);
plot (RFI_Exp_x, RFI_Exp_y, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold off; 

% Steady State Activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Act_Sim = load('testAct_param.txt');
Act_Sim_x = Act_Sim(:,2);
Act_Sim_y = Act_Sim(:,6);

Act_Exp = load('WT_Act.txt');
Act_Exp_x = Act_Exp(:,1);
Act_Exp_y = Act_Exp(:,2);

figure (5);
hold on;
plot (Act_Sim_x, Act_Sim_y, '-k', 'Linewidth', 3);
plot (Act_Exp_x, Act_Exp_y, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold off; 



