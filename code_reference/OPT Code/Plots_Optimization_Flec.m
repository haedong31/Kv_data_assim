
% Steady State Availability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WT_SSA = load('Control_SSA_param.txt');

WT_Sim_C_x = WT_SSA(:,3);
WT_Sim_C_y = WT_SSA(:,4);


WT_SSA2 = load('testSSA_param.txt');

WT_Sim_F_x = WT_SSA2(:,3);
WT_Sim_F_y = WT_SSA2(:,4);


SSA_EXP = load('WT_SSA.txt');
WT_x = SSA_EXP(:, 1);
WT_C = SSA_EXP(:,2);
WT_F = SSA_EXP(:,3);
WT_L = SSA_EXP(:,4);

figure (1);
%axes('YTick', [0 0.25 0.50 0.75 1.0],'XTick',[-125 -100 -75 -50 -25],'LineWidth',6);
hold on;
plot (WT_Sim_C_x, WT_Sim_C_y, '--r', 'Linewidth', 3 );
plot (WT_x, WT_C, 'squarek', 'MarkerSize', 16, 'MarkerFaceColor', 'k');
plot(WT_Sim_F_x, WT_Sim_F_y, '-k', 'Linewidth', 3);
plot (WT_x, WT_F, 'ok', 'MarkerSize', 16, 'MarkerFaceColor', 'k');
hold off;

% Frequency Dependent UDB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDUDB1_Sim = load('testFDUDB1_param.txt');
FDUDB1_Sim_x = FDUDB1_Sim(:,1);
FDUDB1_Sim_y = FDUDB1_Sim(:,2);
FDUDB2_Sim = load('testFDUDB2_param.txt');
FDUDB2_Sim_x = FDUDB2_Sim(:,1);
FDUDB2_Sim_y = FDUDB2_Sim(:,2);

  FDUDB_Exp_x = [1; 2; 5; 10];
  FDUDB_Exp_y = [32.57; 38.08; 50.92; 67.28];  %%% WT with 10uM Flec 
  %FDUDB_Exp_y = [10.28; 19.06; 50.89; 75.83];  %%% WT with 300uM Lido

figure (2);
axes('XScale','log');
hold on;
plot (FDUDB1_Sim_x, FDUDB1_Sim_y, '-k', 'Linewidth', 3);
plot (FDUDB2_Sim_x, FDUDB2_Sim_y, '-k', 'Linewidth', 3);
plot (FDUDB_Exp_x, FDUDB_Exp_y, 'ok', 'MarkerSize', 16, 'MarkerFaceColor', 'k');
hold off; 

% Recovery from UDB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RUDB_Control = load('Control_RUDB100_param.txt');
RUDB_C_x = RUDB_Control(:,1);
RUDB_C_y = RUDB_Control(:,2);

RUDB_Sim = load('testRUDB100_param.txt');
RUDB_F_x = RUDB_Sim(:,1);
RUDB_F_y = RUDB_Sim(:,2);


RUDB_Exp = load('WT_rec_exp.txt');
RUDB_Exp_x = RUDB_Exp(:,1);
RUDB_Exp_C = RUDB_Exp(:,2);
RUDB_Exp_F = RUDB_Exp(:,3);

figure (3);
axes('XScale','log');
hold on;

plot (RUDB_C_x, RUDB_C_y, '--r', 'Linewidth', 3);
plot (RUDB_Exp_x, RUDB_Exp_C, 'squarek', 'MarkerSize', 16, 'MarkerFaceColor', 'k');

plot (RUDB_F_x, RUDB_F_y, '-k', 'Linewidth', 3);
plot (RUDB_Exp_x, RUDB_Exp_F, 'ok', 'MarkerSize', 16, 'MarkerFaceColor', 'k');
hold off;

% UDB Isotherm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UDB_Sim = load('testBlock_param.txt');
UDB_x = UDB_Sim(:,1);
UDB_y = UDB_Sim(:,2);

Drug = [0; 1; 3; 10; 30; 100; 300; 1000];
WTFlec = [100; 87; 79; 54; 25; 13; 1; 1];

figure (4);
axes('XScale','log');
hold on;
plot (UDB_x, UDB_y, '-k', 'Linewidth', 3);
plot (Drug, WTFlec, 'ok', 'MarkerSize', 16, 'MarkerFaceColor', 'k');
hold off; 

% TB Isotherm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TB_x = UDB_Sim(:,1);
TB_y = UDB_Sim(:,3);

Drug = [0; 1; 3; 10; 30; 100; 300];
WT_Flec_TB = [100; 99.3; 95.7; 91.0; 78.7; 33.2; 6.86]; %%%Abriel
WT_Flec_TB2 = [100; 100; 98.9; 94.2; 79.5; 55.3; 31.6]; %%%Nagatomo


figure (5);
axes('XScale','log');
hold on;
plot (TB_x, TB_y, '-k', 'Linewidth', 3);
plot (Drug, WT_Flec_TB, 'ok', 'MarkerSize', 16, 'MarkerFaceColor', 'k');
plot (Drug, WT_Flec_TB2, 'squarek', 'MarkerSize', 16, 'MarkerFaceColor', 'k');
hold off; 

% CELL UV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CELL_Sim = load('Cell_param.txt');
CELL_x = CELL_Sim(:,16);
CELL_y = CELL_Sim(:,17);

figure (6);
hold on;
plot (CELL_x, CELL_y, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold off; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
