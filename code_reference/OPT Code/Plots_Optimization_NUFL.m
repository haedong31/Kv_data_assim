


FDUDB_Sim = load('testFDUDB_NUFL_param.txt');
FDUDB_Sim_y = FDUDB_Sim(:,2);

x = [10, 100];
FDUDB_Exp_y = [23; 56];  %%% WT 10Hz: 10uM NUFL, then 100uM NUFL

figure (1)
hold on;
plot (x, FDUDB_Exp_y, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot (x, FDUDB_Sim_y, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;

%  **********************************
RUDB = load('testRUDB100_NUFL_param.txt');
RUDB_x = RUDB(:,1);
RUDB_y = RUDB(:,2);
t1 = [1.005; 5.005; 10.005; 25.005; 50.005; 100.005; 150.005; 250.005; 500.005; 1000; 2000; 3000; 5000; 7500; 10000];

A1_EXP = 0.5551; T1_EXP = 363.2; A2_EXP = 0.2915; T2_EXP = 35.21; C1_EXP = 0.8643; %% Taken from Common Molecular fig 3A
 
Y_EXP = C1_EXP - A1_EXP*exp(-t1/T1_EXP) - A2_EXP*exp(-t1/T2_EXP);
 
figure (2)
axes('XScale','log');
hold on;
plot (t1, Y_EXP, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot (RUDB_x, RUDB_y,'-k', 'Linewidth', 3);
hold off;

