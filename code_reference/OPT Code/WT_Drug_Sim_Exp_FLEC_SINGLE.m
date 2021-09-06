function Total_Error = WT_Drug_Sim_Exp_FLEC_SINGLE(Inputs)

i = 1;

% %Neutral flecainide protocols
% parfor i = 1:2;
% if i==1  
%     main_RUDB(Inputs); % Make sure portion = 100% neutral, and it's 100uM
% end 
% if i==2 
%     main_FDUDB_NUFL(Inputs);
% end
% end;

 %Charged flecainide protocols
% parfor i = 1:6;
% if i==1  
%     main_SSA(Inputs);
% end 
% if i==2 
%     main_RUDB(Inputs);
% end
% if i==3  
%     main_BLOCK(Inputs);
% end
% if i==4  
%     main_FDUDB1(Inputs);
% end
% if i==5 
%     main_FDUDB2(Inputs);
% end
% if i ==6
%     main_CELL(Inputs);
% end;
% end;

main_SSA(Inputs);

% 100uM Neutral Flecainide RUDB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUDB = load('testRUDB100_NUFL_param.txt');
% RUDB_x = RUDB(:,1);
% RUDB_y = RUDB(:,2);
% t1 = [1.005; 5.005; 10.005; 25.005; 50.005; 100.005; 150.005; 250.005; 500.005; 1000; 2000; 3000; 5000; 7500; 10000];
% 
% A1_EXP = 0.5551; T1_EXP = 363.2; A2_EXP = 0.2915; T2_EXP = 35.21; C1_EXP = 0.8643; %% Taken from Common Molecular fig 3A
%  
% Y_EXP = C1_EXP - A1_EXP*exp(-t1/T1_EXP) - A2_EXP*exp(-t1/T2_EXP);
%  
% Err_RUDB_NUFL = sum((RUDB_y*100 - Y_EXP*100).^2)   / 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% NUFL FDUDB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDUDB_Sim = load('testFDUDB_NUFL_param.txt');
% FDUDB_Sim_y = FDUDB_Sim(:,2);
% 
% FDUDB_Exp_y = [23; 56];  %%% WT 10Hz: 10uM NUFL, then 100uM NUFL
% 
% Err_FDUDB_NUFL = sum((FDUDB_Sim_y - FDUDB_Exp_y).^2) / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R UDB -100mV  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUDB = load('testRUDB100_param.txt');
% RUDB_x = RUDB(:,1);
% RUDB_y = RUDB(:,2);
% t1 = [1.005; 5.005; 10.005; 25.005; 50.005; 100.005; 150.005; 250.005; 500.005; 1000; 2000; 3000; 5000; 7500; 10000];
% 
% %A1_EXP = 0.7283; T1_EXP = 6.805; A2_EXP = .2883; T2_EXP = 291.2; C1_EXP = 0.9826;  %%% WT Drug Free
% %A1_EXP = 0.31; T1_EXP = 103; A2_EXP = 0.17; T2_EXP = 2530; C1_EXP = 0.45;  %%% WT 10uM Flec (Liu, 2002)
% A1_EXP = 0.2547; T1_EXP = 96.9; A2_EXP = 0.1815; T2_EXP = 2162; C1_EXP = 0.443;
% 
% Y_EXP = C1_EXP - A1_EXP*exp(-t1/T1_EXP) - A2_EXP*exp(-t1/T2_EXP);
% 
% Err_RUDB100 = sum((RUDB_y*100 - Y_EXP*100).^2) / 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% SSA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SSA_Sim = load('testSSA_param.txt');
SSA_Sim_y = SSA_Sim(:,4);

V = [-125:5:-30]'; 
%%%Vhalf_EXP = -70.46;   %%% Drug Free
Vhalf_EXP = -70.46 - 2.3; %%% 10uM Flecainide
K_EXP = 6.34;
 
Y_EXP = 1./(1+exp((V - Vhalf_EXP)/K_EXP));
 
Err_SSA = sum((SSA_Sim_y*100 - Y_EXP*100).^2) / 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FDUDB1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FDUDB1 is for 1 Hz...split up 1, 2, 5, 10Hz for faster processing; FDUDB2 = 2, 5, 10 hz
% FDUDB_Sim = load('testFDUDB1_param.txt');
% FDUDB_Sim_y = FDUDB_Sim(:,2);
% 
% %%%FDUDB_Exp_x = [1];
% FDUDB_Exp_y = [32.57];  %%% WT with 10uM Flec (Liu, 2002)    %%32.57 for 1H
% 
% Err_FDUDB1 = sum((FDUDB_Sim_y - FDUDB_Exp_y).^2) / 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FDUDB2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDUDB_Sim = load('testFDUDB2_param.txt');
% FDUDB_Sim_y = FDUDB_Sim(:,2);
% 
% %FDUDB_Exp_x = [2; 5; 10];
% FDUDB_Exp_y = [38.08; 50.92; 67.28];  %%% WT with 10uM Flec (Liu, 2002)
% 
% Err_FDUDB2 = sum((FDUDB_Sim_y - FDUDB_Exp_y).^2) / 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% UDB Isotherm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Block_Sim = load('testBlock_param.txt');
% UDB_Sim_y = Block_Sim(:,2);
% 
% Drug = [0; 10; 100; 1000];
% EC50 = 11.2; n = 1.27; %%% WT 10uM Flec from Liu, 2002
% 
% Y_EXP =( 1./(1 + (Drug/EC50).^n  ))*100;
% 
% Err_UDB = sum((Y_EXP - UDB_Sim_y).^2) / 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TB Isotherm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TB_Sim_y = Block_Sim(:,3);
% 
% Drug = [0; 10; 100; 1000];
% %WTFlec = [100; 99.31; 95.51; 90.51; 78.8; 33.28; 7.18];  %%% Data from Abriel 2000
% EC50 = 59.3; n = 1.529; %%% WT 10uM Flec from Abriel, 2000  59.3, 1.529
% 
% Y_EXP = (1./(1 + (Drug/EC50).^n  ))*100;
% 
% Err_TB = sum((TB_Sim_y - Y_EXP).^2) / 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell_Sim = load('Cell_UV.txt');
% Cell_Sim_y = Cell_Sim(:,6);
% 
% Err_Cell = Cell_Sim_y *20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % NUFL optimization
% Total_Error = Err_FDUDB_NUFL + Err_RUDB_NUFL;

% Charged Rate Optimization
%Total_Error =  Err_RUDB100 + Err_SSA + Err_UDB + Err_TB + Err_FDUDB1 + Err_FDUDB2 + Err_Cell;

Total_Error = Err_SSA;
 
end