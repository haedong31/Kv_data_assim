% Sequential optimization strategy for WT Charged Flecainide Drug Fitting:
% As in the NU FL example, 4 protocols (SSA, RUDB, FDUDB, and BLOCK) were optimized. 
% However, for this simulation, each protocol was optimized independently for 100 iterations, and then 
% an additional protocol was added and the two were optimized together for an additional 100 iterations, etc. 
% until 4 protocols were optimized for a total of 400 iterations. 

parpool('local',5);

%%%% % Fitting Charged Flecainide
mex main_SSA.cpp
mex main_RUDB.cpp
mex main_FDUDB1.cpp
mex main_FDUDB2.cpp
mex main_BLOCK.cpp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Initial conditions
lb = [0; 0; 0; 0; 0; 0; 0; 0];
ub = [Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf];
options = optimset('Display','iter','TolFun',1e-2, 'TolX', 1e-2, 'MaxFunEvals', 500000000, 'MaxIter', 100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% SIMBLOCK 1: Simulate protocol (1)
Inputs = [1; 1; 1; 1; 1; 1; 1; 1];
SIMBLOCK = 1;
Inputs_Final = fminsearchbnd(@WT_Drug_Sim_Exp_FLEC_SEQ, Inputs, lb, ub, options, SIMBLOCK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% SIMBLOCK 2: Use the final inputs from SIMBLOCK 1, and restart simulation
% wtih protocols (1) and (2)
Inputs = Inputs_Final;
SIMBLOCK = 2;
Inputs_Final = fminsearchbnd(@WT_Drug_Sim_Exp_FLEC_SEQ, Inputs, lb, ub, options, SIMBLOCK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%SIMBLOCK 3: As above; restart simulation with final inputs from SIMBLOCK
% 2, and restart with protcols (1), (2), and (3)
Inputs = Inputs_Final;
SIMBLOCK = 3;
Inputs_Final = fminsearchbnd(@WT_Drug_Sim_Exp_FLEC_SEQ, Inputs, lb, ub, options, SIMBLOCK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%SIMBLOCK 4: As above; restart simulation with final inputs from SIMBLOCK
% 3, and restart with protcols (1), (2), (3), and (4)
Inputs = Inputs_Final;
SIMBLOCK = 4;
Inputs_Final = fminsearchbnd(@WT_Drug_Sim_Exp_FLEC_SEQ, Inputs, lb, ub, options, SIMBLOCK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


 
delete(gcp('nocreate'));