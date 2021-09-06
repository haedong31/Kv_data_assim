
% This is the optimization code for the Na channels. To use it, invoke the
% following protocols to optimize over:
% For WT drug free: SSA, Act, RUDB, RFI, Tau.

% For Drug Optimization, include SSA, RUDB, FDUDB, and BLOCK (and CELL) on
% some occiasions

% This code requires the parallel computing toolbox

%Open up a parallel pool of MATLAB workers (ideally equal to the number of
%protocols to be simulated)
parpool('local',5);

%%%%% Fitting Drug Free
mex main_SSA.cpp
mex main_ACT.cpp
mex main_RUDB.cpp
mex main_RFI.cpp
mex main_Tau.cpp

% %%%%% % Fitting NUFL
% mex main_RUDB.cpp
% mex main_FDUDB_NUFL.cpp

%%%% % Fitting Charged Flecainide
% mex main_SSA.cpp
% mex main_RUDB.cpp
% mex main_FDUDB1.cpp
% mex main_FDUDB2.cpp
% mex main_BLOCK.cpp
% mex main_CELL.cpp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL GUESSES HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Wild Type Initial Inputs from STM paper (Moreno et al, 2011)
 InputsInitial_STM = [
4.0;
0.1027;
0.23;
0.1539;
0.1605;
0.1766;
2.624e-3;
24.7;
0.4996;
25.14;
3.4480e+00
2.3070e+01
5.0000e-02
2.2222e-02];

%Initial inputs extracted from experimental data as noted in the
%SUPPLEMENTARY INFORMATION (1 exponential of a11)
Inputs = [
    0.1027;  % 1/(Act1*exp(-V/Act2) ) **Directly from Misyye Noma
    9.3;      % Act2 **Directly from Misyye Noma
    1.0;       % a12 multiplier
    1.0;       % a13 multiplier
    19.832;  % 1/(Deact1*exp(V/Deact2) **Directly from Misyye Noma
    20.3;    % Deact2 **Directly from Misyye Noma
    1.0;       % b12 multiplier
    1.0;       % b13 multiplier
    4.301e-5; %a3_v1
    15.934;   %a3_v2
    0.4996;   %b3_v1
    25.14;    %b3_v2
    3.448;    %a2_v1
    23.07;    %a2_v2
    5.0000e-02; %ax
    2.2222e-02]; %bx

%AFter 500 errror = 48
Inputs_Final =[
   2.1129e-02;
   2.3061e+01;
   1.8491e+00;
   1.1701e+00;
   5.6441e+00;
   6.9120e+00;
   8.1873e-01;
   1.3887e+00;
   3.9970e-05;
   9.8433e+00;
   4.4874e-01;
   3.2013e+01;
   8.7050e+00;
   3.6508e+01;
   2.0150e-02;
   1.8169e-02];

%After 500, error 24
Inputs_Final =[
   8.0040e-03;
   3.2474e+01;
   1.4777e+00;
   1.2295e-01;
   2.4236e+00;
   8.4315e+00;
   1.2144e-03;
   2.7224e+00;
   3.2080e-05;
   9.5932e+00;
   1.1844e+00;
   2.1881e+01;
   1.4971e+01;
   3.9244e+01;
   3.8651e-02;
   2.0173e-02];
%CONVERGED after 235, error = 23.9
Inputs_Final =[
   7.6178e-03;
   3.2764e+01;
   5.8871e-01;
   1.5422e-01;
   2.5898e+00;
   8.5072e+00;
   1.3760e-03;
   2.8880e+00;
   3.2459e-05;
   9.5951e+00;
   1.3771e+00;
   2.1126e+01;
   1.1086e+01;
   4.3725e+01;
   4.1476e-02;
   2.0802e-02];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NU FL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Initial Inputs for the NU FL were set at all 1's (as below):
Inputs = [1; 1; 1; 1; 1; 1; 1; 1];

%Converged after 381, error = 6.2
Inputs_Final = [
   1.6936e+01; %ax2
   1.0714e-01; %a13n
   5.9858e-01; %a_22
   5.6001e-01; %b_33
   2.8919e-07; %a_44
   1.2546e+01; %b_44
   6.6880e-01; %ki_on
   5.2802e-06]; %kc_on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHARGED FLEC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Trial 1: As an initial guess, set all the charged inputs to 1:
Inputs = [1; 1; 1; 1; 1; 1; 1; 1];

%CONVERGED after 256, error 480 (from all 1's)
Inputs_Final_Ones =[
   2.4481e+00; %ax1
   1.7709e-02; %bx1
   3.6791e+00; %a13c
   3.0744e-01; %a22
   5.2879e-01; %b33
   1.5913e+00; %a33
   4.6704e-09; %a44
   1.0091e+00];%b44
%%%%%%%%*******************************%%%%%%%%%%%%%%
% Trial 2: RESULT FROM FACTORIAL DESIGN. The best sequence was SSA, BLOCK, RUDB, FDUDB
% After 400 iterations, (SBRF, error = 286). See 'Global_test_SEQ.m'.
Inputs_Factorial = [
    1.9541E-02;
    4.3132E-04;
    1.6732E+02;
    1.3239E-02;
    2.7360E-02;
    5.6535E-04;
    7.7744E+00;
    1.6783E+00];
%CONVERGED AFTER 311, error 272 (noted as "Sequential Trial" in paper)
Inputs_Final_Factorial =[
   1.5090e-02;
   1.2162e-09;
   5.6531e+02;
   1.4486e-02;
   2.3675e-02;
   5.3763e-04;
   9.9942e-02;
   3.7766e+00];
%%%%%%%%*******************************%%%%%%%%%%%%%%
% Trial 3: Start with STM Final Rates (Moreno, 2011)for CHARGED FLEC
STM_Rates = [
    5.7839e-05;
	1.6689e-08;
	3.6324e-03;
	1.4847e+03;
	1.7352e-06;
	6.7505e-05;
	2.4135e+00;
	4.9001e-02];
%CONVERGED after 432; error 165 (no cell sim)
Inputs_Final =[
   1.0836e-05;
   4.2106e-08;
   2.4824e-03;
   1.2663e+02;
   4.8810e-06;
   1.8309e-04;
   2.5183e+00;
   4.6378e-02];
%%%%%%%%*******************************%%%%%%%%%%%%%%
%Trial 4: Selected inuts 
Inputs = [1; 0.0001; 1; 1; 0.0001; 0.0001; 1; 1];

%Converged after 388, error = 172
Inputs_Final_Selected =[
   1.4439e-01;
   3.6637e-04;
   3.4583e-03;
   2.1129e+00;
   4.8304e-08;
   2.6205e-04;
   2.3428e+00;
   1.0219e-02];

%%%%%%%%%%%%%%%%%%%%%%%%%% ROBUSTNESS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For robustness testing, the following command was used to generate a
% vector of random numbers:
% a = 0.95; b = 1.05; r = a + (b-a).*rand(16,1) (5% deviation)
% Inputs = r.*Inputs_WT_DrugFree (vector multiplier)

Inputs_WT_DrugFree =[
   7.6178e-03;
   3.2764e+01;
   5.8871e-01;
   1.5422e-01;
   2.5898e+00;
   8.5072e+00;
   1.3760e-03;
   2.8880e+00;
   3.2459e-05;
   9.5951e+00;
   1.3771e+00;
   2.1126e+01;
   1.1086e+01;
   4.3725e+01;
   4.1476e-02;
   2.0802e-02];

% 5% Deviation Trial 1:
Inputs = [
   7.8576e-03;
   3.4094e+01;
   5.6675e-01;
   1.6060e-01;
   2.6241e+00;
   8.1648e+00;
   1.3455e-03;
   2.9015e+00;
   3.3944e-05;
   1.0041e+01;
   1.3299e+00;
   2.2120e+01;
   1.1593e+01;
   4.3661e+01;
   4.2721e-02;
   2.0057e-02];
%Converged after 211, error = 17.5501
Inputs_Final =[
   9.1058e-03;
   3.8634e+01;
   5.4867e-01;
   1.6237e-01;
   3.1673e+00;
   1.3324e+01;
   1.3013e-03;
   2.0818e+00;
   3.2845e-05;
   1.2884e+01;
   1.3870e+00;
   1.9185e+01;
   8.1914e+00;
   4.1562e+01;
   4.2109e-02;
   1.8890e-02];

% 5% Deviation Trial 2
Inputs =[
   7.5582e-03;
   3.4126e+01;
   6.0591e-01;
   1.6131e-01;
   2.6301e+00;
   8.1122e+00;
   1.4240e-03;
   3.0133e+00;
   3.3039e-05;
   9.8424e+00;
   1.4106e+00;
   2.0898e+01;
   1.1258e+01;
   4.2287e+01;
   4.2331e-02;
   1.9828e-02];
%Converged after 196, error 17.5505
Inputs_Final =[
   8.7046e-03;
   3.8729e+01;
   5.8737e-01;
   1.6398e-01;
   3.1914e+00;
   1.3353e+01;
   1.3807e-03;
   2.1776e+00;
   3.2307e-05;
   1.2545e+01;
   1.4748e+00;
   1.7910e+01;
   7.9433e+00;
   4.0366e+01;
   4.1839e-02;
   1.8473e-02];

% 5% Deviation Trial 3
Inputs =[
   7.4479e-03;
   3.1277e+01;
   5.6499e-01;
   1.5921e-01;
   2.6403e+00;
   8.3516e+00;
   1.4380e-03;
   2.7535e+00;
   3.2260e-05;
   9.4815e+00;
   1.4137e+00;
   2.1750e+01;
   1.0739e+01;
   4.3680e+01;
   4.1250e-02;
   2.1106e-02];
%Converged after 206, error 17.5502
Inputs_Final =[
   8.3681e-03;
   3.5006e+01;
   5.4958e-01;
   1.6323e-01;
   3.1359e+00;
   1.3190e+01;
   1.4199e-03;
   2.0629e+00;
   3.1053e-05;
   1.1831e+01;
   1.4810e+00;
   1.8874e+01;
   7.9225e+00;
   4.1411e+01;
   4.0555e-02;
   1.9839e-02];

% 10% Deviation Trial 1
Inputs =[
   7.9368e-03;
   3.4433e+01;
   5.6234e-01;
   1.5976e-01;
   2.6701e+00;
   7.9332e+00;
   1.2711e-03;
   2.8871e+00;
   3.5444e-05;
   9.2888e+00;
   1.4006e+00;
   1.9959e+01;
   1.1643e+01;
   4.1583e+01;
   4.1525e-02;
   2.1630e-02];
%Converged after 195, error = 17.5505
Inputs_Final =[
   9.1363e-03;
   3.9017e+01;
   5.4261e-01;
   1.6124e-01;
   3.2386e+00;
   1.3067e+01;
   1.2331e-03;
   2.0874e+00;
   3.4708e-05;
   1.1876e+01;
   1.4651e+00;
   1.7053e+01;
   8.1889e+00;
   3.9568e+01;
   4.1113e-02;
   2.0312e-02];
   
% 10% Deviation Trial 2
   Inputs =[
   8.2134e-03;
   3.5774e+01;
   5.9427e-01;
   1.4307e-01;
   2.4081e+00;
   8.0946e+00;
   1.4698e-03;
   2.7461e+00;
   3.4499e-05;
   9.1029e+00;
   1.4953e+00;
   2.0492e+01;
   1.0413e+01;
   4.1548e+01;
   4.2439e-02;
   2.0691e-02];
%Converged after 252 iterations, 17.5498
Inputs_Final =[
   9.0698e-03;
   4.0539e+01;
   6.5552e-01;
   1.3590e-01;
   2.8457e+00;
   1.3108e+01;
   1.4510e-03;
   1.9860e+00;
   3.2732e-05;
   1.1630e+01;
   1.5591e+00;
   2.1206e+01;
   7.4054e+00;
   3.8170e+01;
   3.9972e-02;
   1.8173e-02];

% 10% Deviation Trial 3
Inputs =[
   7.3918e-03;
   3.4932e+01;
   5.9875e-01;
   1.5575e-01;
   2.8059e+00;
   8.1428e+00;
   1.4468e-03;
   3.0346e+00;
   3.1683e-05;
   9.7253e+00;
   1.2603e+00;
   1.9241e+01;
   1.1154e+01;
   4.6166e+01;
   4.5076e-02;
   1.9262e-02];
%Converged afer 237, error 17.5499
Inputs_Final =[
   8.4841e-03;
   3.9573e+01;
   5.8112e-01;
   1.5797e-01;
   3.3805e+00;
   1.3249e+01;
   1.4048e-03;
   2.1880e+00;
   3.0639e-05;
   1.2362e+01;
   1.3042e+00;
   1.6833e+01;
   7.8827e+00;
   4.3999e+01;
   4.4573e-02;
   1.8114e-02];

% 25% Deviation Trial 1
Inputs =[
   7.8799e-03;
   3.2263e+01;
   4.4504e-01;
   1.4166e-01;
   2.1524e+00;
   9.7590e+00;
   1.2461e-03;
   2.9292e+00;
   2.7033e-05;
   1.0084e+01;
   1.2139e+00;
   2.2754e+01;
   1.2135e+01;
   4.9150e+01;
   4.0450e-02;
   1.6473e-02];
%Converged after 220, error = 17.5499
Inputs_Final =[
   5.4642e-03;
   3.7201e+01;
   5.0519e-01;
   1.3557e-01;
   3.0341e+00;
   1.4019e+01;
   1.2424e-03;
   1.7769e+00;
   2.5635e-05;
   1.2904e+01;
   1.2270e+00;
   2.2582e+01;
   1.0388e+01;
   4.3592e+01;
   3.7157e-02;
   1.3227e-02];

% 25% Deviation Trial 2
Inputs =[
   6.5855e-03;
   3.9535e+01;
   4.8639e-01;
   1.7934e-01;
   2.6394e+00;
   1.0618e+01;
   1.0858e-03;
   2.8052e+00;
   2.6075e-05;
   1.1811e+01;
   1.0360e+00;
   2.4030e+01;
   1.2845e+01;
   5.1786e+01;
   3.2858e-02;
   1.9760e-02];
%Converged after 179, error 17.5504
Inputs_Final =[
   7.3787e-03;
   2.9081e+01;
   4.8999e-01;
   1.7543e-01;
   3.3278e+00;
   1.5554e+01;
   1.2770e-03;
   2.1194e+00;
   2.5898e-05;
   1.1352e+01;
   8.8548e-01;
   2.7818e+01;
   1.0971e+01;
   4.7372e+01;
   3.2472e-02;
   1.7579e-02];

% 25% Deviation Trial 3
Inputs =[
   6.7032e-03;
   3.7680e+01;
   5.6852e-01;
   1.8589e-01;
   2.1778e+00;
   7.5025e+00;
   1.1321e-03;
   2.3625e+00;
   3.8452e-05;
   9.9775e+00;
   1.4114e+00;
   1.7376e+01;
   1.3043e+01;
   4.6393e+01;
   3.8385e-02;
   2.0940e-02];
%Converged after 192, error 17.551
Inputs_Final =[
   7.7542e-03;
   4.2844e+01;
   5.4122e-01;
   1.8406e-01;
   2.7286e+00;
   1.3297e+01;
   1.0709e-03;
   1.6739e+00;
   3.7428e-05;
   1.3609e+01;
   1.5137e+00;
   1.4090e+01;
   8.6616e+00;
   4.3690e+01;
   3.9129e-02;
   1.9191e-02];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Drug Free Fitting
% Define lower (lb), and upper (ub) bounds
lb = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
ub = [Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf];
  
% % Neutral Flecainide Fitting
%  Inputs = [1; 1; 1; 1; 1; 1; 1; 1];
%  lb = [0; 0; 0; 0; 0; 0; 0; 0];
%  ub = [Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf];

%Charged Flecainide Fitting
%Inputs = [1; 1; 1; 1; 1; 1; 1; 1];
%lb = [0; 0; 0; 0; 0; 0; 0; 0];
%ub = [Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf];

%Inputs = STM_Initial;
    
  %Nelder Mead:
  options = optimset('Display','iter','TolFun',1e-2, 'TolX', 1e-2, 'MaxFunEvals', 500000000, 'MaxIter', 500); 
  Inputs_Final = fminsearchbnd(@WT_REDUCED_CHANNEL_Sim_Exp_Norm, Inputs, lb, ub, options);
  
  %options = optimset('Algorithm', 'sqp','Display','iter','MaxIter', 250);
  %Inputs_Final = fmincon(@WT_Drug_Sim_Exp_FLEC, Inputs,[],[],[], [], lb, ub, [], options);
  
  %Patternsearch
  %options = psoptimset('Display','iter', 'MaxIter', 10);
  %Inputs_Final = patternsearch(@WT_Drug_Sim_Exp_FLEC, Inputs,[],[],[], [], lb, ub, [], options);
  
  delete(gcp('nocreate'));
  
  
 
  
