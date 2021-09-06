function Total_Error = WT_REDUCED_CHANNEL_Sim_Exp_Norm(Inputs)

%Parallel loop of protocols to simulate:
% SSA = steady state availability; ACT = steady state activation; RUDB =
% recovery from use-dependent block; TAU = tau to 50% decay; RFI = recovery
% from 1 pulse inactivation

i = 1;
parfor i = 1:5;
if i==1  
    main_SSA(Inputs);
end 
if i==2 
    main_ACT(Inputs);
end
if i==3  
    main_RUDB(Inputs);
end
if i==4  
    main_Tau(Inputs);
end
if i==5 
    main_RFI(Inputs);
end

end;

% PROTOCOLS AND COMPARISION TO EXPERIMENTAL DATA BEGIN HERE
% Data from:
% Liu, H., Tateyama, M., Clancy, C., Abriel, H., and Kass, R., Channel Openings Are Necessary but not Sufficient for 
% Use?dependent Block of Cardiac Na+ Channels by Flecainide Evidence from the Analysis of Disease?linked Mutations. 
% Journal of General Physiology 120 (1), 39 (2002). 
% Rivolta, I., Abriel, H., Tateyama, M., Liu, H., Memmi, M., Vardas, P., Napolitano, C., Priori, S. G., and Kass, R. S., 
% Inherited Brugada and long QT?3 syndrome mutations of a single residue of the cardiac sodium channel confer distinct channel 
% and clinical phenotypes. J Biol Chem 276 (33), 30623 (2001). 

% RFI Peak at -100mV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RFI100 = load('testRFI100_param.txt');
RFI100_y = RFI100(:,2);
t1 = [1.005; 5.005; 10.005; 25.01; 50; 100; 150.01; 250.01; 500.01; 1000.01; 2000.01; 3000; 5000; 7500; 10000];

A1_EXP = 0.65; T1_EXP = 2.68; A2_EXP = 0.31; T2_EXP = 12.47; C1_EXP = 0.983;

Y_100 = C1_EXP - A1_EXP*exp(-t1/T1_EXP) - A2_EXP*exp(-t1/T2_EXP);
Err_RFI100 = sum((RFI100_y*100 - Y_100*100).^2)  / 15;

Err_RFI = Err_RFI100;

% RUDB  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RUDB100 = load('testRUDB100_param.txt');

RUDB100_y = RUDB100(:,2);

t1 = [1.005; 5.005; 10.005; 25.01; 50; 100; 150.01; 250.01; 500.01; 1000.01; 2000.01; 3000; 5000; 7500; 10000];

A1_EXP = 0.7283; T1_EXP = 6.805; A2_EXP = .2883; T2_EXP = 291.2; C1_EXP = 0.9826;
Y_100 = C1_EXP - A1_EXP*exp(-t1/T1_EXP) - A2_EXP*exp(-t1/T2_EXP);
Err_RUDB100 = sum((RUDB100_y*100 - Y_100*100).^2)  / 15;

Err_RUDB = Err_RUDB100;

%%% Activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = (20:-5:-75)';
Act_Sim = load('testAct_param.txt');
Act_Sim_y = Act_Sim(:,6);

Vhalf_EXP = -23.3; 
K_EXP = -7.1;      

Y_EXP = 1./(1 + exp((V - Vhalf_EXP)/K_EXP));

Err_Act = sum((Act_Sim_y*100 - Y_EXP*100).^2)  / 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SSA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = (-125:5:-30)';

SSA_Sim = load('testSSA_param.txt');
SSA_Sim_y = SSA_Sim(:,4); 

Vhalf_EXP = -70.4;  
K_EXP = 6.34; 
 
Y_EXP = 1./(1+exp((V - Vhalf_EXP)/K_EXP));
 
Err_SSA = sum((SSA_Sim_y*100 - Y_EXP*100).^2) / 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Tau H %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tau_Sim = load('testTau_param.txt');
Tau_Sim_y = Tau_Sim(:,7);

%Tau_Exp = load('Tau50.txt');
%Tau_Exp_y = Tau_Exp(:,2);
Tau_Exp_x = [20; 15; 10; 5; 0; -5; -10; -15; -20];
Tau_Exp_y = [0.4775; 0.4955; 0.5174; 0.5493; 0.6050; 0.6488; 0.7284; 0.8438; 0.9831 ]; 

Err_Tau = sum((Tau_Sim_y*10 - Tau_Exp_y*10).^2)   /  9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Err_MOT = abs(2 - 0.19245*(Inputs(6)*exp(40/20.3) + ...
%          (1 + Inputs(13))*(Inputs(11)*exp(-30/Inputs(12)))));


%MAKE SURE TO CHECK THIS DEPENDING ON THE PARAMS FOR THE WT DRUG FREE
%this is verified for MOT for the 1 exponential a11 WT model
Err_MOT = abs(2.0 - 0.19245*(  (1+Inputs(15))*(Inputs(13)*exp(-30/Inputs(14)))+ ...
          Inputs(8)/(Inputs(5)*exp(-30/Inputs(6))) ) );

% This is for MOT = 0.5 at -30mV, so 1/MOT = 2, 0.19245 is the Tfactor, and
% the rest of the equation is 1/ sum of the egress from open states = 1/
% (a2 +b13 + ax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective function to be minimized: sum of the above individual protocol
% errors
Total_Error = Err_MOT + Err_RUDB + Err_Act + Err_SSA + Err_Tau + Err_RFI;

end