clc
close all
clear variables

% 1024-run designs
dgn_model2_Ito = readtable('./dgn-model2-Ito.csv');
dgn_model2_IKslow = readtable('./dgn-model2-IKslow.csv');
dgn_model3_Ito = readtable('./dgn-model3-Ito.csv');
dgn_model3_IKslow = readtable('./dgn-model3-IKslow.csv');

% levels of model2 Ito
low2to = [0, 0, 0, 0.7, 0.018064, 0.003577, ...
          0.03956, 0.006237, 0.0000152, 0.0067083, 0.000095, 0.0051335];
high2to = [60, 17, 67, 14, 0.36128, 0.07154, ...
           0.7912, 0.12474, 0.000304, 0.134166, 0.0019, 0.10267];

% levels of model2 IKslow
low2kslow = [0, 0.77, 0.0493, 0.00629, 0.2058, ...
             0, 0.57, 400, 1, 0, 0.57];
high2kslow = [45, 15.3, 0.986, 0.1258, 4.116, ...
              90.4, 11.4, 3000, 340, 90.4, 11.4];

% % levels of model3 Ito
% low3to = [];
% high3to = [];
% 
% % levels of model3 IKslow
% low3kslow = [];
% high3kslow = [];

% variation of responses of model 2
time_space = cell(1, 3);
time_step = 1;
holdT = 500;
P1T = 25000;
tH = 0:time_step:holdT;
tP1 = (holdT+time_step):time_step:P1T;
tP1_adj = tP1 - tP1(1);
t = [tH, tP1];

time_space{1} = t;
time_space{2} = tH;
time_space{3} = tP1_adj;

num_runs = 1024;
res_ctl_to = zeros(num_runs, 4);
res_ctl_kslow = zeros(num_runs, 4);
for i=1:num_runs
    fprintf('### Exp %i/%i \n', i, num_runs)
    
    dgn_to = table2array(dgn_model2_Ito(i, :));
    dgn_kslow = table2array(dgn_model2_IKslow(i, :));

    ito_param = zeros(1, 12);
    ikslow_param = zeros(1, 11);

    for j=1:12
        if dgn_to(j) == 1
            ito_param(j) = high2to(j);
        else
            ito_param(j) = low2to(j);
        end
    end
    
    for k=1:11
        if dgn_kslow(k) == 1
            ikslow_param(k) = high2kslow(k);
        else
            ikslow_param(k) = low2kslow(k);
        end
    end
    
    param = cell(1, 2);
    param{1} = ito_param;
    param{2} = ikslow_param;
    [res_to, res_kslow] = Ktrace2_doe(param, -70, 50, time_space);

    res_ctl_to(i, :) = res_to;
    res_ctl_kslow(i, :) = res_kslow;
end

res_tbl1 = array2table(res_ctl_to, 'VariableNames',{'peak','time_const','ssa','ssi'});
res_tbl2 = array2table(res_ctl_kslow, 'VariableNames',{'peak','time_const','ssa','ssi'});

writetable(res_tbl1, 'res-model2-Ito.csv')
writetable(res_tbl2, 'res-model2-IKslow.csv')
