%% whole current model
clc
close all
clear variables

dgn_mx = table2array(readtable(fullfile(pwd, 'kcurrent_dgn_matrix.csv')));

% protocol
hold_volt = -70;
volt = 10;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
ek = -91.1; 

time_space = cell(3, 1);
t = 1:ideal_end_time;
time_space{1} = t;
time_space{2} = t(1:ideal_hold_time);
pulse_t = t(ideal_hold_time+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

protocol = cell(4, 1);
protocol{1} = hold_volt;
protocol{2} = volt;
protocol{3} = time_space;
protocol{4} = ek;

% model structure
% field 1
current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikur', 'ikss', 'ik1'};
num_currents = length(current_names);

% field 2
tune_idx1_kto = 1:14;
tune_idx1_kslow1 = 1:10;
tune_idx1_kslow2 = 9:10;
tune_idx1_kur = 9:10;
tune_idx1_kss = 4:6;
tune_idx1_k1 = 1:10;
idx_info1 = {tune_idx1_kto, ...
    tune_idx1_kslow1, ...
    tune_idx1_kslow2, ...
    tune_idx1_kur, ...
    tune_idx1_kss, ...
    tune_idx1_k1};

% field 3
idx_info2 = cell(1, num_currents);
cul_idx_len = 0;

for i = 1: num_currents
    cul_idx_len_old = cul_idx_len;
    cul_idx_len = cul_idx_len + length(idx_info1{i});
    idx_info2{i} = (1+cul_idx_len_old):(cul_idx_len);
end

% create model structure
model_info = [current_names; idx_info1; idx_info2];
field_names = {'name', 'idx1', 'idx2'};
model_struct = cell2struct(model_info, field_names, 1);

% lower and upper bound
[num_runs, num_vars] = size(dgn_mx);
% p0 = [33, 15.5, 20, 8, 7, 0.3956, 0.00095, 0.051335, 0.14067, 0.387, ...
%     22.5, 45.2, 40, 5.7, 2.058, 803, 0.05766, 0.07496, ...
%     5334, 0.05766, ...
%     1050, 0 ...
%     13.17, 0.0428, ...
%     59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];
lb = [-50, -50, -50, -10, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, ...
    -70, -70, -70, eps, eps, eps, eps, eps, eps, eps, ...
    5000, eps, ...
    eps, 500, ...
    eps, eps, eps, ...
    eps, -20, eps, -30, eps, eps, eps, eps, eps, eps
];
ub = [70, 50, 50, 40, 50, 30, 1, 1, 1, 10, 0.005, 0.3, 0.005, 0.5, ...
    50, 50, 50, 10, 50, 50, 1, 100, 1000, 50, ...
    10000, 10000, ...
    500, 2000, ...
    1, 2000, 100, ...
    120, 30, 1000, 30, 3, 1, 2, 0.25, 0.25, 1
];

res = zeros(num_runs, 1);
for i = 1:num_runs
    fprintf("Degisn [%i/%i] \n", i, num_runs)
    
    % set up parameters according to given design
    dgn = dgn_mx(i, :);
    p = zeros(num_vars, 1);
    for j = 1:num_vars
        if dgn(j) == 1
            p(j) = ub(j);
        else
            p(j) = lb(j);
        end
    end

    [yksum, ~] = kcurrent_model(p, model_struct, protocol);
end
