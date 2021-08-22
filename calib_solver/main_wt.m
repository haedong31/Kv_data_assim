clc
close all
clear variables

warning('off', 'all')

% matching table
group_name = 'wt';
matching_table = readtable(fullfile(pwd, 'data', strcat('matching-table-', group_name, '.xlsx')));
[num_files, ~] = size(matching_table);

% file names and capacitance values
file_names = matching_table.trace_file_name_4half;

% exclude row not having 4.5-sec data
loop_idx = [];
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

% regularization
% field 1: selection of currents
current_names = {'ikto', 'ikslow1', 'ikss', 'ik1'};
num_currents = length(current_names);

% field 2: tunning index in individual current models
tune_idx1_kto = [1, 2, 15, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 9, 12, 13];
% tune_idx1_kslow2 = [10, 11];
tune_idx1_kss = [4, 5, 6, 7];
tune_idx1_k1 = [1, 3, 5, 7];
idx_info1 = {tune_idx1_kto, ...
    tune_idx1_kslow1, ...
    tune_idx1_kss, ...
    tune_idx1_k1};

% field 3: tunning index in decision variable, p
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

% optimization options
max_evals = 1e+6;
num_iters = 30;
options = optimoptions('fmincon', 'MaxFunctionEvaluations',max_evals, 'Display','off');

% optimization constraints 
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

p0 = [33, 15.5, 0.2087704319, 0.14067, 0.387, ...
    22.5, 45.2, 40.0, 803.0, 0.05766, 0.07496, ...
    0.0862, 1235.5, 13.17, 0.0428, ...
    59.215, 594.31, 1.02, 0.8];
lb = [-50, -50, eps, eps, eps, ...
    -70, -70, -70, eps, eps, eps, ...
    eps, eps, eps, eps, ...
    eps, eps, eps, eps
];
ub = [70, 50, 1, 1, 1 ...
    50, 50, 50, 2000, 0.5, 0.5, ...
    1, 2000, 100, 0.5, ...
    120, 1000, 3, 2
];

% default values
kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow2_default = [5334, 4912, 0.05766];
kur_default = [270, 1050, 0];
kss_default = [0.0862, 1235.5, 13.17, 0.0428];
k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

% protocol
hold_volt = -70;
volts = -50:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
Ek = -91.1;

volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = Ek;

% main loop
num_files = length(loop_idx);
sols = zeros(num_files, cul_idx_len);
% mkdir(fullfile(pwd, group_name));
for i = loop_idx
    fprintf('[%i/%i] %s \n', i, num_files, file_names{i})

    % read data
    file_path = fullfile(pwd, 'data', strcat(group_name, '-preprocessed2'), file_names{i});
    trace_data = table2array(readtable(file_path));
        
    t = trace_data(:, 1);
    yksum = trace_data(:, 2:end);

    % estimate time points
    [~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
    [~, ideal_end_idx] = min(abs(t - ideal_end_time));

    t = t(1:ideal_end_idx);
    yksum = yksum(1:ideal_end_idx, :);

    % time space
    time_space = cell(1, 3);
    time_space{1} = t;
    time_space{2} = t(1:ideal_hold_idx);
    pulse_t = t(ideal_hold_idx+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    time_space{3} = pulse_t_adj;

    % objective function
     % obj_rmse(p0, model_struct, volt_space, time_space, yksum);
    opt_fun = @(p) obj_rmse(p, model_struct, volt_space, time_space, yksum);

    % run optimization
    rmse_list = zeros(num_iters, 1);
    sol_list = cell(num_iters, 1);

    % first run with p0
    [sol, fval] = fmincon(opt_fun, p0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    sol_list{1} = sol;
    rmse_list(1) = fval;

    for j = 2:num_iters
        fprintf('[%i/%i] \n', j, num_iters)

        % random intialization
        running_p0 = zeros(length(p0), 1);
        for k = 1:length(p0)
            unif_dist = makedist('Uniform', 'lower',lb(k), 'upper',ub(k));
            running_p0(k) = random(unif_dist, 1);
        end
        
        % optimization
        try
            [sol, fval] = fmincon(opt_fun, running_p0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            sol_list{j} = sol;
            rmse_list(j) = fval;
        catch
            rmse_list(j) = 1e+3;
        end
    end

    [~, best_fit_idx] = min(rmse_list);
    sol = sol_list{best_fit_idx};
    sols(i, :) = sol;

    % save calibrated solution
    save_path = fullfile(pwd, group_name, strcat('calib_param_', file_names{i}));

    sol_kto = kto_default;
    sol_kslow1 = kslow1_default;
    sol_kss = kss_default;
    sol_k1 = k1_default;
   
    sol_kto(idx_info1{1}) = sol(idx_info2{1});
    sol_kslow1(idx_info1{2}) = sol(idx_info2{2});
    sol_kss(idx_info1{3}) = sol(idx_info2{3});
    sol_k1(idx_info1{4}) = sol(idx_info2{4});
    
    save_current_names = ["IKto", "IKslow1", "IKss", "IK1"];
    writematrix(save_current_names, save_path, "Sheet","Parameters", "Range","A1:D1");
    writematrix(sol_kto', save_path, "Sheet","Parameters", "Range","A2");
    writematrix(sol_kslow1', save_path, "Sheet","Parameters", "Range","B2");
    writematrix(sol_kss', save_path, "Sheet","Parameters", "Range","C2");
    writematrix(sol_k1', save_path, "sheet","Parameters", "Range","E2");
end

%% viz
clc
close all

[~, ~, ~, ~, yksum_hat] = kcurrent_model(sol, hold_volt, volts(1), time_space, Ek);
plot(t, yksum(:, 1))
for i = 2:num_volts
    [~, ~, ~, ~, yksum_hat] = kcurrent_model(sol, hold_volt, volts(i), time_space, Ek);
    
    hold on
    plot(t, yksum_hat)
    hold off
    axis tight
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    legend('Experimental','Model Prediction')
end
