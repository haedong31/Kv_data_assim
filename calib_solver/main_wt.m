clc
close all
clear variables

warning('off', 'all')

% matching table
matching_table = readtable('./data/matching-table-wt.xlsx');
[num_files, ~] = size(matching_table);

% file names and capacitance values
file_names = matching_table.trace_file_name_4half;
% caps = matching_table.cap;

% exclude row not having 4.5-sec data
loop_idx = [];
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

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

p0 = [33, 15.5, 20, 8, 7, 0.3956, 0.00095, 0.051335, 0.14067, 0.387, ...
    22.5, 45.2, 40, 5.7, 2.058, 803, 0.05766, 0.07496, ...
    5334, 0.05766 ...
    13.17, 0.0428];
lb = [-50, -50, -50, eps, eps, eps, eps, eps, eps, eps, ...
    -70, -70, -70, eps, eps, eps, eps, eps, ...
    5000, eps, ...
    eps, eps
];
ub = [70, 50, 50, 50, 30, 10, 0.005, 0.5, 1, 1, ...
    50, 50, 50, 50, 100, 1000, 0.5, 0.5, ...
    10000, 0.5, ...
    100, 0.5
];

% default values
kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow10 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow20 = [5334, 4912, 0.05766];
kss0 = [0.0862, 1235.5, 13.17, 0.0428];

param_kto = zeros(17,1);
param_kslow1 = zeros(13,1);
param_kslow2 = zeros(3,1);
param_kss = zeros(4,1);

% index for tunning parameters
fixed_kto_idx = [4, 7, 8, 9, 11, 12, 15];
tune_kto_idx = setdiff(1:17, fixed_kto_idx);

fixed_kslow1_idx = [4, 6, 7, 10, 11];
tune_kslow1_idx = setdiff(1:13, fixed_kslow1_idx);

fixed_kslow2_idx = 1;
tune_kslow2_idx = [2, 3];

fixed_kss_idx = [1, 2];
tune_kss_idx = [3, 4];

% protocol
hold_volt = -70;
volts = -50:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
Ek = -91.1;

% main loop
num_files = length(loop_idx);
num_variables = 22;
sols = zeros(num_files, num_variables);
for i = loop_idx
    fprintf('[%i/%i] %s \n', i, num_files, file_names{i})

    % read data
    file_path = fullfile(pwd, 'data', 'wt-preprocessed2', file_names{i});
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
    opt_fun = @(p) obj_rmse(p, hold_volt, volts, time_space, yksum, Ek, false);

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
    save_path = fullfile(pwd, 'wt', strcat('calib_param_', file_names{i}));

    sol_kto = zeros(17, 1);
    sol_kslow1 = zeros(13, 1);
    sol_kslow2 = zeros(3, 1);
    sol_kss = zeros(4, 1);

    sol_kto(fixed_kto_idx) = kto0(fixed_kto_idx);
    sol_kslow1(fixed_kslow1_idx) = kslow10(fixed_kslow1_idx);
    sol_kslow2(fixed_kslow2_idx) = kslow20(fixed_kslow2_idx);
    sol_kss(fixed_kss_idx) = kss0(fixed_kss_idx);
   
    sol_kto(tune_kto_idx) = sol(1:10);
    sol_kslow1(tune_kslow1_idx) = sol(11:18);
    sol_kslow2(tune_kslow2_idx) = sol(19:20);
    sol_kss(tune_kss_idx) = sol(21:22);
    
    current_names = ["IKto", "IKslow1", "IKslow2", "IKss"];
    writematrix(current_names, save_path, "Sheet","Parameters", "Range","A1:D1");            

    writematrix(sol_kto, save_path, "Sheet","Parameters", "Range","A2");
    writematrix(sol_kslow1, save_path, "Sheet","Parameters", "Range","B2");
    writematrix(sol_kslow2, save_path, "Sheet","Parameters", "Range","C2");
    writematrix(sol_kss, save_path, "Sheet","Parameters", "Range","D2");
end

%% three currents (ikto, ikslow ikss; kcurrent_model2)
clc
close all
clear variables

warning('off', 'all')

% matching table
matching_table = readtable('./data/matching-table-wt.xlsx');
[num_files, ~] = size(matching_table);

% file names and capacitance values
file_names = matching_table.trace_file_name_4half;
% caps = matching_table.cap;

% exclude row not having 4.5-sec data
loop_idx = [];
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

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

p0 = [33, 15.5, 20, 8, 7, 0.3956, 0.00095, 0.051335, 0.14067, 0.387, ...
    22.5, 45.2, 40, 5.7, 2.058, 803, 0.05766, 0.07496, ...
    13.17, 0.0428];
lb = [-50, -50, -50, eps, eps, eps, eps, eps, eps, eps, ...
    -70, -70, -70, eps, eps, eps, eps, eps, ...
    eps, eps
];
ub = [70, 50, 50, 50, 30, 10, 0.005, 0.5, 1, 1, ...
    50, 50, 50, 50, 100, 1000, 0.5, 0.5, ...
    100, 0.5
];

% default values
kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kss0 = [0.0862, 1235.5, 13.17, 0.0428];

param_kto = zeros(17,1);
param_kslow = zeros(13,1);
param_kss = zeros(4,1);

% index for tunning parameters
fixed_kto_idx = [4, 7, 8, 9, 11, 12, 15];
tune_kto_idx = setdiff(1:17, fixed_kto_idx);

fixed_kslow_idx = [4, 6, 7, 10, 11];
tune_kslow_idx = setdiff(1:13, fixed_kslow_idx);

fixed_kss_idx = [1, 2];
tune_kss_idx = [3, 4];

% protocol
hold_volt = -70;
volts = -50:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
Ek = -91.1;

% main loop
num_files = length(loop_idx);
num_variables = 20;
sols = zeros(num_files, num_variables);
for i = loop_idx
    fprintf('[%i/%i] %s \n', i, num_files, file_names{i})

    % read data
    file_path = fullfile(pwd, 'data', 'wt-preprocessed2', file_names{i});
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
    opt_fun = @(p) obj_rmse2(p, hold_volt, volts, time_space, yksum, Ek, false);

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
    save_path = fullfile(pwd, 'wt', strcat('calib_param_', file_names{i}));

    sol_kto = zeros(17, 1);
    sol_kslow = zeros(13, 1);
    sol_kss = zeros(4, 1);

    sol_kto(fixed_kto_idx) = kto0(fixed_kto_idx);
    sol_kslow(fixed_kslow_idx) = kslow0(fixed_kslow_idx);
    sol_kss(fixed_kss_idx) = kss0(fixed_kss_idx);
   
    sol_kto(tune_kto_idx) = sol(1:10);
    sol_kslow(tune_kslow_idx) = sol(11:18);
    sol_kss(tune_kss_idx) = sol(19:20);
    
    current_names = ["IKto", "IKslow", "IKss"];
    writematrix(current_names, save_path, "Sheet","Parameters", "Range","A1:C1");            

    writematrix(sol_kto, save_path, "Sheet","Parameters", "Range","A2");
    writematrix(sol_kslow, save_path, "Sheet","Parameters", "Range","B2");
    writematrix(sol_kss, save_path, "Sheet","Parameters", "Range","C2");    
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
