clc
close all
clear variables

warning('off', 'all')

% matching table
matching_table = readtable('./data/matching-table-wt.xlsx');
[num_files, ~] = size(matching_table);

% file names and capacitance values
file_names = matching_table.trace_file_name_4half;
caps = matching_table.cap;

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
    file_path = fullfile(pwd, 'data', 'wt', file_names{i});
    trace_data = table2array(readtable(file_path));

    % downsample
    trace_data = downsample(trace_data, 20);
    
    t = trace_data(:, 1);
    yksum = trace_data(:, 2:end);

    % estimate time points
    [~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
    [~, ideal_end_idx] = min(abs(t - ideal_end_time));

    t = t(1:ideal_end_idx);
    yksum = yksum(1:ideal_end_idx, :);

    % normalize data
    yksum = yksum ./ caps(i);

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
    sols(i, :) = sol_list{best_fit_idx};
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
