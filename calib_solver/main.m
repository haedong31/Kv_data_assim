% 4.5-sec voltage-dependent data 15o29009 (from WT)
clc
close all
clear variables

% load experimental data
trace_data = table2array(readtable('./data/wt-preprocessed/15o29009.xlsx'));
t = trace_data(:,1);

volt_steps = 6:11;
min_volt = -50;
volts = zeros(length(volt_steps), 1);
for i = 1:length(volt_steps)
    volts(i) = min_volt + (volt_steps(i)-1)*10;
end

yksum = trace_data(:, (volt_steps + 1));
[~, num_volts] = size(yksum);

% visualize experimental data
figure(1)
plot(t, yksum(:, 1))
hold on
for i=2:num_volts
    plot(t, yksum(:, i))
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

% estimate holding time
ideal_hold_time = 124;
[~, ideal_hold_idx] = min(abs(t-ideal_hold_time));
init_stable_val = sqrt(var(yksum(ideal_hold_idx,:)));

% generate time space
time_space = cell(1, 3);
time_space{1} = t;
time_space{2} = t(1:ideal_hold_idx);
pulse_t = t(ideal_hold_idx+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;   

% protocol
hold_volt = -70;
Ek = -91.1;

% optimization constraints 
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
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

options = optimoptions('fmincon', 'MaxFunctionEvaluations',1e+6);
opt_fun = @(p) obj_rmse(p, hold_volt, volts, time_space, yksum, Ek, true);

num_iters = 15;
rmse_list = zeros(num_iters, 1);
sol_list = cell(num_iters, 1);

% first run with p0
[sol, fval] = fmincon(opt_fun, p0, A, b, Aeq, beq, lb, ub, nonlcon, options);
sol_list{1} = sol;
rmse_list(1) = fval;

for i = 2:num_iters
    fprintf('[%i/%i] \n', i, num_iters)
    
    % random intialization
    running_p0 = zeros(length(p0), 1);
    for j = 1:length(p0)
        unif_dist = makedist('Uniform', 'lower',lb(j), 'upper',ub(j));
        running_p0(j) = random(unif_dist, 1);
    end
    
    % optimization
    [sol, fval] = fmincon(opt_fun, running_p0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    sol_list{i} = sol;
    rmse_list(i) = fval;
end

[~, best_fit_idx] = min(rmse_list);
sol = sol_list{best_fit_idx};

% visualization calibration result
volt = volts(1);
[~, ~, ~, ~, yksum_hat] = kcurrent_model(sol, hold_volt, volt, time_space, Ek);

figure(2)
plot(t, yksum_hat)
hold on
for i=2:num_volts
    volt = volts(i);
    [~, ~, ~, ~, yksum_hat] = kcurrent_model(sol, hold_volt, volt, time_space, Ek);
    plot(t, yksum_hat)
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

%% plot for each voltage
clc
close all

peak = zeros(num_volts, 1);
for i = 1:num_volts
    volt = volts(i);
    
    [~, ~, ~, ~, yksum_hat] = kcurrent_model(sol, hold_volt, volt, time_space, Ek);
    peak(i) = max(yksum_hat);
    
    figure(i)
    plot(t, yksum(:, i))
    hold on
    plot(t, yksum_hat)
    hold off
    axis tight
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    legend('Experimental','Model Prediction')
end
