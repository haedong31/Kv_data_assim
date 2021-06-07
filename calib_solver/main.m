%% WT
% 4.5-sec voltage-dependent data
clc
close all
clear variables

% load experimental data
trace_data = table2array(readtable('./data/wt-preprocessed/15o29009.xlsx'));
t = trace_data(:,1);

volt_steps = 5:11;
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
ideal_hold_time = 0.120*1000;
[~, ideal_hold_idx] = min(abs(t-ideal_hold_time));
init_stable_val = sqrt(var(yksum(ideal_hold_idx,:)));

hold_idx = zeros(num_volts,1);
for i = 1:num_volts
    trc = yksum(:,i);

    counter = 1;
    statbility_est = abs(trc(ideal_hold_idx+counter) - trc(ideal_hold_idx-counter));
    if statbility_est > 10*init_stable_val
        hold_idx(i) = ideal_hold_idx;
        continue
    end
    
    while true
        % update counter and stable value
        counter = counter + 1;
        stable_val = statbility_est;
        
        % check stability
        statbility_est = abs(trc(ideal_hold_idx+counter) - trc(ideal_hold_idx-counter));
        if statbility_est > 10*stable_val
            hold_idx(i) = ideal_hold_idx + (counter-1);
            break
        end
    end
end

% protocol
hold_volt = -70;
Ek = -91.1;

% optimization constraints 
goal = [0, 0];
weight = [1, 1];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];

p0 = [33,15.5,20,8,7,0.3956,0.00095,0.051335,0.14067,0.387,22.5,45.2,40,2.058,803,0.05766,0.07496,5334,13.17,0.0428];
lb = p0 - 10*p0;
ub = p0 + 10*p0;
lb([5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 19, 20]) = eps;

low = p0 - 2*p0;
high = p0 + 2*p0;
low([5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 19, 20]) = eps;

options = optimoptions('fgoalattain', 'MaxFunctionEvaluations',1e+6);
% obj_multi2(p0, hold_volt, hold_idx, volts, t, yksum, Ek, true);
opt_fun = @(p) [obj_ks(p, hold_volt, hold_idx, volts, t, yksum, Ek, true);
    obj_rmse(p, hold_volt, hold_idx, volts, t, yksum, Ek, true, false)];
% opt_fun = @(p) obj_multi(p, hold_volt, hold_idx, volts, t, yksum, Ek, true, false);
% opt_fun = @(p) obj_rmse(p, hold_volt, hold_idx, volts, t, yksum, Ek, true, false);

num_iters = 1;
rmse_list = zeros(num_iters, 1);
sol_list = cell(num_iters, 1);

% first run with p0
[sol, fval] = fgoalattain(opt_fun, p0, goal, weight, A, b, Aeq, beq, lb, ub, nonlcon, options);
sol_list{1} = sol;
rmse_list(1) = fval;

for i = 2:num_iters
    fprintf('[%i/%i] \n', i, num_iters)
    
    % random intialization
    running_p0 = zeros(length(p0), 1);
    for j = 1:length(p0)
        unif_dist = makedist('Uniform', 'lower',low(j), 'upper',high(j));
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
volt_idx = 1;
volt = volts(volt_idx);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx(volt_idx));
pulse_t = t(hold_idx(volt_idx)+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

[~, ~, ~, ~, yksum_hat] = reduced_model(sol, hold_volt, volt, time_space, Ek);
figure(2)
plot(t, yksum_hat)
hold on
for i=2:num_volts
    volt_idx = i;
    volt = volts(volt_idx);
    
    time_space = cell(1,3);
    time_space{1} = t;
    time_space{2} = t(1:hold_idx(volt_idx));
    pulse_t = t(hold_idx(volt_idx)+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    time_space{3} = pulse_t_adj;

    [~, ~, ~, ~, yksum_hat] = reduced_model(sol, hold_volt, volt, time_space, Ek);
    plot(t, yksum_hat)
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

%% plot for each voltage
clc
close all

for i = 1:num_volts
    volt = volts(i);
    
    time_space = cell(1,3);
    time_space{1} = t;
    time_space{2} = t(1:hold_idx(i));
    pulse_t = t(hold_idx(i)+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    time_space{3} = pulse_t_adj;
    
    [~, ~, ~, ~, yksum_hat] = reduced_model(sol, hold_volt, volt, time_space, Ek);
    
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

%% same routine for Mgat1KO
% 4.5-sec voltage-dependent data
clc
close all
clear variables

% load experimental data
trace_data = table2array(readtable('./data/ko-preprocessed/15n24005.xlsx'));
t = trace_data(:,1);

volt_steps = 5:11;
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
ideal_hold_time = 0.120*1000;
[~, ideal_hold_idx] = min(abs(t-ideal_hold_time));
init_stable_val = sqrt(var(yksum(ideal_hold_idx,:)));

hold_idx = zeros(num_volts,1);
for i = 1:num_volts
    trc = yksum(:,i);

    counter = 1;
    statbility_est = abs(trc(ideal_hold_idx+counter) - trc(ideal_hold_idx-counter));
    if statbility_est > 10*init_stable_val
        hold_idx(i) = ideal_hold_idx;
        continue
    end
    
    while true
        % update counter and stable value
        counter = counter + 1;
        stable_val = statbility_est;
        
        % check stability
        statbility_est = abs(trc(ideal_hold_idx+counter) - trc(ideal_hold_idx-counter));
        if statbility_est > 10*stable_val
            hold_idx(i) = ideal_hold_idx + (counter-1);
            break
        end
    end
end

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

p0 = [33,15.5,20,8,7,0.3956,0.00095,0.051335,0.14067,0.387,22.5,45.2,40,2.058,803,0.05766,0.07496,5334,13.17,0.0428];
lb = p0 - 10*p0;
ub = p0 + 10*p0;
lb([5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 19, 20]) = eps;

low = p0 - 2*p0;
high = p0 + 2*p0;
low([5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 19, 20]) = eps;

options = optimoptions('fmincon', 'MaxFunctionEvaluations',1e+5);
opt_fun = @(p) obj_multi(p, hold_volt, hold_idx, volts, t, yksum, Ek, true, false);
% opt_fun = @(p) obj_rmse(p, hold_volt, hold_idx, volts, t, yksum, Ek, true, false);

num_iters = 1;
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
        unif_dist = makedist('Uniform', 'lower',low(j), 'upper',high(j));
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
volt_idx = 1;
volt = volts(volt_idx);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx(volt_idx));
pulse_t = t(hold_idx(volt_idx)+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

[~, ~, ~, ~, yksum_hat] = reduced_model(sol, hold_volt, volt, time_space, Ek);
figure(2)
plot(t, yksum_hat)
hold on
for i=2:num_volts
    volt_idx = i;
    volt = volts(volt_idx);
    
    time_space = cell(1,3);
    time_space{1} = t;
    time_space{2} = t(1:hold_idx(volt_idx));
    pulse_t = t(hold_idx(volt_idx)+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    time_space{3} = pulse_t_adj;

    [~, ~, ~, ~, yksum_hat] = reduced_model(sol, hold_volt, volt, time_space, Ek);
    plot(t, yksum_hat)
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

%% plot for each voltage
clc
close all

for i = 1:num_volts
    volt = volts(i);
    
    time_space = cell(1,3);
    time_space{1} = t;
    time_space{2} = t(1:hold_idx(i));
    pulse_t = t(hold_idx(i)+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    time_space{3} = pulse_t_adj;
    
    [~, ~, ~, ~, yksum_hat] = reduced_model(sol, hold_volt, volt, time_space, Ek);
    
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
