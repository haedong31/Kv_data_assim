%% WT
% 4.5-sec voltage-dependent data
clc
close all
clear variables

% load experimental data
trace_data = table2array(readtable('./4.5s-avg-wt.csv'));
t = trace_data(:,1);

volts = -50:10:50;
yksum = trace_data(:, 2:end);
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

% optimization
options = optimoptions('fmincon', 'MaxFunctionEvaluations',6e+3);

opt_fun1 = @(p) obj_current_stats(p, hold_volt, hold_idx, volts, t, yksum, Ek);
[p1, ~] = fmincon(opt_fun1, p0, A, b, Aeq, beq, lb, ub, nonlcon);

opt_fun1 = @(p) obj_rmse_over_volt(p, hold_volt, hold_idx, volts, t, yksum, Ek, true);
[sol, ~] = fmincon(opt_fun1, p1, A, b, Aeq, beq, lb, ub, nonlcon, options);

% visualization - warm up
volt_idx = 1;
volt = volts(volt_idx);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx(volt_idx));
pulse_t = t(hold_idx(volt_idx)+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

[~, ~, ~, ~, yksum] = reduced_model(p1, hold_volt, volt, time_space, Ek);
figure(2)
plot(t, yksum)
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

    [~, ~, ~, ~, yksum] = reduced_model(p1, hold_volt, volt, time_space, Ek);
    plot(t, yksum)
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

% visualization - final
volt_idx = 1;
volt = volts(volt_idx);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx(volt_idx));
pulse_t = t(hold_idx(volt_idx)+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

[~, ~, ~, ~, yksum] = reduced_model(sol, hold_volt, volt, time_space, Ek);
figure(3)
plot(t, yksum)
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

    [~, ~, ~, ~, yksum] = reduced_model(sol, hold_volt, volt, time_space, Ek);
    plot(t, yksum)
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

% save parameters
% save('wt_param.mat', 'p2')
% 
% volt_idx = 11;
% time_space = cell(1,3);
% time_space{1} = t;
% time_space{2} = t(1:hold_idx(volt_idx));
% time_space{3} = t(hold_idx(volt_idx)+1:end) - t(hold_idx(volt_idx)+1);
% 
% [ykto, ykslow1, ykslow2, ykss] = reduced_model(p2, hold_volt, volts(volt_idx), time_space, Ek);
% 
% plot(t, ykto)
% hold on
% plot(t, ykslow1)
% plot(t, ykslow2)
% plot(t, ykss)
% hold off
% title('WT with 50 mV')
% xlabel('Time (ms)')
% ylabel('Current (pA/pF)')
% legend("I_{Kto}", "I_{Kslow1}", "I_{Kslow2}", "I_{Kss}");

%% same routine for Mgat1KO
% 4.5-sec voltage-dependent data
clc
close all
clear variables

% load experimental data
trace_data = table2array(readtable('./4.5s-avg-ko.csv'));
t = trace_data(:,1);

volts = -50:10:50;
yksum = trace_data(:, 2:end);
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

% optimization
options = optimoptions('fmincon', 'MaxFunctionEvaluations',6e+3);

opt_fun1 = @(p) obj_current_stats(p, hold_volt, hold_idx, volts, t, yksum, Ek);
[p1, ~] = fmincon(opt_fun1, p0, A, b, Aeq, beq, lb, ub, nonlcon);

opt_fun1 = @(p) obj_rmse_over_volt(p, hold_volt, hold_idx, volts, t, yksum, Ek, true);
[sol, ~] = fmincon(opt_fun1, p1, A, b, Aeq, beq, lb, ub, nonlcon, options);

% visualization - warm up
volt_idx = 1;
volt = volts(volt_idx);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx(volt_idx));
pulse_t = t(hold_idx(volt_idx)+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

[~, ~, ~, ~, yksum] = reduced_model(p1, hold_volt, volt, time_space, Ek);
figure(2)
plot(t, yksum)
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

    [~, ~, ~, ~, yksum] = reduced_model(p1, hold_volt, volt, time_space, Ek);
    plot(t, yksum)
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

% visualization - final
volt_idx = 1;
volt = volts(volt_idx);

time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:hold_idx(volt_idx));
pulse_t = t(hold_idx(volt_idx)+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

[~, ~, ~, ~, yksum] = reduced_model(sol, hold_volt, volt, time_space, Ek);
figure(3)
plot(t, yksum)
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

    [~, ~, ~, ~, yksum] = reduced_model(sol, hold_volt, volt, time_space, Ek);
    plot(t, yksum)
end
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

