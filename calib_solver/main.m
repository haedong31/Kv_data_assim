%% model calibration phase 1 
% for Mgat1KO
% fully parameterized model
clc
close all
clear variables

% data to calculate objectives
exp_ksum = readtable('./ko1.csv');
exp_ksum.current = exp_ksum.current./254.3;
peak_vals = [17.6, 3.1, 3.6, 3.9];

% phase1: peak values
hold_volt = -70;
volt = 50;

time_space = cell(1,3);
hold_pt = 450;
end_pt = 25*1000;

hold_t = 0:hold_pt;
pulse_t = (hold_pt + 1):end_pt;
pulse_t_adj = pulse_t - pulse_t(1);
t = [hold_t, pulse_t];

time_space{1} = t;
time_space{2} = hold_t;
time_space{3} = pulse_t_adj;

Ek = -91.1;

fun1 = @(p) obj_peak(p, hold_volt, volt, time_space, Ek, peak_vals);

% initialization with default values; note shared parameters of Ikslow1, Ikslow2, and Ikss
% 17 + 13 + 2 + 4 = 36 parameters
% p0_kto = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
% p0_kslow1 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
% p0_kslow2 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
% p0_kss = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];
p0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387, ...
    22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496, ...
    5334, 4912, ...
    0.0862, 1235.5, 13.17, 0.0428];

% fix phophorylation rates
pp0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.14067, 0.387, ...
    22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.05766, 0.07496, ...
    5334, 4912, ...
    0.0862, 1235.5, 13.17, 0.0428];

% box constraints
lb = pp0 - pp0*5;
ub = pp0 + pp0*5;

% adjust lower bound for parameters that has to be non-zero
lb(6) = 1;
lb(7:16) = eps;

lb(20:21) = 1;
lb(22:24) = eps;
lb(25) = 10;
lb(26) = 1;
lb(27:28) = eps;

lb(29) = 10;
lb(30) = 1;

lb(31) = eps;
lb(32) = 1;
lb(33:34) = eps;

% other constraints
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

[p1, fval1] = fmincon(fun1, pp0, A, b, Aeq, beq, lb, ub, nonlcon);
disp(fval1)

% simulate with found solution
[ykto, ykslow1, ykslow2, ykss, yksum] = full_model(p1, hold_volt, volt, time_space, Ek);

figure(1)
plot(t, yksum, "--", "Color","black", "LineWidth",1.5)
hold on
plot(exp_ksum.time, exp_ksum.current, "Color","red", "LineWidth",1.5)
title('Phase 1')
legend('Simulated', 'Experimenta')
hold off

figure(2)
plot(t, ykto)
hold on
plot(t, ykslow1)
plot(t, ykslow2)
plot(t, ykss)
hold off
title('Phase 1')
legend("I_{Kto}", "I_{Kslow1}", "I_{Kslow2}", "I_{Kss}");

% phase 2: RMSE
t = exp_ksum.time;
ksum = exp_ksum.current;

[~, peak_idx] = max(ksum);
early_ksum = ksum(1:peak_idx);
stable_val = min(early_ksum);
hold_idx = find(early_ksum == stable_val, 1, 'last');

time_space = cell(1, 3);
time_space{1} = t;

tH = t(1:hold_idx);
time_space{2} = tH;

tP1 = t(hold_idx+1:end);
tP1_adj = tP1 - tP1(1);
time_space{3} = tP1_adj;

obj_rmse(p1, hold_volt, volt, time_space, Ek, ksum);
fun2 = @(p) obj_rmse(p, hold_volt, volt, time_space, Ek, ksum);

[p2, fval2] = fmincon(fun2, p1, A, b, Aeq, beq, lb, ub, nonlcon);
disp(fval2)

% simulate with found solution
[ykto, ykslow1, ykslow2, ykss, yksum] = full_model(p2, hold_volt, volt, time_space, Ek);

figure(3)
plot(t, yksum, "--", "Color","black", "LineWidth",1.5)
hold on
plot(exp_ksum.time, exp_ksum.current, "Color","red", "LineWidth",1.5)
title('Phase 2')
legend('Simulated', 'Experimenta')
hold off

figure(4)
plot(t, ykto)
hold on
plot(t, ykslow1)
plot(t, ykslow2)
plot(t, ykss)
hold off
title('Phase 2')
legend("I_{Kto}", "I_{Kslow1}", "I_{Kslow2}", "I_{Kss}");

%% reduced model
clc
close all
clear variables
