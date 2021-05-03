%% initial value for tri-exponential fitting
clc
close all
clear variables

% amplitude and tau ratios
wt_amp = [24.8, 17.1, 7.3, 3.7];
wt_amp./sum(wt_amp)

ko_amp = [17.6, 3.1, 3.6, 3.9];
ko_amp./sum(ko_amp)

wt_tau = [105.2, 1119.6, 7271.9];
wt_tau./sum(wt_tau)

ko_tau = [111.2, 1115.1, 11266.1];
ko_tau./sum(ko_tau)

% check the experimental data
exp_ksum = table2array(readtable('./4.5s-avg-ko.csv'));
t = exp_ksum(:,1);
yksum = exp_ksum(:, 2:end);
[~, num_volts] = size(yksum);

hold_idx = zeros(num_volts,1);
for i = 1:num_volts
    current_trc = yksum(:,i);
    [~, peak_idx] = max(current_trc);

    early_current_trc = current_trc(1:peak_idx);
    stable_val = min(early_current_trc);
    hold_idx(i) = find(early_current_trc == stable_val, 1, 'last');
end

volt_idx = 10;
[amp, tau] = trace_stat(t, yksum(:,volt_idx), hold_idx(volt_idx));
plot(t, yksum(:,volt_idx))

%% preprocessing for 4.5-2-avg-ko data
clc
close all
clear variables

exp_ksum = readtable('./4.5s-avg-ko.csv');
col_names = exp_ksum.Properties.VariableNames;

exp_ksum = table2array(exp_ksum);

t = exp_ksum(:,1);
[~, num_volts] = size(exp_ksum);
num_volts = num_volts - 1;

% remove sharp negative peaks at the early phase
for i=1:num_volts
    y = exp_ksum(:,i+1);
    y(y < 0) = 0;
    exp_ksum(:,i+1) = y;
end

% visualize experimental data 
plot(t, exp_ksum(:,1+1))
hold on
for i=2:num_volts
    plot(t, exp_ksum(:,i+1))
end
hold off

exp_ksum = array2table(exp_ksum, 'VariableNames',col_names);
writetable(exp_ksum,'./4.5s-avg-ko.csv');

%% test reduced_model.m
clc
close all
clear variables

% set p with default values
kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow10 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow20 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];

fixed_kto_idx = [4, 7, 8, 9, 11, 12, 15];
fixed_kslow1_idx = [4, 5, 6, 7, 10, 11];
fixed_kslow2_idx = [4, 5, 6, 7, 10];
fixed_kss_idx = [3, 4, 5];

tune_kto_idx = setdiff(1:17, fixed_kto_idx);
tune_kslow1_idx = setdiff(1:13, fixed_kslow1_idx);
tune_kslow2_idx = setdiff(1:11, fixed_kslow2_idx);
tune_kss_idx = setdiff(1:7, fixed_kss_idx);

p1 = kto0(tune_kto_idx);
p2 = kslow10(tune_kslow1_idx);
p = [p1, p2, kslow20(9), kss0(6:7)];

% protocol
hold_volt = -70;

volt = 50;

time_space = cell(1,3);
hold_pt = 100;
end_pt = 4.5*1000;
hold_t = 0:hold_pt;
pulse_t = (hold_pt + 1):end_pt;
pulse_t_adj = pulse_t - pulse_t(1);
t = [hold_t, pulse_t];
time_space{1} = t;
time_space{2} = hold_t;
time_space{3} = pulse_t_adj;

Ek = -91.1;

[ykto, ykslow1, ykslow2, ykss, yksum] = reduced_model(p, hold_volt, volt, time_space, Ek);

plot(t, ykto)
hold on
plot(t, ykslow1)
plot(t, ykslow2)
plot(t, ykss)
plot(t, yksum, 'Color','red', 'LineWidth',2)
hold off
legend('ykto','ykslow1','ykslow2','ykss','yksum')
