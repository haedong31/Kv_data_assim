%% check calibration result for general kcurrent model
clc
close all 
clear variables
format shortG

% specify model structure
current_names = {'ikto', 'ikslow1', 'ikss', 'ik1'};
num_currents = length(current_names);

tune_idx1_kto = [1, 2, 15, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 9, 12, 13];
tune_idx1_kss = [4, 5, 6, 7];
tune_idx1_k1 = [1, 3, 5, 7];
idx_info1 = {tune_idx1_kto, ...
    tune_idx1_kslow1, ...
    tune_idx1_kss, ...
    tune_idx1_k1};

idx_info2 = cell(1, num_currents);
cul_idx_len = 0;
for i = 1: num_currents
    cul_idx_len_old = cul_idx_len;
    cul_idx_len = cul_idx_len + length(idx_info1{i});
    idx_info2{i} = (1+cul_idx_len_old):(cul_idx_len);
end

model_info = [current_names; idx_info1; idx_info2];
field_names = {'name', 'idx1', 'idx2'};
model_struct = cell2struct(model_info, field_names, 1);

% specify result file to check
file_group = 'ko';
file_name = '15n03002.xlsx';
calib = table2array(readtable(fullfile(pwd, file_group, strcat('calib_param_', file_name))));

sol = zeros(cul_idx_len, 1);
for i=1:num_currents
    running_calib = calib(:, i);
    running_calib(~isnan(running_calib));
    
    sol(idx_info2{i}) = running_calib(idx_info1{i});
end

% experimental data to be compared
exp_data = table2array(readtable(fullfile(pwd,'data', strcat(file_group, '-preprocessed2'), file_name)));
t = exp_data(:, 1);
yksum = exp_data(:, 2:end);

% protocol
protocol = cell(4, 1);

hold_volt = -70;
volts = -50:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
ek = -91.1;

[~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
[~, ideal_end_idx] = min(abs(t - ideal_end_time));

t = t(1:ideal_end_idx);
yksum = yksum(1:ideal_end_idx, :);

time_space = cell(1, 3);
time_space{1} = t;
time_space{2} = t(1:ideal_hold_idx);
pulse_t = t(ideal_hold_idx+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

protocol{1} = hold_volt;
protocol{3} = time_space;
protocol{4} = ek;

% generate yksum_hat
for i=1:length(volts)
    protocol{2} = volts(i);
    [yksum_hat, ~] = kcurrent_model(sol, model_struct, protocol);

    figure(i)
    plot(t, yksum(:, i))
    hold on
    plot(t, yksum_hat)
    hold off
    legend('Experimental','Model')
end

%% test general kcurrent model
clc
close all
clear variables
format shortG

% structure with three fields
% field 1: current names
current_names = {'ikto', 'ikslow1', 'ikss'};

% field 2: tunning index in individual parameter list
param_info1 = {[1, 2, 3, 5, 6, 10, 13, 14, 16, 17], ...
    [1, 2, 3, 5, 8, 9, 12, 13], ...
    [6, 7]};

% field 3: tunning index in p
len1 = length(param_info1{1});
len2 = length(param_info1{2});
len3 = length(param_info1{3});
param_info2 = {1:len1, ...
    (len1+1):(len1+len2), ...
    (len1+len2+1):(len1+len2+len3)};

model_info = [current_names; param_info1; param_info2];

% create a structure
row_names = {'name', 'idx1', 'idx2'};
model_struct = cell2struct(model_info, row_names, 1);

% protocol information
hold_volt = -70;
volts = -50:10:50;
Ek = -91.1;

time_space = cell(1, 3);
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;

t = 1:ideal_end_time;
[~, ideal_hold_idx] = min(abs(t - ideal_hold_time));

time_space{1} = t;
time_space{2} = t(1:ideal_hold_idx);
pulse_t = t(ideal_hold_idx+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

protocol = cell(4, 1);
protocol{1} = hold_volt;
% protocol{2} = volt;
protocol{3} = time_space;
protocol{4} = Ek;

p0 = [33, 15.5, 20, 8, 7, 0.3956, 0.00095, 0.051335, 0.14067, 0.387, ...
    22.5, 45.2, 40, 5.7, 2.058, 803, 0.05766, 0.07496, ...
    13.17, 0.0428];

hold on
for i = 1:length(volts)
    protocol{2} = volts(i);
    [yksum, ~] = kcurrent_model(p0, model_struct, protocol);
    
    plot(t, yksum)
end
hold off

% % toy example
% kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
%     0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
% kslow10 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
% kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];
% 
% param_kto = zeros(length(kto0), 1);
% param_kslow1 = zeros(length(kslow10), 1);
% param_kss0 = zeros(length(kss0), 1);
% 
% p = [100, 100, 200, 200, 300, 300];
% 
% current_names = model_struct.name;
% 
% tune_kto_idx1 = model_struct(1).idx1;
% tune_kslow1_idx1 = model_struct(2).idx1;
% tune_kss_idx1 = model_struct(3).idx1;
% 
% tune_kto_idx2 = model_struct(1).idx2;
% tune_kslow1_idx2 = model_struct(2).idx2;
% tune_kss_idx2 = model_struct(3).idx2;
% 
% fixed_kto_idx = setdiff(1:length(kto0), tune_kto_idx1);
% fixed_kslow1_idx = setdiff(1:length(kslow10), tune_kslow1_idx1);
% fixed_kss_idx = setdiff(1:length(kss0), tune_kss_idx1);
% 
% param_kto(fixed_kto_idx) = kto0(fixed_kto_idx);
% param_kto(tune_kto_idx1) = p(tune_kto_idx2);
% disp(param_kto)
% 
% param_kslow1(fixed_kslow1_idx) = kslow10(fixed_kslow1_idx);
% param_kslow1(tune_kslow1_idx1) = p(tune_kslow1_idx2);
% disp(param_kslow1)
% 
% param_kss(fixed_kss_idx) = kss0(fixed_kss_idx);
% param_kss(tune_kss_idx1) = p(tune_kss_idx2);
% disp(param_kss)

%% tunnning parameters
fixed_kto_idx = [4, 7, 8, 9, 11, 12, 15];
tune_kto_idx = setdiff(1:17, fixed_kto_idx);

fixed_kslow1_idx = [4, 6, 7, 10, 11];
tune_kslow1_idx = setdiff(1:13, fixed_kslow1_idx);

shared_ikslow2_idx = 1:8;
fixed_kslow2_idx = 9;
tune_kslow2_idx = setdiff(1:11, [shared_ikslow2_idx, fixed_kslow2_idx]);

shared_ikss_idx = 1:3;
fixed_kss_idx = [4, 5];
tune_kss_idx = setdiff(1:7, [shared_ikss_idx, fixed_kss_idx]);

%% fix excel files of exp4 calibration result of wt
clc
close all
clear variables

group = 'wt';
file_names = dir(fullfile(pwd, 'wt', '*.xlsx'));

for i=1:length(file_names)
    sol = readtable(fullfile(file_names(i).folder, file_names(i).name));
    sol = sol(:, [1, 2, 5, 6]);
    sol.Properties.VariableNames = {'ikto', 'ikslow1', 'ikss', 'ik1'};
    writetable(sol, 'test.xlsx')
end

%% transform calibrated solutions to fit into xlsx file
clc
close all
clear variables

matching_table = readtable('./data/matching-table-ko.xlsx');
file_names = matching_table.trace_file_name_4half;
load('sols_ko.mat')

% exclude row not having 4.5-sec data
loop_idx = [];
[num_files, ~] = size(matching_table);
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

% remove rows with all zeros
% sols = sols(any(sols, 2), :);

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

for i = loop_idx
    save_path = fullfile(pwd, strcat('calib_param_', file_names{i}));
    sol = sols(i, :);

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

    writematrix(sol_kto', save_path, "Sheet","Parameters", "Range","A2");
    writematrix(sol_kslow1', save_path, "Sheet","Parameters", "Range","B2");
    writematrix(sol_kslow2', save_path, "Sheet","Parameters", "Range","C2");
    writematrix(sol_kss', save_path, "Sheet","Parameters", "Range","D2");
end

%% check calibration result for kcurrent3
clc
close all
clear variables

% protocol
hold_volt = -70;
volts = -50:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
Ek = -91.1;

% experimental data
file_group = 'ko';
file_name = '15o27002.xlsx';

trace_data = table2array(readtable(fullfile( ...
    pwd, 'data', strcat(file_group, '-preprocessed2'), file_name)));
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

% calibrated parameter
calib_param = readtable(fullfile(pwd, file_group, ...
    strcat('calib_param_', file_name)));

% default values
kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow10 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow20 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
kur0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 270, 1050, 0];
kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];

param_kto = calib_param.IKto;
param_kslow1 = calib_param.IKslow1;
param_kslow1 = param_kslow1(~isnan(param_kslow1));
param_kslow2 = zeros(11,1);
param_kur = zeros(11,1);
param_kss = zeros(7,1);
param_k1 = calib_param.IK1;
param_k1 = param_k1(~isnan(param_k1));

% ikslow2
shared_ikslow2_idx = 1:8;
param_kslow2(shared_ikslow2_idx) = param_kslow1(shared_ikslow2_idx);
param_kslow2(9:11) = calib_param.IKslow2(~isnan(calib_param.IKslow2));

% ikur
shared_ikur_idx = 1:8;
param_kur(shared_ikur_idx) = param_kslow1(shared_ikur_idx);
param_kur(9:11) = calib_param.IKur(~isnan(calib_param.IKur));

% ikss
shared_ikss_idx = 1:3;
param_kss(shared_ikss_idx) = param_kslow1([1, 3, 4]);
param_kss(4:7) = calib_param.IKss(~isnan(calib_param.IKss));

num_volts = length(volts);
for i = 1:num_volts
    % generate K+ currents
    ykto = ikto(param_kto, hold_volt, volts(i), time_space, Ek);
    ykslow1 = ikslow1(param_kslow1, hold_volt, volts(i), time_space, Ek);
    ykslow2 = ikslow2(param_kslow2, hold_volt, volts(i), time_space, Ek);
    ykur = ikur(param_kur, hold_volt, volts(i), time_space, Ek);
    ykss = ikss(param_kss, hold_volt, volts(i), time_space, Ek);
    yk1 = ik1(param_k1, hold_volt, volts(i), time_space, Ek);
    yksum_hat = ykto + ykslow1 + ykslow2 + ykur + ykss + yk1;

    figure(i)
    plot(t, yksum(:, i))
    hold on
    plot(t, yksum_hat)
    hold off
end

%% test ik1.m & ikur.m
clc
close all
clear variables

% default values
p0_k1 = [59.215, 5.476, 594.31, 4.753];
p0_kur = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 270, 1050, 0.05766];

% current model arguments
hold_volt = -70;
volts = -50:10:50;
time_space = cell(1,3);
Ek = -91.1;

% number of voltages
num_volts = length(volts);

% time space
hold_pt = 100;
end_pt = 4.5*1000;
hold_t = 0:hold_pt;
pulse_t = (hold_pt + 1):end_pt;
pulse_t_adj = pulse_t - pulse_t(1);
t = [hold_t, pulse_t];

time_space{1} = t;
time_space{2} = hold_t;
time_space{3} = pulse_t_adj;

figure(1)
hold on
for i = 1:num_volts
    yk1 = ik1(p0_k1, hold_volt, volts(i), time_space, Ek);
    plot(t, yk1)    
end
hold off

figure(2)
hold on
for i = 1:num_volts
    ykur = ikur(p0_kur, hold_volt, volts(i), time_space, Ek);
    plot(t, ykur)    
end
hold off

%% test kcurrent_model2
clc
close all
clear variables

% default values
kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];

% current model arguments
hold_volt = -70;
volts = -50:10:50;
time_space = cell(1,3);
Ek = -91.1;

% number of voltages
num_volts = length(volts);

% time space
hold_pt = 100;
end_pt = 4.5*1000;
hold_t = 0:hold_pt;
pulse_t = (hold_pt + 1):end_pt;
pulse_t_adj = pulse_t - pulse_t(1);
t = [hold_t, pulse_t];

time_space{1} = t;
time_space{2} = hold_t;
time_space{3} = pulse_t_adj;

p0 = [33, 15.5, 20, 8, 7, 0.3956, 0.00095, 0.051335, 0.14067, 0.387, ...
    22.5, 45.2, 40, 5.7, 2.058, 803, 0.05766, 0.07496, ...
    13.17, 0.0428];

hold on
for i = 1:num_volts
    [ykto, ykslow, ykss, yksum] = kcurrent_model2(p0, hold_volt, volts(i), time_space, Ek);
    plot(t, ykto)    
end
hold off

%% test kcurrent_model
clc
close all
clear variables

% default values
kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow10 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow20 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];

% current model arguments
hold_volt = -70;
volts = -50:10:50;
time_space = cell(1,3);
Ek = -91.1;

% number of voltages
num_volts = length(volts);

% time space
hold_pt = 100;
end_pt = 4.5*1000;
hold_t = 0:hold_pt;
pulse_t = (hold_pt + 1):end_pt;
pulse_t_adj = pulse_t - pulse_t(1);
t = [hold_t, pulse_t];

time_space{1} = t;
time_space{2} = hold_t;
time_space{3} = pulse_t_adj;

% calculate variances of early phases
var_df = zeros(num_volts, 5);
for i = 1:num_volts
    ykto = ikto(kto0, hold_volt, volts(i), time_space, Ek);
    ykslow1 = ikslow1(kslow10, hold_volt, volts(i), time_space, Ek);
    ykslow2 = ikslow2(kslow20, hold_volt, volts(i), time_space, Ek);
    ykss = ikss(kss0, hold_volt, volts(i), time_space, Ek);
    yksum = ykto + ykslow1 + ykslow2 + ykss;

    var_df(i, 1) = var(ykto(1:hold_pt));
    var_df(i, 2) = var(ykslow1(1:hold_pt));
    var_df(i, 3) = var(ykslow2(1:hold_pt));
    var_df(i, 4) = var(ykss(1:hold_pt));
    var_df(i, 5) = var(yksum(1:hold_pt));
end

disp(mean(var_df, 1))

%% initial value for tri-exponential fitting
clc
close all
clear variables

% amplitude and tau ratios
wt_amp = [24.8, 17.1, 7.3, 3.7];
disp(wt_amp./sum(wt_amp))

ko_amp = [17.6, 3.1, 3.6, 3.9];
disp(ko_amp./sum(ko_amp))

wt_tau = [105.2, 1119.6, 7271.9];
disp(wt_tau./sum(wt_tau))

ko_tau = [111.2, 1115.1, 11266.1];
disp(ko_tau./sum(ko_tau))

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
[amp, time_stat] = trace_stat(t, yksum(:,volt_idx), hold_idx(volt_idx));
plot(t, yksum(:,volt_idx))

for i = 1:num_volts
    t(hold_idx(i))
end

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
