%% run kcurrent_model2 with default parameter values
clc
close all
clearvars

current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% tunning index in individual current models
tune_idx1_kto = 1:13;
tune_idx1_kslow1 = 1:11;
tune_idx1_kslow2 = 1:3;
tune_idx1_kss = 1:4;

% protocol
hold_volt = -70;
min_volt = -50;
volt = 0;
ek = -91.1;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;

idx_info1 = cell(1, num_currents);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        idx_info1{i} = tune_idx1_kto;
    case 'ikslow1'
        idx_info1{i} = tune_idx1_kslow1;
    case 'ikslow2'
        idx_info1{i} = tune_idx1_kslow2;
    case 'ikss'
        idx_info1{i} = tune_idx1_kss;
    case 'ikur'
        idx_info1{i} = tune_idx1_kur;
    case 'ik1'
        idx_info1{i} = tune_idx1_k1;
    end
end

% tunning index in decision variable, p
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

% default values
kto_default = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.3846];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, 6.1, 18.0, 2.058, 803.0, 0.16];
kslow2_default = [4912, 5334, 0.16];
kss_default = [0.0862, 1235.5, 13.17, 0.0611];
kur_default = [270, 1050, 0];
k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volt;
volt_space{3} = ek;

p0 = zeros(1, cul_idx_len);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        p0(idx_info2{i}) = kto_default(tune_idx1_kto);
    case 'ikslow1'
        p0(idx_info2{i}) = kslow1_default(tune_idx1_kslow1);
    case 'ikslow2'
        p0(idx_info2{i}) = kslow2_default(tune_idx1_kslow2);
    case 'ikss'
        p0(idx_info2{i}) = kss_default(tune_idx1_kss);
    case 'ikur'
        p0(idx_info2{i}) = kur_default(tune_idx1_kur);
    case 'ik1'
        p0(idx_info2{i}) = k1_default(tune_idx1_kur);
    end    
end

time_space = cell(1, 3);
t = 1:ideal_end_time;
time_space{1} = t;
time_space{2} = t(1:ideal_hold_time);
pulse_t = t(ideal_hold_time+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

protocol = cell(4, 1);
protocol{1} = volt_space{1};
protocol{2} = volt_space{2};
protocol{3} = time_space;
protocol{4} = volt_space{3};

[~, comp_currents] = kcurrent_model2(p0, model_struct, protocol);

figure('Color','w', 'Position',[100,100,550,500])
subplot(2,2,1)
plot(t,comp_currents{1},'LineWidth',2,'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

subplot(2,2,2)
plot(t,comp_currents{2}, 'LineWidth',2, 'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

subplot(2,2,3)
plot(t,comp_currents{3},'LineWidth',2,'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

subplot(2,2,4)
plot(t,comp_currents{4}, 'LineWidth',2, 'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

% in one figure
iksum = comp_currents{1} + comp_currents{2} + comp_currents{3} + comp_currents{4};

figure('Color','w')
plot(t,comp_currents{1}, '-', 'Color','black', 'LineWidth',1.5)
hold on
plot(t,comp_currents{2}, '--','Color','black', 'LineWidth',1.5)
plot(t,comp_currents{3}, ':','Color','black', 'LineWidth',1.5)
plot(t,comp_currents{4}, '-.','Color','black', 'LineWidth',1.5)
plot(t,iksum, 'Color','red', 'LineWidth',1.5)
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
legend('I_{Kto,f}','I_{Kslow1}','I_{Kslow2}','I_{Kss}','I_{Ksum}')
set(gca,'FontName','Arial','FontSize',11','FontWeight','bold')

%% run kcurrent_model1 with default parameter values
clc
close all
clear variables

current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% tunning index in individual current models
tune_idx1_kto = [1, 2, 6, 10, 13, 14, 15, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 4, 5, 8, 9, 11, 12, 13];
tune_idx1_kslow2 = [1, 3];
tune_idx1_kss = [3, 4];
tune_idx1_kur = [1, 3];
tune_idx1_k1 = [1, 3, 5, 7];

% protocol
hold_volt = -70;
min_volt = -50;
volt = 0;
ek = -91.1;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;

idx_info1 = cell(1, num_currents);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        idx_info1{i} = tune_idx1_kto;
    case 'ikslow1'
        idx_info1{i} = tune_idx1_kslow1;
    case 'ikslow2'
        idx_info1{i} = tune_idx1_kslow2;
    case 'ikss'
        idx_info1{i} = tune_idx1_kss;
    case 'ikur'
        idx_info1{i} = tune_idx1_kur;
    case 'ik1'
        idx_info1{i} = tune_idx1_k1;
    end
end

% tunning index in decision variable, p
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

% default values
kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow2_default = [5334, 4912, 0.05766];
kss_default = [0.0862, 1235.5, 13.17, 0.0428];
kur_default = [270, 1050, 0];
k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volt;
volt_space{3} = ek;

p0 = zeros(1, cul_idx_len);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        p0(idx_info2{i}) = kto_default(tune_idx1_kto);
    case 'ikslow1'
        p0(idx_info2{i}) = kslow1_default(tune_idx1_kslow1);
    case 'ikslow2'
        p0(idx_info2{i}) = kslow2_default(tune_idx1_kslow2);
    case 'ikss'
        p0(idx_info2{i}) = kss_default(tune_idx1_kss);
    case 'ikur'
        p0(idx_info2{i}) = kur_default(tune_idx1_kur);
    case 'ik1'
        p0(idx_info2{i}) = k1_default(tune_idx1_kur);
    end    
end

time_space = cell(1, 3);
t = 1:ideal_end_time;
time_space{1} = t;
time_space{2} = t(1:ideal_hold_time);
pulse_t = t(ideal_hold_time+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

protocol = cell(4, 1);
protocol{1} = volt_space{1};
protocol{2} = volt_space{2};
protocol{3} = time_space;
protocol{4} = volt_space{3};

[yksum_hat, comp_currents] = kcurrent_model1(p0, model_struct, protocol);

figure('Color','w', 'Position',[100,100,550,500])
subplot(2,2,1)
plot(t,comp_currents{1},'LineWidth',2,'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

subplot(2,2,2)
plot(t,comp_currents{2}, 'LineWidth',2, 'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

subplot(2,2,3)
plot(t,comp_currents{3},'LineWidth',2,'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

subplot(2,2,4)
plot(t,comp_currents{4}, 'LineWidth',2, 'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',10','FontWeight','bold')
box off

% in one figure
iksum = comp_currents{1} + comp_currents{2} + comp_currents{3} + comp_currents{4};

figure('Color','w')
plot(t,comp_currents{1}, '-', 'Color','black', 'LineWidth',1.5)
hold on
plot(t,comp_currents{2}, '--','Color','black', 'LineWidth',1.5)
plot(t,comp_currents{3}, ':','Color','black', 'LineWidth',1.5)
plot(t,comp_currents{4}, '-.','Color','black', 'LineWidth',1.5)
plot(t,iksum, 'Color','red', 'LineWidth',1.5)
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
legend('I_{Kto,f}','I_{Kslow1}','I_{Kslow2}','I_{Kss}','I_{Ksum}')
set(gca,'FontName','Arial','FontSize',11','FontWeight','bold')

%% test trust-region-reflective
% options = optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true);
options = optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true);

fun = @rosenbrockwithgrad;
x0 = [-1,2];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-2,-2];
ub = [2,2];
nonlcon = [];
[x, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
disp(fval)

%% test obj_rmse_grad
clc
close all
clear variables

warning('off', 'all')

% code arguments for calibration
group_name = 'wt';

% selection of currents
current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% voltages
volt_range = 3:11;

% tunning index in individual current models
tune_idx1_kto = [1, 2, 6, 10, 13, 14, 15, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 4, 5, 8, 9, 11, 12, 13];
tune_idx1_kslow2 = [1, 3];
tune_idx1_kss = [3, 4];
tune_idx1_kur = [1, 3];
tune_idx1_k1 = [1, 3, 5, 7];

% protocol
hold_volt = -70;
min_volt = -50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
ek = -91.1;

volts = NaN(length(volt_range), 1);
for i = 1:length(volt_range)
    volts(i) = min_volt + (volt_range(i)-1)*10;
end

idx_info1 = cell(1, num_currents);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        idx_info1{i} = tune_idx1_kto;
    case 'ikslow1'
        idx_info1{i} = tune_idx1_kslow1;
    case 'ikslow2'
        idx_info1{i} = tune_idx1_kslow2;
    case 'ikss'
        idx_info1{i} = tune_idx1_kss;
    case 'ikur'
        idx_info1{i} = tune_idx1_kur;
    case 'ik1'
        idx_info1{i} = tune_idx1_k1;
    end
end

% tunning index in decision variable, p
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

% matching table
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
len_loop_idx = length(loop_idx);

% default values
kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow2_default = [5334, 4912, 0.05766];
kss_default = [0.0862, 1235.5, 13.17, 0.0428];
kur_default = [270, 1050, 0];
k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

global lb_kto
global lb_kslow1
global lb_kslow2
global lb_kss
global lb_kur
global lb_k1

lb_kto = [-70, -70, -70, -10, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps];
lb_kslow1 = [-70, -70, -70, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps];
lb_kslow2 = [5000, eps, eps];
lb_kss = [eps, eps, eps, eps];
lb_kur = [eps, 500, eps];
lb_k1 = [eps, -20, eps, -30, eps, eps, eps, eps, eps, eps];

global ub_kto
global ub_kslow1
global ub_kslow2
global ub_kss
global ub_kur
global ub_k1

ub_kto = [70, 70, 70, 50, 50, 30, 1, 1, 1, 10, 0.005, 0.3, 0.005, 0.5, 1, 1, 1];
ub_kslow1 = [50, 50, 50, 10, 50, 50, 1, 100, 1000, 50, 1, 1, 1];
ub_kslow2 = [10000, 5000, 1];
ub_kss = [1, 2000, 100, 1];
ub_kur = [500, 2000, 1];
ub_k1 = [120, 30, 1000, 30, 3, 1, 2, 0.25, 0.25, 1];

% voltage info
volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = ek;

% optimization constraints 
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

max_param_len = 0;
p0 = zeros(1, cul_idx_len);
lb = zeros(1, cul_idx_len);
ub = zeros(1, cul_idx_len);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        p0(idx_info2{i}) = kto_default(tune_idx1_kto);
        lb(idx_info2{i}) = lb_kto(tune_idx1_kto);
        ub(idx_info2{i}) = ub_kto(tune_idx1_kto);
        if max_param_len < length(kto_default)
            max_param_len = length(kto_default);
        end
    case 'ikslow1'
        p0(idx_info2{i}) = kslow1_default(tune_idx1_kslow1);
        lb(idx_info2{i}) = lb_kslow1(tune_idx1_kslow1);
        ub(idx_info2{i}) = ub_kslow1(tune_idx1_kslow1);
        if max_param_len < length(kslow1_default)
            max_param_len = length(kslow1_default);
        end
    case 'ikslow2'
        p0(idx_info2{i}) = kslow2_default(tune_idx1_kslow2);
        lb(idx_info2{i}) = lb_kslow2(tune_idx1_kslow2);
        ub(idx_info2{i}) = ub_kslow2(tune_idx1_kslow2);
        if max_param_len < length(kslow2_default)
            max_param_len = length(kslow2_default);
        end    
    case 'ikss'
        p0(idx_info2{i}) = kss_default(tune_idx1_kss);
        lb(idx_info2{i}) = lb_kss(tune_idx1_kss);
        ub(idx_info2{i}) = ub_kss(tune_idx1_kss);
        if max_param_len < length(kss_default)
            max_param_len = length(kss_default);
        end    
    case 'ikur'
        p0(idx_info2{i}) = kur_default(tune_idx1_kur);
        lb(idx_info2{i}) = lb_kur(tune_idx1_kur);
        ub(idx_info2{i}) = ub_kur(tune_idx1_kur);
        if max_param_len < length(kur_default)
            max_param_len = length(kur_default);
        end
    case 'ik1'
        p0(idx_info2{i}) = k1_default(tune_idx1_kur);
        lb(idx_info2{i}) = lb_k1(tune_idx1_k1);
        ub(idx_info2{i}) = ub_k1(tune_idx1_k1);
        if max_param_len < length(k1_default)
            max_param_len = length(k1_default);
        end
    end    
end

% main loop
for l = 1:len_loop_idx
    i = loop_idx(l);

    % read data
    file_path = fullfile(pwd, 'data', strcat(group_name, '-preprocessed'), file_names{i});
    trace_data = table2array(readtable(file_path));

    t = trace_data(:, 1);
    yksum = trace_data(:, (volt_range+1));

    % estimate time points
    [~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
    [~, ideal_end_idx] = min(abs(t - ideal_end_time));

    t = t(1:ideal_end_idx);
    yksum = yksum(1:ideal_end_idx, :);

    % time space
    time_space = cell(1, 4);
    time_space{1} = t;
    time_space{2} = t(1:ideal_hold_idx);
    pulse_t = t(ideal_hold_idx+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    time_space{3} = pulse_t_adj;
    time_space{4} = ideal_hold_idx;
    time_space{5} = ideal_end_idx;
    
    % objective function
    obj_rmse_grad(p0, model_struct, volt_space, time_space, yksum);
end

%% calculate RMSEs from calibration results
clc
close all
clear variables

group_name = 'wt';
exp_name = strcat('calib_exp28_', group_name);

% selection of currents
current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% voltages
volt_range = 3:11;

% tunning index in individual current models
tune_idx1_kto = [1, 2, 4, 5, 7, 11, 13];
tune_idx1_kslow1 = [1, 2, 4, 5, 9, 10, 11];
tune_idx1_kslow2 = [2, 3];
tune_idx1_kss = [1, 2, 3, 4];
tune_idx1_kur = [1, 3];
tune_idx1_k1 = [1, 3, 5, 7];

% protocol
hold_volt = -70;
min_volt = -50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
ek = -91.1;

% voltage info
volts = NaN(length(volt_range), 1);
for i = 1:length(volt_range)
    volts(i) = min_volt + (volt_range(i)-1)*10;
end

volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = ek;

% create model structure
idx_info1 = cell(1, num_currents);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        idx_info1{i} = tune_idx1_kto;
    case 'ikslow1'
        idx_info1{i} = tune_idx1_kslow1;
    case 'ikslow2'
        idx_info1{i} = tune_idx1_kslow2;
    case 'ikss'
        idx_info1{i} = tune_idx1_kss;
    case 'ikur'
        idx_info1{i} = tune_idx1_kur;
    case 'ik1'
        idx_info1{i} = tune_idx1_k1;
    end
end

% tunning index in decision variable, p
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

% matching table
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

for l = 1:length(loop_idx)
    i = loop_idx(l);
    trace_data = table2array(readtable(fullfile(pwd, 'data', strcat(group_name, '-preprocessed'), file_names{i})));
    sol_mx = table2array(readtable(fullfile(pwd, exp_name, file_names{i})));

    t = trace_data(:, 1);
    yksum = trace_data(:, (volt_range+1));

    % estimate time points
    [~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
    [~, ideal_end_idx] = min(abs(t - ideal_end_time));

    t = t(1:ideal_end_idx);
    yksum = yksum(1:ideal_end_idx, :);

    % time space
    time_space = cell(1, 4);
    time_space{1} = t;
    time_space{2} = t(1:ideal_hold_idx);
    pulse_t = t(ideal_hold_idx+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    time_space{3} = pulse_t_adj;
    time_space{4} = ideal_hold_idx;
    time_space{5} = ideal_end_idx;
    
    % solution vector
    sol = NaN(1, cul_idx_len);
    for j = 1:length(model_struct)
        running_sol = sol_mx(:, j);
        running_sol = running_sol(~isnan(running_sol));
        sol(model_struct(j).idx2) = running_sol(model_struct(j).idx1);
    end
    
    r = obj_rmse(sol, @kcurrent_model2, model_struct, volt_space, time_space, yksum);
    fprintf('[File %i/%i] %s Min RMSE: %f \n', l, length(loop_idx), file_names{i}, r)
end

%% check calibration result for general kcurrent model
clc
close all 
clear variables
format shortG

% specify model structure
current_names = {'ikto', 'ikslow1', 'ikss'};
num_currents = length(current_names);

tune_idx1_kto = [1, 2, 3, 5, 6, 15, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 9, 12, 13];
tune_idx1_kss = [6, 7];
% tune_idx1_k1 = [1, 3, 5, 7];
idx_info1 = {tune_idx1_kto, ...
    tune_idx1_kslow1, ...
    tune_idx1_kss};
    
%     tune_idx1_k1};

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
file_name = '15o27002.xlsx';
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
hold_idx = length(time_space{2});

protocol{1} = hold_volt;
protocol{3} = time_space;
protocol{4} = ek;

% generate yksum_hat
for i=1:length(volts)
    protocol{2} = volts(i);
    yksum_i = yksum(:, i);
    [yksum_hat, comp] = kcurrent_model(sol, model_struct, protocol);

    figure(i)
    plot(t, yksum_i)
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
clear variables

% [1, 2, 3, 5, 6, 10, 13, 14, 16, 17];
fixed_kto_idx = [4, 7, 8, 9, 11, 12, 15];
tune_kto_idx = setdiff(1:17, fixed_kto_idx);

% [1, 2, 3, 5, 8, 9, 12, 13];
fixed_kslow1_idx = [4, 6, 7, 10, 11];
tune_kslow1_idx = setdiff(1:13, fixed_kslow1_idx);

% 
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
    writetable(sol, fullfile(pwd, 'wt-re', file_names(i).name))
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
