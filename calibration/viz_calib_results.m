%% multiple bar graph
clc
clearvars
close all

exp_num1 = "exp35";
exp_num2 = "exp36";
exp_num3 = "exp37";
group = "wt";

rmse_df1 = readtable(fullfile(pwd,"log",strcat(exp_num1,"_",group,".csv")));
rmse_df2 = readtable(fullfile(pwd,"log",strcat(exp_num2,"_",group,".csv")));
rmse_df3 = readtable(fullfile(pwd,"log",strcat(exp_num3,"_",group,".csv")));
rmse_df_comb = join(rmse_df1,rmse_df2, 'Keys','File');
rmse_df_comb = join(rmse_df_comb,rmse_df3, 'Keys','File');

file_names = categorical(rmse_df_comb.File);
rmse = table2array(rmse_df_comb(:,2:end));

figure('Color','w', 'Position',[50,50,900,770])
barh(file_names,rmse)
legend("Interior point","SQP","Active set")
set(gca, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')

%% bar graph of RMSE
clc
close all

file_group = "wt";
exp_num = "exp35";

rmse_df = readtable(fullfile(pwd,"log",strcat(exp_num,"_",file_group,".csv")));
file_names = categorical(rmse_df.File);
rmse = rmse_df.RMSE;

figure('Color','w', 'Position',[100,100,670,550])
barh(file_names,rmse, 'blue')
set(gca, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')

%% compare with experimental data (25-sec data)
clc
close all
clearvars

group = "ko";
exp_num = "exp41";
file_name = "15o27006.xlsx";

data_path = fullfile(pwd,"data",strcat(group,"-preprocessed-25s"),file_name);
calib_path = fullfile(pwd,strcat("calib_",exp_num,"_",group),file_name);

current_names = {'iktof', 'ikslow1', 'ikslow2', 'ikss'};
tune_idx1{1} = [1, 2, 4, 5, 7, 11, 13];
tune_idx1{2} = [1, 2, 4, 5, 9, 10, 11];
tune_idx1{3} = [2, 3];
tune_idx1{4} = [1, 2, 3, 4];

% import the experimental data and calibration results
exp_mx = table2array(readtable(data_path));
sol_mx = table2array(readtable(calib_path));
t = exp_mx(:,1);
yksum = exp_mx(:,2);

% protocol
ideal_hold_time = 470;
ideal_end_time = 25*1000;
hold_volt = -70;
volt = 50;
ek = -91.1;

% estimate the critical time points
[~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
[~, ideal_end_idx] = min(abs(t - ideal_end_time));
t = t(1:ideal_end_idx);
yksum = yksum(1:ideal_end_idx);

% time space
time_space = cell(1,3);
time_space{1} = t;
time_space{2} = t(1:ideal_hold_idx);
pulse_t = t(ideal_hold_idx+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

% protocol
protocol = cell(1,4);
protocol{1} = hold_volt;
protocol{2} = ek;
protocol{3} = volt;
protocol{4} = time_space;

% prediction from the calibration resukt
[model_info, psize] = gen_mdl_struct(current_names, tune_idx1);
sol = gen_sol_vec(sol_mx, model_info, psize);
[z,yksum_hat,comp_hat] = rmsey(sol, yksum, model_info, protocol);

% figure 1: yksum
figure('Color','w', 'Position',[100,100,560,420])
plot(t,yksum, 'Color','black')
hold on
plot(t, yksum_hat{1}, '--', 'Color','red');
hold off
set(gca, 'XLimSpec','tight')
legend("Experimental Data","Model Prediction")

% figure 2: individual K+ currents
comp_hat = comp_hat{1};
figure('Color','w', 'Position',[600,100,560,420])
plot(t,comp_hat{1})
hold on
plot(t,comp_hat{2})
plot(t,comp_hat{3})
plot(t,comp_hat{4})
hold off
set(gca, 'XLimSpec','tight')
legend("I_{Kto,f}","I_{Kslow1}","I_{Kslow2}","I_{Kss}")

%% compare with experimental data (kcurrent_model2)
clc
close all 
clearvars

file_group = 'ko';
exp_num = 'exp27';

% specify result files and voltages to check
file_name = '15o29024.xlsx';
% file_name = file_names(1);
save_dir = strcat('calib_', exp_num);

% voltages info
min_volt = -50;
volt_range = 3:11;
volts = NaN(length(volt_range), 1);
for i = 1:length(volt_range)
    volts(i) = min_volt + (volt_range(i)-1)*10;
end

% import calibration results
calib1 = table2array(readtable(fullfile(pwd, strcat(save_dir, '_', file_group), file_name)));

% specify model structure
% field 1: selection of currents
current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% field 2: tunning index in individual current models
tune_idx1{1} = [1, 2, 4, 5, 7, 11, 13];
tune_idx1{2} = [1, 2, 4, 5, 9, 10, 11];
tune_idx1{3} = [2, 3];
tune_idx1{4} = [3, 4];
tune_idx1{5} = [1, 3];
tune_idx1{6} = [1, 3, 5, 7];

[model_struct1, psize1] = gen_mdl_struct(current_names, tune_idx1);

sol1 = gen_sol_vec(calib1, model_struct1, psize1);

% experimental data to be compared
exp_data = table2array(readtable(fullfile(pwd,'data', strcat(file_group, '-preprocessed'), file_name)));
t = exp_data(:, 1);
yksum = exp_data(:, (volt_range+1));

% protocol
protocol = cell(4, 1);

hold_volt = -70;
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
time_space{4} = ideal_hold_idx;
time_space{5} = ideal_end_idx;

volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = ek;

protocol{1} = hold_volt;
protocol{3} = time_space;
protocol{4} = ek;

% generate yksum_hat
r = obj_rmse(sol1, @kcurrent_model2, model_struct1, volt_space, time_space, yksum);
fprintf('File %s Min RMSE: %f \n', file_name, r)

figure('Position',[50,50,800,800])
for i=1:length(volts)
    protocol{2} = volts(i);
    yksum_i = yksum(:, i);

    [yksum_hat1, ~] = kcurrent_model2(sol1, model_struct1, protocol);
    
    subplot(3, 3, i)
    plot(t, yksum_i, 'Color','red')
    hold on
    plot(t, yksum_hat1, '--', 'Color','green', 'LineWidth',2)
    hold off
    axis tight
    grid on
%    legend('Experimental Data','Model Prediction 1', 'Model Prediction 2', 'Location','best')
    title(strcat('Clamp Voltage: ', string(volts(i)), ' mV'), 'FontSize',10, 'FontWeight','bold');
    xlabel('Time (ms)', 'FontSize',10, 'FontWeight','bold');
    ylabel('Current (pA/pF)', 'FontSize',10, 'FontWeight','bold');
    get(gcf,'CurrentAxes');
    set(gca,'YDir','normal');
    set(gca,'LineWidth',2, 'FontSize',10, 'FontWeight','bold');
    set(gca,'GridLineStyle','--')
end

figure(2)
plot(t, yksum(:, 1))
hold on
for i = 2:length(volts)
    plot(t, yksum(:, i))
end
hold off

%% compare with experimental data (kcurrent_model1)
% clc
% close all 

file_group = 'wt';
exp_num = 'exp16';

% specify result files and voltages to check
file_name = '15o29024.xlsx';
% file_name = file_names(1);
save_dir = strcat('calib_', exp_num);

% voltages info
min_volt = -50;
volt_range = 3:11;
volts = NaN(length(volt_range), 1);
for i = 1:length(volt_range)
    volts(i) = min_volt + (volt_range(i)-1)*10;
end

% import calibration results
calib1 = table2array(readtable(fullfile(pwd, strcat(save_dir, '_', file_group), file_name)));

% specify model structure
% field 1: selection of currents
current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% field 2: tunning index in individual current models
tune_idx1 = cell(1, 6);
tune_idx1{1} = [1, 2, 6, 10, 13, 14, 15, 16, 17];
tune_idx1{2} = [1, 2, 3, 4, 5, 8, 9, 11, 12, 13];
tune_idx1{3} = [1, 3];
tune_idx1{4} = [3, 4];
tune_idx1{5} = [1, 3];
tune_idx1{6} = [1, 3, 5, 7];

% tune_idx1{1} = [1, 2, 6, 7, 9, 13, 14, 15, 16, 17];
% tune_idx1{2} = [1, 2, 4, 5, 8, 9, 11, 12, 13];
% tune_idx1{3} = [1, 3];
% tune_idx1{4} = [1, 2, 3, 4];
% tune_idx1{5} = [1, 3];
% tune_idx1{6} = [1, 3, 5, 7];

[model_struct1, psize1] = gen_mdl_struct(current_names, tune_idx1);

sol1 = gen_sol_vec(calib1, model_struct1, psize1);

% experimental data to be compared
exp_data = table2array(readtable(fullfile(pwd,'data', strcat(file_group, '-preprocessed'), file_name)));
t = exp_data(:, 1);
yksum = exp_data(:, (volt_range+1));

% protocol
protocol = cell(4, 1);

hold_volt = -70;
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
time_space{4} = ideal_hold_idx;
time_space{5} = ideal_end_idx;

volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = ek;

protocol{1} = hold_volt;
protocol{3} = time_space;
protocol{4} = ek;

% generate yksum_hat
r = obj_rmse(sol1, @kcurrent_model1, model_struct1, volt_space, time_space, yksum);
fprintf('File %s Min RMSE: %f \n', file_name, r)

figure('Position',[50,50,800,800])
for i=1:length(volts)
    protocol{2} = volts(i);
    yksum_i = yksum(:, i);

    [yksum_hat1, ~] = kcurrent_model1(sol1, model_struct1, protocol);
    
    subplot(3, 3, i)
    plot(t, yksum_i, 'Color','red')
    hold on
    plot(t, yksum_hat1, '--', 'Color','green', 'LineWidth',2)
    hold off
    axis tight
    grid on
%    legend('Experimental Data','Model Prediction 1', 'Model Prediction 2', 'Location','best')
    title(strcat('Clamp Voltage: ', string(volts(i)), ' mV'), 'FontSize',10, 'FontWeight','bold');
    xlabel('Time (ms)', 'FontSize',10, 'FontWeight','bold');
    ylabel('Current (pA/pF)', 'FontSize',10, 'FontWeight','bold');
    get(gcf,'CurrentAxes');
    set(gca,'YDir','normal');
    set(gca,'LineWidth',2, 'FontSize',10, 'FontWeight','bold');
    set(gca,'GridLineStyle','--')
end

figure(2)
plot(t, yksum(:, 1))
hold on
for i = 2:length(volts)
    plot(t, yksum(:, i))
end
hold off

%% custom functions
function [model_struct, psize] = gen_mdl_struct(current_names, tune_idx)
    num_currents = length(current_names);
    
    % current model strcuture index 1
    idx_info1 = cell(1, num_currents);
    for i = 1:num_currents
        switch current_names{i}
            case 'iktof'
                idx_info1{i} = tune_idx{1};
            case 'ikslow1'
                idx_info1{i} = tune_idx{2};
            case 'ikslow2'
                idx_info1{i} = tune_idx{3};
            case 'ikss'
                idx_info1{i} = tune_idx{4};
        end
    end

    % current model strcuture index 2
    idx_info2 = cell(1, num_currents);
    psize = 0;
    for i = 1: num_currents
        psize_old = psize;
        psize = psize + length(idx_info1{i});
        idx_info2{i} = (1+psize_old):(psize);
    end

    model_info = [current_names; idx_info1; idx_info2];
    field_names = {'name', 'idx1', 'idx2'};
    model_struct = cell2struct(model_info, field_names, 1);
end

function sol = gen_sol_vec(sol_mx, model_struct, psize)
    sol = zeros(1, psize);
    for i = 1:length(model_struct)
        running_sol = sol_mx(:, i);
        running_sol = running_sol(~isnan(running_sol));
        sol(model_struct(i).idx2) = running_sol(model_struct(i).idx1);
    end
end

function [z, yksum_set, ycomp_set] = rmsey(p, y, model_info, protocol_info)
    time_space = protocol_info{4};
    hold_idx = length(time_space{2});
    
    volts = protocol_info{3};
    num_volts = length(volts);

    rmse_list = NaN(num_volts,1);
    yksum_set = cell(num_volts,1);
    ycomp_set = cell(num_volts,1);
    for i = 1:num_volts
        y_i = y(:, i);
        protocol_info{3} = volts(i);

        [yhat, comp_hat] = kcurrent_model_temp(p, model_info, protocol_info);
        rmse_list(i) = sqrt(mean((y_i((hold_idx + 1):end) - yhat((hold_idx + 1):end)).^2));
        yksum_set{i} = yhat;
        ycomp_set{i} = comp_hat;
    end
    z = sum(rmse_list);
end

function [yksum, comp_currents] = kcurrent_model_temp(p, model_struct, protocol_info)
    global param_kslow1
    
    num_currents = length(model_struct);
    current_names = cell(num_currents, 1);
    
    for i = 1:num_currents
        current_names{i} = model_struct(i).name;
    end
    
    % declare shared parameters of ikslow1 as global variable
    matching_idx = strcmp(current_names, 'ikslow1');
    if any(matching_idx)    
        num_kslow1_param = 11;
        kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, ...
            6.1, 18.0, 2.058, 803.0, 0.16];
        
        tune_idx1 = model_struct(matching_idx).idx1;
        tune_idx2 = model_struct(matching_idx).idx2;
        fixed_idx = setdiff(1:num_kslow1_param, tune_idx1);

        param_kslow1 = zeros(num_kslow1_param, 1);
        param_kslow1(tune_idx1) = p(tune_idx2);
        param_kslow1(fixed_idx) = kslow1_default(fixed_idx);
    end

    time_space = protocol_info{4};
    yksum = zeros(length(time_space{1}), 1);
    comp_currents = cell(num_currents, 1);
    for i = 1:num_currents
        comp_currents{i} = gen_matching_current(p, model_struct(i), protocol_info);
        yksum = yksum + comp_currents{i};
    end
end

function current_trace = gen_matching_current(p, model_info, protocol_info)
    global param_kslow1

    % current model info
    current_name = model_info.name;
    tune_idx1 = model_info.idx1;
    tune_idx2 = model_info.idx2;

    % protocol info
    hold_volt = protocol_info{1};
    ek = protocol_info{2};
    volt = protocol_info{3};
    time_space = protocol_info{4};

    % generate current
    switch current_name
    case 'iktof'
        % generate ikto
        num_param = 13;

        kto_default = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
            0.000152, 0.067083, 0.00095, 0.051335, 0.4067];
        fixed_idx = setdiff(1:num_param, tune_idx1);
        
        param = zeros(num_param, 1);
        param(tune_idx1) = p(tune_idx2);
        param(fixed_idx) = kto_default(fixed_idx);
        
        current_trace = ikto(param, hold_volt, volt, time_space, ek);
    case 'ikslow1'
        % generate ikslow1
        current_trace = ikslow1(param_kslow1, hold_volt, volt, time_space, ek);
    case 'ikslow2'
        % generate ikslow2
        num_param = 11;
        kslow2_default = [4912, 5334, 0.16];
        shared_idx = [1:7, 9];
        uniq_idx = setdiff(1:num_param, shared_idx);

        param = zeros(num_param, 1);
        param(shared_idx) = param_kslow1(shared_idx);

        uniq_param = kslow2_default;
        uniq_param(tune_idx1) = p(tune_idx2);
        param(uniq_idx) = uniq_param;

        current_trace = ikslow2(param, hold_volt, volt, time_space, ek);
    case 'ikss'
        % generate ikss
        num_param = 7;
        kss_default = [0.0862, 1235.5, 13.17, 0.0611];
        shared_idx1 = 1:3;
        shared_idx2 = [1, 3, 4];
        uniq_idx = setdiff(1:num_param, shared_idx1);

        param = zeros(num_param, 1);
        param(shared_idx1) = param_kslow1(shared_idx2);

        uniq_param = kss_default;
        uniq_param(tune_idx1) = p(tune_idx2);
        param(uniq_idx) = uniq_param;

        current_trace = ikss(param, hold_volt, volt, time_space, ek);
    end
end

function current_trc = ikto(p, hold_volt, volt, time_space, ek)    
    % constants & initial values
    gmax = p(13);
    act0 = 0.4139033547E-02;
    inact0 = 0.9999623535E+00;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding 
    gv_hold = ikto_gating_variables(p, hold_volt);

    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));      

    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(hold_volt-ek);

    % current equation at pulse voltage
    gv_pulse = ikto_gating_variables(p, volt);

    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));

    current_trc((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-ek);
end

function gv = ikto_gating_variables(p, V)
    gv = zeros(4,1);
    
    % for Ikto
    alpha1 = p(7).*exp(p(5).*(V+p(1)));
    beta1 = p(8).*exp(-p(6).*(V+p(1)));
    
    alpha2_temp1 = p(9).*exp((V+p(2))./(-1.0.*p(4)));
    alpha2_temp2 = p(10).*exp((V+p(2)+p(3))./(-1.0.*p(4)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(11).*exp((V+p(2)+p(3))./p(4));
    beta2_temp2 = p(12).*exp((V+p(2)+p(3))./p(4));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    gv(1) = alpha1./(alpha1+beta1);
    gv(2) = alpha2./(alpha2+beta2);
    gv(3) = 1./(alpha1+beta1);
    gv(4) = 1./(alpha2+beta2);
end

function current_trc = ikslow1(p, hold_volt, volt, time_space, ek)
    % constants & initial values
    gmax = p(11);
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;
    
    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding 
    gv_hold = ikslow1_gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));        
    current_trc(1:hold_idx) = (gmax).*(act_hold).*(inact_hold).*(hold_volt-ek);

    % current equation at pulse voltage
    gv_pulse = ikslow1_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = (gmax).*(act_pulse).*(inact_pulse).*(volt-ek);
end

function gv = ikslow1_gating_variables(p, V)
    % gv(1:3) = gv(1:3) in Ikslow2
    % gv(1:3) = gv(1:3) in Ikur
    % gv(1) = gv(1) in Ikss

    gv = zeros(4, 1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    gv(3) = p(7)./(exp(p(6)*(V+p(3))) + exp(-p(6)*(V+p(3))))+p(9); % taua
    gv(4) = p(10) - p(8)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function current_trc = ikslow2(p, hold_volt, volt, time_space, Ek)
    % constants & initial values
    gmax = p(11); % 0.05766
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % current equation at holding
    gv_hold = ikslow2_gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(3));
    inact_hold = hh_model(hold_t, inact0, gv_hold(2), gv_hold(4));
    current_trc(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = ikslow2_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(3));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(2), gv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - Ek);
end

function gv = ikslow2_gating_variables(p, V)
    % gv(1:3) = gv(1:3) in Ikslow1
    % gv(4) = p(9) - p(1)*[gv(2) in Ikslow1]
    % {p(8): p1, p(10): p2}

    gv = zeros(4,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    gv(3) = p(7)./(exp(p(6)*(V+p(3))) + exp(-p(6)*(V+p(3))))+p(9); % taua
    gv(4) = p(10) - p(8)./(1.0+exp((p(2)+V)./p(5))); % taui    
end

function current_trc = ikss(p, hold_volt, volt, time_space, Ek)
    % 7 parameters; {p(7): gmax}

    % constants & initial values
    gmax = p(7); % 0.0428
    act0 = 0.5091689794e-03;

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);
    
    current_trc = zeros(length(t), 1);

    % current equation at holding
    gv_hold = ikss_gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(2));
    current_trc(1:hold_idx) = gmax.*act_hold.*(hold_volt - Ek);

    % current equation at pulse voltage
    gv_pulse = ikss_gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(2));
    current_trc((hold_idx + 1):end) = gmax.*act_pulse.*(volt - Ek);
end

function gv = ikss_gating_variables(p, V)
    % gv(1) = gv(1) in Ikslow1

    gv = zeros(2,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    gv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
