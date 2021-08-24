%% comparison with experimental data
clc
close all 
clear variables

% specify result file to check
file_group = 'wt';
file_name = '15o26002.xlsx';
calib = table2array(readtable(fullfile(pwd, file_group, strcat('calib_param_', file_name))));

% specify model structure
% field 1: selection of currents
current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% field 2: tunning index in individual current models
tune_idx1_kto = [1, 2, 3, 5, 6, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 9, 12, 13];
tune_idx1_kslow2 = [9, 11];
tune_idx1_kss = [6, 7];
% tune_idx1_k1 = [1, 3, 5, 7];
idx_info1 = {tune_idx1_kto, ...
    tune_idx1_kslow1, ...
    tune_idx1_kslow2, ...
    tune_idx1_kss};

% field 3: tunning index in decision variable, p
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

sol = zeros(cul_idx_len, 1);
for i=1:num_currents
    running_calib = calib(:, i);
    running_calib(~isnan(running_calib));
    sol(idx_info2{i}) = running_calib(idx_info1{i});
    % for older code
%     if i==3
%         sol(idx_info2{i}) = running_calib((idx_info1{i}-3));
%     else
%         sol(idx_info2{i}) = running_calib(idx_info1{i});
%     end
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

%% box plots of calibration parameters
clc
close all
clear variables

load('sols_ko.mat')
load('sols_wt.mat')

% remove rows with all zeros
sols_wt = sols_wt(any(sols_wt, 2), :);
sols_ko = sols_ko(any(sols_ko, 2), :);

v = 3;
g = 2;
n = 33;

boxplot_mat = zeros(v, g, n);

for i = 1:v
    for j = 1:g
        for k = 1:n
            if j == 1
                boxplot_mat(i, j, k) = sols_wt(k, i);
            else
                boxplot_mat(i, j, k) = sols_ko(k, i);
            end
        end
    end
end

bp = boxplot2(boxplot_mat, 1:v);
cmap = get(0, 'defaultaxescolororder');
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), bp);
end
legend('WT','Mgat1KO', 'Location','best')
