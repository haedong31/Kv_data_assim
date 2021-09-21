%% comparison with experimental data
clc
close all 
clear variables

% specify result files and voltages to check
file_group = 'ko';
cell_name = '15o27002';
file_name = strcat(cell_name, '.xlsx');

save_dir1 = 'calib_exp16';
save_dir2 = 'calib_exp17';
save_dir3 = 'calib_exp18';
save_dir4 = 'calib_exp19';

% voltages info
min_volt = -50;
volt_range = 3:11;
volts = NaN(length(volt_range), 1);
for i = 1:length(volt_range)
    volts(i) = min_volt + (volt_range(i)-1)*10;
end

% import calibration results
calib1 = table2array(readtable(fullfile(pwd, strcat(save_dir1, '_', file_group), file_name)));
calib2 = table2array(readtable(fullfile(pwd, strcat(save_dir2, '_', file_group), file_name)));
calib3 = table2array(readtable(fullfile(pwd, strcat(save_dir3, '_', file_group), file_name)));
calib4 = table2array(readtable(fullfile(pwd, strcat(save_dir4, '_', file_group), file_name)));

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

[model_struct1, psize1] = gen_mdl_struct(current_names, tune_idx1);

sol1 = gen_sol_vec(calib1, model_struct1, psize1);
sol2 = gen_sol_vec(calib2, model_struct1, psize1);
sol3 = gen_sol_vec(calib3, model_struct1, psize1);
sol4 = gen_sol_vec(calib4, model_struct1, psize1);

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

protocol{1} = hold_volt;
protocol{3} = time_space;
protocol{4} = ek;

% generate yksum_hat
figure
for i=1:length(volts)
    protocol{2} = volts(i);
    yksum_i = yksum(:, i);

    [yksum_hat1, ~] = kcurrent_model1(sol1, model_struct1, protocol);
    [yksum_hat2, ~] = kcurrent_model1(sol2, model_struct1, protocol);
    [yksum_hat3, ~] = kcurrent_model1(sol3, model_struct1, protocol);
    [yksum_hat4, ~] = kcurrent_model1(sol4, model_struct1, protocol); 
    
    subplot(3, 3, i)
    plot(t, yksum_i, 'Color','red')
    hold on
    plot(t, yksum_hat1, '--', 'Color','green', 'LineWidth',2)
    plot(t, yksum_hat2, '--', 'Color','blue', 'LineWidth',2)
    plot(t, yksum_hat4, '--', 'Color','cyan', 'LineWidth',2)
    plot(t, yksum_hat3, '--', 'Color','black', 'LineWidth',2)
    hold off
    axis tight
    grid on
    legend('Experimental Data','Interior Point', 'SQP', 'Active Set', 'GA', 'Location','best')
    title(strcat(cell_name, '(', file_group, ')', ' Clamp Voltage: ', string(volts(i)), ' mV'), 'FontSize',10, 'FontWeight','bold');
    xlabel('Time (ms)', 'FontSize',10, 'FontWeight','bold');
    ylabel('Current (pA/pF)', 'FontSize',10, 'FontWeight','bold');
    get(gcf,'CurrentAxes');
    set(gca,'YDir','normal');
    set(gca,'LineWidth',2, 'FontSize',10, 'FontWeight','bold');
    set(gca,'GridLineStyle','--')
end

function [model_struct, psize] = gen_mdl_struct(current_names, tune_idx)
    num_currents = length(current_names);
    
    % current model strcuture index 1
    idx_info1 = cell(1, num_currents);
    for i = 1:num_currents
        switch current_names{i}
        case 'ikto'
            idx_info1{i} = tune_idx{1};
        case 'ikslow1'
            idx_info1{i} = tune_idx{2};
        case 'ikslow2'
            idx_info1{i} = tune_idx{3};
        case 'ikss'
            idx_info1{i} = tune_idx{4};
        case 'ikur'
            idx_info1{i} = tune_idx{5};
        case 'ik1'
            idx_info1{i} = tune_idx{6};
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
