clc
close all
clear variables

warning('off', 'all')

%
% input arguments for calibration
% 
group_name = 'ko';
save_dir = strcat('calib_exp13_', group_name);

% selection of currents
current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% tunning index in individual current models

tune_idx_phase1 = cell(6, 1);
tune_idx_phase1{1} = [2, 4, 7, 10];
tune_idx_phase1{2} = [2, 5, 9];
tune_idx_phase1{3} = 1;
tune_idx_phase1{4} = 3;
tune_idx_phase1{5} = [1, 3];
tune_idx_phase1{6} = [1, 3, 5, 7];

tune_idx_phase2 = cell(6, 1);
tune_idx_phase2{1} = 16;
tune_idx_phase2{2} = 11;
tune_idx_phase2{3} = 3;
tune_idx_phase2{4} = 4;
tune_idx_phase2{5} = [1, 3];
tune_idx_phase2{6} = [1, 3, 5, 7];

tune_idx_phase3 = cell(6, 1);
tune_idx_phase3{1} = [1, 2, 4, 6, 7, 10, 16];
tune_idx_phase3{2} = [1, 2, 3, 4, 5, 8, 9, 11];
tune_idx_phase3{3} = [1, 3];
tune_idx_phase3{4} = [3, 4];
tune_idx_phase3{5} = [1, 3];
tune_idx_phase3{6} = [1, 3, 5, 7];

% model structure
field_names = {'name', 'idx1', 'idx2'};
[idx_info1, idx_info21, psize1] = gen_idx_struct(current_names, tune_idx_phase1);
[idx_info12, idx_info22, psize2] = gen_idx_struct(current_names, tune_idx_phase2);
[idx_info13, idx_info23, psize3] = gen_idx_struct(current_names, tune_idx_phase3);

model_info1 = [current_names; idx_info1; idx_info21];
model_info2 = [current_names; idx_info12; idx_info22];
model_info3 = [current_names; idx_info13; idx_info23];

model_struct1 = cell2struct(model_info1, field_names, 1);
model_struct2 = cell2struct(model_info2, field_names, 1);
model_struct3 = cell2struct(model_info3, field_names, 1);

[p01, lb1, ub1, max_param_len] = gen_optim_bound(model_struct1);
[p02, lb2, ub2, ~] = gen_optim_bound(model_struct2);
[p03, lb3, ub3, ~] = gen_optim_bound(model_struct3);

% voltages info
min_volt = -50;
volt_range = 3:11;
volts = NaN(length(volt_range), 1);
for i = 1:length(volt_range)
    volts(i) = min_volt + (volt_range(i)-1)*10;
end

%
% optimization set up
%
% protocol
hold_volt = -70;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
ek = -91.1;

% voltage info
volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = ek;

% optimization options
max_evals = 1e+6;
options = optimoptions(@fmincon, 'MaxFunctionEvaluations',max_evals, ...
    'Display','off');

% optimization constraints 
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

% matching table
matching_table = readtable(fullfile(pwd, 'data', strcat('matching-table-', group_name, '.xlsx')));
[num_files, ~] = size(matching_table);
mkdir(fullfile(pwd, save_dir));

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

%
% main loop
%
for l = 1:len_loop_idx
    i = loop_idx(l);

    % read data
    file_path = fullfile(pwd, 'data', strcat(group_name, '-preprocessed2'), file_names{i});
    trace_data = table2array(readtable(file_path));

    t = trace_data(:, 1);
    yksum = trace_data(:, (volt_range+1));

    % estimate time points
    [~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
    [~, ideal_end_idx] = min(abs(t - ideal_end_time));

    t = t(1:ideal_end_idx);
    yksum = yksum(1:ideal_end_idx, :);

    hold_t = t(1:ideal_hold_idx);
    pulse_t = t(ideal_hold_idx+1:end);
    pulse_t_adj = pulse_t - pulse_t(1);
    
    % time space
    time_space = cell(5, 1);
    time_space{1} = t;
    time_space{2} = hold_t;
    time_space{3} = pulse_t_adj;
    time_space{4} = ideal_hold_idx;
    time_space{5} = ideal_end_idx;

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{2} = 30;
    protocol{3} = time_space;
    protocol{4} = volt_space{3};

    %
    % phase 1
    % 
%     obj_rmse(p01, 'tail', @kcurrent_model2, model_struct1, volt_space, time_space, yksum);
    opt_fun1 = @(p) obj_rmse(p, 'tail', @kcurrent_model2, model_struct1, volt_space, time_space, yksum);
    [sol1, fval] = fmincon(opt_fun1, p01, A, b, Aeq, beq, lb1, ub1, nonlcon, options);
    fprintf('[File %i/%i] %s Phase 1: %f \n', l, len_loop_idx, file_names{i}, fval)
%     obj_rmse(sol1, 'tail', @kcurrent_model2, model_struct1, volt_space, time_space, yksum)

    pretrn_sol1 = cell(1, num_currents);
    for j = 1:num_currents
        pretrn_sol1{j} = sol1(model_struct1(j).idx2);
    end
    pretrn_mdl1 = cell2struct([model_info1; pretrn_sol1], {'name', 'idx1', 'idx2', 'sol'});

    %
    % phase 2
    %
%     obj_pretrn_rmse(p02, 'early', @kcurrent_model3, model_struct2, volt_space, time_space, yksum, pretrn_mdl1);
    opt_fun2 = @(p) obj_pretrn_rmse(p, 'early', @kcurrent_model3, model_struct2, volt_space, time_space, yksum, pretrn_mdl1);
    [sol2, fval] = fmincon(opt_fun2, p02, A, b, Aeq, beq, lb2, ub2, nonlcon, options);
    fprintf('[File %i/%i] %s Phase 2: %f \n', l, len_loop_idx, file_names{i}, fval)
%     obj_rmse(sol2, 'early', @kcurrent_model2, model_struct1, volt_space, time_space, yksum)
    
    for j = 1:num_currents
        pretrn1_idx1 = model_struct1(j).idx1;
        pretrn1_idx2 = model_struct1(j).idx2;

        pretrn2_idx1 = model_struct2(j).idx1;
        pretrn2_idx2 = model_struct2(j).idx2;
        switch model_struct3(j).name
            case 'ikto'
                kto_default = [33.0, 40.0, 33.0, 15.5, 20.0, 7.7, 5.7, ...
                    0.03577, 0.06237, 7.0, 0.18064, 0.3956, 0.0226585, 0.00095, 0.051335, 0.14067];
                kto_default(pretrn1_idx1) = sol1(pretrn1_idx2);
                kto_default(pretrn2_idx1) = sol2(pretrn2_idx2);
                p03(model_struct3(j).idx2) = kto_default(model_struct3(j).idx1);
            case 'ikslow1'
                kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.05766];
                kslow1_default(pretrn1_idx1) = sol1(pretrn1_idx2);
                kslow1_default(pretrn2_idx1) = sol2(pretrn2_idx2);
                p03(model_struct3(j).idx2) = kslow1_default(model_struct3(j).idx1);
            case 'ikslow2'
                kslow2_default = [5334, 4912, 0.05766];
                kslow2_default(pretrn1_idx1) = sol1(pretrn1_idx2);
                kslow2_default(pretrn2_idx1) = sol2(pretrn2_idx2);
                p03(model_struct3(j).idx2) = kslow2_default(model_struct3(j).idx1);
            case 'ikss'
                kss_default = [0.0862, 1235.5, 13.17, 0.0428];
                kss_default(pretrn1_idx1) = sol1(pretrn1_idx2);
                kss_default(pretrn2_idx1) = sol2(pretrn2_idx2);
                p03(model_struct3(j).idx2) = kss_default(model_struct3(j).idx1);
            case 'ikur'
                kur_default = [270, 1050, 0];
                kur_default(pretrn1_idx1) = sol1(pretrn1_idx2);
                kur_default(pretrn2_idx1) = sol2(pretrn2_idx2);
                p03(model_struct3(j).idx2) = kur_default(model_struct3(j).idx1);
            case 'ik1'
                k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];
                k1_default(pretrn1_idx1) = sol1(pretrn1_idx2);
                k1_default(pretrn2_idx1) = sol2(pretrn2_idx2);
                p03(model_struct3(j).idx2) = k1_default(model_struct3(j).idx1);
        end
    end
    %
    % phase 3
    %
    opt_fun3 = @(p) obj_rmse(p, 'all', @kcurrent_model2, model_struct3, volt_space, time_space, yksum);
    [sol3, fval] = fmincon(opt_fun3, p03, A, b, Aeq, beq, lb3, ub3, nonlcon, options);
    fprintf('[File %i/%i] %s Phase 3: %f \n', l, len_loop_idx, file_names{i}, fval)
%     obj_rmse(sol3, 'all', @kcurrent_model2, model_struct3, volt_space, time_space, yksum)

    save_path = fullfile(pwd, save_dir, file_names{i});
    sol_mx = zeros(max_param_len, num_currents);
    for j = 1:num_currents
        switch current_names{j}
        case 'ikto'
            sol_kto = kto_default;
            sol_kto(model_struct3(j).idx1) = sol3(model_struct3(j).idx2);
            sol_kto = [sol_kto, NaN(1, (max_param_len-length(sol_kto)))];
            sol_mx(:, j) = sol_kto;
        case 'ikslow1'
            sol_kslow1 = kslow1_default;
            sol_kslow1(model_struct3(j).idx1) = sol3(model_struct3(j).idx2);
            sol_kslow1 = [sol_kslow1, NaN(1, (max_param_len-length(sol_kslow1)))];
            sol_mx(:, j) = sol_kslow1;
        case 'ikslow2'
            sol_kslow2 = kslow2_default;
            sol_kslow2(model_struct3(j).idx1) = sol3(model_struct3(j).idx2);
            sol_kslow2 = [sol_kslow2, NaN(1, (max_param_len-length(sol_kslow2)))];
            sol_mx(:, j) = sol_kslow2;
        case 'ikss'
            sol_kss = kss_default;
            sol_kss(model_struct3(j).idx1) = sol3(model_struct3(j).idx2);
            sol_kss = [sol_kss, NaN(1, (max_param_len-length(sol_kss)))];
            sol_mx(:, j) = sol_kss;
        case 'ikur'
            sol_kur = kur_default;
            sol_kur(model_struct3(j).idx1) = sol3(model_struct3(j).idx2);
            sol_kur = [sol_kur, NaN(1, (max_param_len-length(sol_kur)))];
            sol_mx(:, j) = sol_kur;
        case 'ik1'
            sol_k1 = k1_default;
            sol_k1(model_struct3(j).idx1) = sol3(model_struct3(j).idx2);
            sol_k1 = [sol_k1, NaN(1, (max_param_len-length(sol_k1)))];
            sol_mx(:, j) = sol_k1;
        end
    end
    writematrix(string(current_names) , save_path, "Sheet","Parameters", "Range","A1");
    writematrix(sol_mx, save_path, "Sheet","Parameters", "Range","A2");
end

function [idx_info1, idx_info2, psize] = gen_idx_struct(current_names, tune_idx)
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
end

function [init_pt, lb, ub, max_param_len] = gen_optim_bound(model_struct)
    kto_default = [33.0, 40.0, 33.0, 15.5, 20.0, 7.7, 5.7, 0.03577, 0.06237, 7.0, 0.18064, 0.3956, 0.0226585, 0.00095, 0.051335, 0.14067];
    kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.05766];
    kslow2_default = [5334, 4912, 0.05766];
    kss_default = [0.0862, 1235.5, 13.17, 0.0428];
    kur_default = [270, 1050, 0];
    k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

    lb_kto = [-70, -70, -70, -50, -20, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps];
    lb_kslow1 = [-70, -70, -70, eps, eps, eps, eps, eps, eps, eps, eps];
    lb_kslow2 = [5000, eps, eps];
    lb_kss = [eps, eps, eps, eps];
    lb_kur = [eps, 500, eps];
    lb_k1 = [eps, -20, eps, -30, eps, eps, eps, eps, eps, eps];

    ub_kto = [70, 70, 70, 50, 20, 20, 20, 1, 1, 20, 1, 1, 1, 1, 1, 1];
    ub_kslow1 = [70, 70, 70, 10, 70, 70, 1, 100, 1000, 70, 0.5];
    ub_kslow2 = [10000, 5000, 0.5];
    ub_kss = [1, 2000, 100, 0.5];
    ub_kur = [500, 2000, 1];
    ub_k1 = [120, 30, 1000, 30, 3, 1, 2, 0.25, 0.25, 1];

    num_currents = length(model_struct);
    cul_idx_len = model_struct(end).idx2(end);

    max_param_len = 0;
    init_pt = zeros(1, cul_idx_len);
    lb = zeros(1, cul_idx_len);
    ub = zeros(1, cul_idx_len);
    for i = 1:num_currents
        model_info = model_struct(i);
        current_name = model_info.name;
        idx_info1 = model_info.idx1;
        idx_info2 = model_info.idx2;        
        
        switch current_name
            case 'ikto'
                init_pt(idx_info2) = kto_default(idx_info1);
                lb(idx_info2) = lb_kto(idx_info1);
                ub(idx_info2) = ub_kto(idx_info1);
                if max_param_len < length(kto_default)
                    max_param_len = length(kto_default);
                end
            case 'ikslow1'
                init_pt(idx_info2) = kslow1_default(idx_info1);
                lb(idx_info2) = lb_kslow1(idx_info1);
                ub(idx_info2) = ub_kslow1(idx_info1);
                if max_param_len < length(kslow1_default)
                    max_param_len = length(kslow1_default);
                end
            case 'ikslow2'
                init_pt(idx_info2) = kslow2_default(idx_info1);
                lb(idx_info2) = lb_kslow2(idx_info1);
                ub(idx_info2) = ub_kslow2(idx_info1);
                if max_param_len < length(kslow2_default)
                    max_param_len = length(kslow2_default);
                end    
            case 'ikss'
                init_pt(idx_info2) = kss_default(idx_info1);
                lb(idx_info2) = lb_kss(idx_info1);
                ub(idx_info2) = ub_kss(idx_info1);
                if max_param_len < length(kss_default)
                    max_param_len = length(kss_default);
                end    
            case 'ikur'
                init_pt(idx_info2) = kur_default(idx_info1);
                lb(idx_info2) = lb_kur(idx_info1);
                ub(idx_info2) = ub_kur(idx_info1);
                if max_param_len < length(kur_default)
                    max_param_len = length(kur_default);
                end
            case 'ik1'
                init_pt(idx_info2) = k1_default(idx_info1);
                lb(idx_info2) = lb_k1(idx_info1);
                ub(idx_info2) = ub_k1(idx_info1);
                if max_param_len < length(k1_default)
                    max_param_len = length(k1_default);
                end
        end    
    end
end
