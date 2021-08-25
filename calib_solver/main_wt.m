clc
close all
clear variables

warning('off', 'all')

% code arguments for calibration
group_name = 'wt';

% selection of currents
current_names = {'ikto', 'ikslow1', 'ikss'};
num_currents = length(current_names);

% tunning index in individual current models
tune_idx1_kto = [1, 2, 3, 5, 6, 15, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 9, 12, 13];
tune_idx1_kslow2 = [1, 3];
tune_idx1_kss = [3, 4];
tune_idx1_kur = [1, 3];
tune_idx1_k1 = [1, 3, 5, 7];

% optimization options
max_evals = 1e+6;
num_iters = 1;
options = optimoptions('fmincon', 'MaxFunctionEvaluations',max_evals, 'Display','off');

% protocol
hold_volt = -70;
volts = -50:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;
ek = -91.1;

idx_info1 = cell(1, num_currents);
for i = 1:num_currents
    switch current_names{i}
    case "ikto"
        idx_info1{i} = tune_idx1_kto;
    case "ikslow1"
        idx_info1{i} = tune_idx1_kslow1;
    case "ikslow2"
        idx_info1{i} = tune_idx1_kslow2;
    case "ikss"
        idx_info1{i} = tune_idx1_kss;
    case "ikur"
        idx_info1{i} = tune_idx1_kur;
    case "ik1"
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
mkdir(fullfile(pwd, strcat('calib_result_', group_name)));

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

lb_kto_default = [-50, -50, -50, -10, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps];
lb_kslow1_default = [-70, -70, -70, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps];
lb_kslow2_default = [5000, eps, eps];
lb_kss_default = [eps, eps, eps, eps];
lb_kur_default = [eps, 500, eps];
lb_k1_default = [eps, -20, eps, -30, eps, eps, eps, eps, eps, eps];

ub_kto_default = [70, 50, 50, 40, 50, 30, 1, 1, 1, 10, 0.005, 0.3, 0.005, 0.5, 1, 1, 1];
ub_kslow1_default = [50, 50, 50, 10, 50, 50, 1, 100, 1000, 50, 1, 0.5, 0.5];
ub_kslow2_default = [10000, 5000, 0.5];
ub_kss_default = [1, 2000, 100, 0.5];
ub_kur_default = [500, 2000, 1];
ub_k1_default = [120, 30, 1000, 30, 3, 1, 2, 0.25, 0.25, 1];

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
    case "ikto"
        p0(idx_info2{i}) = kto_default(tune_idx1_kto);
        lb(idx_info2{i}) = lb_kto_default(tune_idx1_kto);
        ub(idx_info2{i}) = ub_kto_default(tune_idx1_kto);
        if max_param_len < length(kto_default)
            max_param_len = length(kto_default);
        end
    case "ikslow1"
        p0(idx_info2{i}) = kslow1_default(tune_idx1_kslow1);
        lb(idx_info2{i}) = lb_kslow1_default(tune_idx1_kslow1);
        ub(idx_info2{i}) = ub_kslow1_default(tune_idx1_kslow1);
        if max_param_len < length(kslow1_default)
            max_param_len = length(kslow1_default);
        end
    case "ikslow2"
        p0(idx_info2{i}) = kslow2_default(tune_idx1_kslow2);
        lb(idx_info2{i}) = lb_kslow2_default(tune_idx1_kslow2);
        ub(idx_info2{i}) = ub_kslow2_default(tune_idx1_kslow2);
        if max_param_len < length(kslow2_default)
            max_param_len = length(kslow2_default);
        end    
    case "ikss"
        p0(idx_info2{i}) = kss_default(tune_idx1_kss);
        lb(idx_info2{i}) = lb_kss_default(tune_idx1_kss);
        ub(idx_info2{i}) = ub_kss_default(tune_idx1_kss);
        if max_param_len < length(kss_default)
            max_param_len = length(kss_default);
        end    
    case "ikur"
        p0(idx_info2{i}) = kur_default(tune_idx1_kur);
        lb(idx_info2{i}) = lb_kur_default(tune_idx1_kur);
        ub(idx_info2{i}) = ub_kur_default(tune_idx1_kur);
        if max_param_len < length(kur_default)
            max_param_len = length(kur_default);
        end
    case "ik1"
        p0(idx_info2{i}) = k1_default(tune_idx1_kur);
        lb(idx_info2{i}) = lb_k1_default(tune_idx1_k1);
        ub(idx_info2{i}) = ub_k1_default(tune_idx1_k1);
        if max_param_len < length(k1_default)
            max_param_len = length(k1_default);
        end
    end    
end

% main loop
num_files = length(loop_idx);
for l = 1:floor(len_loop_idx/2)
    i = loop_idx(l);
    fprintf('[%i/%i] %s \n', i, num_files, file_names{i})

    % read data
    file_path = fullfile(pwd, 'data', strcat(group_name, '-preprocessed2'), file_names{i});
    trace_data = table2array(readtable(file_path));

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

    % objective function
%     disp(obj_rmse(p0, model_struct, volt_space, time_space, yksum))
    opt_fun = @(p) obj_rmse(p, model_struct, volt_space, time_space, yksum);

    % run optimization
    rmse_list = zeros(num_iters, 1);
    sol_list = cell(num_iters, 1);

    % first run with p0
    [sol, fval] = fmincon(opt_fun, p0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    sol_list{1} = sol;
    rmse_list(1) = fval;

    for j = 2:num_iters
        fprintf('[%i/%i] \n', j, num_iters)

        % random intialization
        running_p0 = zeros(length(p0), 1);
        for k = 1:length(p0)
            unif_dist = makedist('Uniform', 'lower',lb(k), 'upper',ub(k));
            running_p0(k) = random(unif_dist, 1);
        end
        
        % optimization
        try
            [sol, fval] = fmincon(opt_fun, running_p0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            sol_list{j} = sol;
            rmse_list(j) = fval;
        catch
            rmse_list(j) = 1e+3;
        end
    end

    [~, best_fit_idx] = min(rmse_list);
    best_sol = sol_list{best_fit_idx};

    % save calibrated solution
    save_path = fullfile(pwd, strcat('calib_result_', group_name), file_names{i});
    
    sol_mx = zeros(max_param_len, num_currents);
    for j = 1:num_currents
        switch current_names{j}
        case 'ikto'
            sol_kto = kto_default;
            sol_kto(tune_idx1_kto) = best_sol(idx_info2{j});
            sol_kto = [sol_kto, NaN(1, (max_param_len-length(sol_kto)))];
            sol_mx(:, j) = sol_kto;
        case 'ikslow1'
            sol_kslow1 = kslow1_default;
            sol_kslow1(tune_idx1_kslow1) = best_sol(idx_info2{j});
            sol_kslow1 = [sol_kslow1, NaN(1, (max_param_len-length(sol_kslow1)))];
            sol_mx(:, j) = sol_kslow1;
        case 'ikslow2'
            sol_kslow2 = kslow2_default;
            sol_kslow2(tune_idx1_kslow2) = best_sol(idx_info2{j});
            sol_kslow2 = [sol_kslow2, NaN(1, (max_param_len-length(sol_kslow2)))];
            sol_mx(:, j) = sol_kslow2;
        case 'ikss'
            sol_kss = kss_default;
            sol_kss(tune_idx1_kss) = best_sol(idx_info2{j});
            sol_kss = [sol_kss, NaN(1, (max_param_len-length(sol_kss)))];
            sol_mx(:, j) = sol_kss;
        case 'ikur'
            sol_kur = kur_default;
            sol_kur(tune_idx1_kur) = best_sol(idx_info2{j});
            sol_kur = [sol_kur, NaN(1, (max_param_len-length(sol_kur)))];
            sol_mx(:, j) = sol_kur;
        case 'ik1'
            sol_k1 = k1_default;
            sol_k1(tune_idx1_k1) = best_sol(idx_info2{j});
            sol_k1 = [sol_k1, NaN(1, (max_param_len-length(sol_k1)))];
            sol_mx(:, j) = sol_k1;
        end
    end
    writematrix(string(current_names) , save_path, "Sheet","Parameters", "Range","A1");
    writematrix(sol_mx, save_path, "Sheet","Parameters", "Range","A2");
end
