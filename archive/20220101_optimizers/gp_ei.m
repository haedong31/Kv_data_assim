clc
close all
clear variables
format shortG

%%%%%
% code arguments for calibration
%%%%%
group_name = 'wt';
save_dir = strcat('calib_result1_', group_name);

% selection of currents
current_names = {'ikto', 'ikslow1', 'ikss'};

% tunning index in individual current models
tune_idx1_kto = [1, 2, 3, 5, 6, 15, 16, 17];
tune_idx1_kslow1 = [1, 2, 3, 9, 12, 13];
tune_idx1_kslow2 = [1, 3];
tune_idx1_kss = [3, 4];
tune_idx1_kur = [1, 3];
tune_idx1_k1 = [1, 3, 5, 7];

% protocol
hold_volt = -70;
volts = -50:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1e+3;
ek = -91.1;

% optimization options
reps = 30;
ninit = length(current_names)*10;
out_size = 100;

%%%%%
% main
%%%%%
% matching table
matching_table = readtable(fullfile(pwd, 'data', strcat('matching-table-', group_name, '.xlsx')));
[num_files, ~] = size(matching_table);
mkdir(fullfile(pwd, save_dir));
file_names = matching_table.trace_file_name_4half;

% exclude row not having 4.5-sec data
loop_idx = [];
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end
num_files = length(loop_idx);

% current model
global model_struct
model_struct = gen_model_struct(current_names, tune_idx1_kto, tune_idx1_kslow1, ...
    tune_idx1_kslow2, tune_idx1_kss, tune_idx1_kur, tune_idx1_k1);
num_var = model_struct(end).idx2(end);

% voltage
global volt_space
volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = ek;

global time_space
global yksum

% default values for initial starting point
global kto_default
global kslow1_default
global kslow2_default
global kss_default
global kur_default
global k1_default

kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow2_default = [5334, 4912, 0.05766];
kss_default = [0.0862, 1235.5, 13.17, 0.0428];
kur_default = [270, 1050, 0];
k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

p0 = NaN(1, num_var);
max_param_len = 0;
for i = 1:length(current_names)
    switch current_names{i}
    case 'ikto'
        p0(model_struct(i).idx2) = kto_default(tune_idx1_kto);
        if max_param_len < length(kto_default)
            max_param_len = length(kto_default);
        end
    case 'ikslow1'
        p0(model_struct(i).idx2) = kslow1_default(tune_idx1_kslow1);
        if max_param_len < length(kslow1_default)
            max_param_len = length(kslow1_default);
        end
    case 'ikslow2'
        p0(model_struct(i).idx2) = kslow2_default(tune_idx1_kslow2);
        if max_param_len < length(kslow2_default)
            max_param_len = length(kslow2_default);
        end    
    case 'ikss'
        p0(model_struct(i).idx2) = kss_default(tune_idx1_kss);
        if max_param_len < length(kss_default)
            max_param_len = length(kss_default);
        end    
    case 'ikur'
        p0(model_struct(i).idx2) = kur_default(tune_idx1_kur);
        if max_param_len < length(kur_default)
            max_param_len = length(kur_default);
        end
    case 'ik1'
        p0(model_struct(i).idx2) = k1_default(tune_idx1_kur);
        if max_param_len < length(k1_default)
            max_param_len = length(k1_default);
        end
    end        
end

for i = 1:floor(num_files/2)
    disp('##### RUN')

    % read data
    file_path = fullfile(pwd, 'data', strcat(group_name, '-preprocessed2'), file_names{loop_idx(i)});
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
    
    % initial point
    initx = lhsdesign(ninit, num_var);
    initx = scale_param(initx);
    initx(1, :) = p0;
    inity = obj_rmse(initx);

    % MC-style repetitions
    prog = NaN(reps, out_size);
    for r = 1:reps
        tic
        [sol, fval, maxei] = optim_ei(@obj_rmse, initx, inity, out_size);
        toc
        
        running_prog = bov(fval);
        prog(r, :) = running_prog';

        fprintf('[File %i/%i] %s [Reps %i/%i] Min RMSE: %f \n', i, num_files, file_names{loop_idx(i)}, r, reps, min(running_prog));
        [~, min_idx] = min(running_prog);
        bsol = sol(min_idx, :);
        bsol = scale_param(bsol);

        save_path = fullfile(pwd, save_dir, file_names{i});
        save_sol(bsol, save_path, max_param_len);
    end
end

%%%%%
% problem specific functions
%%%%%
function rmse_list = obj_rmse(p)
    % RMSE objective function; main optimization objective
    global model_struct
    global volt_space
    global time_space
    global yksum

    [num_pt, ~] = size(p);

    hold_idx = length(time_space{2});
    volts = volt_space{2};
    num_volts = length(volts);

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{3} = time_space;
    protocol{4} = volt_space{3};

    rmse_list = NaN(num_pt, 1);
    for i = 1:num_pt
        running_rmse = NaN(num_volts, 1);
        for j = 1:num_volts
            yksum_i = yksum(:, j);
            protocol{2} = volts(j);
    
            [yksum_hat, ~] = kcurrent_model(p(i, :), model_struct, protocol);
            
            % check validity
            [~, peak_idx] = max(yksum_hat);
            check_pt1 = any(isnan(yksum_hat));
            check_pt2 = any(yksum_hat < 0);
            check_pt3 = var(yksum_hat(1:hold_idx)) > 0.4812e-4; % not stable at hold_volt
            check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse
    
            if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
                running_rmse(j) = 1e+3; % arbitrary big number
            else
                running_rmse(j) = sqrt(mean((yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2));
            end
        end
        rmse_list(i) = sum(running_rmse);
    end
end

function scaledp = scale_param(unitp)
    global model_struct

    lb_kto = [-50, -50, -50, -10, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps];
    lb_kslow1 = [-70, -70, -70, eps, eps, eps, eps, eps, eps, eps, eps, eps, eps];
    lb_kslow2 = [5000, eps, eps];
    lb_kss = [eps, eps, eps, eps];
    lb_kur = [eps, 500, eps];
    lb_k1 = [eps, -20, eps, -30, eps, eps, eps, eps, eps, eps];
    
    ub_kto = [70, 50, 50, 40, 50, 30, 1, 1, 1, 10, 0.005, 0.3, 0.005, 0.5, 1, 1, 1];
    ub_kslow1 = [50, 50, 50, 10, 50, 50, 1, 100, 1000, 50, 1, 0.5, 0.5];
    ub_kslow2 = [10000, 5000, 0.5];
    ub_kss = [1, 2000, 100, 0.5];
    ub_kur = [500, 2000, 1];
    ub_k1 = [120, 30, 1000, 30, 3, 1, 2, 0.25, 0.25, 1];

    num_currents = length(model_struct);
    [num_pt, num_var] = size(unitp);
    scaledp = NaN(num_pt, num_var);

    for i = 1:num_pt
        row = NaN(1, num_var);
        for j = 1:num_currents
            switch model_struct(j).name
            case 'ikto'
                lb = lb_kto(model_struct(j).idx1);
                ub = ub_kto(model_struct(j).idx1);
                row(model_struct(j).idx2) = unitp(i, model_struct(j).idx2).*(ub-lb) + lb;
            case 'ikslow1'
                lb = lb_kslow1(model_struct(j).idx1);
                ub = ub_kslow1(model_struct(j).idx1);
                row(model_struct(j).idx2) = unitp(i, model_struct(j).idx2).*(ub-lb) + lb;
            case 'ikslow2'
                lb = lb_kslow2(model_struct(j).idx1);
                ub = ub_kslow2(model_struct(j).idx1);
                row(model_struct(j).idx2) = unitp(i, model_struct(j).idx2).*(ub-lb) + lb;
            case 'ikss'
                lb = lb_kss(model_struct(j).idx1);
                ub = ub_kss(model_struct(j).idx1);
                row(model_struct(j).idx2) = unitp(i, model_struct(j).idx2).*(ub-lb) + lb;
            case 'ikur'
                lb = lb_kur(model_struct(j).idx1);
                ub = ub_kur(model_struct(j).idx1);
                row(model_struct(j).idx2) = unitp(i, model_struct(j).idx2).*(ub-lb) + lb;
            case 'ik1'
                lb = lb_k1(model_struct(j).idx1);
                ub = ub_k1(model_struct(j).idx1);
                row(model_struct(j).idx2) = unitp(i, model_struct(j).idx2).*(ub-lb) + lb;
            end
        end
        scaledp(i, :) = row;
    end
end

function model_struct = gen_model_struct(current_names, tune_idx1_kto, tune_idx1_kslow1, ...
    tune_idx1_kslow2, tune_idx1_kss, tune_idx1_kur, tune_idx1_k1)
    num_currents = length(current_names);
    
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
end

function save_sol(sol, save_path, max_param_len)
    global model_struct
    global kto_default
    global kslow1_default
    global kslow2_default
    global kss_default
    global kur_default
    global k1_default

    num_currents = length(model_struct);
    current_names = cell(1, num_currents);
    sol_mx = NaN(max_param_len, num_currents);

    for i = 1:num_currents
        current_name = model_struct(i).name;
        current_names{i} = current_name;

        switch current_name
        case 'ikto'
            sol_kto = kto_default;
            sol_kto(model_struct(i).idx1) = sol(model_struct(i).idx2);
            sol_kto = [sol_kto, NaN(1, (max_param_len-length(sol_kto)))];
            sol_mx(:, i) = sol_kto;
        case 'ikslow1'
            sol_kslow1 = kslow1_default;
            sol_kslow1(model_struct(i).idx1) = sol(model_struct(i).idx2);
            sol_kslow1 = [sol_kslow1, NaN(1, (max_param_len-length(sol_kslow1)))];
            sol_mx(:, i) = sol_kslow1;
        case 'ikslow2'
            sol_kslow2 = kslow2_default;
            sol_kslow2(model_struct(i).idx1) = sol(model_struct(i).idx2);
            sol_kslow2 = [sol_kslow2, NaN(1, (max_param_len-length(sol_kslow2)))];
            sol_mx(:, i) = sol_kslow2;
        case 'ikss'
            sol_kss = kss_default;
            sol_kss(model_struct(i).idx1) = sol(model_struct(i).idx2);
            sol_kss = [sol_kss, NaN(1, (max_param_len-length(sol_kss)))];
            sol_mx(:, i) = sol_kss;
        case 'ikur'
            sol_kur = kur_default;
            sol_kur(model_struct(i).idx1) = sol(model_struct(i).idx2);
            sol_kur = [sol_kur, NaN(1, (max_param_len-length(sol_kur)))];
            sol_mx(:, i) = sol_kur;
        case 'ik1'
            sol_k1 = k1_default;
            sol_k1(model_struct(i).idx1) = sol(model_struct(i).idx2);
            sol_k1 = [sol_k1, NaN(1, (max_param_len-length(sol_k1)))];
            sol_mx(:, i) = sol_k1;
        end        
    end
    writematrix(string(current_names) , save_path, "Sheet","Parameters", "Range","A1");
    writematrix(sol_mx, save_path, "Sheet","Parameters", "Range","A2");
end

%%%%%
% surrogate-assisted optimization functions
%%%%%
function [sol, fval, maxei] = optim_ei(f, initx, inity, out_size)    
    [num_init_pt, num_var] = size(initx);
    
    % ops = cell(1,2);
    % ops{1} = 'ardsquaredexponential';
    % ops{2} = 'fic';
    gpmdl = fitrgp(initx, inity, 'KernelFunction','ardsquaredexponential', 'FitMethod','fic');

    % sequential acquisition
    sol = NaN(out_size, num_var);
    fval = NaN(out_size, 1);
    maxei = NaN(out_size-num_init_pt, 1);
    
    sol(1:num_init_pt, :) = initx;
    fval(1:num_init_pt) = inity;

    for i = (num_init_pt+1):out_size
        % new point
        [xnew, running_fval] = ei_acquisition(sol(1:(i-1), :), fval(1:(i-1)), gpmdl, 5, eps);
        [running_maxei, max_idx] = max(running_fval);
        maxei(i-num_init_pt) = running_maxei;

        % augment new point and update GP model
        xnew = xnew(max_idx, :);
        sol(i, :) = xnew;
        fval(i) = f(xnew);

        gpmdl = fitrgp(sol(1:i, :), fval(1:i), 'KernelFunction','ardsquaredexponential', 'FitMethod','fic');
    end
end

function [xnew, fval] = ei_acquisition(x, y, gpmdl, num_multi_start, tol)
    % EI-based acquisition function
    [fmin, min_idx] = min(y);
    [~, num_var] = size(x);

    % starting points of search
    if num_multi_start > 1
        start_pts = NaN(num_multi_start, num_var);
        start_pts(1, :) = x(min_idx, :);
        start_pts(2:end, :) = lhsdesign(num_multi_start-1, num_var);
    else
        start_pts = x(min_idx, :);
    end

    % optimization for searching next acquisition point
    xnew = NaN(num_multi_start, num_var);
    fval = NaN(num_multi_start, 1);
    obj_fun = @(x) obj_ei(x, fmin, gpmdl);
    optim_opts = optimoptions('fminunc', 'Display','off');

    for i = 1:num_multi_start
        [opt_out, running_fval] = fminunc(obj_fun, start_pts(i, :), optim_opts);
        xnew(i, 1:end) = opt_out;
        fval(i) = -running_fval;
    end

    valid_idx = fval > tol;
    xnew = xnew(valid_idx, :);
    fval = fval(valid_idx);
end

function z = obj_ei(x, fmin, gpmdl)
    % EI objective function in subroutine; - to formulate as minimization problem
    z = -ei(x, fmin, gpmdl);
end

function z = ei(x, fmin, gpmdl)
    % expected improvement
    [ypred, ysd, ~] = predict(gpmdl, x);
    d = fmin - ypred;
    dsig = d./ysd;
    z = d.*normcdf(dsig) + ysd.*normpdf(dsig);
end

function prog = bov(y)
    % best objective value over course of procedure
    ylen = length(y);
    prog = y;

    for i = 2:ylen
        if (isnan(prog(i)) || prog(i) > prog(i-1)) 
            prog(i) = prog(i-1);
        end
    end
end
