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
ideal_end_time = 4.6*1000;
ek = -91.1;

% optimization options
reps = 100;
ninit = 12;
out_size = 50;

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

global model_struct
model_struct = gen_model_struct(current_names, tune_idx1_kto, tune_idx1_kslow1, ...
    tune_idx1_kslow2, tune_idx1_kss, tune_idx1_kur, tune_idx1_k1);

global volt_space
volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volts;
volt_space{3} = ek;

global time_space
global yksum

for i = 1:floor(num_files/2)
    fprintf('##### RUN')

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
    
    % MC-style repetitions
    prog = NaN(reps, out_size);
    for r = 1:reps
        [sol, fval, maxei] = optim_ei(@obj_rmse, ninit, 2, out_size);
        running_prog = bov(fval);
        prog(r, :) = running_prog';

        fprintf('[File %i/%i] %s [Reps %i/%i] Min RMSE: %f \n', i, file_names{loop_idx(i)}, r, reps, min(running_prog));
    end
end

% plot bov
plot(mean(prog, 1))
hold on
xline(ninit, '--', 'Color',[0.5,0.5,0.5])
hold off

%%%%%
% problem specific functions
%%%%%
function z = obj_rmse(p)
    % RMSE objective function; main optimization objective
    global model_struct
    global volt_space
    global time_space
    global yksum

    hold_idx = length(time_space{2});
    volts = volt_space{2};
    num_volts = length(volts);

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{3} = time_space;
    protocol{4} = volt_space{3};

    rmse_list = zeros(num_volts, 1);    
    for i = 1:num_volts
        yksum_i = yksum(:, i);
        protocol{2} = volts(i);

        [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol);
        
        % check validity
        [~, peak_idx] = max(yksum_hat);
        check_pt1 = any(isnan(yksum_hat));
        check_pt2 = any(yksum_hat < 0);
        check_pt3 = var(yksum_hat(1:hold_idx)) > 0.4812e-4; % not stable at hold_volt
        check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse

        if(check_pt1 || check_pt2 || check_pt3 || check_pt4)
            rmse_list(i) = 1e+3; % arbitrary big number
        else
            rmse_list(i) = sqrt(mean((yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2));
        end
    end
    z = sum(rmse_list);
end

function scaledp = scale_param(p)
end

function model_struct = gen_model_struct(current_names)
    num_currents = length(current_names);
    
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
end

%%%%%
% surrogate-assisted optimization functions
%%%%%
function [sol, fval, maxei] = optim_ei(f, num_init_pt, num_var, out_size)
    % initialization
    initx = lhsdesign(num_init_pt, num_var);
    inity = f(initx);
    
    ops = cell(1,2);
    ops{1} = 'ardsquaredexponential';
    ops{2} = 'fic';
    gpmdl = fitrgp(initx, inity, 'KernelFunction',ops{1}, 'FitMethod',ops{2});

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

        gpmdl = fitrgp(sol(1:i, :), fval(1:i), 'KernelFunction','ardsquaredexponential');
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
