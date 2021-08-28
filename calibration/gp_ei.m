clc
close all
clear variables

%
% code arguments for calibration
%
group_name = 'wt';
save_dir = strcat('calib_result1_', group_name);

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


reps = 100;
ninit = 12;
out_size = 50;
prog_ei = NaN(reps, out_size);
for r = 1:reps
    [sol, fval, maxei] = optim_ei(@gpfun, ninit, 2, out_size);
    running_prog_ei = bov(fval);
    prog_ei(r, :) = running_prog_ei';
end

plot(mean(prog_ei, 1))
hold on
xline(ninit, '--', 'Color',[0.5,0.5,0.5])
hold off

% RMSE objective function; main optimization objective
function z = obj_rmse(p)
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

% EI-based acquisition function
function [xnew, fval] = ei_acquisition(x, y, gpmdl, num_multi_start, tol)
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

% EI objective function in subroutine; - to formulate as minimization problem
function z = obj_ei(x, fmin, gpmdl)
    z = -ei(x, fmin, gpmdl);
end

% expected improvement
function z = ei(x, fmin, gpmdl)
    [ypred, ysd, ~] = predict(gpmdl, x);
    d = fmin - ypred;
    dsig = d./ysd;
    z = d.*normcdf(dsig) + ysd.*normpdf(dsig);
end

% best objective value over course of procedure
function prog = bov(y)
    ylen = length(y);
    prog = y;

    for i = 2:ylen
        if (isnan(prog(i)) || prog(i) > prog(i-1)) 
            prog(i) = prog(i-1);
        end
    end
end
