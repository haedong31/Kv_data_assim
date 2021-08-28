clc
close all
clear variables

tic
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
toc

% optim_ei(@obj_rmse, 100);
function [sol, fval, maxei] = optim_ei(f, num_init_pt, num_var, out_size)
    % initialization
    initx = lhsdesign(num_init_pt, num_var);
    inity = f(initx);
    
    gpmld = fitrgp(initx, inity, 'KernelFunction','ardsquaredexponential');

    % sequential acquisition
    sol = NaN(out_size, num_var);
    fval = NaN(out_size, 1);
    maxei = NaN(out_size-num_init_pt, 1);
    
    sol(1:num_init_pt, :) = initx;
    fval(1:num_init_pt) = inity;

    for i = (num_init_pt+1):out_size
        % new point
        [xnew, running_fval] = ei_acquisition(sol(1:(i-1), :), fval(1:(i-1)), gpmld, 5, eps);
        [running_maxei, max_idx] = max(running_fval);
        maxei(i-num_init_pt) = running_maxei;

        % augment new point and update GP model
        xnew = xnew(max_idx, :);
        sol(i, :) = xnew;
        fval(i) = f(xnew);

        gpmld = fitrgp(sol(1:i, :), fval(1:i), 'KernelFunction','ardsquaredexponential');
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

% objective function; - to formulate as minimization problem
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

function prog = bov(y)
    ylen = length(y);
    prog = y;

    for i = 2:ylen
        if (isnan(prog(i)) || prog(i) > prog(i-1)) 
            prog(i) = prog(i-1);
        end
    end
end

% Goldstein-Price function for test
function y = gpfun(x)
    m = 8.6928;
    s = 2.4269;
    x1 = 4.*x(:, 1) - 2;
    x2 = 4.*x(:, 2) - 2;

    a = 1 + (x1 + x2 + 1).^2.*(19 - 14.*x1 + 3.*x1.^2 - 14.*x2 + 6.*x1.*x2 + 3.*x2.^2);
    b = 30 + (2.*x1 - 3.*x2).^2.*(18 - 32.*x1 + 12.*x1.^2 + 48.*x2 - 36.*x1.*x2 + 27.*x2.^2);

    y = log(a.*b);
    y = (y - m)./s;
end
