clc
close all
clear variables

function z = optim_ei(x)
    disp(x)
end

% EI-based acquisition function
function xnew = ei_acquisition(x, y, gpmdl, num_multi_start, tol)
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
    xnew = NaN(num_multi_start, num_var+1);
    obj_fun = @(x) obj_ei(x, fmin, gpmdl);
    for i = 1:num_multi_start
        [opt_out, fval] = fminunc(obj_fun, start_pts(i, :));
        xnew(i, 1:end-1) = opt_out;
        xnew(i, end) = fval;
    end
end

% objective function; - for mnimization
function z = obj_ei(x, fmin, gpmdl)
    z = -ei(x, fmin, gpmdl);
end

% expected improvement
function z = ei(x, fmin, gpmdl)
    [ypred, ysd, ~] = predict(gpmdl, x);
    d = fmin - ypred;
    dsig = d./ysd;
    z = d*normcdf(dsig) + ysd*normpdf(dsig);
end
