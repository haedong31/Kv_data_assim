function [amp, tau, s] = tri_exp_fit(t, yksum, hold_idx)
    % objective function
    t_trunc = t(hold_idx:end);
    t_trunc = t_trunc - t_trunc(1);
    yksum_trunc = yksum(hold_idx:end);
    obj = @(x) rmse_obj(x, yksum_trunc, t_trunc);

    [amp_sum, time_stats] = trace_stat(t, yksum, hold_idx);
    
    % initial value
    x0 = zeros(7,1);
    x0(1:3) = amp_sum*[0.55, 0.15, 0.15];
    x0(4:6) = time_stats(1:3);
    x0(7) = amp_sum*0.15;

    % constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = [];

    % lower bound (positive)
    lb = zeros(7,1);

    % upper bound
    ub = [ones(3,1)*amp_sum; ones(3,1)*t_trunc(end); amp_sum];

    [sol, min_avl] = fmincon(obj, x0, A, b, Aeq, beq, lb, ub, nonlcon);
    disp(min_avl)

    amp = sol(1:3);
    tau = sol(4:6);
    s = sol(7);
end

function [z] = rmse_obj(x, y, t)
    yfit = tri_exp_fun(t, x(1:3), x(4:6), x(7));
    z = sqrt(mean((y-yfit).^2));
end

function [y] = tri_exp_fun(t, amp, tau, s)
    y = amp(1)*exp(-t./tau(1)) + amp(2)*exp(-t./tau(2)) + amp(3)*exp(-t./tau(3)) + s;
end
