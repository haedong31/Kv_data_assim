function [amp, tau, ss] = tri_exp_fit(t, yksum)
    % objective function
    obj = @(x) rmse_obj(x, yksum, t);

    % initial value
    A0 = [1];
    tau0 = [1];
    C0 = [1];

    % constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = [];

    [sol, min_avl] = fmincon(obj, x0, A, b, Aeq, beq, lb, ub, nonlcon);

end

function [z] = rmse_obj(x, y, t)
    yfit = tri_exp_fun(t, x(1:3), x(4:6), x(7));
    z = sqrt(mean((y-yfit).^2));
end

function [y] = tri_exp_fun(t, amp, tau, ss)
    y = amp(1)*exp(-t./tau(1)) + amp(2)*exp(-t./tau(2)) + amp(3)*exp(-t./tau(3)) + ss;
end
