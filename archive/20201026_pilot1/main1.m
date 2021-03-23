clc
close all
clear variables
format shortG

load ds_Ktrace_ko
cap_ko = 254.3;

t = ds_Ktrace_ko.time;
y = ds_Ktrace_ko.KO/cap_ko;

[bfs, bchroms] = fit1(t, y, 50, 10, 4, 50000);

holdV = -70;
[~, hold_idx] = max(y);
P1 = 50;
A = Ktrace1(bchroms(end, :), holdV, hold_idx, P1, t);

plot(t, A(:, 1))
hold on
    plot(t, A(:, 2))
    plot(t, A(:, 3))
    plot(t, A(:, 4))
hold off
legend('Ito', 'IKslow1', 'IKslow2', 'Iss')
