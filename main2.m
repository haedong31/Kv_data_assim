clc
close all
clear variables
format shortG

load ds_Ktrace_ko
cap_ko = 254.3;

t = ds_Ktrace_ko.time;
y = ds_Ktrace_ko.KO/cap_ko;

% time space
time_space = cell(1, 3);
time_space{1} = t;
hold_idx = 98; % manually found hold_idx
tH = t(1:hold_idx);
time_space{2} = tH;
tP1 = t(hold_idx+1:end);
tP1_adj = tP1 - tP1(1);
time_space{3} = tP1_adj;

[bfs, bchroms] = fit2(time_space, y, 50, 10, 4, 100);

param = cell(1, 5);
param{1} = bchroms(end, 1:11);
param{2} = bchroms(end, 12:19);
param{3} = bchroms(end, 20:27);
param{4} = bchroms(end, 28);
param{5} = bchroms(end, 29:31);
A = Ktrace2(param, -70, 50, time_space);

plot(t, A(:, 1))
hold on
    plot(t, A(:, 2))
    plot(t, A(:, 3))
    plot(t, A(:, 4))
hold off
legend('Ito', 'IKslow1', 'IKslow2', 'Iss')
