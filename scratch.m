%% GUI with interactive response-plot updates
% https://www.mathworks.com/help/control/ug/build-app-with-interactive-plot-updates.html
% display the step response of a second-order dynamic system
% of fixed natural frequency

clc
close all
clear variables

% initial values
zeta = 0.5; % damping ratio
wn = 2; % natural frequency
sys = tf(wn^2, [1, 2*zeta*wn, wn^2]);

f = figure;
ax = axes('Parent',f, 'position',[0.13, 0.39, 0.77, 0.54]);
h = stepplot(ax, sys);
setoptions(h, 'Xlim',[0, 10], 'Ylim',[0, 2]);

b = uicontrol('Parent',f, 'Style','slider', 'Position',[81, 54, 419, 23], ...
    'value',zeta, 'min',0, 'max',1);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f, 'Style','text', 'Position',[50, 54, 23, 23], ...
    'String','0', 'BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f, 'Style','text', 'Position',[500, 54, 23, 23], ...
    'String','1', 'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f, 'Style','text', 'Position',[240, 25, 100, 23], ...
    'String','Damping Ratio', 'BackgroundColor',bgcolor);

b.Callback = @(es, ed) updateSystem(h, tf(wn^2, [1, 2*(es.Value)*wn, wn^2]));

%% check validity of Ktrace2.m
clc
close all
clear variables

param = cell(1, 5);
param{1} = [30.0, 13.5, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335];
param{2} = [22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 1200.0, 170.0, 45.2, 5.7];
param{3} = [22.5, 7.7, 0.493, 0.0629, 2.058, 45.2, 5.7, 1200.0, 170.0, 45.2, 5.7];
param{4} = 2.5;
param{5} = [0.4067, 0.16, 0.16];

time_space = cell(1, 3);
time_step = 0.1;
holdT = 10;
P1T = 5000;
tH = 0:time_step:holdT;
tP1 = (holdT+time_step):time_step:P1T;
tP1_adj = tP1 - tP1(1);
t = [tH, tP1];
time_space{1} = t;
time_space{2} = tH;
time_space{3} = tP1_adj;

k = Ktrace2(param, -80, 50, time_space);
plot(t, k(:, 1))
hold on
    plot(t, k(:, 2))
    plot(t, k(:, 4))
hold off
axis tight
legend('I_{to}', 'I_{Kslow}', 'I_{ss}')
xlabel('Time (ms)')

%% check validity of Ktrace1.m
clc
close all
clear variables

holdV = -80;
holdT = 500; %msec
P1 = 50; %mV
P1T = 5000; % msec

param = [1, 2];

% generate time space
time_step = 0.001;
t = 0:time_step:P1T;
[~, hold_idx] = min(abs(t - holdT));

A = Ktrace1(param, holdV, hold_idx, P1, t);

%% pilot1 - fitting
clc
close all
clear variables

load ds_Ktrace_ko
cap_ko = 254.3;

t = ds_Ktrace_ko.time;
y = ds_Ktrace_ko.KO/cap_ko;

[~, hold_idx] = max(y);

%% pilot 1
clc
close all
clear variables

holdV = -70; %mV
holdT = 450; %ms
P1 = 50; %mV
sim_time = 25*1000; % ms
N0 = 50;

rng(7981);
init_gen = unifrnd(0, 0.1, N0, 27);
init_gen(:, 9) = unifrnd(0, 1, N0, 1);
init_gen(:, 18) = unifrnd(0, 1, N0, 1);
init_gen(:, 27) = unifrnd(0, 1, N0, 1);

[t, A] = Ktrace1(init_gen(7,:), holdV, holdT, P1, sim_time);

plot(t, A(:, 1))
hold on
    plot(t, A(:, 2))
    plot(t, A(:, 3))
hold off
legend('Ito', 'IKslow1', 'IKslow2')

%% compare the original Rasmusson model and the analytical model
clc
close all
clear variables

holding_p = -80; %mV
holding_t = 10; %msec
P1s = -70:10:50; %mV
P1_t = 25000; % msec

for i=1:length(P1s)
    [t, ~, A, ~] = RasmussonUnparam(holding_p, holding_t, P1s(i), P1_t, P1s(i), P1_t);
    
    [t1, A1] = analytic_form(holding_p, holding_t, P1s(i), P1_t);
    
    figure(1)
    hold on
        plot(t, A(:,61))
    hold off
    title('I_{to} - Original')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
    
    figure(2)
    hold on
        plot(t1, A1(:,1))
    hold off
    title('I_{to} - Analytical')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
    
    figure(3)
    hold on
        plot(t, A(:,65)) 
    hold off
    title('I_{Kslow} - Original')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
    
    figure(4)
    hold on
        plot(t1, A1(:,2))
    hold off
    title('I_{Kslow} - Analytical')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
end

%% only for one clamp voltage
clc
close all
clear variables

holdV = -80;
holdT = 500; %msec
P1 = 50; %mV
P1T = 5000; % msec

[t, ~, A, ~] = RasmussonUnparam(holdV, holdT, P1, P1T, P1, P1T);
[t1, A1] = analytic_form(holdV, holdT, P1, P1T);

figure(1)
hold on
    plot(t, A(:,61))
hold off
title('I_{to} - Original')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
axis tight

figure(2)
hold on
    plot(t1, A1(:,1))
hold off
title('I_{to} - Analytical')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
axis tight

figure(3)
hold on
    plot(t, A(:,65)) 
hold off
title('I_{Kslow} - Original')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
axis tight

figure(4)
hold on
    plot(t1, A1(:,2))
hold off
title('I_{Kslow} - Analytical')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
axis tight
