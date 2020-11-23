%% compare the original Rasmusson model and the analytical model
clc
close all
clear variables

holding_p = -80; %mV
holding_t = 10; %msec
P1s = -70:10:50; %mV
P1_t = 5000; % msec

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

%% pilot1 - fitting
clc
close all
clear variables

load ds_Ktrace_ko
cap_ko = 254.3;

t = ds_Ktrace_ko.time;
y = ds_Ktrace_ko.KO/cap_ko;

[~, hold_idx] = max(y);
