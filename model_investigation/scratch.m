%% line objects
x = 0:0.1:2*pi;
h(1) = plot(x, sin(x));
hold on
h(2) = plot(x, cos(x));
hold off

legend(h(2), "hello world");
set(h(2), "Color","red", "LineWidth",2);

%% check IKss
clc
close all
clear variables

time_space = cell(1, 3);
time_step = 1;
holdT = 100;
P1T = 5000;
tH = 0:time_step:holdT;
tP1 = (holdT+time_step):time_step:P1T;
tP1_adj = tP1 - tP1(1);
t = [tH, tP1];

time_space{1} = t;
time_space{2} = tH;
time_space{3} = tP1_adj;

p0 = [22.5, 7.7, 39.3, 0.0862, 13.17, 0.05, -91.1];
volt = -90:10:50;

for i = 1:length(volt)
    current_trc = IKss(p0, -100, volt(i), time_space);
    hold on
        plot(t, current_trc)
    hold off
end

%% under holding potential
clc
close all
clear variables

holding_p = -80; %mV
holding_t = 500; %msec
P1s = -70:10:50; %mV
P1_t = 5000; % msec

TR_ori = cell(1, length(P1s));
states = cell(1, length(P1s));
TR_anal = cell(1, length(P1s)); 
for i=1:length(P1s)
    [t, S, A, ~] = RasmussonUnparam(holding_p, holding_t, P1s(i), P1_t, P1s(i), P1_t);
    states{i} = S;
    
    TR = zeros(8, 2);
    TR(1, 1) = A(1, 5);
    TR(1, 2) = A(end, 5);
    TR(2, 1) = A(1, 15);
    TR(2, 2) = A(end, 15);
    TR(3, 1) = A(1, 6);
    TR(3, 2) = A(end, 6);
    TR(4, 1) = A(1, 16);
    TR(4, 2) = A(end, 16);
    
    TR(5, 1) = A(1, 7);
    TR(5, 2) = A(end, 7);
    TR(6, 1) = A(1, 20);
    TR(6, 2) = A(end, 20);
    TR(7, 1) = A(1, 8);
    TR(7, 2) = A(end, 8);
    TR(8, 1) = A(1, 21);
    TR(8, 2) = A(end, 21);
    TR_ori{i} = TR;

    [t1, I1, TR1] = analytic_form(holding_p, holding_t, P1s(i), P1_t);
    TR_anal{i} = TR1;
    
    [~, hold_idx] = min(abs(t - holding_t));
    [~, hold_idx1] = min(abs(t1 - holding_t));
    
    figure(1)
    hold on
        plot(t(1:hold_idx-1), A(1:hold_idx-1,61))
    hold off
    title('I_{to} - Original')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
    
    figure(2)
    hold on
        plot(t1(1:hold_idx1), I1(1:hold_idx1,1))
    hold off
    title('I_{to} - Analytical')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
    
    figure(3)
    hold on
        plot(t(1:hold_idx-1), A(1:hold_idx-1,65)) 
    hold off
    title('I_{Kslow} - Original')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
    
    figure(4)
    hold on
        plot(t1(1:hold_idx1), I1(1:hold_idx1,2))
    hold off
    title('I_{Kslow} - Analytical')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
end


%% check analytical solutions
clc
close all
clear variables

alphaa = 0.0302049585228441;
betaa = 8.94547949889369;

alphai = 0.330343529000106;
betai = 1.23805046797418e-06;
% alphai = 0.217886158136666;
% betai = 5.16495883881682e-06;

ass10 = 0.265563e-2;
iss10 = 0.999977;
ass20 = 0.417069e-3;
iss20 = 0.998543;

ass1 = alphaa/(alphaa + betaa);
atau1 = 1/(alphaa + betaa);
iss1 = alphai/(alphai + betai);
itau1 = 1/(alphai + betai);

ass2 = 0.000571010111115601;
atau2 = 77.6049181019847;
iss2 = 0.997773872537225;
itau2 = 1030.37844166867;

t = 0:0.01:500;
ag1 = ass1 - (ass1 - ass10).*exp(-t./atau1);
ig1 = iss1 - (iss1 - iss10).*exp(-t./itau1);

ag2 = ass2 - (ass2 - ass20).*exp(-t./atau2);
ig2 = iss2 - (iss2 - iss20).*exp(-t./itau2);

plot(t, ag1.*ig1)
plot(t, ag2.*ig2)

%% long holding potential
clc
close all
clear variables

holding_p = -80; %mV
holding_t = 10; %msec
P1s = -70:10:50; %mV
P1_t = 5000; % msec

TR_ori = cell(1, length(P1s));
states = cell(1, length(P1s));
TR_anal = cell(1, length(P1s)); 
for i=1:length(P1s)
    [t, S, A, ~] = RasmussonUnparam(holding_p, holding_t, P1s(i), P1_t, P1s(i), P1_t);
    states{i} = S;
    
    TR = zeros(8, 2);
    TR(1, 1) = A(1, 5);
    TR(1, 2) = A(end, 5);
    TR(2, 1) = A(1, 15);
    TR(2, 2) = A(end, 15);
    TR(3, 1) = A(1, 6);
    TR(3, 2) = A(end, 6);
    TR(4, 1) = A(1, 16);
    TR(4, 2) = A(end, 16);
    
    TR(5, 1) = A(1, 7);
    TR(5, 2) = A(end, 7);
    TR(6, 1) = A(1, 20);
    TR(6, 2) = A(end, 20);
    TR(7, 1) = A(1, 8);
    TR(7, 2) = A(end, 8);
    TR(8, 1) = A(1, 21);
    TR(8, 2) = A(end, 21);
    TR_ori{i} = TR;

    [t1, I1, TR1] = analytic_form(holding_p, holding_t, P1s(i), P1_t);
    TR_anal{i} = TR1;
        
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
        plot(t1, I1(:,1))
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
        plot(t1, I1(:,2))
    hold off
    title('I_{Kslow} - Analytical')
    xlabel('Time (ms)')
    ylabel('Current (pA/pF)')
    axis tight
end

%% for each current
clc
close all
clear variables

holding_p = -80; %mV
holding_t = 500; %msec
P1 = 30; %mV
P1_t = 5000; % msec

[t, I, TR] = analytic_form(holding_p, holding_t, P1, P1_t);

plot(t, I(:,1))
title('I_{to} - Analytical')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
axis tight
