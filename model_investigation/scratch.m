%% DOE for IKr
clc
close all
clear variables

p0 = [0.022348, 0.01176, 0.047002, 0.0631, 0.013733, 0.038198, 6.89e-05, 0.04178, 5, 0.090821, 0.023391, 0.006497, 0.03268];
cv = IKr(p0, 50);

kf = 0.023761;  % kf; Rate constant for rapid delayed-rectifier K+ current:ms^-1
kb = 0.036778;  % kb; Rate constant for rapid delayed-rectifier K+ current:ms^-1

syms ck1(t) ck2(t) ok(t) ik(t)
ode1 = diff(ck1) == cv(1) - (cv(1) + cv(2) + kf)*ck1(t) + (-cv(1) + kb)*ck2(t) - cv(1)*ok(t) - cv(1)*ik(t);
ode2 = diff(ck2) == kf*ck1(t) - (kb + cv(3))*ck2(t) + cv(4)*ok(t);
ode3 = diff(ok) == cv(3)*ck2(t) - (cv(4) + cv(5))*ok(t) + cv(6)*ik(t);
ode4 = diff(ik) == cv(5)*ok(t) - cv(6)*ik(t);
odes = [ode1; ode2; ode3; ode4];

cond1 = ck1(0) == 0.992513e-3; % CK1; mERG channel closed state
cond2 = ck2(0) == 0.641229e-3; % CK2; mERG channel closed state
cond3 = ok(0) == 0.175298e-3; % OK; mERG channel open state
cond4 = ik(0) == 0.319129e-4; % IK; mERG channel inactivated state
conds = [cond1; cond2; cond3; cond4];

[ck1_sol(t), ck2_sol(t), ok_sol(t), ik_sol(t)] = dsolve(odes, conds);

figure(1)
fplot(ck1_sol)

double(ck1_sol(45))

tt = 1:1000;
y = zeros(length(tt), 1);
for i = 1:length(tt)
    disp(i)
   y(i) = double(ck1_sol(tt(i))); 
end

% solutions have root of z and vpa() forces to find roots
fn = matlabFunction(vpa(ck1_sol));
plot(tt, fn(tt))

%% check IKs
clc
close all
clear variables

time_step = 1;
holdT = 100;
P1T = 5000;

tH = 0:time_step:holdT;
tP1 = (holdT+time_step):time_step:P1T;
tP1_adj = tP1 - tP1(1);
t = [tH, tP1];

time_space = cell(1, 3);
time_space{1} = t;
time_space{2} = tH;
time_space{3} = tP1_adj;

p0 = [26.5, 0.128, 0.038, 4.81333e-06, 9.53333e-05, 0.00575, -91.1];
holdV = -80;
volt = -70:10:50;

for i = 1:length(volt)
    yKs = IKs(p0, holdV, volt(i), time_space);
    hold on
    plot(t, yKs)
    hold off
end

%% check IKI
clc
close all
clear variables

extra_Kconcent = 5400;
Ek = -91.1;
param = [0.2938, 210, 0.0896];
volt = -150:10:0;

yKI = zeros(1, length(volt));
for i = 1:length(volt)
    yKI(i) = param(3).*(extra_Kconcent./(extra_Kconcent + param(2))).*((volt(i) - Ek)./(1 + exp(param(1).*(volt(i) - Ek))));
end

plot(volt, yKI, "-o", "LineWidth",2)

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
