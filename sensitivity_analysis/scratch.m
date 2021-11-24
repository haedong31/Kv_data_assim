clc
clearvars
close all

kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];

kslow2_default = zeros(1,11);
kslow2_default(1:8) = kslow1_default(1:8);
kslow2_default(9:11) = [5334, 4912, 0.05766];

kss_default = zeros(1,7);
kss_default(1:3) = kslow1_default([1,3,4]);
kss_default(4:7) = [0.0862, 1235.5, 13.17, 0.0428];

volts = -30:10:50;
num_volts = length(volts);

kto_gv = NaN(num_volts,8);
kslow1_gv = NaN(num_volts,4);
kslow2_gv = NaN(num_volts,4);
kss_gv = NaN(num_volts,2);
for i=1:length(volts)
    kto_gv(i,:) = ikto_gating_variables(kto_default, volts(i));
    kslow1_gv(i,:) = ikslow1_gating_variables(kslow1_default, volts(i));
    kslow2_gv(i,:) = ikslow2_gating_variables(kslow2_default, volts(i));
    kss_gv(i,:) = ikss_gating_variables(kss_default, volts(i));
end

figure('Name','ikto','Color','w')
subplot(2,4,1)
plot(volts, kto_gv(:,1), '-o', 'LineWidth',1.5,'Color','red')
xlabel('Voltages (mV)')
ylabel('a_{ss}')
grid on
axis tight

subplot(2,4,2)
plot(volts, kto_gv(:,2), '-o', 'LineWidth',1.5,'Color','blue')
xlabel('Voltages (mV)')
ylabel('i_{ss}')
grid on
axis tight

subplot(2,4,5)
plot(volts, kto_gv(:,3), '-o', 'LineWidth',1.5,'Color','red')
xlabel('Voltages (mV)')
ylabel('tau_{a}')
grid on
axis tight

subplot(2,4,6)
plot(volts, kto_gv(:,3), '-o', 'LineWidth',1.5,'Color','blue')
xlabel('Voltages (mV)')
ylabel('tau_{i}')
grid on
axis tight

subplot(2,4,3)
plot(volts, kto_gv(:,5), '-o', 'LineWidth',1.5,'Color','red')
xlabel('Voltages (mV)')
ylabel('ap_{ss}')
grid on
axis tight

subplot(2,4,4)
plot(volts, kto_gv(:,6), '-o', 'LineWidth',1.5,'Color','blue')
xlabel('Voltages (mV)')
ylabel('ip_{ss}')
grid on
axis tight

subplot(2,4,7)
plot(volts, kto_gv(:,7), '-o', 'LineWidth',1.5,'Color','red')
xlabel('Voltages (mV)')
ylabel('taup_{a}')
grid on
axis tight

subplot(2,4,8)
plot(volts, kto_gv(:,8), '-o', 'LineWidth',1.5,'Color','blue')
xlabel('Voltages (mV)')
ylabel('taup_{i}')
grid on
axis tight

figure('Name','ikslow1','Color','w')
subplot(2,2,1)
plot(volts, kslow1_gv(:,1), '-o', 'LineWidth',1.5,'Color','red')
xlabel('Voltages (mV)')
ylabel('a_{ss}')
grid on
axis tight

subplot(2,2,2)
plot(volts, kslow1_gv(:,2), '-o', 'LineWidth',1.5,'Color','blue')
xlabel('Voltages (mV)')
ylabel('i_{ss}')
grid on
axis tight

subplot(2,2,3)
plot(volts, kslow1_gv(:,3), '-o', 'LineWidth',1.5,'Color','red')
xlabel('Voltages (mV)')
ylabel('tau_{a}')
grid on
axis tight

subplot(2,2,4)
plot(volts, kslow1_gv(:,4), '-o', 'LineWidth',1.5,'Color','blue')
xlabel('Voltages (mV)')
ylabel('tau_{i}')
grid on
axis tight

function gv = ikto_gating_variables(p, V)
    gv = zeros(8,1);
    
    % for Ikto
    alpha1 = p(9).*exp(p(7).*(V+p(1)));
    beta1 = p(10).*exp(-p(8).*(V+p(1)));
    
    alpha2_temp1 = p(11).*exp((V+p(2))./(-1.0.*p(6)));
    alpha2_temp2 = p(12).*exp((V+p(2)+p(3))./(-1.0.*p(6)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(13).*exp((V+p(2)+p(3))./p(6));
    beta2_temp2 = p(14).*exp((V+p(2)+p(3))./p(6));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    gv(1) = alpha1./(alpha1+beta1);
    gv(2) = alpha2./(alpha2+beta2);
    gv(3) = 1./(alpha1+beta1);
    gv(4) = 1./(alpha2+beta2);

    % for Ikto phosphorylated
    alpha1p = p(9).*exp(p(7).*(V+p(1)-p(4)));
    beta1p = p(10).*exp(-p(8).*(V+p(1)-p(4)));
    
    alpha2_temp1p = p(11).*exp((V+p(2)-p(5))./(-1.0.*p(6)));
    alpha2_temp2p = p(12).*exp((V+p(2)+p(3)-p(5))./(-1.0.*p(6)));
    alpha2p = alpha2_temp1p./(1.0+alpha2_temp2p);
    
    beta2_temp1p = p(13).*exp((V+p(2)+p(3)-p(5))./p(6));
    beta2_temp2p = p(14).*exp((V+p(2)+p(3)-p(5))./p(6));
    beta2p = beta2_temp1p./(1.0+beta2_temp2p);

    gv(5) = alpha1p./(alpha1p+beta1p);
    gv(6) = alpha2p./(alpha2p+beta2p);
    gv(7) = 1./(alpha1p+beta1p);
    gv(8) = 1./(alpha2p+beta2p);
end

function gv = ikslow1_gating_variables(p, V)
    gv = zeros(4, 1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    gv(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    gv(4) = p(9)-p(10)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function gv = ikslow2_gating_variables(p, V)
    gv = zeros(4,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    gv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    gv(3) = p(6)./(exp(p(7)*(V+p(3))) + exp(-p(7)*(V+p(3))))+p(8); % taua
    gv(4) = p(9) - p(10)./(1.0+exp((p(2)+V)./p(5)));
end

function gv = ikss_gating_variables(p, V)
    gv = zeros(2,1);
    gv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    gv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end