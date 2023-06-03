%% run kcurrent_model2 with default parameter values
clc
close all
clearvars

current_names = {'ikto', 'ikslow1', 'ikslow2', 'ikss'};
num_currents = length(current_names);

% tunning index in individual current models
tune_idx1_kto = 1:13;
tune_idx1_kslow1 = 1:11;
tune_idx1_kslow2 = 1:3;
tune_idx1_kss = 1:4;

% protocol
hold_volt = -70;
min_volt = -50;
volt = 0;
ek = -91.1;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;

idx_info1 = cell(1, num_currents);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        idx_info1{i} = tune_idx1_kto;
    case 'ikslow1'
        idx_info1{i} = tune_idx1_kslow1;
    case 'ikslow2'
        idx_info1{i} = tune_idx1_kslow2;
    case 'ikss'
        idx_info1{i} = tune_idx1_kss;
    case 'ikur'
        idx_info1{i} = tune_idx1_kur;
    case 'ik1'
        idx_info1{i} = tune_idx1_k1;
    end
end

% tunning index in decision variable, p
idx_info2 = cell(1, num_currents);
cul_idx_len = 0;
for i = 1: num_currents
    cul_idx_len_old = cul_idx_len;
    cul_idx_len = cul_idx_len + length(idx_info1{i});
    idx_info2{i} = (1+cul_idx_len_old):(cul_idx_len);
end

% create model structure
model_info = [current_names; idx_info1; idx_info2];
field_names = {'name', 'idx1', 'idx2'};
model_struct = cell2struct(model_info, field_names, 1);

% default values
kto_default = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.3846];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, 6.1, 18.0, 2.058, 803.0, 0.16];
kslow2_default = [4912, 5334, 0.16];
kss_default = [0.0862, 1235.5, 13.17, 0.0611];
kur_default = [270, 1050, 0];
k1_default = [59.215, 5.476, 594.31, 4.753, 1.02, 0.2385, 0.8, 0.08032, 0.06175, 0.5143];

volt_space = cell(3, 1);
volt_space{1} = hold_volt;
volt_space{2} = volt;
volt_space{3} = ek;

p0 = zeros(1, cul_idx_len);
for i = 1:num_currents
    switch current_names{i}
    case 'ikto'
        p0(idx_info2{i}) = kto_default(tune_idx1_kto);
    case 'ikslow1'
        p0(idx_info2{i}) = kslow1_default(tune_idx1_kslow1);
    case 'ikslow2'
        p0(idx_info2{i}) = kslow2_default(tune_idx1_kslow2);
    case 'ikss'
        p0(idx_info2{i}) = kss_default(tune_idx1_kss);
    case 'ikur'
        p0(idx_info2{i}) = kur_default(tune_idx1_kur);
    case 'ik1'
        p0(idx_info2{i}) = k1_default(tune_idx1_kur);
    end    
end

time_space = cell(1, 3);
t = 1:ideal_end_time;
time_space{1} = t;
time_space{2} = t(1:ideal_hold_time);
pulse_t = t(ideal_hold_time+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

protocol = cell(4, 1);
protocol{1} = volt_space{1};
protocol{2} = volt_space{2};
protocol{3} = time_space;
protocol{4} = volt_space{3};

[~, comp_currents] = kcurrent_model2(p0, model_struct, protocol);

figure('Color','w', 'Position',[100,100,520,400])
subplot(2,2,1)
plot(t,comp_currents{1},'LineWidth',2,'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',11','FontWeight','bold')
box off

subplot(2,2,2)
plot(t,comp_currents{2}, 'LineWidth',2, 'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',11','FontWeight','bold')
box off

subplot(2,2,3)
plot(t,comp_currents{3},'LineWidth',2,'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',11','FontWeight','bold')
box off

subplot(2,2,4)
plot(t,comp_currents{4}, 'LineWidth',2, 'Color','red')
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
set(gca,'FontName','Arial','FontSize',11','FontWeight','bold')
box off

% in one figure
iksum = comp_currents{1} + comp_currents{2} + comp_currents{3} + comp_currents{4};

figure('Color','w')
plot(t,comp_currents{1}, '-', 'Color','black', 'LineWidth',1.5)
hold on
plot(t,comp_currents{2}, '--','Color','black', 'LineWidth',1.5)
plot(t,comp_currents{3}, ':','Color','black', 'LineWidth',1.5)
plot(t,comp_currents{4}, '-.','Color','black', 'LineWidth',1.5)
plot(t,iksum, 'Color','red', 'LineWidth',1.5)
hold off
axis tight
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
legend('I_{Kto}','I_{Kslow1}','I_{Kslow2}','I_{Kss}','I_{Ksum}')
set(gca,'FontName','Arial','FontSize',11','FontWeight','bold')
