clc
clearvars
close all

% import design matrices
dgn_kto = table2array(readtable('ikto_dgn.csv'));
dgn_kslow1 = table2array(readtable('ikslow1_dgn.csv'));
dgn_kslow2 = table2array(readtable('ikslow2_dgn.csv'));
dgn_kss = table2array(readtable('ikss_dgn.csv'));

% bounds
kto_default = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, 0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
kslow1_default = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
kslow2_default = [5334, 4912, 0.05766];
kss_default = [0.0862, 1235.5, 13.17, 0.0428];

lb_kto = [-70, eps, eps, eps, eps, 1, eps, eps, eps, eps, eps, eps, eps, eps];
lb_kto(7:14) = kto_default(7:14)*0.05;
lb_kslow1 = [-70, -70, -70, 1, 1, eps, eps, eps, 50, eps];

lb_kslow2 = [5000, eps];
lb_kss = [eps, eps, eps];

ub_kto = [70, 30, 40, 35, 20, 30, 1, 1, 1, 10, 0.005, 0.3, 0.005, 0.5];
ub_kto(7:14) = kto_default(7:14)*1.95;
ub_kslow1 = [50, 50, 50, 10, 50, 50, 1, 100, 1000, 50];
ub_kslow2 = [10000, 5000];
ub_kss = [1, 2000, 100];

% voltages
volt_range = 3:11;
min_volt = -50;
volts = NaN(length(volt_range), 1);
for i = 1:length(volt_range)
    volts(i) = min_volt + (volt_range(i)-1)*10;
end

% time space
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;

time_space = cell(1,3); % t, holdt, pulset
t = 1:ideal_end_time;
time_space{1} = t;

time_space{2} = t(1:ideal_hold_time);

pulse_t = t(ideal_hold_time+1:end);
pulse_t_adj = pulse_t - pulse_t(1);
time_space{3} = pulse_t_adj;

protocol_info = cell(4, 1);
protocol_info{1} = -70;
protocol_info{3} = -91.1;
protocol_info{4} = time_space;

%% ikto
tic
p = NaN(17,1);
p(15:17) = kto_default(15:17);
[num_dgn, num_param] = size(dgn_kto);
response_kto = NaN(num_dgn,6);
for i=1:num_dgn
    dgn = dgn_kto(i,:);
    
    for j=1:num_param
        if dgn(j) == 1
            p(j) = ub_kto(j);
        else
            p(j) = lb_kto(j);
        end
    end
    
    o = NaN(length(volts),6);
    for j=1:length(volts)
        protocol_info{2} = volts(j);
         o(j,:) = ikto_biomarkers(p, protocol_info);
    end
    response_kto(i,:) = sum(o,1);
end
toc

%% ikslow1
tic
p = NaN(13,1);
p(11:13) = kslow1_default(11:13);
[num_dgn, num_param] = size(dgn_kslow1);
response_kslow1 = NaN(num_dgn,6);
for i=1:num_dgn
    dgn = dgn_kslow1(i,:);
    
    for j=1:num_param
        if dgn(j) == 1
            p(j) = ub_kslow1(j);
        else
            p(j) = lb_kslow1(j);
        end
    end
    
    o = NaN(length(volts),6);
    for j=1:length(volts)
        protocol_info{2} = volts(j);
         o(j,:) = ikslow1_biomarkers(p, protocol_info);
    end
    response_kslow1(i,:) = sum(o,1);
end
toc

%% ikslow2
tic
p = NaN(11,1);
p(1:8) = kslow1_default(1:8);
p(11) = kslow2_default(3);
[num_dgn, num_param] = size(dgn_kslow2);
response_kslow2 = NaN(num_dgn,6);
for i=1:num_dgn
    dgn = dgn_kslow2(i,:);
    
    for j=1:num_param
        if dgn(j) == 1
            p(j+8) = ub_kslow2(j);
        else
            p(j+8) = lb_kslow2(j);
        end
    end
    
    o = NaN(length(volts),6);
    for j=1:length(volts)
        protocol_info{2} = volts(j);
         o(j,:) = ikslow2_biomarkers(p, protocol_info);
    end
    response_kslow2(i,:) = sum(o,1);
end
toc

%% ikss
tic
p = NaN(7,1);
p(1:3) = kslow1_default([1,3,4]);
p(7) = kss_default(4);
[num_dgn, num_param] = size(dgn_kss);
response_kss = NaN(num_dgn,6);
for i=1:num_dgn
    dgn = dgn_kss(i,:);
    
    for j=1:num_param
        if dgn(j) == 1
            p(j+3) = ub_kss(j);
        else
            p(j+3) = lb_kss(j);
        end
    end
    
    o = NaN(length(volts),6);
    for j=1:length(volts)
        protocol_info{2} = volts(j);
         o(j,:) = ikss_biomarkers(p, protocol_info);
    end
    response_kss(i,:) = sum(o,1);
end
toc

%% analysis ikto
clc
clearvars
close all

load('response_kto.mat')
dgn_kto = table2array(readtable('ikto_dgn.csv'));

fig = figure('Color','w','Position',[100,100,900,500]);
orient(fig,'landscape')

subplot(2,3,1)
[theta,idx,phiinv] = doe_plot(response_kto(:,1),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 1', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,2)
[theta,idx,phiinv] = doe_plot(response_kto(:,2),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 2', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,3)
[theta,idx,phiinv] = doe_plot(response_kto(:,3),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 3', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,4)
[theta,idx,phiinv] = doe_plot(response_kto(:,4),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.8,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 4', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,5)
[theta,idx,phiinv] = doe_plot(response_kto(:,5),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 5', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,6)
[theta,idx,phiinv] = doe_plot(response_kto(:,6),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+500,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 6', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

%% ikslow1
clc
clearvars
close all

load('response_kslow1.mat')
dgn_kslow1 = table2array(readtable('ikslow1_dgn.csv'));

fig = figure('Color','w','Position',[100,100,900,500]);
orient(fig,'landscape')

subplot(2,3,1)
[theta,idx,phiinv] = doe_plot(response_kslow1(:,1),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 1', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,2)
[theta,idx,phiinv] = doe_plot(response_kslow1(:,2),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1.5,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 2', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,3)
[theta,idx,phiinv] = doe_plot(response_kslow1(:,3),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 3', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,4)
[theta,idx,phiinv] = doe_plot(response_kslow1(:,4),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 4', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,5)
[theta,idx,phiinv] = doe_plot(response_kslow1(:,5),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1.8,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 5', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,6)
[theta,idx,phiinv] = doe_plot(response_kslow1(:,6),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+500,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 6', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

%% ikslow2
clc
clearvars
close all

load('response_kslow2.mat')
dgn_kslow2 = table2array(readtable('ikslow2_dgn.csv'));

fig = figure('Color','w','Position',[100,100,900,500]);
orient(fig,'landscape')

subplot(2,3,1)
[theta,idx,phiinv] = doe_plot(response_kslow2(:,1),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.001,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 1', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,2)
[theta,idx,phiinv] = doe_plot(response_kslow2(:,2),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.2,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 2', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,3)
[theta,idx,phiinv] = doe_plot(response_kslow2(:,3),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.3,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 3', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,4)
[theta,idx,phiinv] = doe_plot(response_kslow2(:,4),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.3,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 4', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,5)
[theta,idx,phiinv] = doe_plot(response_kslow2(:,5),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.002,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 5', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,6)
[theta,idx,phiinv] = doe_plot(response_kslow2(:,6),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 6', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

%% ikss
clc
clearvars
close all

load('response_kss.mat')
dgn_kss = table2array(readtable('ikss_dgn.csv'));

fig = figure('Color','w','Position',[100,100,900,500]);
orient(fig,'landscape')

subplot(2,3,1)
[theta,idx,phiinv] = doe_plot(response_kss(:,1),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.8,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 1', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,2)
[theta,idx,phiinv] = doe_plot(response_kss(:,2),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.3,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2)-1,'Biomarker 2', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,3)
[theta,idx,phiinv] = doe_plot(response_kss(:,3),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.07,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 3', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,4)
[theta,idx,phiinv] = doe_plot(response_kss(:,4),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.03,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2)-0.2,'Biomarker 4', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,5)
[theta,idx,phiinv] = doe_plot(response_kss(:,5),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.5,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2)-2,'Biomarker 5', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

subplot(2,3,6)
[theta,idx,phiinv] = doe_plot(response_kss(:,6),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+100,int2str(idx(i)),'FontSize',10);
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),'Biomarker 6', 'FontSize',11,'FontWeight','bold')
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontSize',11,'FontWeight','bold')

%% consolidate half-normal plots
clc
clearvars
close all

load('response_kto.mat')
load('response_kslow1.mat')
load('response_kslow2.mat')
load('response_kss.mat')

dgn_kto = table2array(readtable('ikto_dgn.csv'));
dgn_kslow1 = table2array(readtable('ikslow1_dgn.csv'));
dgn_kslow2 = table2array(readtable('ikslow2_dgn.csv'));
dgn_kss = table2array(readtable('ikss_dgn.csv'));

c = parula(5);

fig = figure('Color','w','Position',[100,100,900,1200]);
orient(fig,'landscape')

% biomarker 1
% kto
subplot(3,2,1)
[theta,idx,phiinv] = doe_plot(response_kto(:,1),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(response_kslow1(:,1),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end

% kslow2
[theta,idx,phiinv] = doe_plot(response_kslow2(:,1),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.001,int2str(idx(i)),'FontSize',10);
end

% kss
[theta,idx,phiinv] = doe_plot(response_kslow2(:,1),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.2,int2str(idx(i)),'FontSize',10);
end
hold off
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"])
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 2
% kto
subplot(3,2,2)
[theta,idx,phiinv] = doe_plot(response_kto(:,2),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(response_kslow1(:,2),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end

% kslow2
[theta,idx,phiinv] = doe_plot(response_kslow2(:,2),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.001,int2str(idx(i)),'FontSize',10);
end

% kss
[theta,idx,phiinv] = doe_plot(response_kslow2(:,2),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.2,int2str(idx(i)),'FontSize',10);
end
hold off
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"])
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 3
% kto
subplot(3,2,3)
[theta,idx,phiinv] = doe_plot(response_kto(:,3),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(response_kslow1(:,3),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end

% kslow2
[theta,idx,phiinv] = doe_plot(response_kslow2(:,3),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.001,int2str(idx(i)),'FontSize',10);
end

% kss
[theta,idx,phiinv] = doe_plot(response_kslow2(:,3),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.2,int2str(idx(i)),'FontSize',10);
end
hold off
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"])
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 4
% kto
subplot(3,2,4)
[theta,idx,phiinv] = doe_plot(response_kto(:,4),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(response_kslow1(:,4),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end

% kslow2
[theta,idx,phiinv] = doe_plot(response_kslow2(:,4),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.001,int2str(idx(i)),'FontSize',10);
end

% kss
[theta,idx,phiinv] = doe_plot(response_kslow2(:,4),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.2,int2str(idx(i)),'FontSize',10);
end
hold off
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"])
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 5
% kto
subplot(3,2,5)
[theta,idx,phiinv] = doe_plot(response_kto(:,5),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(response_kslow1(:,5),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end

% kslow2
[theta,idx,phiinv] = doe_plot(response_kslow2(:,5),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.001,int2str(idx(i)),'FontSize',10);
end

% kss
[theta,idx,phiinv] = doe_plot(response_kslow2(:,5),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.2,int2str(idx(i)),'FontSize',10);
end
hold off
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"])
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 6
% kto
subplot(3,2,6)
[theta,idx,phiinv] = doe_plot(response_kto(:,6),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+2,int2str(idx(i)),'FontSize',10);
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(response_kslow1(:,6),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+1,int2str(idx(i)),'FontSize',10);
end

% kslow2
[theta,idx,phiinv] = doe_plot(response_kslow2(:,6),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.001,int2str(idx(i)),'FontSize',10);
end

% kss
[theta,idx,phiinv] = doe_plot(response_kslow2(:,6),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i)+0.2,int2str(idx(i)),'FontSize',10);
end
hold off
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"])
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

%% consolidate half-normal plots - normalized
clc
clearvars
close all

load('response_kto.mat')
load('response_kslow1.mat')
load('response_kslow2.mat')
load('response_kss.mat')

dgn_kto = table2array(readtable('ikto_dgn.csv'));
dgn_kslow1 = table2array(readtable('ikslow1_dgn.csv'));
dgn_kslow2 = table2array(readtable('ikslow2_dgn.csv'));
dgn_kss = table2array(readtable('ikss_dgn.csv'));

c = parula(5);

fig = figure('Color','w','Position',[100,100,900,1200]);
orient(fig,'landscape')

% biomarker 1
% kto
subplot(3,2,1)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,1),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,1),'range'),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kslow2
[theta,idx,phiinv] = doe_plot(normalize(response_kslow2(:,1),'range'),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,1),'range'),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end
hold off
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"],'Location','northwest')
title("Anchor Point 1")
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 2
% kto
subplot(3,2,2)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,2),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,2),'range'),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kslow2
[theta,idx,phiinv] = doe_plot(normalize(response_kslow2(:,2),'range'),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,2),'range'),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end
hold off
title("Anchor Point 2")
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 3
% kto
subplot(3,2,3)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,3),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,3),'range'),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kslow2
[theta,idx,phiinv] = doe_plot(normalize(response_kslow2(:,3),'range'),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,3),'range'),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end
hold off
title("Anchor Point 3")
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 4
% kto
subplot(3,2,4)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,4),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,4),'range'),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kslow2
[theta,idx,phiinv] = doe_plot(normalize(response_kslow2(:,4),'range'),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,4),'range'),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end
hold off
title("Anchor Point 4")
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 5
% kto
subplot(3,2,5)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,5),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,5),'range'),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kslow2
[theta,idx,phiinv] = doe_plot(normalize(response_kslow2(:,5),'range'),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,5),'range'),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end
hold off
title("Anchor Point 5")
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% biomarker 6
% kto
subplot(3,2,6)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,6),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,6),'range'),dgn_kslow1);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kslow2
[theta,idx,phiinv] = doe_plot(normalize(response_kslow2(:,6),'range'),dgn_kslow2);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,6),'range'),dgn_kss);
plot(phiinv,theta,'+','MarkerSize',10,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','bottom');
end
hold off
title("Anchor Point 6")
xlabel('Half-normal Quantiles')
ylabel('Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

%% custrom functions
function [theta,idx,phiinv] = doe_plot(y,dgn)
    stat_beta = regress(y,dgn);
    beta = 2*stat_beta;
    beta = abs(beta);

    [theta,idx] = sort(beta);
    p = 0.5 + (0.5*(1:length(beta))-0.5)/length(beta);
    phiinv = norminv(p,0,1);
end