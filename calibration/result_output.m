%% Model fitness - nonlinear R-squared measure
clc
close all
clearvars

exp_num = "exp51";
calib_dir = strcat("calib_",exp_num);

pa_wt = readtable(fullfile(calib_dir,"r2l2_measures_wt.csv"));
pa_ko = readtable(fullfile(calib_dir,"r2l2_measures_mgat1ko.csv"));

files_wt = strings(length(pa_wt.file_name),1);
files_ko = strings(length(pa_ko.file_name),1);

for i=1:length(pa_wt.file_name)
    [~,files_wt(i),~] = fileparts(pa_wt.file_name{i});
end

for i=1:length(pa_ko.file_name)
    [~,files_ko(i),~] = fileparts(pa_ko.file_name{i});
end

r2_wt = pa_wt.mean_R2;
l2_wt = pa_wt.mean_L2;
r2_ko = pa_ko.mean_R2;
l2_ko = pa_ko.mean_L2;

%----- Without the worst fitting point in MGAT1KO -----%
[~,rmv_idx1] = mink(r2_wt,2);
[~,rmv_idx2] = min(r2_ko);

files_wt(rmv_idx1) = [];
% files_wt = categorical(files_wt);
r2_wt(rmv_idx1) = [];
files_ko(rmv_idx2) = [];
% files_ko = categorical(files_ko);
r2_ko(rmv_idx2) = [];

fig = figure('Color','w','Position',[0,0,750,600]);
orient(fig,'landscape')
subplot(6,4,[1,5,9])
barh(categorical(1:length(files_wt)), r2_wt,'b');
xlim([0,1])
grid on
% title("WT")
% xlabel("Nonlinear R^{2}")
set(gca,'FontWeight','bold')

subplot(6,4,[13,17,21])
barh(categorical(1:length(files_ko)),r2_ko,'r')
xlim([0,1])
grid on
% title("MGAT1KO")
% xlabel("Nonlinear R^{2}")
set(gca,'FontWeight','bold')

%----- Representative example of showing actual fitness -----%
volts = -30:10:50;
[~,max_wt] = max(r2_wt);
exp_data = table2array(...
    readtable(fullfile("mgat1ko_data/wt-preprocessed",pa_wt.file_name{max_wt})));
yksum = table2array(...
    readtable(fullfile(calib_dir,"wt_yhat",pa_wt.file_name{max_wt})));
t = exp_data(:,1);
pa_row = pa_wt(string(pa_wt.file_name) == pa_wt.file_name{max_wt},:);

for i=1:size(yksum,2)
    yi = exp_data(:,i+1);
    yhati = yksum(:,i);
    r2 = pa_row{1,2+(i-1)*2};
    l2 = pa_row{1,3+(i-1)*2};

    if i==1
        subplot(6,4,1+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
%         xlabel("Time (ms)")
        ylabel(" ")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')
    end

    if ismember(i,2:3)
        subplot(6,4,1+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
%         xlabel("Time (ms)")
%         ylabel("Current (pA/pF)")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')
    end

    if ismember(i,4:6)
        subplot(6,4,2+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
%         xlabel("Time (ms)")
%         ylabel("Current (pA/pF)")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')        
    end

    if ismember(i,7:9)
        subplot(6,4,3+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
%         xlabel("Time (ms)")
%         ylabel("Current (pA/pF)")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')        

        if i==9
            legend(["Experimental I_{K}", "Model Prediction"])
            legend Box off
        end
    end
end

[~,max_ko] = max(r2_ko);
exp_data = table2array(...
    readtable(fullfile("mgat1ko_data/mgat1ko-preprocessed",pa_ko.file_name{max_ko})));
yksum = table2array(...
    readtable(fullfile(calib_dir,"mgat1ko_yhat",pa_ko.file_name{max_ko})));
t = exp_data(:,1);
pa_row = pa_ko(string(pa_ko.file_name) == pa_ko.file_name{max_ko},:);

for i=1:size(yksum,2)
    yi = exp_data(:,i+1);
    yhati = yksum(:,i);
    r2 = pa_row{1,2+(i-1)*2};
    l2 = pa_row{1,3+(i-1)*2};

    if ismember(i,1:3)
        subplot(6,4,13+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
%         xlabel("Time (ms)")
%         ylabel("Current (pA/pF)")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')
    end

    if ismember(i,4:6)
        subplot(6,4,14+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
%         xlabel("Time (ms)")
%         ylabel("Current (pA/pF)")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')        
    end

    if ismember(i,7:9)
        subplot(6,4,15+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
%         xlabel("Time (ms)")
%         ylabel("Current (pA/pF)")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')        
    end
end

%% Model fitness - RMSE
%----- protocol -----%
clc
clearvars
close all
load("file_names_4half4.mat")

exp_nums = ["exp45","exp46","exp47"];
current_names = {'iktof','ikslow1','ikslow2','ikss'};
tune_idx = cell(5,1);
tune_idx{1} = [1, 2, 4, 5, 7, 11, 13];
tune_idx{2} = [2, 3];
tune_idx{3} = [1, 2, 4, 5, 9, 10 11];
tune_idx{4} = [2, 3];
tune_idx{5} = [1, 2, 3, 4];
[mdl_struct,psize] = gen_mdl_struct(current_names,tune_idx);

pdefault = cell(5,1);
pdefault{1} = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.3846];
pdefault{2} = [-1050, 270, 0.0629];
pdefault{3} = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, 6.1, 18, 2.058, 803.0, 0.16];
pdefault{4} = [4912, 5334, 0.16];
pdefault{5} = [0.0862, 1235.5, 13.17, 0.0611];

volts = -100:10:50;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;

protocol = cell(6,1);
protocol{1} = -70; % hold_volt
protocol{2} = -91.1; % ek

%----- RMSE bar graph and current fitting -----%
% plot 1
fig = figure('Color','w','Position',[50,50,900,350]);
orient(fig,'landscape')
subplot(2,4,[1,5])
sub_rmse_df = rmse_df_wt(~rmv_idx_wt,:);
rmse = table2array(sub_rmse_df(:,2:end));
b = barh(categorical(1:length(sub_rmse_df.File)),rmse);
c = parula(3);
for i=1:3
    b(i).FaceColor = c(i,:);
end
grid on
title("WT")
xlabel("\Sigma RMSE")
ylabel("Data Index")
set(gca,'FontWeight','bold')

% plot 2
c = parula(3);
fidx = 26;
dpath = fullfile(pwd,"mgat1ko_data","wt-preprocessed",strcat(files_wt(fidx),".xlsx"));
exp_data = table2array(readtable(dpath));
t = exp_data(:,1);
yksum = exp_data(:,2:end);

[~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
[~, ideal_end_idx] = min(abs(t - ideal_end_time));
t = t(1:ideal_end_idx);
yksum = yksum(1:ideal_end_idx,:);

protocol{4} = t;
protocol{5} = t(1:ideal_hold_idx);
pulse_t = t(ideal_hold_idx+1:end);
protocol{6} = pulse_t - pulse_t(1);

subplot(2,4,2)
plot(t,yksum(:,3),'-','Color','red','LineWidth',1)
hold on
for i=1:length(exp_nums)
    fpath = fullfile(pwd,strcat("calib_",exp_nums(i)),"wt",strcat(files_wt(fidx),".xlsx"));
    sol_mx = table2array(readtable(fpath));
    sol = gen_sol_vec(sol_mx, mdl_struct, psize);
    protocol{3} = volts(3);
    [ymx,~] = kv_model(sol,mdl_struct,pdefault,protocol);
    yhat = sum(ymx,2);
    plot(t,yhat,'-','Color',c(i,:),'LineWidth',1.5)
 end
hold off
ylim([0,1.3])
grid on
title(strcat("Data Index ",string(fidx)," at ", string(volts(3))," mV"))
xlabel("Time (ms)")
ylabel("Current (pA/pF)")
set(gca,'FontWeight','bold')

% plot 3
subplot(2,4,6)
plot(t,yksum(:,end),'-','Color','red','LineWidth',2)
hold on
for i=1:length(exp_nums)
    fpath = fullfile(pwd,strcat("calib_",exp_nums(i)),"wt",strcat(files_wt(fidx),".xlsx"));
    sol_mx = table2array(readtable(fpath));
    sol = gen_sol_vec(sol_mx, mdl_struct, psize);
    protocol{3} = volts(end);
    [ymx,~] = kv_model(sol,mdl_struct,pdefault,protocol);
    yhat = sum(ymx,2);
    plot(t,yhat,'-','Color',c(i,:),'LineWidth',1.5)
end
hold off
ylim([0,25])
grid on
title(strcat("Data Index ",string(fidx)," at ", string(volts(end))," mV"))
% legend(["","BFGS","SQP","Active Set"])
xlabel("Time (ms)")
ylabel("Current (pA/pF)")
set(gca,'FontWeight','bold')

% plot 4
subplot(2,4,[3,7])
sub_rmse_df = rmse_df_ko(~rmv_idx_ko,:);
rmse = table2array(sub_rmse_df(:,2:end));
b = barh(categorical(1:length(sub_rmse_df.File)),rmse);
for i=1:3
    b(i).FaceColor = c(i,:);
end
xlim([0,6])
grid on
title("MGAT1KO")
xlabel("\Sigma RMSE")
ylabel("Data Index")
legend(["BFGS","SQP","Active Set"],'Location','southeast')
set(gca,'FontWeight','bold')

% plot 5
fidx = 4;
dpath = fullfile(pwd,"mgat1ko_data","ko-preprocessed",strcat(files_ko(fidx),".xlsx"));
exp_data = table2array(readtable(dpath));
t = exp_data(:,1);
yksum = exp_data(:,2:end);

[~, ideal_hold_idx] = min(abs(t - ideal_hold_time));
[~, ideal_end_idx] = min(abs(t - ideal_end_time));
t = t(1:ideal_end_idx);
yksum = yksum(1:ideal_end_idx,:);

protocol{4} = t;
protocol{5} = t(1:ideal_hold_idx);
pulse_t = t(ideal_hold_idx+1:end);
protocol{6} = pulse_t - pulse_t(1);

subplot(2,4,4)
plot(t,yksum(:,3),'Color','red','LineWidth',1)
hold on
for i=1:length(exp_nums)
    fpath = fullfile(pwd,strcat("calib_",exp_nums(i)),"mgat1ko",strcat(files_ko(fidx),".xlsx"));
    sol_mx = table2array(readtable(fpath));
    sol = gen_sol_vec(sol_mx, mdl_struct, psize);
    protocol{3} = volts(3);
    [ymx,~] = kv_model(sol,mdl_struct,pdefault,protocol);
    yhat = sum(ymx,2);
    plot(t,yhat,'-','Color',c(i,:),'LineWidth',1.5)    
end
hold off
ylim([0,1.3])
grid on
title(strcat("Data Index ",string(fidx)," at ", string(volts(3))," mV"))
xlabel("Time (ms)")
ylabel("Current (pA/pF)")
set(gca,'FontWeight','bold')

% plot 6
subplot(2,4,8)
plot(t,yksum(:,end),'Color','red','LineWidth',2)
hold on
for i=1:length(exp_nums)
    fpath = fullfile(pwd,strcat("calib_",exp_nums(i)),"mgat1ko",strcat(files_ko(fidx),".xlsx"));
    sol_mx = table2array(readtable(fpath));
    sol = gen_sol_vec(sol_mx, mdl_struct, psize);
    protocol{3} = volts(end);
    [ymx,~] = kv_model(sol,mdl_struct,pdefault,protocol);
    yhat = sum(ymx,2);
    plot(t,yhat,'-','Color',c(i,:),'LineWidth',1.5)    
end
hold off
ylim([0,25])
grid on
title(strcat("Data Index ",string(fidx)," at ", string(volts(end))," mV"))
legend(["Experimental","BFGS","SQP","Active Set"])
xlabel("Time (ms)")
ylabel("Current (pA/pF)")
set(gca,'FontWeight','bold')

%% Kinetics modeling - Estimation
clc
close all

%----- protocol -----%
load("file_names_4half4.mat")

exp_num = "exp51";
current_names = {'iktof','ikslow1','ikslow2','ikss'};
tune_idx = cell(5,1);
tune_idx{1} = [1, 2, 3, 4, 5, 7, 11, 13];
tune_idx{2} = [2, 3];
tune_idx{3} = [1, 2, 4, 5, 9, 10 11];
tune_idx{4} = [2, 3];
tune_idx{5} = [1, 2, 3, 4];
[mdl_struct,psize] = gen_mdl_struct(current_names,tune_idx);

pdefault = cell(5,1);
pdefault{1} = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.3846];
pdefault{2} = [-1050, 270, 0.0629];
pdefault{3} = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, 6.1, 18, 2.058, 803.0, 0.16];
pdefault{4} = [4912, 5334, 0.16];
pdefault{5} = [0.0862, 1235.5, 13.17, 0.0611];

volts = -100:10:70;
ideal_hold_time = 120;
ideal_end_time = 4.6*1000;

protocol = cell(6,1);
protocol{1} = -70; % hold_volt
protocol{2} = -91.1; % ek
sol_dir = fullfile(pwd,strcat("calib_",exp_num));

atof_wt = NaN(length(files_wt),length(volts));
itof_wt = NaN(length(files_wt),length(volts));
taua_tof_wt = NaN(length(files_wt),length(volts));
taui_tof_wt = NaN(length(files_wt),length(volts));
akslow1_wt = NaN(length(files_wt),length(volts));
ikslow1_wt = NaN(length(files_wt),length(volts));
taua_kslow1_wt = NaN(length(files_wt),length(volts));
taui_kslow1_wt = NaN(length(files_wt),length(volts));
taui_kslow2_wt = NaN(length(files_wt),length(volts));
taua_kss_wt = NaN(length(files_wt),length(volts));

atof_ko = NaN(length(files_ko),length(volts));
itof_ko = NaN(length(files_ko),length(volts));
taua_tof_ko = NaN(length(files_ko),length(volts));
taui_tof_ko = NaN(length(files_ko),length(volts));
akslow1_ko = NaN(length(files_ko),length(volts));
ikslow1_ko = NaN(length(files_ko),length(volts));
taua_kslow1_ko = NaN(length(files_ko),length(volts));
taui_kslow1_ko = NaN(length(files_ko),length(volts));
taui_kslow2_ko = NaN(length(files_ko),length(volts));
taua_kss_ko = NaN(length(files_ko),length(volts));
for i=1:length(files_wt)
    sol_mx = table2array(readtable(fullfile(sol_dir,"wt",strcat(files_wt(i),".xlsx"))));
    sol = gen_sol_vec(sol_mx, mdl_struct, psize);
    
    for j=1:length(volts)
        protocol{3} = volts(j);
        [~,kv] = kv_model(sol,mdl_struct,pdefault,protocol);
        atof_wt(i,j) = kv{1}(1);
        itof_wt(i,j) = kv{1}(2);
        taua_tof_wt(i,j) = kv{1}(3);
        taui_tof_wt(i,j) = kv{1}(4);
        akslow1_wt(i,j) = kv{2}(1);
        ikslow1_wt(i,j) = kv{2}(2);
        taua_kslow1_wt(i,j) = kv{2}(3);
        taui_kslow1_wt(i,j) = kv{2}(4);
        taui_kslow2_wt(i,j) = kv{3}(4);
        taua_kss_wt(i,j) = kv{4}(2);        
    end
end

for i=1:length(files_ko)
    sol_mx = table2array(readtable(fullfile(sol_dir,"mgat1ko",strcat(files_ko(i),".xlsx"))));
    sol = gen_sol_vec(sol_mx, mdl_struct, psize);

    for j=1:length(volts)
        protocol{3} = volts(j);
        [~,kv] = kv_model(sol,mdl_struct,pdefault,protocol);
        atof_ko(i,j) = kv{1}(1);
        itof_ko(i,j) = kv{1}(2);
        taua_tof_ko(i,j) = kv{1}(3);
        taui_tof_ko(i,j) = kv{1}(4);
        akslow1_ko(i,j) = kv{2}(1);
        ikslow1_ko(i,j) = kv{2}(2);
        taua_kslow1_ko(i,j) = kv{2}(3);
        taui_kslow1_ko(i,j) = kv{2}(4);
        taui_kslow2_ko(i,j) = kv{3}(4);
        taua_kss_ko(i,j) = kv{4}(2);        
    end
end

%% Kinetics modeling visualization - IKto
clc
close all

figure('Color','w','Position',[50,50,600,500]);
t = tiledlayout(2,2);
xlabel(t,"Voltage (mV)",'FontWeight','bold')

% plot1
nexttile
r = 3:16;
atof_sem = std(atof_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(atof_wt(:,r),1),atof_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
atof_sem = std(atof_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(atof_ko(:,r),1),atof_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
set(gca,'XLimSpec','tight')
ylim([0,1])
xticks(volts(r))
ylabel("a_{ss}^{(1)}")
legend(["WT","MGAT1KO"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)

% plot2
nexttile
r = 3:12;
itof_sem = std(itof_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(itof_wt(:,r),1),itof_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
itof_sem = std(itof_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(itof_ko(:,r),1),itof_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
set(gca,'XLimSpec','tight')
ylim([0,1])
xticks(volts(r))
ylabel("i_{ss}^{(1)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

% plot3
nexttile
r = 3:16;
taua_tof_sem = std(taua_tof_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taua_tof_wt(:,r),1),taua_tof_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taua_tof_sem = std(taua_tof_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taua_tof_ko(:,r),1),taua_tof_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
ylabel("\tau_{a}^{(1)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

% plot4
nexttile
r = 3:12;
taui_tof_sem = std(taui_tof_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taui_tof_wt(:,r),1),taui_tof_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taui_tof_sem = std(taui_tof_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taui_tof_ko(:,r),1),taui_tof_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
ylabel("\tau_{i}^{(1)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

%% Kinetics modeling visualization - Rectifier currents (2x3)
clc
close all

fig = figure('Color','w','Position',[50,50,900,500]);
orient(fig,'landscape')

subplot(2,3,1)
r = 3:16;
akslow1_sem = std(akslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(akslow1_wt(:,r),1),akslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
akslow1_sem = std(akslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(akslow1_ko(:,r),1),akslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
ylim([0,1])
set(gca,'XLimSpec','tight')
xticks(volts(r))
xlabel("Voltage (mV)")
ylabel("a_{ss}^{(2)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

subplot(2,3,2)
r = 4:11;
ikslow1_sem = std(ikslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(ikslow1_wt(:,r),1),ikslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
ikslow1_sem = std(ikslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(ikslow1_ko(:,r),1),ikslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
ylim([0,1])
set(gca,'XLimSpec','tight')
xticks(volts(r))
xlabel("Voltage (mV)")
ylabel("i_{ss}^{(2)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

subplot(2,3,3)
r = 3:16;
taua_kslow1_sem = std(taua_kslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taua_kslow1_wt(:,r),1),taua_kslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taua_kslow1_sem = std(taua_kslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taua_kslow1_ko(:,r),1),taua_kslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
xlabel("Voltage (mV)")
ylabel("\tau_{a}^{(2)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

subplot(2,3,4)
r = 3:16;
taui_kslow1_sem = std(taui_kslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taui_kslow1_wt(:,r),1),taui_kslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taui_kslow1_sem = std(taui_kslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taui_kslow1_ko(:,r),1),taui_kslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
xlabel("Voltage (mV)")
ylabel("\tau_{i}^{(2)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

subplot(2,3,5)
taui_kslow2_sem = std(taui_kslow2_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taui_kslow2_wt(:,r),1),taui_kslow2_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taui_kslow2_sem = std(taui_kslow2_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taui_kslow2_ko(:,r),1),taui_kslow2_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
xlabel("Voltage (mV)")
ylabel("\tau_{i}^{(3)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

subplot(2,3,6)
taua_kss_sem = std(taua_kss_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taua_kss_wt(:,r),1),taua_kss_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taua_kss_sem = std(taua_kss_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taua_kss_ko(:,r),1),taua_kss_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
xlabel("Voltage (mV)")
ylabel("\tau_{a}^{(4)}")
set(gca,'FontWeight','bold','LineWidth',1.5)
legend(["WT","MGAT1KO"],'Location','northeast')

%% Kinetics modeling visualization - Rectifier currents (3x2)
clc
close all

figure('Color','w','Position',[50,50,600,750]);
t = tiledlayout(3,2);
xlabel(t,"Voltage (mV)",'FontWeight','bold')

nexttile
r = 3:16;
akslow1_sem = std(akslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(akslow1_wt(:,r),1),akslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
akslow1_sem = std(akslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(akslow1_ko(:,r),1),akslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
ylim([0,1])
set(gca,'XLimSpec','tight')
xticks(volts(r))
ylabel("a_{ss}^{(2)}")
legend(["WT","MGAT1KO"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)

nexttile
r = 4:11;
ikslow1_sem = std(ikslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(ikslow1_wt(:,r),1),ikslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
ikslow1_sem = std(ikslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(ikslow1_ko(:,r),1),ikslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
ylim([0,1])
set(gca,'XLimSpec','tight')
xticks(volts(r))
ylabel("i_{ss}^{(2)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

nexttile
r = 3:16;
taua_kslow1_sem = std(taua_kslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taua_kslow1_wt(:,r),1),taua_kslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taua_kslow1_sem = std(taua_kslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taua_kslow1_ko(:,r),1),taua_kslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
ylabel("\tau_{a}^{(2)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

nexttile
r = 3:16;
taui_kslow1_sem = std(taui_kslow1_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taui_kslow1_wt(:,r),1),taui_kslow1_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taui_kslow1_sem = std(taui_kslow1_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taui_kslow1_ko(:,r),1),taui_kslow1_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
ylabel("\tau_{i}^{(2)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

nexttile
taui_kslow2_sem = std(taui_kslow2_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taui_kslow2_wt(:,r),1),taui_kslow2_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taui_kslow2_sem = std(taui_kslow2_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taui_kslow2_ko(:,r),1),taui_kslow2_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
ylabel("\tau_{i}^{(3)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

nexttile
taua_kss_sem = std(taua_kss_wt(:,r),0,1)/sqrt(length(files_wt));
errorbar(volts(r),mean(taua_kss_wt(:,r),1),taua_kss_sem,'-o','Color','blue','MarkerFaceColor','blue','LineWidth',1.5)
hold on
taua_kss_sem = std(taua_kss_ko(:,r),0,1)/sqrt(length(files_ko));
errorbar(volts(r),mean(taua_kss_ko(:,r),1),taua_kss_sem,'--s','Color','red','MarkerFaceColor','red','LineWidth',1.5)
hold off
grid on
axis tight
xticks(volts(r))
ylabel("\tau_{a}^{(4)}")
set(gca,'FontWeight','bold','LineWidth',1.5)

%% Parameter distribution
clc
clearvars
close all
load('r2_file_names.mat')

exp_num = "exp51";
base_dir = fullfile(pwd,strcat("calib_",exp_num));

% pktof
pktof = readtable(fullfile(base_dir,"pktof.csv"));
pidx = unique(pktof.param);

fig = figure('Color','w','Position',[50,50,830,320]);
orient(fig,'landscape')
tl = tiledlayout(2,4);
ylabel(tl,"Frequency",'FontWeight','bold')
for i=1:length(pidx)
    psub = pktof(pktof.param==pidx(i),:);
    psub_wt = psub(string(psub.Group)=="WT",:);
    psub_ko = psub(string(psub.Group)=="Mgat1KO",:);

%     [f1,xi1] = ksdensity(psub_wt.value);
%     [f2,xi2] = ksdensity(psub_ko.value);
    
%     subplot(2,4,i)
%     plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    nexttile
    h1 = histogram(psub_wt.value,15,'FaceColor','b','FaceAlpha',0.5);
    bin_edges1 = h1.BinEdges;
    hold on
%     plot(xi2,f2, 'Color','red', 'LineWidth',2)
    h2 = histogram(psub_ko.value,'BinEdges',bin_edges1,'FaceColor','r','FaceAlpha',0.5);
    hold off
    axis tight
    grid on
    xlabel(strcat("p_{", num2str(pidx(i)), "} Range"))
%     ylabel('Density')
%     ylabel("Frequency")
    set(gca, 'FontName','Arial','FontWeight','bold','LineWidth',1.5)

    if i==length(pidx)
        legend(["WT","Mgat1KO"])
        legend box off
        xlabel("G_{Kto} Range")
    end
end
nrow = size(pktof,1)/length(pidx);
pktof2 = array2table(NaN(nrow,length(pidx)+2));

pktof2.Properties.VariableNames(1:2) = {'Group','File'};
pktof2.Group = pktof(pktof.param==pidx(1),:).Group;
pktof2.File = pktof(pktof.param==pidx(1),:).file;

for i=1:length(pidx)
    pktof2.Properties.VariableNames(i+2) = {strcat('P',num2str(pidx(i)))};
    pktof2(:,i+2)= pktof(pktof.param==pidx(i),'value');
end

% pkslow1
pkslow1 = readtable(fullfile(base_dir,"pkslow1.csv"));
pidx = unique(pkslow1.param);

fig = figure('Color','w','Position',[100,100,830,640]);
orient(fig,'landscape')
for i=1:length(pidx)
    psub = pkslow1(pkslow1.param==pidx(i),:);
    psub_wt = psub(string(psub.Group)=="WT",:);
    psub_ko = psub(string(psub.Group)=="Mgat1KO",:);

%     [f1,xi1] = ksdensity(psub_wt.value);
%     [f2,xi2] = ksdensity(psub_ko.value);
    
    subplot(4,4,i)
%     plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    h1 = histogram(psub_wt.value,15,'FaceColor','b','FaceAlpha',0.5);
    bin_edges1 = h1.BinEdges;
    hold on
%     plot(xi2,f2, 'Color','red', 'LineWidth',2)
    histogram(psub_ko.value,'BinEdges',bin_edges1,'FaceColor','r','FaceAlpha',0.5)
    hold off
    axis tight
    grid on
    xlabel(strcat("p_{", num2str(pidx(i)), "} Range"))
%     ylabel('Density')
%     ylabel("Frequency")
    set(gca, 'FontName','Arial','FontWeight','bold','LineWidth',1.5)
    
    if i==length(pidx)
        xlabel("G_{Kslow1} Range")
        legend('WT','Mgat1KO', 'Location','best')
        legend box off
    end
end
nrow = size(pkslow1,1)/length(pidx);
pkslow12 = array2table(NaN(nrow,length(pidx)+1));
pkslow12.Properties.VariableNames(1) = {'Group'};
pkslow12.Group = pkslow1(pkslow1.param==pidx(1),:).Group;

for i=1:length(pidx)
    pkslow12.Properties.VariableNames(i+1) = {strcat('P',num2str(pidx(i)))};
    pkslow12(:,i+1)= pkslow1(pkslow1.param==pidx(i),'value');
end

% pkslow2
% figure('Color','w','Position',[150,150,610,320]);
pkslow2 = readtable(fullfile(base_dir,"pkslow2.csv"));
pidx = unique(pkslow2.param);

for i=1:length(pidx)
    psub = pkslow2(pkslow2.param==pidx(i),:);
    psub_wt = psub(string(psub.Group)=="WT",:);
    psub_ko = psub(string(psub.Group)=="Mgat1KO",:);

%     [f1,xi1] = ksdensity(psub_wt.value);
%     [f2,xi2] = ksdensity(psub_ko.value);
    

    subplot(4,4,i+8)
%     plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    h1 = histogram(psub_wt.value,15,'FaceColor','b','FaceAlpha',0.5);
    bin_edges1 = h1.BinEdges;
    hold on
%     plot(xi2,f2, 'Color','red', 'LineWidth',2)
    histogram(psub_ko.value,'BinEdges',bin_edges1,'FaceColor','r','FaceAlpha',0.5)
    hold off
    axis tight
    grid on
    xlabel(strcat("p_{", num2str(pidx(i)), "} Range"))
%     ylabel('Density')
%     ylabel("Frequency")
    set(gca, 'FontName','Arial','FontWeight','bold','LineWidth',1.5)
    
    if i==length(pidx)
        xlabel("G_{kslow2}")
    end
end
nrow = size(pkslow2,1)/length(pidx);
g = categorical(pkslow2(pkslow2.param==pidx(1),:).Group);
v1 = pkslow2(pkslow2.param==pidx(1),:).value;
v2 = pkslow2(pkslow2.param==pidx(2),:).value;
pkslow22 = table(g,v1,v2);


% pkss
pkss = readtable(fullfile(base_dir,"pkss.csv"));
pidx = unique(pkss.param);
for i=1:length(pidx)
    psub = pkss(pkss.param==pidx(i),:);
    psub_wt = psub(string(psub.Group)=="WT",:);
    psub_ko = psub(string(psub.Group)=="Mgat1KO",:);

%     [f1,xi1] = ksdensity(psub_wt.value);
%     [f2,xi2] = ksdensity(psub_ko.value);
    
    subplot(4,4,i+12)
%     plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    h1 = histogram(psub_wt.value,15,'FaceColor','b','FaceAlpha',0.5);
    bin_edges1 = h1.BinEdges;
    hold on
%     plot(xi2,f2, 'Color','red', 'LineWidth',2)
    histogram(psub_ko.value,'BinEdges',bin_edges1,'FaceColor','r','FaceAlpha',0.5)
    hold off
    axis tight
    grid on
    xlabel(strcat("p_{", num2str(pidx(i)), "} Range"))
%     ylabel('Density')
%     ylabel("Frequency")
    set(gca, 'FontName','Arial','FontWeight','bold','LineWidth',1.5)

    if i==length(pidx)
        xlabel("G_{Kss} Range")
    end
end
pkss2 = array2table(NaN(nrow,length(pidx)+1));
pkss2.Properties.VariableNames(1) = {'Group'};
pkss2.Group = pkss(pkss.param==pidx(1),:).Group;

for i=1:length(pidx)
    pkss2.Properties.VariableNames(i+1) = {strcat('P',num2str(pidx(i)))};
    pkss2(:,i+1)= pkss(pkss.param==pidx(i),'value');
end

%% Low-dimensional embedding
clc
close all

% g = NaN(size(pktof2,1),1);
% g(pktof2.Group=="Mgat1KO") = 1;
% g(pktof2.Group=="WT") = 2;
% c = brewermap(length(unique(g)),'Set1');
c = NaN(size(pktof2,1),3);
idx1 = pktof2.Group=="WT";
idx2 = pktof2.Group=="Mgat1KO";

c(idx1,:) = repmat([0 0 1],sum(idx1),1);
c(idx2,:) = repmat([1 0 0],sum(idx2),1);

pmx = [table2array(pktof2(:,3:end)),table2array(pkslow12(:,2:end)),table2array(pkslow22(:,2:end)),table2array(pkss2(:,2:end))];

rng(2118)
embd = tsne(normalize(pmx),'Algorithm','exact','Distance','minkowski','NumDimensions',3);

figure('Color','w','Position',[100,100,560,510])
s1 = scatter3(embd(idx1,1),embd(idx1,2),embd(idx1,3),'o','CData',c(pktof2.Group=="WT",:));
s1.SizeData = 90;
s1.LineWidth = 1.5;
hold on
s2 = scatter3(embd(idx2,1),embd(idx2,2),embd(idx2,3),'^','CData',c(pktof2.Group=="Mgat1KO",:));
s2.SizeData = 90;
s2.LineWidth = 1.5;
hold off
xlabel("Embedding Dim1")
ylabel("Embedding Dim2")
zlabel("Embedding Dim3")
axis tight
grid off
legend(["WT","MGAT1KO"],'location','best')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',1.5)


%% custom functions
function [mdl_struct, psize] = gen_mdl_struct(current_names, tune_idx)
    num_currents = length(current_names);
    
    % index 1
    idx_info1 = cell(1,num_currents);
    for i = 1:num_currents
        switch current_names{i}
            case 'iktof'
                idx_info1{i} = tune_idx{1};
            case 'iktos'
                idx_info1{i} = tune_idx{2};
            case 'ikslow1'
                idx_info1{i} = tune_idx{3};
            case 'ikslow2'
                idx_info1{i} = tune_idx{4};
            case 'ikss'
                idx_info1{i} = tune_idx{5};
        end
    end

    % index 2
    idx_info2 = cell(1,num_currents);
    psize = 0;
    for i = 1: num_currents
        psize_old = psize;
        psize = psize + length(idx_info1{i});
        idx_info2{i} = (1+psize_old):(psize);
    end

    model_info = [current_names; idx_info1; idx_info2];
    field_names = ["name","idx1","idx2"];
    mdl_struct = cell2struct(model_info, field_names, 1);    
end

function sol = gen_sol_vec(sol_mx, model_struct, psize)
    sol = zeros(1, psize);
    for i = 1:length(model_struct)
        running_sol = sol_mx(:, i);
        running_sol = running_sol(~isnan(running_sol));
        sol(model_struct(i).idx2) = running_sol(model_struct(i).idx1);
    end
end

function [ymx, kv] = kv_model(p,mdl_struct,pdefault,protocol)
    current_names = {mdl_struct.name};
    idx_info1 = {mdl_struct.idx1};
    idx_info2 = {mdl_struct.idx2};
    
    % define param_kslow1 first just in case it would be used in advance
    matching_idx = strcmp(current_names, 'ikslow1');
    param_kslow1 = pdefault{3};
    if any(matching_idx)
        param_kslow1(idx_info1{matching_idx}) = p(idx_info2{matching_idx});
    end

    ymx = zeros(length(protocol{4}),length(current_names));
    kv = cell(length(current_names),1);
    for i = 1:length(current_names)
        switch current_names{i}
            case 'iktof'
                param = pdefault{1};
                param(idx_info1{i}) = p(idx_info2{i});
                [ymx(:,i),kv{i}] = iktof(param, protocol);
            case 'iktos'
                num_param = 11;
                shared_idx = [1:7,9];
                uniq_idx = setdiff(1:num_param, shared_idx);

                param = NaN(num_param,1);
                param(shared_idx) = param_kslow1(shared_idx);
                uniq_param = pdefault{2};
                uniq_param(idx_info1{i}) = p(idx_info2{i});
                param(uniq_idx) = uniq_param;
                [ymx(:,i),kv{i}] = iktos(param, protocol);
            case 'ikslow1'
                [ymx(:,i),kv{i}] = ikslow1(param_kslow1, protocol);
            case 'ikslow2'
                num_param = 11;
                shared_idx = [1:7,9];
                uniq_idx = setdiff(1:num_param, shared_idx);

                param = NaN(num_param,1);
                param(shared_idx) = param_kslow1(shared_idx);
                uniq_param = pdefault{4};
                uniq_param(idx_info1{i}) = p(idx_info2{i});
                param(uniq_idx) = uniq_param;
                [ymx(:,i),kv{i}] = ikslow2(param, protocol);
            case 'ikss'
                num_param = 7;
                shared_idx1 = 1:3;
                shared_idx2 = [1,3,4];
                uniq_idx = setdiff(1:num_param, shared_idx1);

                param = NaN(num_param,1);
                param(shared_idx1) = param_kslow1(shared_idx2);
                uniq_param = pdefault{5};
                uniq_param(idx_info1{i}) = p(idx_info2{i});
                param(uniq_idx) = uniq_param;
                [ymx(:,i),kv{i}] = ikss(param, protocol);
        end
    end    
end

function [y,kv_pulse] = iktof(p, protocol)
    gmax = p(13);
    act0 = 0.4139033547e-02;
    inact0 = 0.9999623535e+00;
    
    holdv = protocol{1};
    ek = protocol{2};
    volt = protocol{3};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding volt
    kv_hold = iktof_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));
    y(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(holdv-ek);

    % current equation at pulse volt
    kv_pulse = iktof_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt-ek);    
end

function [y,kv_pulse] = iktos(p,protocol)
    gmax = p(11);
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689e+00;

    holdv = protocol{1};
    ek = protocol{2};
    volt = protocol{3};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));
    y(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(holdv - ek);
    
    % current equation at pulse voltage
    kv_pulse = rectifier_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - ek);    
end

function [y,kv_pulse] = ikslow1(p,protocol)
    gmax = p(11);
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;
    
    holdv = protocol{1};
    ek = protocol{2};
    volt = protocol{3};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};

    y = zeros(length(t),1);

    % current equation at holding 
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));        
    y(1:hold_idx) = (gmax).*(act_hold).*(inact_hold).*(holdv-ek);

    % current equation at pulse voltage
    kv_pulse = rectifier_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = (gmax).*(act_pulse).*(inact_pulse).*(volt-ek);
end

function [y,kv_pulse] = ikslow2(p,protocol)
    gmax = p(11); 
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;

    holdv = protocol{1};
    ek = protocol{2};
    volt = protocol{3};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(3));
    inact_hold = hh_model(holdt, inact0, kv_hold(2), kv_hold(4));
    y(1:hold_idx) = gmax.*(act_hold).*(inact_hold).*(holdv - ek);

    % current equation at pulse voltage
    kv_pulse = rectifier_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(2), kv_pulse(4));
    y((hold_idx+1):end) = gmax.*(act_pulse).*(inact_pulse).*(volt - ek);
end

function [y,kv_pulse] = ikss(p,protocol)
    gmax = p(7);
    act0 = 0.5091689794e-03;

    holdv = protocol{1};
    ek = protocol{2};
    volt = protocol{3};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(length(t),1);

    % current equation at holding
    kv_hold = ikss_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(1), kv_hold(2));
    y(1:hold_idx) = gmax.*act_hold.*(holdv - ek);

    % current equation at pulse voltage
    kv_pulse = ikss_kinetic_variables(p, volt);
    act_pulse = hh_model(pulset, act0, kv_pulse(1), kv_pulse(2));
    y((hold_idx+1):end) = gmax.*act_pulse.*(volt - ek);
end

function kv = iktof_kinetic_variables(p, V)
    kv = NaN(8,1);
    
    alpha1 = p(7).*exp(p(5).*(V+p(1)));
    beta1 = p(8).*exp(-p(6).*(V+p(1)));
    
    alpha2_temp1 = p(9).*exp((V+p(2))./(-1.0.*p(4)));
    alpha2_temp2 = p(10).*exp((V+p(2)+p(3))./(-1.0.*p(4)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(11).*exp((V+p(2)+p(3))./p(4));
    beta2_temp2 = p(12).*exp((V+p(2)+p(3))./p(4));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    kv(1) = alpha1./(alpha1+beta1);
    kv(2) = alpha2./(alpha2+beta2);
    kv(3) = 1./(alpha1+beta1);
    kv(4) = 1./(alpha2+beta2);
    kv(5) = alpha1;
    kv(6) = beta1;
    kv(7) = alpha2;
    kv(8) = beta2;
end

function kv = rectifier_kinetic_variables(p, V)
    kv = NaN(4,1);
    kv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(4))); % ass
    kv(2) = 1.0./(1.0+exp((p(2)+V)./p(5))); % iss
    kv(3) = p(7)./(exp(p(6)*(V+p(3))) + exp(-p(6)*(V+p(3))))+p(9); % taua
    kv(4) = p(10) - p(8)./(1.0+exp((p(2)+V)./p(5))); % taui
end

function kv = ikss_kinetic_variables(p, V)
    kv = NaN(2,1);
    kv(1) = 1.0./(1.0+exp(-(p(1)+V)./p(3)));
    kv(2) = p(5)./(exp(p(4)*(V+p(2))) + exp(-p(4)*(V+p(2)))) + p(6);
end

function y = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end
