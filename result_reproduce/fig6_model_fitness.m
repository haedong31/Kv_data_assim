%% Figure 6 Model fitness
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
r2_ko = pa_ko.mean_R2;

[~,rmv_idx1] = mink(r2_wt,2);
[~,rmv_idx2] = min(r2_ko);

files_wt(rmv_idx1) = [];
% files_wt = categorical(files_wt);
r2_wt(rmv_idx1) = [];
files_ko(rmv_idx2) = [];
% files_ko = categorical(files_ko);
r2_ko(rmv_idx2) = [];

%----- Bar graphs -----%
figure('Color','w','Position',[0,0,600,750]);
subplot(6,3,[1,4,7])
barh(categorical(1:length(files_wt)), r2_wt,'b');
xlim([0,1])
grid on
title("WT")
% xlabel("Nonlinear R^{2}")
ylabel("Cell Index")
set(gca,'FontWeight','bold')

subplot(6,3,[10,13,16])
barh(categorical(1:length(files_ko)),r2_ko,'r')
xlim([0,1])
grid on
title("MGAT1KO")
xlabel("Nonlinear R^{2}")
ylabel("Cell Index")
set(gca,'FontWeight','bold')

%----- Representative example of showing actual fitness -----%
volts = [-30,-10,0,20,30,50];
vidx = [1,3,4,6,7,9];
[~,max_wt] = max(r2_wt);
exp_data = table2array(...
    readtable(fullfile("mgat1ko_data/wt-preprocessed",pa_wt.file_name{max_wt})));
yksum = table2array(...
    readtable(fullfile(calib_dir,"wt_yhat",pa_wt.file_name{max_wt})));
t = exp_data(:,1);
pa_row = pa_wt(string(pa_wt.file_name) == pa_wt.file_name{max_wt},:);

for i=1:length(volts)
    i2 = vidx(i);
    yi = exp_data(:,i2+1);
    yhati = yksum(:,i2);
    r2 = pa_row{1,2+(i2-1)*2};
    l2 = pa_row{1,3+(i2-1)*2};

    if i==1
        subplot(6,3,1+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
        ylabel(" ")
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--');
    end

    if i==2
        subplot(6,3,1+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));

        legend(["Experimental I_{K}", "Model Prediction"])
        legend Box off
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--');
    end

    if ismember(i,3:4)
        subplot(6,3,2+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--');
    end

    if ismember(i,5:6)
        subplot(6,3,3+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--');
    end
end

[~,max_ko] = max(r2_ko);
exp_data = table2array(...
    readtable(fullfile("mgat1ko_data/mgat1ko-preprocessed",pa_ko.file_name{max_ko})));
yksum = table2array(...
    readtable(fullfile(calib_dir,"mgat1ko_yhat",pa_ko.file_name{max_ko})));
t = exp_data(:,1);
pa_row = pa_ko(string(pa_ko.file_name) == pa_ko.file_name{max_ko},:);

for i=1:length(volts)
    i2 = vidx(i);
    yi = exp_data(:,i2+1);
    yhati = yksum(:,i2);
    r2 = pa_row{1,2+(i2-1)*2};
    l2 = pa_row{1,3+(i2-1)*2};

    if ismember(i,1:2)
        subplot(6,3,10+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')
    end

    if ismember(i,3:4)
        subplot(6,3,11+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')        
    end

    if ismember(i,5:6)
        subplot(6,3,12+i)
        plot(t,yi,'Color','magenta')
        hold on
        plot(t,yhati,'--','Color','green','LineWidth',1.5)
        hold off
        axis tight
        grid on
        title(strcat(...
            string(volts(i))," mV / ",num2str(round(r2,4))));
        set(gca,'LineWidth',1.5,'FontWeight','bold');
        set(gca,'GridLineStyle','--')        
    end
end
