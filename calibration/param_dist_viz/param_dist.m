clc
clearvars
close all

%% visualize distributions of ikto parameters
pkto = readtable('pkto.csv');
param = unique(pkto.param);
psub = pkto(pkto.param==param(1),:);
psub_wt = psub(string(psub.Group)=='WT',:);
psub_ko = psub(string(psub.Group)=='Mgat1KO',:);

[f1,xi1] = ksdensity(psub_wt.value);
[f2,xi2] = ksdensity(psub_ko.value);

figure('Color','w', 'Position',[100, 100, 600, 500])
subplot(3,3,1)
plot(xi1,f1, 'Color','blue', 'LineWidth',2)
hold on
plot(xi2,f2, 'Color','red', 'LineWidth',2)
hold off
axis tight
xlabel(strcat('Parameter', num2str(param(1))))
ylabel('Density')
set(gca, 'LineWidth',2, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')
legend('WT','Mgat1KO', 'Location','best')
legend box off

for i=2:length(param)
    psub = pkto(pkto.param==param(i),:);
    psub_wt = psub(string(psub.Group)=='WT',:);
    psub_ko = psub(string(psub.Group)=='Mgat1KO',:);

    [f1,xi1] = ksdensity(psub_wt.value);
    [f2,xi2] = ksdensity(psub_ko.value);
    
    subplot(3,3,i)
    plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    hold on
    plot(xi2,f2, 'Color','red', 'LineWidth',2)
    hold off
    axis tight
    xlabel(strcat('Parameter', num2str(param(i))))
    ylabel('Density')
    set(gca, 'LineWidth',2, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')
end

%% visualize distributions of ikslow1 parameters