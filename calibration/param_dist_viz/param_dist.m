%% visualize distributions of ikto parameters
clc
clearvars
close all

pkto = readtable('pkto.csv');
param = unique(pkto.param);

figure('Color','w', 'Position',[100, 100, 600, 510])
for i=1:length(param)
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
    set(gca, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')

    if i==4
        legend('WT','Mgat1KO', 'Location','best')
        legend box off
    end
end

%% visualize distributions of ikslow1 parameters
clc
clearvars
close all

pkslow1 = readtable('pkslow1.csv');
param = unique(pkslow1.param);

figure('Color','w', 'Position',[100, 100, 600, 680])
for i=1:length(param)
    psub = pkslow1(pkslow1.param==param(i),:);
    psub_wt = psub(string(psub.Group)=='WT',:);
    psub_ko = psub(string(psub.Group)=='Mgat1KO',:);

    [f1,xi1] = ksdensity(psub_wt.value);
    [f2,xi2] = ksdensity(psub_ko.value);
    
    subplot(4,3,i)
    plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    hold on
    plot(xi2,f2, 'Color','red', 'LineWidth',2)
    hold off
    axis tight
    xlabel(strcat('Parameter', num2str(param(i))))
    ylabel('Density')
    set(gca, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')

    if i==4
        legend('WT','Mgat1KO', 'Location','best')
        legend box off
    end
end

%% ikslow2
clc
clearvars
close all

pkslow2 = readtable('pkslow2.csv');
param = unique(pkslow2.param);

figure('Color','w', 'Position',[100, 100, 400, 170])
for i=1:length(param)
    psub = pkslow2(pkslow2.param==param(i),:);
    psub_wt = psub(string(psub.Group)=='WT',:);
    psub_ko = psub(string(psub.Group)=='Mgat1KO',:);

    [f1,xi1] = ksdensity(psub_wt.value);
    [f2,xi2] = ksdensity(psub_ko.value);
    
    subplot(1,2,i)
    plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    hold on
    plot(xi2,f2, 'Color','red', 'LineWidth',2)
    hold off
    axis tight
    xlabel(strcat('Parameter', num2str(param(i))))
    ylabel('Density')
    set(gca, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')

    if i==1
        legend('WT','Mgat1KO', 'Location','best')
        legend box off
    end
end

%% ikss
clc
clearvars
close all

pkss = readtable('pkss.csv');
param = unique(pkss.param);

figure('Color','w', 'Position',[100, 100, 400, 170])
for i=1:length(param)
    psub = pkss(pkss.param==param(i),:);
    psub_wt = psub(string(psub.Group)=='WT',:);
    psub_ko = psub(string(psub.Group)=='Mgat1KO',:);

    [f1,xi1] = ksdensity(psub_wt.value);
    [f2,xi2] = ksdensity(psub_ko.value);
    
    subplot(1,2,i)
    plot(xi1,f1, 'Color','blue', 'LineWidth',2)
    hold on
    plot(xi2,f2, 'Color','red', 'LineWidth',2)
    hold off
    axis tight
    xlabel(strcat('Parameter', num2str(param(i))))
    ylabel('Density')
    set(gca, 'FontName','Arial', 'FontSize',11, 'FontWeight','bold')

    if i==1
        legend('WT','Mgat1KO', 'Location','best')
        legend box off
    end
end