%% consolidate half-normal plots - normalized (2x3 layout)
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

fig = figure('Color','w','Position',[100,100,770,420]);
orient(fig,'landscape')

% marker 1
% kto
subplot(2,3,1)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,1),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',8,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,1),'range'),dgn_kslow1);
plot(phiinv,theta,'x','MarkerSize',8,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kslow2
idx_kslow2 = [2,1];
[theta,~,phiinv] = doe_plot(normalize(response_kslow2(:,1),'range'),dgn_kslow2);
plot(phiinv,theta,'*','MarkerSize',8,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(i),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,1),'range'),dgn_kss);
plot(phiinv,theta,'o','MarkerSize',8,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end
hold off
title("Marker a")
legend(["I_{Kto}","I_{Kslow1}","I_{Kslow2}","I_{Kss}"],...
    'Location','northeast')
xlabel('Half-normal Quantiles')
ylabel('Normalized Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% marker 2
% kto
subplot(2,3,2)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,2),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',8,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,2),'range'),dgn_kslow1);
plot(phiinv,theta,'x','MarkerSize',8,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kslow2
[theta,~,phiinv] = doe_plot(normalize(response_kslow2(:,2),'range'),dgn_kslow2);
plot(phiinv,theta,'*','MarkerSize',8,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(i),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,2),'range'),dgn_kss);
plot(phiinv,theta,'o','MarkerSize',8,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end
hold off
title("Marker b")
xlabel('Half-normal Quantiles')
ylabel('Normalized Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% marker 3
% kto
subplot(2,3,3)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,3),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',8,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,3),'range'),dgn_kslow1);
plot(phiinv,theta,'x','MarkerSize',8,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kslow2
[theta,~,phiinv] = doe_plot(normalize(response_kslow2(:,3),'range'),dgn_kslow2);
plot(phiinv,theta,'*','MarkerSize',8,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(i),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,3),'range'),dgn_kss);
plot(phiinv,theta,'o','MarkerSize',8,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end
hold off
title("Marker c")
xlabel('Half-normal Quantiles')
ylabel('Normalized Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% marker 4
% kto
subplot(2,3,4)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,4),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',8,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,4),'range'),dgn_kslow1);
plot(phiinv,theta,'x','MarkerSize',8,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kslow2
[theta,~,phiinv] = doe_plot(normalize(response_kslow2(:,4),'range'),dgn_kslow2);
plot(phiinv,theta,'*','MarkerSize',8,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(i),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,4),'range'),dgn_kss);
plot(phiinv,theta,'o','MarkerSize',8,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end
hold off
title("Marker d")
xlabel('Half-normal Quantiles')
ylabel('Normalized Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% marker 5
% kto
subplot(2,3,5)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,5),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',8,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,5),'range'),dgn_kslow1);
plot(phiinv,theta,'x','MarkerSize',8,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kslow2
[theta,~,phiinv] = doe_plot(normalize(response_kslow2(:,5),'range'),dgn_kslow2);
plot(phiinv,theta,'*','MarkerSize',8,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(i),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,5),'range'),dgn_kss);
plot(phiinv,theta,'o','MarkerSize',8,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end
hold off
title("Marker e")
xlabel('Half-normal Quantiles')
ylabel('Normalized Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

% marker 6
% kto
subplot(2,3,6)
[theta,idx,phiinv] = doe_plot(normalize(response_kto(:,6),'range'),dgn_kto);
plot(phiinv,theta,'+','MarkerSize',8,'Color',c(1,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

hold on
% kslow1
[theta,idx,phiinv] = doe_plot(normalize(response_kslow1(:,6),'range'),dgn_kslow1);
plot(phiinv,theta,'x','MarkerSize',8,'Color',c(2,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kslow2
[theta,~,phiinv] = doe_plot(normalize(response_kslow2(:,6),'range'),dgn_kslow2);
plot(phiinv,theta,'*','MarkerSize',8,'Color',c(3,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(i),'HorizontalAlignment','right','VerticalAlignment','top');
end

% kss
[theta,idx,phiinv] = doe_plot(normalize(response_kss(:,6),'range'),dgn_kss);
plot(phiinv,theta,'o','MarkerSize',8,'Color',c(4,:),'LineWidth',1.5)
for i=1:length(theta)
    text(phiinv(i),theta(i),int2str(idx(i)),'HorizontalAlignment','right','VerticalAlignment','top');
end
hold off
title("Marker f")
xlabel('Half-normal Quantiles')
ylabel('Normalized Absolute Effects')
set(gca,'FontWeight','bold','LineWidth',1.5)

%% custom functions
function [theta,idx,phiinv] = doe_plot(y,dgn)
    stat_beta = regress(y,dgn);
    beta = 2*stat_beta;
    beta = abs(beta);

    [theta,idx] = sort(beta);
    p = 0.5 + (0.5*(1:length(beta))-0.5)/length(beta);
    phiinv = norminv(p,0,1);
end
