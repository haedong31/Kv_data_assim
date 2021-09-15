clc
close all
clear variables

rmse_df = readtable('bar_graph_df2.csv');
rmse_df = table2array(rmse_df(:,2:3));

avg1 = mean(rmse_df(:,1));
avg2 = mean(rmse_df(:,2));

b = bar(rmse_df);
b(1).FaceColor = 'green';
b(2).FaceColor = 'blue';
hold on
yline(avg1, 'Color', 'green', 'LineWidth',2)
yline(avg2, 'Color', 'blue', 'LineWidth',2)
hold off
axis tight
grid on
title('WT')
xlabel('Data Index', 'FontSize',10, 'FontWeight','bold');
ylabel('Sum of RMSEs (-30 mV ~ 50 mV)', 'FontSize',10, 'FontWeight','bold');
legend('Calibration 1', 'Calibration 2')
get(gcf,'CurrentAxes');
set(gca,'YDir','normal');
set(gca,'LineWidth',2, 'FontSize',10, 'FontWeight','bold');
set(gca,'GridLineStyle','--')
