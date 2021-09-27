clc
close all
clear variables

rmse_df = readtable('bar_graph_exp20-21-22-23_ko.csv');
rmse_df = table2array(rmse_df(:,2:end));
[num_files, ~] = size(rmse_df);

avg1 = mean(rmse_df(:,1));
avg2 = mean(rmse_df(:,2));
avg3 = mean(rmse_df(:,3));
avg4 = mean(rmse_df(:,4));

figure('Position',[100,100,1000,800])
b = bar(rmse_df);
b(1).FaceColor = 'green';
b(2).FaceColor = 'blue';
n(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = 'cyan';

hold on
yline(avg1, 'Color', 'green', 'LineWidth',2)
yline(avg2, 'Color', 'blue', 'LineWidth',2)
yline(avg3, 'Color', [0.9290 0.6940 0.1250], 'LineWidth',2)
yline(avg4, 'Color', 'cyan', 'LineWidth',2)
hold off
axis tight
grid on
title('Mgat1KO')
xlabel('Data Index', 'FontSize',10, 'FontWeight','bold');
ylabel('Sum of RMSEs (-30 mV ~ 50 mV)', 'FontSize',10, 'FontWeight','bold');
legend('Interior Point', 'Trust Region', 'SQP', 'Active Set', 'Location','best')
xticks(1:num_files)
get(gcf,'CurrentAxes');
set(gca,'YDir','normal');
set(gca,'LineWidth',2, 'FontSize',10, 'FontWeight','bold');
set(gca,'GridLineStyle','--')
