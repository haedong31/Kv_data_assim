%% box plots of calibration parameters
clc
close all
clear variables

load('sols_ko.mat')
load('sols_wt.mat')

% remove rows with all zeros
sols_wt = sols_wt(any(sols_wt, 2), :);
sols_ko = sols_ko(any(sols_ko, 2), :);

v = 3;
g = 2;
n = 33;

boxplot_mat = zeros(v, g, n);

for i = 1:v
    for j = 1:g
        for k = 1:n
            if j == 1
                boxplot_mat(i, j, k) = sols_wt(k, i);
            else
                boxplot_mat(i, j, k) = sols_ko(k, i);
            end
        end
    end
end

bp = boxplot2(boxplot_mat, 1:v);
cmap = get(0, 'defaultaxescolororder');
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), bp);
end
legend('WT','Mgat1KO', 'Location','best')
