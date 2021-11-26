% Dr. Hui Yang
% prepared for  ESI 6247 Statistical Design Models
% Deaprtment of Industrial and Management Systems Engineering
% University of South Florida
% Email: huiyang@usf.edu

% Chapter 5 Leaf Spring Experiment - Fractional Factorial Experiment

clear all
close all
clc

x = xlsread('LeafSpring.xlsx','LeafSpring','A1:E16');
height = xlsread('LeafSpring.xlsx','LeafSpring','F1:H16');
ym = xlsread('LeafSpring.xlsx','LeafSpring','I1:I16');
lns2 = xlsread('LeafSpring.xlsx','LeafSpring','K1:K16');

%% factorial effects, adapted epitaxial layer growth experiment

% 2 factor interaction model matrix
x(:,6) = x(:,1).*x(:,5);
x(:,7) = x(:,2).*x(:,5);
x(:,8) = x(:,3).*x(:,5);
x(:,9) = x(:,4).*x(:,5);
x(:,10) = x(:,1).*x(:,2);
x(:,11) = x(:,1).*x(:,3);
x(:,12) = x(:,2).*x(:,3);

% 3 factor interaction model matrix
x(:,13) = x(:,1).*x(:,2).*x(:,5);
x(:,14) = x(:,1).*x(:,3).*x(:,5);
x(:,15) = x(:,1).*x(:,4).*x(:,5);

stats = regress(ym,x);
ymbeta = stats.*2;
stats = regress(lns2,x);
lns2beta = stats.*2;


%% command window outputs
fprintf('\n')
fprintf('Factorial effects, leaf spring experiment');
fprintf('\n\n')

effectname = {'B','C','D','E','Q','BQ','CQ','DQ','EQ','BC','BD','CD','BCQ','BDQ','BEQ'}';

fprintf('%15s','Effect','Y-mean','lns2');
fprintf('\n');
for i = 1:length(ymbeta)
    fprintf('%15s',char(effectname(i)));
    fprintf('%14.3f',ymbeta(i),lns2beta(i));
    fprintf('\n');
end


%% Half-normal plot of location effects
absymbeta = abs(ymbeta);
[theta,idx] = sort(absymbeta);
p = 0.5+(0.5*(1:length(ymbeta))-0.5)/length(ymbeta);
phiinv = norminv(p,0,1);
figure('color','w')
plot(phiinv,theta,'.','MarkerSize',10);
for i = 1:length(ymbeta)
    text(phiinv(i),theta(i)-0.008,effectname(idx(i)));
end
xlabel('Half-normal quantiles');
ylabel('absolute effects');
title('Half-normal plot of location effects');
axis([-0.2 2.5 -0.01 0.28]);

%% Half-normal plot of dispersion effects
abslns2beta = abs(lns2beta);
[theta,idx] = sort(abslns2beta);
p = 0.5+(0.5*(1:length(lns2beta))-0.5)/length(lns2beta);
phiinv = norminv(p,0,1);
figure('color','w')
plot(phiinv,theta,'.','MarkerSize',10);
for i = 1:length(lns2beta)
    text(phiinv(i),theta(i)-0.08,effectname(idx(i)));
end
xlabel('Half-normal quantiles');
ylabel('absolute effects');
title('Half-normal plot of dispersion effects');
axis([-0.2 2.5 -0.01 2]);

%% Lenth's method, |tPSE| values
absymbeta = abs(ymbeta);
abslns2beta = abs(lns2beta);
s0ym = 1.5*median(absymbeta);
s0lns2 = 1.5*median(abslns2beta);
PSEym = 1.5*median(absymbeta(find(absymbeta<2.5*s0ym)));
PSElns2 = 1.5*median(abslns2beta(find(abslns2beta<2.5*s0lns2)));

t_pse_ym = ymbeta./PSEym;
t_pse_lns2 = lns2beta./PSElns2;

%% command window outputs
fprintf('\n')
fprintf('Lenth method, |tPSE| values, leaf spring experiment');
fprintf('\n\n')

fprintf('%15s','Effect','Y-mean','lns2');
fprintf('\n');
for i = 1:length(ymbeta)
    fprintf('%15s',char(effectname(i)));
    fprintf('%14.3f',abs(t_pse_ym(i)),abs(t_pse_lns2(i)));
    fprintf('\n');
end
