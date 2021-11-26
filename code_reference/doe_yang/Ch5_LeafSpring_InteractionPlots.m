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

effectname = {'B','C','D','E','Q','BQ','CQ','DQ','EQ','BC','BD','CD','BCQ','BDQ','BEQ'}';

%% B-against-Q   - location effects
zBnQn = mean(ym(find(x(:,1)==-1 & x(:,5)==-1)));
zBpQn = mean(ym(find(x(:,1)==1 & x(:,5)==-1)));
zBnQp = mean(ym(find(x(:,1)==-1 & x(:,5)==1)));
zBpQp = mean(ym(find(x(:,1)==1 & x(:,5)==1)));

figure('color','w');
plot(1:2,[zBnQn,zBnQp],'.-','LineWidth',2);hold on;
plot(1:2,[zBpQn,zBpQp],'.--','LineWidth',2);hold off;
text(2.2,zBnQp,'B-');
text(2.2,zBpQp,'B+');
set(gca,'XTick',1:2,'XTickLabel',{'-';'+'});
xlim([0.5 2.5]);
xlabel('Q','FontSize',10,'FontWeight','bold');
title('BQ - location effects','FontSize',10,'FontWeight','bold');
set(gca,'LineWidth',1,'FontSize',10,'FontWeight','bold');

%% C-against-Q  - location effects
zCnQn = mean(ym(find(x(:,2)==-1 & x(:,5)==-1)));
zCpQn = mean(ym(find(x(:,2)==1 & x(:,5)==-1)));
zCnQp = mean(ym(find(x(:,2)==-1 & x(:,5)==1)));
zCpQp = mean(ym(find(x(:,2)==1 & x(:,5)==1)));

figure('color','w');
plot(1:2,[zCnQn,zCnQp],'.-','LineWidth',2);hold on;
plot(1:2,[zCpQn,zCpQp],'.--','LineWidth',2);hold off;
text(2.2,zCnQp,'C-');
text(2.2,zCpQp,'C+');
set(gca,'XTick',1:2,'XTickLabel',{'-';'+'});
xlim([0.5 2.5]);
xlabel('Q','FontSize',10,'FontWeight','bold');
title('CQ - location effects','FontSize',10,'FontWeight','bold');
set(gca,'LineWidth',1,'FontSize',10,'FontWeight','bold');


%% D-against-Q   - dispersion effects
zDnQn = mean(lns2(find(x(:,3)==-1 & x(:,5)==-1)));
zDpQn = mean(lns2(find(x(:,3)==1 & x(:,5)==-1)));
zDnQp = mean(lns2(find(x(:,3)==-1 & x(:,5)==1)));
zDpQp = mean(lns2(find(x(:,3)==1 & x(:,5)==1)));

figure('color','w');
plot(1:2,[zDnQn,zDnQp],'.-','LineWidth',2);hold on;
plot(1:2,[zDpQn,zDpQp],'.--','LineWidth',2);hold off;
text(2.2,zDnQp,'D-');
text(2.2,zDpQp,'D+');
set(gca,'XTick',1:2,'XTickLabel',{'-';'+'});
xlim([0.5 2.5]);
xlabel('Q','FontSize',10,'FontWeight','bold');
title('DQ - dispersion effects','FontSize',10,'FontWeight','bold');
set(gca,'LineWidth',1,'FontSize',10,'FontWeight','bold');


%% BCQ   - dispersion effects
zBnCnQn = mean(lns2(find(x(:,1)==-1 & x(:,2)==-1 & x(:,5)==-1)));
zBnCnQp = mean(lns2(find(x(:,1)==-1 & x(:,2)==-1 & x(:,5)==1)));

zBnCpQn = mean(lns2(find(x(:,1)==-1 & x(:,2)==1 & x(:,5)==-1)));
zBnCpQp = mean(lns2(find(x(:,1)==-1 & x(:,2)==1 & x(:,5)==1)));

zBpCnQn = mean(lns2(find(x(:,1)==1 & x(:,2)==-1 & x(:,5)==-1)));
zBpCnQp = mean(lns2(find(x(:,1)==1 & x(:,2)==-1 & x(:,5)==1)));

zBpCpQn = mean(lns2(find(x(:,1)==1 & x(:,2)==1 & x(:,5)==-1)));
zBpCpQp = mean(lns2(find(x(:,1)==1 & x(:,2)==1 & x(:,5)==1)));

figure('color','w');
plot(1:2,[zBnCpQn,zBnCpQp],'.-','LineWidth',2);hold on;
plot(1:2,[zBpCpQn,zBpCpQp],'.--','LineWidth',2);
plot(1:2,[zBpCnQn,zBpCnQp],'.-.','LineWidth',2);
plot(1:2,[zBnCnQn,zBnCnQp],'.:','LineWidth',2);
hold off;
text(2.2,zBnCnQp,'B-C-');
text(2.2,zBnCpQp,'B-C+');
text(2.2,zBpCnQp,'B+C-');
text(2.2,zBpCpQp+0.2,'B+C+');

set(gca,'XTick',1:2,'XTickLabel',{'-';'+'});
xlim([0.5 2.5]);
xlabel('Q','FontSize',10,'FontWeight','bold');
title('BCQ - dispersion effects','FontSize',10,'FontWeight','bold');
set(gca,'LineWidth',1,'FontSize',10,'FontWeight','bold');

