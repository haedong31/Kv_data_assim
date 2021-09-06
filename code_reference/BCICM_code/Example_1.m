clear all
close all
clc

Color_List = [[0 0.5 0];
    [0.8 0.4 0.4];
    [.8 .6 .2];
    [119/255 158/255 203/255];]
x = [0:0.05:0.8]';
n = length(x);

sig = 1/50;
y = x.*(4+sin(5*x))+normrnd(0,sig,n,1);
s = 1;
X = [0:0.005:1]';
N = length(X);
[~,d] = ismember(round(x/0.005),round(X/0.005));

q = X.*(4+sin(5*X));
dq = (4+sin(5*X))+5*X.*cos(5*X);

tv = [2.1:.001:4.10];
psi = 1/2;
Rv = corr_matern3(psi,X,X);
rv = corr_matern3(psi,x,X);

QR = figure(1)
set(QR,'position',[0 0 900 250])

L1 =@(t) mean((q - t*X).^2);
L2 =@(t) mean((q - t*X).^2) +mean((dq - t).^2);
C = chol(Rv)
L3 =@(t) sum(((q - t*X)'/C').^2);

for j = 1:length(tv)
    L1v(j) = L1(tv(j));
    L2v(j) = L2(tv(j));
    L3v(j) = L3(tv(j));
end

R = subplot(1,2,1)
plot(tv,(L1v-min(L1v))/(max(L1v)-min(L1v)),'-','linewidth',3,'color',[0 .5 0])
text(max(tv)+0.15,(L1v(end)-min(L1v))/(max(L1v)-min(L1v))-0.05,'$L_{L^2}$','interpreter','latex','HorizontalAlignment','center','fontsize',16,'color','k')
hold on
plot(tv,(L2v-min(L2v))/(max(L2v)-min(L2v)),'-','linewidth',3,'color',[1 .5 .5])
text(max(tv)+0.15,(L2v(end)-min(L2v))/(max(L2v)-min(L2v))-0.05,'$L_{W_1^2}$','interpreter','latex','HorizontalAlignment','center','fontsize',16,'color','k')

xlim([min(tv) max(tv)])

hold on
plot(tv,(L3v-min(L3v))/(max(L3v)-min(L3v)),'k--','linewidth',3)
text(max(tv)+0.15,(L3v(end)-min(L3v))/(max(L3v)-min(L3v))-0.05,'$L_{KO}$','interpreter','latex','HorizontalAlignment','center','fontsize',16,'color','k');
ylabel({'normalized loss'})
ylim([0,1.05])

set(R,'ytick',[])
xlabel({'$t$'},'FontSize',18,'interpreter','latex');


Lab_val = {'$L_{L^2}$OGP','$L_{W_1^2}$OGP','KO','No bias'};
theta_0 = (X'*X)\(X'*q);
theta_1 = (X'*X+length(X))\(X'*q +ones(1,length(X))*dq);
theta_3 = (X'*Rv^(-1)*X)\(X'*Rv^(-1)*q);
for k = 1:4
if k ==1 
h =  (exp(-2*X)/2 + (3*exp(2*X - 2))/2 + 2).*X + (3*exp(-2*X))/4 - (13*exp(2*X - 2))/4;
C = 5.17932849197041*corr_matern3(psi,X,X)-23.516174328494*h*h';
C = 0.5^2*1/2*(C+C');
colorv = Color_List(1,:);
type = '-';
end
if k ==2 
h =  (2 - exp(2*X - 2)/2 - (3*exp(-2*X))/2).*X - exp(-2*X)/4 - exp(2*X - 2)/4;
  C = 1.273090329685822754*corr_matern3(psi,X,X)-0.7828374331437503133170342*h*h';
C = 0.5^2*1/2*(C+C');
colorv = Color_List(2,:);
type = '-';
end
if k ==3
C =0.5^2*corr_matern3(psi,X,X);
colorv = 'k';
type = '--';
end
if k ==4
C =zeros(N,N);
colorv = 'k';
type = ':';
end

Cp = C(d,d)+sig^2*eye(length(d));

mu_T = (x'*Cp^(-1)*x + s)^(-1)*x'*Cp^(-1)*y;
S_T = (x'*Cp^(-1)*x +s)^(-1);
mu_z =  C(:,d)*Cp^(-1)*(y-mu_T*x);
s_z = S_T*(C(:,d)*Cp^(-1)*x)*(C(:,d)*Cp^(-1)*x)'+ (C-C(:,d)*Cp^(-1)*C(d,:));

subplot(1,2,2)
hold on
plot(tv,normpdf(tv,mu_T,sqrt(S_T))/normpdf(mu_T,mu_T,sqrt(S_T)),type,'linewidth',3,'color',colorv)

text(mu_T,1+0.05,Lab_val{k},'interpreter','latex','fontsize',12,'VerticalAlignment','middle','HorizontalAlignment','center');
end

ylim([0,1.1])


subplot(1,2,2)



plot(theta_0,0,'o','color',Color_List(1,:),'markersize',15,'linewidth',3,'markerfacecolor',Color_List(1,:))
text(theta_0,-0.1,'$\theta, L_{L^2}$','interpreter','latex','HorizontalAlignment','center','fontsize',12,'color',Color_List(1,:));
hold on
plot(theta_1,0,'rs','markersize',15,'linewidth',3,'markerfacecolor',Color_List(2,:),'color',Color_List(2,:))
text(theta_1,-0.1,'$\theta, L_{W_1^2}$','interpreter','latex','HorizontalAlignment','center','fontsize',12,'color','r');
ylabel('normalized posterior density')

xlim([min(tv) max(tv)])
set(gca,'XTick',[min(tv),max(tv)])
set(gca,'yTick',[])
set(gca,'yaxislocation','left')


xlabel({'$\theta$'},'FontSize',18,'interpreter','latex');
tightfig