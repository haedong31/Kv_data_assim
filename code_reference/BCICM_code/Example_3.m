clear all, close all, clc
load('heartdata')


Color_List = [[0 0.5 0];
    [0.8 0.4 0.8];
    [0.8 .5 .2];
    [119/255 158/255 203/255];];

cV =@(x1,x2) corr_matern3(1,x1,x2);
dcV =@(x1,x2) corr_dmatern3(1,x1,x2);
ddcV =@(x1,x2) corr_ddmatern3(1,x1,x2);
c = cV(X,X);

r =@(w) -postKO(w,X,Y,c);
[xs,~,~,~,~,W] = fminunc(r,[1 1 1.5]');
xs;
s1 = MH_me(r,W^(-1/2),1,xs,10000);

c1 = cV(X,Xd);
C1 = cV(Xd,Xd);
c2 = dcV(X,Xd);
C2 = ddcV(Xd,Xd);
C12 = dcV(Xd,Xd);

r2 =@(w) -postOGPL2(w,X,Y,Xd,c,c1,C1);
[xs2,~,~,~,~,W2] = fminunc(r2,[1 1 1.5]');
s2 = MH_me(r2,W2^(-1/2),1,xs2,10000);

heart_f(X,[1 1 1.5]')
r3 =@(w) -postOGPW1(w,X,Y,Xd,c,c1,C1,c2,C2,C12);
[xs3,~,~,~,~,W3] = fminunc(r3,[1 1 1.5]');
s3 = MH_me(r3,W3^(-1/2),1,xs3,10000);

QR = figure(1);
set(QR,'Position',[0 0 900 250])

subplot(1,3,1)
YL = 0.05;
for i = size(s1,1)-500:100:size(s1,1)
    theta = s1(i,:)';
    fKO = heart_f(Xd,theta);  
    fKO(100)
    
    subplot(2,3,1)
    plot(Xd,fKO,'color',[.4 .4 .4]);
    hold on
    
    Yh = heart_f(X,theta);
    [~, dYhd,~,dYhdx] = heart_f(Xd,theta);
    
      
    R = c;
    Rd  = c1;
    Rd2  = C1;
    
    S = 0.02^2*R+0.001^2*eye(length(X));
    s = 0.02^2*Rd;
    
    zh = (s'/S)*(Y-Yh);
    Sh = 0.02^2*Rd2 - (s'/S)*s+10^(-9)*eye(length(Xd));
    Sh = 1/2*(Sh + Sh'); 
    
    subplot(2,3,4)
    hold on
    Z = mvnrnd(zh,Sh);
    plot(Xd,Z,'color',[.4 .4 .4]);
    xlim([min(X) max(X)])
    YL = max(YL,max(abs(Z)));
    ylim([-YL,YL])
end
ylabel('bias')
xlabel('log(time)')
hold on
    subplot(2,3,1)
title('KO','interpreter','latex')
plot(X,Y,'o','markerfacecolor',[0 0 0],'color',[0 0 0],'markersize',3)
ylim([0,.18])
xlim([min(X) max(X)])
set(gca,'xticklabel',{})
ylabel('normalized current')

for i = size(s1,1)-500:100:size(s1,1)
    theta = s2(i,:)';
    fKO = heart_f(Xd,theta);
    subplot(2,3,2)
    hold on
    plot(Xd,fKO,'color',Color_List(1,:));
    Yh = heart_f(X,theta);
    [~, dYhd,~,dYhdx] = heart_f(Xd,theta);
    
    h = c1*dYhd';
    hd =  C1*dYhd';
    H = dYhd*C1*dYhd';
    
    
    R  = c -  h*H^(-1)*(h');%
    WC = mean(diag(R));
    
    R= R/WC;
    Rd  = (c1 -  h*H^(-1)*(hd'))/WC;
    Rd2  = (C1 -  hd*H^(-1)*(hd'))/WC;
    
    S = 0.02^2*R+0.001^2*eye(length(X));
    s = 0.02^2*Rd;
        
    zh = (s'/S)*(Y-Yh);
    Sh = 0.02^2*Rd2 - (s'/S)*s+10^(-9)*eye(length(Xd));
    Sh = 1/2*(Sh + Sh'); 
    
    subplot(2,3,5)
    hold on
    plot(Xd,mvnrnd(zh,Sh),'color',Color_List(1,:));
    xlim([min(X) max(X)])
    ylim([-YL,YL])
    
end
hold on
xlabel('log(time)')
    subplot(2,3,2)
title('OGP, $L_{L^2}$','interpreter','latex')
plot(X,Y,'o','markerfacecolor',[0 0 0],'color',[0 0 0],'markersize',3)
ylim([0,.18])
set(gca,'ytick',[]);
xlim([min(X) max(X)])
set(gca,'xticklabel',{})



for i = size(s1,1)-500:100:size(s1,1)
    theta = s3(i,:)';
    fKO = heart_f(Xd,theta);

subplot(2,3,3)
hold on
    plot(Xd,fKO,'color',Color_List(2,:));
    
    
    Yh = heart_f(X,theta);
    [~, dYhd,~,dYhdx] = heart_f(Xd,theta);
    subplot(2,3,3)
    h = c1*dYhd'+c2*dYhdx';
    hd =  C1*dYhd'+C12*dYhdx';
    H = dYhd*C1*(dYhd') + dYhd*C12*(dYhdx')+ (dYhd*C12*(dYhdx'))'+(dYhdx*C2*(dYhdx'))+10^(-6)*eye(3);
    
    
    R  = c -  h*H^(-1)*(h');
    WC = mean(diag(R));
    
    R= R/WC;
    Rd  = (c1 -  h*H^(-1)*(hd'))/WC;
    Rd2  = (C1 -  hd*H^(-1)*(hd'))/WC;
    
    S = 0.02^2*R+0.001^2*eye(length(X));
    s = 0.02^2*Rd;
        
    zh = (s'/S)*(Y-Yh);
    Sh = 0.02^2*Rd2 - (s'/S)*s+10^(-9)*eye(length(Xd));
    Sh = 1/2*(Sh + Sh'); 
    
    subplot(2,3,6)
    hold on
    plot(Xd,mvnrnd(zh,Sh),'color',Color_List(2,:));
    hold on
    xlim([min(X) max(X)])
    ylim([-YL,YL])
end
hold on
xlabel('log(time)')
    figure(1)
    subplot(2,3,3)
title('OGP, $L_{W_1^2}$','interpreter','latex')
plot(X,Y,'o','markerfacecolor',[0 0 0],'color',[0 0 0],'markersize',3)
ylim([0,.18])
set(gca,'ytick',[]);
xlim([min(X) max(X)])
set(gca,'xticklabel',{})
tightfig
