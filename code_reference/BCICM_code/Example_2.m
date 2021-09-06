close all, clear all
y =@(x)  (5*log(-50/49*tanh(atanh(2^(1/2)/10) + 2^(1/2)*x).^2 + 50/49))/2+ 8
G =@(x) [ x;
    x.^2];

Q = figure(8)
set(Q, 'Position', [0, 0, 600, 250]);
G =@(x) [ x;
    x.^2];
phi = 1;

Y = [0.1 0.15 0.2 0.4 0.7 0.75 0.95 1.00 rand(1,1001)]';
[X,I] = sort(Y);
unsorted = 1:length(Y);
I(I) = unsorted;
d = I(1:7);
n= length(d);
Gv = G(X(d)')';
Gt = G(X')';

e = normrnd(0,1,length(d),1);
for k = 1:4
    (Gt'*Gt)^(-1)*Gt'*(y(X)-4);
    (Gt(:,2)'*Gt(:,2))^(-1)*Gt(:,2)'*(y(X)-4);
    
    
    
   
    if k == 1
         sig = 0.05;
    elseif k==2
         sig = 0.1;
    elseif k ==3
         sig = 0.25;
    elseif k ==4
         sig = 0.5;
    end
    RKO = corr_matern3(phi,X(d),X(d));
    S = corr_matern3(phi,X,X);
    
    Y = y(X(d))+sig*e;
    
    
    s0 = 0.1^2/mean(diag(RKO )) ;
    s = s0;
    for i = 1:2000
        So = s*RKO+ sig^2*eye(n);
        
        mu_KO = (Gv'*So^(-1)*Gv)^(-1)*Gv'*So^(-1)*(Y-8);
        S_KO = (Gv'*So^(-1)*Gv)^(-1);
        S_KO = .5*(S_KO+S_KO');
        
        thetasKO(:,i) = mvnrnd(mu_KO,S_KO)';
        
        mu_z =  s*RKO*So^(-1)*(Y-8-Gv*thetasKO(:,i));
        S_z = s*RKO - s^2*RKO*So^(-1)*RKO;
        S_z = .5*(S_z+S_z');
        zs = mvnrnd(mu_z,S_z)';
        s = 1/gamrnd(0.01+n/2,2/(0.01+zs'*RKO^(-1)*zs)) ;
    end
    
    ROGP  = S(d,d) -  S(d,:)*Gt*(Gt'*S*Gt)^(-1)*(S(d,:)*Gt)';
    s0 = 0.1^2/mean(diag(ROGP )) ;
    s = s0;
    for i = 1:2000
        So = s*ROGP+ sig^2*eye(n);
        
        mu_OGP = (Gv'*So^(-1)*Gv)^(-1)*Gv'*So^(-1)*(Y-8);
        S_OGP = (Gv'*So^(-1)*Gv)^(-1);
        S_OGP = .5*(S_OGP+S_OGP');
        
        thetasOGP(:,i) = mvnrnd(mu_OGP,S_OGP)';
        
        mu_z =  s*ROGP*So^(-1)*(Y-8-Gv*thetasOGP(:,i));
        S_z = s*ROGP - s^2*ROGP*So^(-1)*ROGP;
        S_z = .5*(S_z+S_z');
        zs = mvnrnd(mu_z,S_z)';
        s = 1/gamrnd(0.01+n/2,2/(0.01+zs'*ROGP^(-1)*zs)) ;
    end
    
    figure(8)
    
    subplot(2,4,k)
    plot(X,8+Gt*thetasKO(:,1500:50:end),'color',[.6 .6 .6]);
    hold on
    plot(X(d),Y,'o','markerfacecolor',[0 0 0],'color',[0 0 0],'markersize',5)
    hold on
    plot(X,y(X),'k')
    ylim([3,8])
    subplot(2,4,k+4)
    plot(X,8+Gt*thetasOGP(:,1500:50:end),'color',[.6 .6 .6]);
    hold on
    plot(X(d),Y,'o','markerfacecolor',[0 0 0],'color',[0 0 0],'markersize',5)
    hold on
    plot(X,y(X),'k')
    ylim([3,8])
    if k == 1
        L = subplot(2,4,k)
        set(L,'xtick',[]);
        ylabel('KO','interpreter','latex')
        %ylabel('Kennedy and O''Hagan')
        L = subplot(2,4,k+4)
        ylabel('OGP','interpreter','latex')
        %ylabel('Orthgonal Gaussian process')
        xlabel('$x$','interpreter','latex')
    elseif k>=2
        L = subplot(2,4,k)
        set(L,'xtick',[]);
        set(L,'ytick',[]);
        L = subplot(2,4,k+4)
        set(L,'ytick',[]);
        xlabel('$x$','interpreter','latex')
    end
    if k == 1
        L = subplot(2,4,k)
        title('$v = 0.05^2$','interpreter','latex')
    elseif k==2
        L = subplot(2,4,k)
        title('$v = 0.10^2$','interpreter','latex')
    elseif k ==3
        L = subplot(2,4,k)
        title('$v = 0.25^2$','interpreter','latex')
    elseif k ==4
        L = subplot(2,4,k)
        title('$v = 0.50^2$','interpreter','latex')
    end
    
end
Q = figure(8)
set(Q, 'Position', [0, 0, 600, 250]);
tightfig;