close all, clear all
y =@(x)  (5*log(-50/49*tanh(atanh(2^(1/2)/10) + 2^(1/2)*x).^2 + 50/49))/2+ 8;
 dy =@(x) -(250*2^(1/2)*tanh(2^(1/2)*x + 5129623148137269/36028797018963968)...
     .*(tanh(2^(1/2)*x + 5129623148137269/36028797018963968).^2 - 1))...
     ./(49*((50*tanh(2^(1/2)*x + 5129623148137269/36028797018963968).^2)/49 - 50/49));
ddy =@(x) (500.*(tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2 - 1).^2)./(49.*((50.*tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2)./49 - 50./49)) + (1000.*tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2.*(tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2 - 1))./(49.*((50.*tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2)./49 - 50./49)) - (50000.*tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2.*(tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2 - 1).^2)./(2401.*((50.*tanh(2.^(1./2).*x + 5129623148137269./36028797018963968).^2)./49 - 50./49).^2);

Color_List = [[0 0.5 0];
    [0.8 0.4 0.8];
    [0.8 .5 .2];
    [119/255 158/255 203/255];]
    
G =@(x) [ x,-x.^2/2];
% 
dG=@(x)[ones(size(x)),-x];
ddG=@(x)[zeros(size(x)),-ones(size(x))];

%Yv =@(xe) exp(xe/10).*sin(x);
X = [0:.001:1]';

%for k = 1:1000
QR = figure(1) 
set(QR,'position',[0 0 800 235])
subplot(1,3,1)
d = [101 151  301 401 651 701 851 951 1001]';
n= length(d)
sig = 0.25;
e = normrnd(0,sig,length(d),1);
Y = y(X(d))+e;
Gv = G(X(d));
Gt = G(X);
figure(1)
hold on
plot(X,y(X),'k','linewidth',4);
hold on
plot(X,8+Gt*([-1.5,6]'),'--','color',[.6 .6 .6],'linewidth',4);
hold on
plot(X,8+Gt*([-1,10]'),':','color',[.6 .6 .6],'linewidth',4);
ylabel('Vertical height')
xlabel('Time')

subplot(1,3,2)
Gv = dG(X(d));
Gt = dG(X);
figure(1)
hold on
plot(X,dy(X),'k','linewidth',4);
hold on
plot(X,Gt*([-1.5,6]'),'--','color',[.6 .6 .6],'linewidth',4);
hold on
plot(X,Gt*([-1,10]'),':','color',[.6 .6 .6],'linewidth',4);
ylabel('Derivative of vertical height')
xlabel('Time')

subplot(1,3,3)
Gv = ddG(X(d));
Gt = ddG(X);
figure(1)
hold on
plot(X,ddy(X),'k','linewidth',4);
ylabel('2^n^d derivative of vertical height')
hold on
plot(X,Gt*([-1.5,6]'),'--','color',[.6 .6 .6],'linewidth',4);
hold on
plot(X,Gt*([-1,10]'),':','color',[.6 .6 .6],'linewidth',4);
xlabel('Time')
ylim([-11,0])
tightfig
