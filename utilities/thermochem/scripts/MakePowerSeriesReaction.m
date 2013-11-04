clear all
close all
clc

Tp = 400;

k0 = 2e11;
Ta0 = 2980;
Ta0lt = 20000;

Tht = 2000:-1:Tp;
Tlt = Tp:-1:300;
xht = 1./Tht;
xlt = 1./Tlt;

logR0ht = log(k0) - Ta0.*xht;
logR0lt = logR0ht(end) - Ta0lt.*(xlt - 1/Tp);

x = [xht xlt];
logR0 = [logR0ht logR0lt];

p = polyfit(x,logR0,4);
pOld = [-1e11 0 0 0 log(5.89e9)];

figure;
hold all
plot(xht,logR0ht)
plot(xlt,logR0lt)
plot(x,polyval(p,x),'--r')
plot(x,polyval(pOld,x),'--c')
T = 300:1:2000;

figure;
hold all
plot(T, 2e11*exp(-2980./T))
plot(1./x, exp(polyval(p,x)),'--r')
plot(1./x, exp(polyval(pOld,x)),'--c')
%xlim([300 500])

fprintf('A = %e\n',exp(p(5)));
fprintf('FIT1/ %e %e %e %e/\n',p(4),p(3),p(2),p(1));