clear all
close all
clc


% Reaction rate shaper

T = 300:1:1000;

i = 1;

r(i).A = 2e11;
r(i).B = 0;
r(i).b = [-5900/1.985 0 0 0];
i = i + 1;


r(i).B = 0;
r(i).b = [0 0 0 -1e11];
r(i).A = r(1).A*exp(r(1).b(1)/400)/(400^r(2).B*exp(r(2).b(1)./400+ r(2).b(2)./400.^2 + r(2).b(3)./400.^3 + r(2).b(4)./400.^4));
i = i + 1;

for i = 1:length(r)
    r(i).r = r(i).A.*T.^r(i).B.*exp(r(i).b(1)./T + r(i).b(2)./T.^2 + r(i).b(3)./T.^3 + r(i).b(4)./T.^4);
end

figure;
hold all

for i = 1:length(r)
    plot(T,r(i).r);
end
%ylim([0 2e8])


figure;
hold all

for i = 1:length(r)
    plot(1000./T,r(i).r);
end
set(gca,'YScale','log')
%ylim([.1 1e20])
%ylim([0 2e8])


% figure;
% hold all
% ratio = r(1).r ./ r(2).r;
% plot(T,ratio)
% 
% % for i = 1:length(r)
% %     plot(T,r(i).r);
% % end
% set(gca,'YScale','log')
% %ylim([.1 1e20])
% %ylim([0 2e8])
% 
% 
% figure;
% plot(1./T, log(1./T))
