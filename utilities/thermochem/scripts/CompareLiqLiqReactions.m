clear all
close all
clc


db = ReadDb();

test(1).A = 1e10;
test(1).B = 0;
test(1).Ea = 5900;

test(2).A = 2e10;
test(2).B = 0;
test(2).Ea = 5900;

test(3).A = 1e11;
test(3).B = 0;
test(3).Ea = 5900;

for i = 1:length(test)
    [t,T] = LiqLiqReactionRateTester(db,test(i).A,test(i).B,test(i).Ea);
    test(i).t = t;
    test(i).T = T;
end

figure;
hold all
h = zeros(size(test));
for i = 1:length(test)
    h(i)=plot(test(i).t.*1e6, test(i).T);
end
set(gca,'XScale','log')
lh = legend(h,'Location','NorthWest');
ylim([300 500])