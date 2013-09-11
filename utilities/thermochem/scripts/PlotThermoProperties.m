clear all
close all
clc


% Load thermodynamic database
db = ReadDb();

species = {'SiO2'};

T = 250:5:3000;

figure;
set(gcf,'Units','inches','Position',[2 2 12 5]);
subplot(1,3,1)
hold all
for i = 1:length(species)
    plot(T,Cp(db.(species{i}),T));
end

subplot(1,3,2)
hold all
for i = 1:length(species)
    plot(T,H(db.(species{i}),T));
end


subplot(1,3,3)
hold all
h = zeros(size(species));
for i = 1:length(species)
    h(i) = plot(T,S(db.(species{i}),T),'DisplayName',species{i});
end
legend(h)


