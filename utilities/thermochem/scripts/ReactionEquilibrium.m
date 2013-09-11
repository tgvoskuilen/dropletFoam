clear all
close all
clc


% Load thermodynamic database
db = ReadDb();

% Define reaction
%rxn.name = 'NO2+OH <=> HNO3';

% rxn.A = 2.41e13; %mol,cm3,s,K units
% rxn.B = 0;
% rxn.Ea = 0; %cal/mol
% rxn.low.A=6.42e32;
% rxn.low.B=-5.49;
% rxn.low.Ea=2349;
x = 0:0.05:1;
dH400 = zeros(size(x));
Nr = zeros(size(x));

for i = 1:length(x)

    rxn.reactants(1).name = 'CH3NHNH2L';
    rxn.reactants(1).e = 1;
    rxn.reactants(2).name = 'HNO3L';
    rxn.reactants(2).e = 4;

    rxn.products(1).name = 'MMHN0x2D3HNO3';
    rxn.products(1).e = 1;

    T = 280:5:1000;
    Kp = CalcKp(rxn,T,db);
    dH = CalcHeatRelease(rxn,T,db); %J/mol
    dH400(i) = dH(T==400);
    Nr(i) = 2*x(i)/(sum([rxn.products.e]));
end

figure;
subplot(1,2,1)
plot(T,dH./1000)

subplot(1,2,2)
plot(T,1-1e5./(8314.*T)./Kp)

figure;
hold all
plot(T,H(db.(rxn.reactants(1).name),T)./1000,'-b')
plot(T,H(db.(rxn.reactants(2).name),T)./1000,'-k')
plot(T,H(db.(rxn.products(1).name),T)./1000,'-r')


Hn = db.MMHN.high + 3*db.HNO3L.high;
Ln = db.MMHN.low + 3*db.HNO3L.low;
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    2\n',Hn(1),Hn(2),Hn(3),Hn(4),Hn(5));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    3\n',Hn(6),Hn(7),Ln(1),Ln(2),Ln(3));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E                   4\n\n',Ln(4),Ln(5),Ln(6),Ln(7));