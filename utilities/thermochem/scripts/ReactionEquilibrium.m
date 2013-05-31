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

rxn.reactants(1).name = 'CH3N3';
rxn.reactants(1).e = 1;
rxn.reactants(2).name = 'HONO';
rxn.reactants(2).e = 1;

rxn.products(1).name = 'CH3OH';
rxn.products(1).e = 1;
rxn.products(2).name = 'N2O';
rxn.products(2).e = 2;

T = 280:5:1000;
Kp = CalcKp(rxn,T,db);
dH = CalcHeatRelease(rxn,T,db); %J/mol


figure;
hold all
h = [];
for r = 1:length(rxn.reactants)
    hN=plot(T,H(db.(rxn.reactants(r).name),T)./1000 );
    set(hN,'DisplayName',rxn.reactants(r).name);
    h = [h; hN];
end
ylabel('dH (kJ/mol)')

for p = 1:length(rxn.products)
    hN=plot(T,H(db.(rxn.products(p).name),T)./1000 );
    set(hN,'DisplayName',rxn.products(p).name);
    h = [h; hN];
end
lh = legend(h);


figure;
plot(T,dH./1000) %kJ/mol

figure;
semilogy(T,Kp)