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

    rxn.reactants(1).name = 'CH3NHNH2';
    rxn.reactants(1).e = 1;
    rxn.reactants(2).name = 'HNO3';
    rxn.reactants(2).e = 3;


    rxn.products(1).name = 'CH3ONO2';
    rxn.products(1).e = 1;
    rxn.products(2).name = 'HONO';
    rxn.products(2).e = 2*(1-x(i));
    rxn.products(3).name = 'OH';
    rxn.products(3).e = 2*x(i);
    rxn.products(4).name = 'NO';
    rxn.products(4).e = 2*x(i);
    rxn.products(5).name = 'H2O';
    rxn.products(5).e = 2;
    rxn.products(6).name = 'N2';
    rxn.products(6).e = 1;


    T = 280:5:1000;
    Kp = CalcKp(rxn,T,db);
    dH = CalcHeatRelease(rxn,T,db); %J/mol
    dH400(i) = dH(T==400);
    Nr(i) = 2*x(i)/(sum([rxn.products.e])-1-2*(1-x(i))-2);
end

figure
subplot(1,2,1)
plot(x,dH400./1000)
hold on
plot([0,1],[0,0],'--k')
xlabel('x')
xlim([0 1])
ylabel('dH (kJ/mol)')

subplot(1,2,2)
plot(x,Nr)
xlabel('x')
ylabel('N_{rad} / N_{prod}')

% 
% figure;
% hold all
% h = [];
% for r = 1:length(rxn.reactants)
%     hN=plot(T,H(db.(rxn.reactants(r).name),T)./1000 );
%     set(hN,'DisplayName',rxn.reactants(r).name);
%     h = [h; hN];
% end
% ylabel('dH (kJ/mol)')
% 
% for p = 1:length(rxn.products)
%     hN=plot(T,H(db.(rxn.products(p).name),T)./1000 );
%     set(hN,'DisplayName',rxn.products(p).name);
%     h = [h; hN];
% end
% lh = legend(h);
% 
% 
% figure;
% plot(T,dH./1000) %kJ/mol
% 
% figure;
% semilogy(T,Kp)