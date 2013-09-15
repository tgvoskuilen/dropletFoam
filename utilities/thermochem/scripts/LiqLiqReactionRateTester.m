function [t,T] = LiqLiqReactionRateTester(db,A,B,Ea)
% clear all
% close all
% clc



% Load thermodynamic database
%db = ReadDb();

% Define reaction
rxn.name = 'CH3NHNH2L+2HNO3L=>MMH2N';

rxn.A = A;%1e12; %mol,cm3,s,K units
rxn.B = B;
rxn.Ea = Ea;%8000; %cal/mol

rxn.low = [];
%rhog = 1e5/(8.314*400)/100^3;

rxn.reactants(1).name = 'CH3NHNH2L';
rxn.reactants(1).e = 1;
rxn.reactants(1).rho = 870/46*1000/100^3; %mol/cm3
rxn.reactants(2).name = 'HNO3L';
rxn.reactants(2).e = 2;
rxn.reactants(2).rho = 1550/57*1000/100^3; %mol/cm3

rxn.products(1).name = 'MMH2N';
rxn.products(1).e = 1;
% rxn.products(2).name = 'CH2O';
% rxn.products(2).e = 1;
% rxn.products(3).name = 'NH3';
% rxn.products(3).e = 1;
% rxn.products(4).name = 'HONO';
% rxn.products(4).e = 2;


% rxn.reactants(1).name = 'CH3NHNH2L';
% rxn.reactants(1).e = 1;
% rxn.reactants(1).rho = 870/46*1000/100^3; %mol/cm3
% rxn.reactants(2).name = 'NO2L';
% rxn.reactants(2).e = 1;
% rxn.reactants(2).rho = 1550/57*1000/100^3; %mol/cm3
% 
% rxn.products(1).name = 'CH3NHNH';
% rxn.products(1).e = 1;
% rxn.products(2).name = 'HONO';
% rxn.products(2).e = 1;


% Create energy equation ODE handles
omega = @(T) CalcMaxLiquidRate(rxn,T); %mol/cm3/s
dH = @(T) CalcRxndH(rxn,T,db); %J/mol
rhoCp = @(T) 0.5*rxn.reactants(1).rho*Cp(db.(rxn.reactants(1).name),T) ...
    + 0.5*rxn.reactants(2).rho*Cp(db.(rxn.reactants(2).name),T); %J/cm3/K

dTdt = @(t,T) omega(T).*dH(T)./rhoCp(T);

%TdH = 300:5:700;
% figure;
% plot(TdH, dH(TdH)./1000)
% ylabel('kJ/mol')

% solve
tend = 3e-3; %2 ms
Ti = 300; %K
[t,T] = ode15s(dTdt,[0 tend],Ti);

% find the 350K point
% Tb = 360;
% ib = find(T>Tb,1,'first');
% if isempty(ib)
%     ib = length(T);
% end
% tb = interp1(T(ib-2:ib),t(ib-2:ib),Tb);
% fprintf('tb = %f us\n',tb*1e6);
% 
% %plot
% figure;
% plot(t*1e6, T)
% hold on
% plot(tb*1e6,Tb,'ok','MarkerFaceColor','r')
% set(gca,'XScale','log')
% xlabel('t (\mus)')
% ylabel('T (K)')
% ylim([300 500])
