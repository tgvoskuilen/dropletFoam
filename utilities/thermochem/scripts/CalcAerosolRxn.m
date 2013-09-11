clear all
close all
clc


% Load thermodynamic database
db = ReadDb();

T = 280:5:1000;

Ha = -0;
Sa = Ha/590;

dH = H(db.MMHN,T) - H(db.HNO3,T) - H(db.CH3NHNH2,T) + Ha;
dS = S(db.MMHN,T) - S(db.HNO3,T) - S(db.CH3NHNH2,T) + Sa;
dG = dH - T.*dS;

dS2 = S(db.CH4,T) + S(db.N2,T) + S(db.H2,T) - S(db.CH3NHNH2,T);

figure
plot(T,dS2)

Kp = exp(-dG./(8.314.*T));
R = 1 - 1.01325e5./(8314.*T)./Kp;

k0 = 1.42e7;
beta = 1.861;
Ea = 5200;

kf = k0.*T.^beta.*exp(-Ea./(1.9872041*T));
kfi = 1.42e5.*T.^1.861;

K = kf.*R;
Ki = kfi.*R;

figure;
subplot(1,3,1)
plot(T,dG./1000)
subplot(1,3,2)
plot(T,dH./1000)
subplot(1,3,3)
plot(T,R)
ylim([-1 1])

figure;
subplot(1,2,1)
plot(T,dS)
subplot(1,2,2)
plot(T,dH./1000)


figure;
hold on
h(1) = plot(T,K,'LineWidth',2,'DisplayName','k_f - k_r');
h(2) = plot(T,kf,'--r','LineWidth',2,'DisplayName','k_f');
h(3) = plot(T,kfi,'--g','LineWidth',2,'DisplayName','Old k_f');
lh = legend(h,'Location','NorthWest');
set(lh,'Box','off');
ylim([0 max([max(K) max(Ki)])*1.2])
set(gca,'Box','on')
xlabel('T (K)')
ylabel('k_f - k_r')



db.MMHNcustom = db.MMHN;

db.MMHNcustom.low(6) = db.MMHNcustom.low(6) + Ha/8.314;
db.MMHNcustom.low(7) = db.MMHNcustom.low(7) + Sa/8.314;
db.MMHNcustom.high(6) = db.MMHNcustom.high(6) + Ha/8.314;
db.MMHNcustom.high(7) = db.MMHNcustom.high(7) + Sa/8.314;



rxn.reactants(1).name = 'CH3NHNH2';
rxn.reactants(1).e = 1;
rxn.reactants(2).name = 'HNO3';
rxn.reactants(2).e = 1;


rxn.products(1).name = 'MMHNcustom';
rxn.products(1).e = 1;

Kp = CalcKp(rxn,T,db);
dH = CalcHeatRelease(rxn,T,db); %J/mol

% figure;
% subplot(1,2,1)
% plot(T,dH./1000)
% 
% subplot(1,2,2)
% plot(T,1-1.01325e5./(8314.*T)./Kp)
% ylim([-1 1])


Ln = db.MMHNcustom.low;
Hn = db.MMHNcustom.high;


fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    2\n',Hn(1),Hn(2),Hn(3),Hn(4),Hn(5));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    3\n',Hn(6),Hn(7),Ln(1),Ln(2),Ln(3));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E                   4\n\n',Ln(4),Ln(5),Ln(6),Ln(7));

