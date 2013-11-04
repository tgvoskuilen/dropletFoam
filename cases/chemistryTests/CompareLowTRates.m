clear all
close all
clc

% Load and compare reaction mechanisms from Chemkin files

Ns = 41;   %select 25, 29, or 41

chemfile = strcat('chemkin',filesep,'chem',num2str(Ns,'%2d'),'.inp');
thermofile = strcat('chemkin',filesep,'therm',num2str(Ns,'%2d'),'.dat');

T = 200:5:2500;
reactions = ReadCHEMKINReactions(chemfile);
species = ReadCHEMKINThermo(thermofile,T);

Ns =  length(fieldnames(species));
Nr = length(reactions);

T = 280:5:500;
R = 1.9858775; %cal/mol-K
Rsi = 8.314; %J/mol-K
P = 101325; %Pa
M = P./(Rsi.*T)./100^3;%mol/cm3

for r = 1:Nr
    kInf = reactions(r).A .* T.^reactions(r).B .* exp(-reactions(r).Ea./(R.*T));
    
    if ~isempty(reactions(r).low)
        kLow =  reactions(r).low.A .* T.^reactions(r).low.B .* exp(-reactions(r).low.Ea./(R.*T));

        Pr = kLow/kInf*M;
        
        if ~isempty(reactions(r).troe)
            tr = reactions(r).troe;
            Fcent = (1-tr.alpha).*exp(-T./tr.T3)+tr.alpha.*exp(-T./tr.T1)+exp(-tr.T2./T);
            d = 0.14;
            n = 0.75-1.27*log(Fcent);
            c = -0.4-0.67*log(Fcent);
            F = exp(log(Fcent)/(1+((log(Pr)+c)./(n-d.*(log(Pr)+c))).^2));
        else
            F = 1;
        end
        
        reactions(r).kf = kInf*(Pr/(1+Pr))*F;
    else
        reactions(r).kf = kInf;
    end
    
end

n = 0;
figure;
hold all
for r = 1:Nr
    rxn = reactions(r);

    pn = {rxn.products.name};
    rn = {rxn.reactants.name};

    if any(strcmp('HONO', pn)) || any(strcmp('HONO', rn))
        semilogy(T,reactions(r).kf)
        n = n + 1;
        fprintf('%d: %s  ->  %5.4e\n',n,reactions(r).name, max(reactions(r).kf));
    end
end
set(gca,'YScale','log')

fprintf('reduced from %d to %d reactions\n',Nr,n)