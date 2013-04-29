clear all
close all
clc


T = 280:5:1000;

% Load thermodynamic database
db = ReadDb();

% Read reactions from CHEMKIN input file
reactions = ReadCHEMKINReactions('chem41.inp');

for r = 1:length(reactions)
    reactions(r).omega = CalcMaxRate(reactions(r),T,1);
    reactions(r).Sh = CalcHeatRelease(reactions(r),T,db);
end


figure;
hold all
h = [];
for r = 1:length(reactions)
    rs = {reactions(r).reactants.name};
    ps = {reactions(r).products.name};
    
    if any(strcmpi(rs,'HONO')) || any(strcmpi(ps,'HONO'))
        hN = plot(T,reactions(r).Sh);
        h = [h; hN];
        set(hN,'DisplayName',reactions(r).name)
    end
end
lh = legend(h,'Location','BestOutside');
%set(gca,'YScale','log')

return


dH1 = H(db.CH3NHNH2,T) + H(db.HNO3,T) - H(db.HNCO,T) - 2*H(db.H2O,T) - H(db.N2H2,T);
dH2 = H(db.HNCO,T) + H(db.HNO3,T) - H(db.CO2,T) - H(db.H2O,T) - H(db.N2O,T);
dH3 = H(db.CH3NHNH2,T) + 2*H(db.HNO3,T) - H(db.CH3ONO2,T) - 0.5*H(db.N2O,T) - 2.5*H(db.H2O,T) - H(db.N2,T);
dH4 = H(db.CH3NHNH2,T) + 2*H(db.HNO3,T) - H(db.CO2,T) - 3*H(db.H2O,T) - H(db.N2O,T) - H(db.N2H2,T);
dH5 = H(db.CH3NHNH2,T) + 3*H(db.HNO3,T) - 2*H(db.H2O,T) - 2*H(db.HONO,T) - H(db.CH3ONO2,T) - H(db.N2,T);
dH6 = H(db.CH3NHNH2,T) + 3*H(db.HNO3,T) - 2*H(db.H2O,T) - 2*H(db.HONO,T) - H(db.CH3O,T) - H(db.NO2,T) - H(db.N2,T);
dH7 = H(db.CH3,T) + H(db.HONO,T) - H(db.CH4,T) - H(db.NO2,T);

figure;
hold all
plot(T,dH1,'-r')
plot(T,dH2,'-k')
plot(T,dH3,'-b')
plot(T,dH4,'-m')
plot(T,dH5,'--c')
plot(T,dH6,'--k')

figure;
hold all
plot(T,dH7,'r')