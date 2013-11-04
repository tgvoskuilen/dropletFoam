clear all
close all
clc


T = 280:5:1000;

% Load thermodynamic database
db = ReadDb();

% Read reactions from CHEMKIN input file
reactions = ReadCHEMKINReactions('../chem/chem41.inp');

for r = 1:length(reactions)
    reactions(r).omega = CalcMaxRate(reactions(r),T,1);
    reactions(r).dH = CalcHeatRelease(reactions(r),T,db);
    reactions(r).Sh = reactions(r).omega.*reactions(r).dH;
    reactions(r).Kp = CalcKp(reactions(r),T,db);
end


figure;
hold all
h = [];
for r = 1:length(reactions)
    rs = {reactions(r).reactants.name};
    ps = {reactions(r).products.name};
    
    if (any(strcmpi(rs,'NO')) || any(strcmpi(ps,'NO'))) && max(reactions(r).Kp) < 100
        hN = plot(T,reactions(r).Kp);
        h = [h; hN];
        set(hN,'DisplayName',reactions(r).name)
    end
end
lh = legend(h,'Location','BestOutside');
set(gca,'YScale','log')
