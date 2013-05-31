function Kp = CalcKp(rxn,T,thermo)

    Rsi = 8.314; %J/mol/K
    
    dG = zeros(size(T));

    Nr = length(rxn.reactants);
    Np = length(rxn.products);
    
    for r = 1:Nr
        dG = dG + G(thermo.(rxn.reactants(r).name),T)*rxn.reactants(r).e;
    end
    
    for p = 1:Np
        dG = dG - G(thermo.(rxn.products(p).name),T)*rxn.products(p).e;
    end
    
    Kp = exp(-dG./(Rsi.*T));