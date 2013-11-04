function dG = CalcDG(rxn,T,thermo)

    dG = zeros(size(T));

    Nr = length(rxn.reactants);
    Np = length(rxn.products);
    
    for r = 1:Nr
        dG = dG - G(thermo.(rxn.reactants(r).name),T)*rxn.reactants(r).e;
    end
    
    for p = 1:Np
        dG = dG + G(thermo.(rxn.products(p).name),T)*rxn.products(p).e;
    end
    