function omega = CalcMaxLiquidRate(rxn,T)
    R = 1.9858775; %cal/mol-K
    
    kf = rxn.A .* T.^rxn.B .* exp(-rxn.Ea./(R.*T));
        
    Nr = length(rxn.reactants);
    
    %assume there is an equal number of moles of all reactants present
    cp = (rxn.reactants(1).rho/Nr) .^ rxn.reactants(1).e;
    
    for r = 2:Nr
        cp = cp .* (rxn.reactants(1).rho/Nr) .^ rxn.reactants(r).e;
    end
    
    omega = kf .* cp; %mol/cm3/s