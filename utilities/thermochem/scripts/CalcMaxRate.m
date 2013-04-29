function omega = CalcMaxRate(rxn,T,p)
    %p in bar
    R = 1.9858775; %cal/mol-K
    Rsi = 8.314; %J/mol-K
    
    kInf = rxn.A .* T.^rxn.B .* exp(-rxn.Ea./(R.*T));
    
    M = p*1e5./(Rsi.*T)/100^3; %mol/cm3
    
    %ignore third body efficiencies for now
    if ~isempty(rxn.low)
        kLow =  rxn.low.A .* T.^rxn.low.B .* exp(-rxn.low.Ea./(R.*T));

        Pr = kLow/kInf*M;
        
        if ~isempty(rxn.troe)
            tr = rxn.troe;
            Fcent = (1-tr.alpha).*exp(-T./tr.T3)+tr.alpha.*exp(-T./tr.T1)+exp(-tr.T2./T);
            d = 0.14;
            n = 0.75-1.27*log(Fcent);
            c = -0.4-0.67*log(Fcent);
            F = exp(log(Fcent)/(1+((log(Pr)+c)./(n-d.*(log(Pr)+c))).^2));
        else
            F = 1;
        end
        
        kf = kInf*(Pr/(1+Pr))*F;
    else
        kf = kInf;
    end
    
    Nr = length(rxn.reactants);
    
    %assume there is an equal number of moles of all reactants present
    cp = (M./Nr) .^ rxn.reactants(1).e;
    
    for r = 2:Nr
        cp = cp .* (M./Nr) .^ rxn.reactants(r).e;
    end
    
    omega = kf .* cp; %mol/cm3/s