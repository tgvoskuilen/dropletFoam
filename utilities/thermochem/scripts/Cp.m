function cp = Cp(entry,T)

    R = 8.314; %J/mol-K
    
    cp = zeros(size(T));
    Tmid = entry.Tmid;

    cp(T<=Tmid) = polyval(entry.low(5:-1:1),T(T<=Tmid))*R;
    cp(T>Tmid) = polyval(entry.high(5:-1:1),T(T>Tmid))*R;