function g = G(entry,T)

    R = 8.314; %J/mol-K

    h = zeros(size(T));
    Tmid = entry.Tmid;

    plc = [entry.low(5:-1:1)./(5:-1:1), entry.low(6)];
    phc = [entry.high(5:-1:1)./(5:-1:1), entry.high(6)];

    h(T<=Tmid) = polyval(plc,T(T<=Tmid))*R;
    h(T>Tmid) = polyval(phc,T(T>Tmid))*R;
 
    s = zeros(size(T));

    plc = [entry.low(5:-1:2)./(4:-1:1), entry.low(7)];
    phc = [entry.high(5:-1:2)./(4:-1:1), entry.high(7)];

    s(T<=Tmid) = (entry.low(1)*log(T(T<=Tmid)) + ...
        polyval(plc,T(T<=Tmid)))*R;
    s(T>Tmid) = (entry.high(1)*log(T(T>Tmid)) + ...
        polyval(phc,T(T>Tmid)))*R;
    
    g = h - T.*s;