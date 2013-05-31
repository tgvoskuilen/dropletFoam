function s = S(entry,T)

    R = 8.314; %J/mol-K

    Tmid = entry.Tmid;

    s = zeros(size(T));

    plc = [entry.low(5:-1:2)./(4:-1:1), entry.low(7)];
    phc = [entry.high(5:-1:2)./(4:-1:1), entry.high(7)];

    s(T<=Tmid) = (entry.low(1)*log(T(T<=Tmid)) + ...
        polyval(plc,T(T<=Tmid)))*R;
    s(T>Tmid) = (entry.high(1)*log(T(T>Tmid)) + ...
        polyval(phc,T(T>Tmid)))*R;
