function h = H(entry,T)

    R = 8314;% / entry.W;

    h = zeros(size(T));
    Tmid = entry.Tmid;

    plc = [entry.low(5:-1:1)./(5:-1:1), entry.low(6)];
    phc = [entry.high(5:-1:1)./(5:-1:1), entry.high(6)];

    h(T<=Tmid) = polyval(plc,T(T<=Tmid))*R;
    h(T>Tmid) = polyval(phc,T(T>Tmid))*R;
