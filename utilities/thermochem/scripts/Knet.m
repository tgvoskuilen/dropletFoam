function RR = Knet(R,P,T)

    pR_h = [R(5:-1:1)./(5:-1:1), R(6)];
    pP_h = [P(5:-1:1)./(5:-1:1), P(6)];
    
    pR_s = [R(5:-1:2)./(4:-1:1), R(7)];
    pP_s = [P(5:-1:2)./(4:-1:1), P(7)];

    
    dH = (polyval(pP_h,T) - polyval(pR_h,T))*8314;
    dS = (pP_s(1)*log(T) + polyval(pP_s,T) - ...
          pR_s(1)*log(T) - polyval(pR_s,T))*8314;
      
    dG = dH - T.*dS;
    
    Kp = exp(-dG./(8314.*T));
    
    Kc = Kp .* (P0/(8314.*T)).^-1;
    
    RR = 1 - 1./Kc;