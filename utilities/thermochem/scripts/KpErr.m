function err = KpErr(Ls,L_MMH,L_HNO3,Ts,Kps,dH_target)

    R = 8.314; %J/mol-K
    
    L = [Ls(1) 0 0 0 0 Ls(2) Ls(3)];
    
    pMMH_h = [L_MMH(5:-1:1)./(5:-1:1), L_MMH(6)];
    pNA_h = [L_HNO3(5:-1:1)./(5:-1:1), L_HNO3(6)];
    pP_h = [L(5:-1:1)./(5:-1:1), L(6)];
    
    pMMH_s = [L_MMH(5:-1:2)./(4:-1:1), L_MMH(7)];
    pNA_s = [L_HNO3(5:-1:2)./(4:-1:1), L_HNO3(7)];
    pP_s = [L(5:-1:2)./(4:-1:1), L(7)];

    
    dH = (polyval(pP_h,Ts) - polyval(pMMH_h,Ts) - polyval(pNA_h,Ts))*R;
    
    dH400 = (polyval(pP_h,400) - polyval(pMMH_h,400) - polyval(pNA_h,400))*R;
    
    dS = (pP_s(1)*log(Ts) + polyval(pP_s,Ts) - ...
          pMMH_s(1)*log(Ts) - polyval(pMMH_s,Ts) - ...
          pNA_s(1)*log(Ts) - polyval(pNA_s,Ts))*R;
    
    dG = dH - Ts.*dS;
    
    Kp = exp(-dG./(R.*Ts));
    
    err = sum(abs(Kp-Kps)) + abs(dH_target + dH400./1000)/100;
    %fprintf('Errs = %f, %f\n',sum(abs(Kp-Kps)),abs(dH_target + dH400./1000)/100);