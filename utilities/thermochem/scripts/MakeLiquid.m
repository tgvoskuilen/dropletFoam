function MakeLiquid(name, gas, dH, dS, Tb, dCp)

    R = 8.314;
    
    Lo = gas.low;
    Ho = gas.high;
    
    Ln = Lo;
    Hn = Ho;
    
    %adjust low
    Ln(1) = Lo(1) - dCp/R;
    Ln(6) = Lo(6) - dH*1000/R + dCp/R*Tb;
    if Tb > 0
        Ln(7) = Lo(7) - dS/R + dCp/R*log(Tb);
    else
        Ln(7) = Lo(7) - dS/R; 
    end

    liq.Tlow = gas.Tlow;
    liq.Tmid = gas.Tmid;
    liq.Thigh = gas.Thigh;
    
    %adjust high
    Hn(1) = Ho(1) - dCp/R;
    Hn(6) = Ho(6) - dH*1000/R + dCp/R*Tb;
    if Tb > 0
        Hn(7) = Ho(7) - dS/R + dCp/R*log(Tb);
    else
        Hn(7) = Ho(7) - dS/R;
    end
    
    liq.low = Ln;
    liq.high = Hn;
    T = 300:5:2000;
    
    figure;
    set(gcf,'Units','inches','OuterPosition',[2 2 8 4])
    subplot(1,3,1)
    hold all
    plot(T,Cp(gas,T),'--r')
    plot(T,Cp(liq,T),'-b')
    
    subplot(1,3,2)
    hold all
    plot(T,H(gas,T)./1000,'--r')
    plot(T,H(liq,T)./1000,'-b')
    
    subplot(1,3,3)
    hold all
    plot(T,S(gas,T),'--r')
    plot(T,S(liq,T),'-b')
    
    
    fprintf('CHEMKIN Polynomial Coefficients for %s\n',name);
    fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    2\n',Hn(1),Hn(2),Hn(3),Hn(4),Hn(5));
    fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    3\n',Hn(6),Hn(7),Ln(1),Ln(2),Ln(3));
    fprintf('%+10.8E%+10.8E%+10.8E%+10.8E                   4\n\n',Ln(4),Ln(5),Ln(6),Ln(7));