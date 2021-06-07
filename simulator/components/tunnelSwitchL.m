function [con] = tunnelSwitchL(d, phi, A, Coff, Con)
    C0 = 10.19;
    J1 = 0.0000471307;
   
    tun = 2/A*d/phi^0.5.*exp(C0.*d*phi^2)/J1;

    con = 1./(tun+1./Con) + Coff;
    
end
