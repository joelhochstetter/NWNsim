function [conductance] = tunnelSwitchL(d, phi, A, Coff, Con)
%{
    Applying Simmons low voltage tunnelling formula
    (https://doi.org/10.1063/1.1702682) to calculate conductance
    of junctions.

    Inputs:
        d: filament gap distance in nm
      phi: filament barrier height in V
        A: amplitude constant
     Coff: off-conductance of junction (called residual in paper)
      Con: on-conductance of junction (called filament in paper)

    Outputs:
        conductance (in units of S): conductance

    Written by Joel Hochstetter

%}


    C0 = 10.19;
    J1 = 0.0000471307;
   
    tun = 2/A*d/phi^0.5.*exp(C0.*d*phi^2)/J1;

    conductance = 1./(tun+1./Con) + Coff;
    
end
