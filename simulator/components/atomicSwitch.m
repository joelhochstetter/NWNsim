function conductance = atomicSwitch(lambda, lambdaCrit, Con, Coff)
    conductance = zeros(size(lambda));
    conductance(lambda < lambdaCrit) = Coff; 
    conductance(lambda >= lambdaCrit) = Con;
    

end