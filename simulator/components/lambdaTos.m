function s = lambdaTos(lambda, critFlux, smax)
% Calculates filament gap (s)
%Inputs:
%   lambda (filament state of junction)
%   critFlux (lambda crit in model)
%   smax (maximum filament gap distance i
%Written by Joel Hochstetter.

    s = (critFlux-abs(lambda))*smax/critFlux;
    s(s < 0) = 0;

end