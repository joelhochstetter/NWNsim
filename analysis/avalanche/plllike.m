function loglike = plllike(xi, alpha, xmin, xmax)
%{
    Log-likehihood of power law:
        xi: data to fit to power law
     alpha: power law exponent
      xmin: lower cut-off
      xmax: upper cut-off
%}
    N = numel(xi);
    x = xmin:xmax;
    loglike = -log(sum(x.^(-alpha))) - alpha/N*sum(log(xi));
end