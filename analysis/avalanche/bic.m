function bicval = bic(k, logL, n)
%{
    Bayersian Information Criterion:
        k: number of parameters in model
        logL: log likelihood of model
        n:   Number of datapoints
%}
   bicval = k*log(n)-2*logL;
end