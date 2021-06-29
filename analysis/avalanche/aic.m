function aicval = aic(logL, k)
%{
    Aikike Information Criterion:
        k: number of parameters in model
        logL: log likelihood of model
%}
   aicval = 2*(k-logL);
end