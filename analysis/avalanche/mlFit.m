function Fit = mlFit(x, fitTrun)
%{
    Uses a maximum likelihood method to compare the fits between:
    For functional form see Marshall et. all 2016
    - power law (PL): uses doubly truncated if PL is true
    - Exponential (E)
    - Log-normal (LN)
    - Weibull (W)
    
    Inputs:
        x: data 
 fitTrunc: true (fit truncated power law)

    Requires:
        - MATLAB Stats Toolbox
        - Custom functions: aic, bic, plllike
        - NCC Toolbox (Marshall (2016) doi: 10.3389/fphys.2016.00250)

    Written by: Joel Hochstetter
%}

    Fit     = struct();
    N = numel(x);
    
    %% Power law
    if fitTrun %truncated
        [tau, xmin, xmax, dtau, p, pcrit, ks, fullResults] = plparams(x); 
        x = x((x >= xmin) & (x <= xmax));
        N = numel(x);
        llike = plllike(x, tau, xmin, xmax);
        AIC  = aic(llike, 1);
        BIC  = bic(llike, 1, N);       
        Fit.PL  = struct('tau', tau, 'xmin', xmin, 'xmax', xmax, ...
            'llike', llike, 'aic', AIC, 'bic', BIC, 'dtau', dtau, 'p', p, ...
            'pcrit', pcrit, 'ks', ks, 'fullResults', fullResults);  
    else %not truncated
        tau  = plmle(x);
        xmin = min(x);
        xmax = max(x);
        llike = plllike(x, tau, xmin, xmax);
        AIC  = aic(llike, 1);
        BIC  = bic(llike, 1, N);    
        Fit.PL  = struct('tau', tau, 'xmin', xmin, 'xmax', xmax, ...
            'llike', llike, 'aic', AIC, 'bic', BIC, 'dtau', 0.0, ...
            'fullResults', struct());
    end
        
    %% Exponential
    [lambda, lambdaci] = expfit(x);
    dlambda = (lambdaci(2) - lambdaci(1))/2;
    llike   = -explike(lambda,x);
    AIC  = aic(llike, 1);
    BIC  = bic(llike, 1, N);     
    Fit.E   = struct('lambda', lambda, 'dlambda', dlambda, 'llike', llike, ...
        'aic', AIC, 'bic', BIC);   
    
    %% Lognormal
    [parmhat, parmci] = lognfit(x, [],[],[], statset('MaxIter',10000, 'MaxFunEvals', 10000));
    llike   = -lognlike(parmhat, x);
    AIC  = aic(llike, 2);
    BIC  = bic(llike, 2, N);     
    Fit.LN  = struct('mu', parmhat(1), 'sigma', parmhat(2), 'dmu', ...
        parmci(1,:),  'dsigma', parmci(2,:),'llike', llike, 'aic', AIC, 'bic', BIC);    
    
    %% Weibull
    [parmhat, parmci] = wblfit(x, [],[],[], statset('MaxIter',10000, 'MaxFunEvals', 10000));
    llike   = -wbllike(parmhat, x);
    AIC  = aic(llike, 2);
    BIC  = bic(llike, 2, N);     
    Fit.WB  = struct('nu', parmhat(1), 'gamma', parmhat(2), 'dnu', parmci(1,:), ...
        'dgamma', parmci(2,:),'llike', llike, 'aic', AIC, 'bic', BIC);    
    
end