function [gamma_m_1, dgamma_m_1, mSize, mLife] = plotAvalancheAveSize(sizeAv, lifeAv, fitP)
%{
    Plots the avalanche average size as a function of lifetime

    Inputs:
    sizeAv: Avalanche sizes
    
     fitP: A struct containing parameters to fit
        fitP.lc:    Lower cut-off of IEI
        fitP.uc:    Upper cut-off of IEI
        fitP.cLevel: confidence level to use for errors


    Option to fit if we provide cut-offs


    Outputs:
        gamma_m_1: exponent to <S> ~ T^gamma obtained by linear fitting on
                log-log plot
       dgamma_m_1: uncertainty in gamma_m_1
            mSize: average avalanche size for a given life-tims
            mLife: life-times of valanches found

    Written by Joel Hochstetter

%}
    %gamma_m_1
    [mSize, mLife] = avalancheAvSize(sizeAv, lifeAv);
    gamma_m_1  = nan;
    dgamma_m_1 = nan;

    if nargin == 3 
        fitPL = true;
    else
        fitPL = false;
    end
    
    if fitPL
        %add defaults for cut-offs for PL
        if ~isfield(fitP, 'lc')
            fitP.lc = min(mLife);
        end

        if ~isfield(fitP, 'uc')
            fitP.uc = max(mLife);
        end
        
        if ~isfield(fitP, 'cLevel')
            fitP.cLevel = 0.5;
        end           
        
    end

    if min(mLife) > fitP.lc
        fitP.lc = min(mLife);
    end
    
    if max(mLife) < fitP.uc
        fitP.uc = max(mLife);
    end


    
    loglog(mLife, mSize, 'bx')
    hold on;
    
    if fitPL
        %only include bins within include range to fit
        fitLives = mLife((mLife >= fitP.lc) & (mLife <= fitP.uc));
        fitSizes = mSize((mLife >= fitP.lc) & (mLife <= fitP.uc));         

        if fitP.lc >= fitP.uc || numel(fitLives) <= 1 || numel(fitSizes) <= 1
            return 
        end
        
        %fit power law
        [beta, dbeta] = fitPowerLawLinearLogLog(fitLives, fitSizes);    
        x =  fitP.lc:fitP.lc/5: fitP.uc;
        A = fitSizes(find(fitLives >= fitP.lc, 1));
        y = A*(x/min(x)).^(-beta);
        loglog(x, y, 'r--');        
        
        text(fitLives(1), fitSizes(1)/3, strcat('T^{', num2str(beta,3),'}'), 'Color','b')
  

        gamma_m_1 = -beta;
        dgamma_m_1 =dbeta;


    end

    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('T (bins)')
    ylabel('<S>(T)')



end