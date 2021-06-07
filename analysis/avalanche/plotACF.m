function [alpha, dalph] = plotACF(G, dt, useAbs, fitP)
%{
    Plots the distribution of DeltaG
    Inputs:
        G: Conductance time-series
       dt: sampling time of G in s
   useAbs: Use absolute value of dG and the ACF for the distribution plot

     fitP: A struct containing parameters to fit
        fitP.lc:    Lower cut-off of ACF
        fitP.uc:    Upper cut-off of ACF
        fitP.lcn:    Lower cut-off of ACF for neg event fit
        fitP.ucn:    Upper cut-off of ACF for neg event fit
        fitP.toInc:  Vector of same length as G. Tells which time-points to
                        include

    Option to fit if we provide cut-offs

%}

    alpha = nan;
    dalph = nan;

    if nargin == 3
        fitPL = 0;
    else
        fitPL = 1;
    end
    
    dG = [diff(G), 0];
    dG(isnan(dG)) = 0;    
    
    if fitPL
        %if we exclude points then we do it here
        if isfield(fitP, 'toInc')
            dG = dG(fitP.toInc);
        end
        
        %add defaults for cut-offs for PL
        if ~isfield(fitP, 'lc')
            fitP.lc = 0;
        end

        if ~isfield(fitP, 'uc')
            fitP.uc = Inf;
        end    
        
        if ~isfield(fitP, 'lcn')
            fitP.lcn = 0;
        end

        if ~isfield(fitP, 'ucn')
            fitP.ucn = Inf;
        end   
        
        if ~isfield(fitP, 'cLevel')
            fitP.cLevel = 0.95;
        end     
    end


    
    if useAbs
       
        [acf,lags,bounds] = autocorr(abs(dG));
        acf = abs(acf);
        lags = lags*dt;
        %loglog(lags,acf, 'bx','HandleVisibility','off');
        hold on;

        
        if fitPL
            %only include bins within include range to fit
            fitLags  = lags((lags > fitP.lc) & (lags <= fitP.uc));
            cutFront = numel(lags(lags <= fitP.lc));
            cutEnd   = numel(lags(lags > fitP.uc));
            fitACF   = acf(1 + cutFront : end - cutEnd);    
            [fitresult, xData, yData, gof] = fitPowerLaw(fitLags, fitACF);
            plot(fitresult, xData, yData, 'gx')
            text(fitLags(1) , fitACF(1)/3 , strcat('t^{-', num2str(-fitresult.b ,3),'}'), 'Color','b') 
        end
        
        %shuffle data and plot
        shuff_dG = shuffle_data(abs(dG));
        [sacf,lags, bounds] = autocorr(abs(shuff_dG));
        sacf = abs(sacf);
        lags = lags*dt;
        
        if fitPL
            %only include bins within include range to fit
            fitLags1 = lags((lags > fitP.lc) & (lags <= fitP.uc));
            cutFront = numel(lags(lags <= fitP.lc));
            cutEnd   = numel(lags(lags > fitP.uc));
            fitACF   = sacf(1 + cutFront : end - cutEnd);    
            [fitresult, xData, yData, gof] = fitPowerLaw(fitLags1, fitACF);
            plot(fitresult, xData, yData, 'mx')
            text(fitLags1(1) , fitACF(1)/3 , strcat('t^{-', num2str(-fitresult.b ,3),'}'), 'Color','r') 

        else
            legend('data', 'shuffled')

        end
        
        xlabel('Time (dt)')
        ylabel('P(t)')
        title('Auto-correlation function')
       
    else %~pn
        
        [acf,lags,bounds] = autocorr((dG));
        lags = lags*dt;
        
        %loglog(lags,acf, 'bx','HandleVisibility','off');
        hold on;

        
        if fitPL
            %only include bins within include range to fit
            fitLags  = lags((lags > fitP.lc) & (lags <= fitP.uc));
            cutFront = numel(lags(lags <= fitP.lc));
            cutEnd   = numel(lags(lags > fitP.uc));
            fitACF   = acf(1 + cutFront : end - cutEnd);    
            [fitresult, xData, yData, gof] = fitPowerLaw(fitLags, fitACF);
            plot(fitresult, xData, yData, 'gx')
            text(fitLags(1) , fitACF(1)/3 , strcat('t^{-', num2str(-fitresult.b ,3),'}'), 'Color','b') 
        end
        
        %shuffle data and plot
        shuff_dG = shuffle_data(dG);
        [sacf,lags,bounds] = autocorr((shuff_dG));
        lags = lags*dt;
        
        if fitPL
            %only include bins within include range to fit
            fitLags1 = lags((lags > fitP.lc) & (lags <= fitP.uc));
            cutFront = numel(lags(lags <= fitP.lc));
            cutEnd   = numel(lags(lags > fitP.uc));
            fitACF   = sacf(1 + cutFront : end - cutEnd);    
            [fitresult, xData, yData, gof] = fitPowerLaw(fitLags1, fitACF);
            plot(fitresult, xData, yData, 'mx')
            text(fitLags1(1) , fitACF(1)/3 , strcat('t^{-', num2str(-fitresult.b ,3),'}'), 'Color','r') 
            alpha = -fitresult.b;
            CI  = confint(fitresult, fitP.cLevel);
            tCI = CI(:,2);
            dalph = (tCI(2) - tCI(1))/2;
        else
            legend('data', 'shuffled')

        end
        
        xlabel('Time (dt)')
        ylabel('P(t)')
        title('Auto-correlation function')
        
    end
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')


end


