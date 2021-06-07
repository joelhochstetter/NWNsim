function plotIEI(G, thr, t, pn, fitP)
%{
    Plots the distribution of inter-event interval
    Inputs:
        G: Conductance time-series
      thr: Threshold in conductance to define an event
       t: time-vector
       pn: if 1 plots +ve and -ve separately on the same axis. Else we plot
            them together

     fitP: A struct containing parameters to fit
        fitP.lc:    Lower cut-off of IEI
        fitP.uc:    Upper cut-off of IEI
        fitP.lcn:    Lower cut-off of IEI for neg event fit
        fitP.ucn:    Upper cut-off of IEI for neg event fit
        fitP.toInc:  Vector of same length as G. Tells which time-points to
                        include
        fitP.logBins: binary depending on whether or not log bins are used

    Option to fit if we provide cut-offs

%}
    dt = t(2) - t(1);
    
    if nargin == 4
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
        
        if ~isfield(fitP, 'logBins')
            fitP.logBins = true;
        end
    end

    
    
    if pn
       
        ddG = dG > thr;
        [~, ieiDat] = IEI(ddG, 1, t);
        if fitP.logBins
           [~, Niei, edgesiei] = LogBin(ieiDat);
        else
            [Niei,edgesiei] = histcounts(ieiDat, 'Normalization', 'probability');            
        end
        
        loglog((edgesiei(1:end-1) + edgesiei(2:end))/2,Niei, 'bx')
        hold on;
        
        %negative fluctuations        
        ddG = dG < thr;
        [~, ieiDat] = IEI(ddG, 1, t);
        [Niei1,edgesiei1] = histcounts(ieiDat, 'Normalization', 'probability');
        loglog((edgesiei1(1:end-1) + edgesiei1(2:end))/2, Niei1, 'rx')        
        legend('\Delta G > 0', '\Delta G < 0')

        
        if fitPL
            %only include bins within include range to fit
            fitEdges = edgesiei((edgesiei >= fitP.lc) & (edgesiei <= fitP.uc));
            cutFront = numel(edgesiei(edgesiei < fitP.lc));
            cutEnd   = numel(edgesiei(edgesiei > fitP.uc));
            edgeCen  = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
            fitN     = Niei(1 + cutFront : end - cutEnd);
         
            %fit power law
            [fitresult, xData, yData, gof] = fitPowerLaw(edgeCen , fitN );    
            plot(fitresult, 'b--', xData, yData, 'gx')
            
            %repeat for dG < 0
            fitEdges = edgesiei1((edgesiei1 >= fitP.lcn) & (edgesiei1 <= fitP.ucn));
            cutFront = numel(edgesiei1(edgesiei1 < fitP.lcn));
            cutEnd   = numel(edgesiei1(edgesiei1 > fitP.ucn));
            edgeCen1 = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
            fitN1    = Niei1(1 + cutFront : end - cutEnd);

            [fitresult1, xData, yData, gof] = fitPowerLaw(edgeCen1, fitN1);
            plot(fitresult1, 'r--', xData, yData, 'mx')
            legend('\Delta G > 0 not fit', '\Delta G < 0 not fit', ...
                '\Delta G > 0 inc fit', '\Delta G > 0 fit', ...
                '\Delta G < 0 inc fit', '\Delta G < 0 fit')
            text(edgeCen(1) , fitN(1)/3 , strcat('t^{-', num2str(-fitresult.b ,3),'}'), 'Color','b') 
            text(edgeCen1(1), fitN1(1)/3, strcat('t^{-', num2str(-fitresult1.b,3),'}'), 'Color','r')

        else
            legend('\Delta G > 0', '\Delta G < 0')
            title('IEI distribution')

        end
       
    else %~pn
        
        ddG = abs(dG) > thr;
        [~, ieiDat] = IEI(ddG, 1, t);
        if fitP.logBins
           [~, Niei, edgesiei] = LogBin(ieiDat);
        else
            [Niei,edgesiei] = histcounts(ieiDat, 'Normalization', 'probability');            
        end
        loglog((edgesiei(1:end-1) + edgesiei(2:end))/2,Niei, 'bx')
        hold on;
        
        if fitPL
            %only include bins within include range to fit
            fitEdges = edgesiei((edgesiei >= fitP.lc) & (edgesiei <= fitP.uc));
            cutFront = numel(edgesiei(edgesiei < fitP.lc));
            cutEnd   = numel(edgesiei(edgesiei > fitP.uc));
            edgeCen  = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
            fitN     = Niei(1 + cutFront : end - cutEnd);         
            %fit power law
            [fitresult, xData, yData, gof] = fitPowerLaw(edgeCen , fitN );    
            plot(fitresult, 'b--', xData, yData, 'gx')

            text(edgeCen(1), fitN(1)/3, strcat('t^{-', num2str(-fitresult.b,3),'}'), 'Color','b')
            legend('not fit', 'inc fit', 'fit')        
        end
        
    end
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('IEI (s)')
    ylabel('P(IEI)') 


end