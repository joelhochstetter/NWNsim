function dGfit = plotDeltaG(G, pn, fitP, joinperiod)
%{
    Plots the distribution of DeltaG
    Inputs:
        G: Conductance time-series
       pn: if 1 plots +ve and -ve separately on the same axis. Else we plot
            them together
     fitP: A struct containing parameters to fit
        fitP.lc:    Lower cut-off of dG in units of G
        fitP.uc:    Upper cut-off of dG in units of G
        fitP.lcn:    Lower cut-off of dG in units of G for dG < 0 fit
        fitP.ucn:    Upper cut-off of dG in units of G for dG < 0 fit
        fitP.toInc:  Vector of same length as G. Tells which time-points to
                        include
        fitP.useML:  Uses Clauset 2007 maximum likelihood to get power law
        exponent of tail.
        joinperiod: for time series which join an ensemble of
        different simulations stores periodicity so ignores events calculated 
        between adjacent simulations

    Option to fit if we provide cut-offs

%}

    if nargin < 4
        joinperiod = -1;
    end


    dGfit = struct();

    alpha = nan;
    dalph = 0.0;
    lc = 0;

    if nargin < 3
        fitPL = 0;
    else
        fitPL = 1;
    end
    
    dG = [diff(G), 0];
    dG(isnan(dG)) = 0;        
    
    if numel(joinperiod) == 1
        if joinperiod > 0
            dG(joinperiod:joinperiod:numel(dG)) = 0;
        end
    else
        dG(joinperiod) = 0;
    end
    
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
            fitP.uc = 10*max(abs(dG));
        end    
        
        if ~isfield(fitP, 'lcn')
            fitP.lcn = 0;
        end

        if ~isfield(fitP, 'ucn')
            fitP.ucn = 10*max(abs(dG));
        end    
        
        if ~isfield(fitP, 'cLevel')
            fitP.cLevel = 0.95;
        end    
        
        if ~isfield(fitP, 'useML')
            fitP.useML = false;
        end             
    end


    if pn
        [N,edges] = histcounts(abs(dG(dG > 0)), 'Normalization', 'probability');
        loglog((edges(1:end-1) + edges(2:end))/2,N, 'bx')
        hold on;
        [N1,edges1] = histcounts(abs(dG(dG < 0)), 'Normalization', 'probability');
        loglog((edges1(1:end-1) + edges1(2:end))/2,N1, 'rx')
        legend('\Delta G > 0', '\Delta G < 0')

        
        if fitPL
            %only include bins within include range to fit
            fitEdges = edges((edges >= fitP.lc) & (edges <= fitP.uc));
            cutFront = numel(edges(edges < fitP.lc));
            cutEnd   = numel(edges(edges > fitP.uc));
            edgeCen  = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
            fitN     = N(1 + cutFront : end - cutEnd);
         
            %fit power law
            [fitresult, xData, yData, gof] = fitPowerLaw(edgeCen , fitN );    
            plot(fitresult, 'b--', xData, yData, 'gx')
            
            %repeat for dG < 0
            fitEdges = edges1((edges1 >= fitP.lcn) & (edges1 <= fitP.ucn));
            cutFront = numel(edges1(edges1 < fitP.lcn));
            cutEnd   = numel(edges1(edges1 > fitP.ucn));
            edgeCen1 = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
            fitN1    = N1(1 + cutFront : end - cutEnd);

            [fitresult1, xData, yData, gof] = fitPowerLaw(edgeCen1, fitN1);
            plot(fitresult1, 'r--', xData, yData, 'mx')
            legend('\Delta G > 0 not fit', '\Delta G < 0 not fit', ...
                '\Delta G > 0 inc fit', '\Delta G > 0 fit', ...
                '\Delta G < 0 inc fit', '\Delta G < 0 fit')
            text(edgeCen(1) , fitN(1)/3 , strcat('\Delta G^{-', num2str(-fitresult.b ,3),'}'), 'Color','b') 
            text(edgeCen1(1), fitN1(1)/3, strcat('\Delta G^{-', num2str(-fitresult1.b,3),'}'), 'Color','r')
%             title(strcat('\Delta G > 0: \Delta G^{-', ...
%                 num2str(-fitresult.b,3),'}', ...
%                 ',  \Delta G < 0: \Delta G^{-', num2str(-fitresult1.b,3),'}'))
            alpha = -fitresult.b;
            CI  = confint(fitresult, fitP.cLevel);
            tCI = CI(:,2);
            dalph = (tCI(2) - tCI(1))/2;

        else
            legend('\Delta G > 0', '\Delta G < 0')
            title('\Delta G distribution')

        end
       
    else %~pn
        [N,edges] = histcounts(abs(dG), 'Normalization', 'probability');
        loglog((edges(1:end-1) + edges(2:end))/2,N, 'bx')
        hold on;
        title('\Delta G distribution')
        xlabel('\Delta G')
        ylabel('P(\Delta G)') 
        
        dG = abs(dG);
        if fitPL
            
            if fitP.useML 
                [alpha,   lc, ~] = plfit(dG);
                [dalph, dlc, ~]  = plvar(dG,'xmin', lc);
                
                x = lc:lc/1000:max(dG);
                A = N(find(edges <= lc, 1));
                x1 = dG(find(edges <= lc, 1)); 
                y = A*(x/x1).^(-alpha);
                loglog(x, y, 'r--');
                text(x(2), y(2)*2, strcat('\Delta G^{-', num2str(alpha,3),'}'), 'Color','r')
                
                            %put goodness of fit in as well

                
            else
                %only include bins within include range to fit
                fitEdges = edges((edges >= fitP.lc) & (edges <= fitP.uc));
                cutFront = numel(edges(edges < fitP.lc));
                cutEnd   = numel(edges(edges > fitP.uc));
                edgeCen  = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
                fitN     = N(1 + cutFront : end - cutEnd);
                
                %fit power law
               [alpha, dalph] = fitPowerLawLinearLogLog(edgeCen, fitN);    
                
                if numel(edgeCen) == 0
                    disp('dG not enough points')
                    return 
                end
               
                lc = edgeCen(1); 
                
                x = lc:lc/100:max(edgeCen);
                A = fitN(find(edges >= lc, 1));
                y = A*(x/x(1)).^(-alpha);
                loglog(x, y, 'r--');
                text(x(1), y(1)*2, strcat('\Delta G^{-', num2str(alpha,3),'}'), 'Color','r')

            end
            
        end


        
    end
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('\Delta G')
    ylabel('P(\Delta G)') 
    
    bins = (edges(1:end-1) + edges(2:end))/2;
    
    dGfit.alpha  = alpha;
    dGfit.dalph  = dalph;
    dGfit.lc     = lc;
    dGfit.bins = bins;
    dGfit.prob = N;
    
end