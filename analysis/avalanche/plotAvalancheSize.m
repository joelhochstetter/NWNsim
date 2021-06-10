function [tau, dta, xmin, xmax, p, pcrit, ks, bins, prob, MLcompare] = plotAvalancheSize(sizeAv, fitP)
%{
    Plots the avalanche size distribution
    Inputs:
    sizeAv: Avalanche sizes
     fitP: A struct containing parameters to fit
        fitP.lc:    Lower cut-off of IEI
        fitP.uc:    Upper cut-off of IEI


    Option to fit if we provide cut-offs

%}
    
    MLcompare = struct();

    if nargin == 1
        fitPL = 0;
    else
        fitPL = 1;
    end
    

    if fitPL
        %add defaults for cut-offs for PL
        if ~isfield(fitP, 'lc')
            fitP.lc = 0;
        end

        if ~isfield(fitP, 'uc')
            fitP.uc = Inf;
        end
        
        if ~isfield(fitP, 'cLevel')
            fitP.cLevel = 0.5;
        end
        
         if ~isfield(fitP, 'useML')
            fitP.useML = false;
         end       
 
         if ~isfield(fitP, 'logBin')
            fitP.logBin = true;
         end             
        
         if ~isfield(fitP, 'fitTrun')
            fitP.fitTrun = true;
         end
         
         if ~isfield(fitP, 'minBinEvents')
            fitP.minBinEvents = 3;
         end         
    end

    

    
    tau = 0.0;
    dta = 0.0;
    xmin = 0.0;
    xmax = 0.0;
    p    = 0.0;
    pcrit = 0.0;
    ks = 0.0;
    
    if fitP.logBin      
        [bins, N, edges] = LogBin(sizeAv);
    else 
        [N,edges] = histcounts(sizeAv, 'Normalization', 'probability');
        bins = (edges(1:end-1) + edges(2:end))/2;
    end    
       
    %%

    loglog(bins, N, 'bx')
    hold on;

    if fitPL
        
        if fitP.useML
            if numel(unique(sizeAv(sizeAv <= fitP.uc))) > 2
                MLcompare = mlFit(sizeAv(sizeAv <= fitP.uc), fitP.fitTrun);
                tau   = MLcompare.PL.tau;
                xmin = MLcompare.PL.xmin;
                xmax = MLcompare.PL.xmax;
                if isfield(MLcompare.PL, 'dtau')
                    dta    = MLcompare.PL.dtau; 
                end
                if isfield(MLcompare.PL, 'p')
                    p    = MLcompare.PL.p; 
                end
                if isfield(MLcompare.PL, 'pcrit')
                    pcrit   = MLcompare.PL.pcrit; 
                end
                if isfield(MLcompare.PL, 'ks')
                    ks     = MLcompare.PL.ks; 
                end                
                x = xmin:0.01:xmax;
                A = N(find(edges <= xmin, 1));
                y = A*x.^(-tau);
                loglog(x, y, 'r--');
                text(x(2), y(2)/3, strcat('S^{-', num2str(tau,3),'}'), 'Color','r')
            end
        else
            %only include bins within include range to fit
            fitEdges = edges((edges >= fitP.lc) & (edges <= fitP.uc));
            cutFront = numel(edges(edges < fitP.lc));
            cutEnd   = numel(edges(edges > fitP.uc));
            edgeCen  = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
            fitN     = N(1 + cutFront : end - cutEnd);         
            
            %fit power law
            [tau, dta] = fitPowerLawLinearLogLog(edgeCen, fitN);         
            
            x = min(edgeCen):min(edgeCen)/100:max(edgeCen);
            A = N(find(edges >= min(fitEdges), 1));
            y = A*(x/min(x)).^(-tau);
            loglog(x, y, 'r--');
            text(x(1), y(1)/3, strcat('S^{-', num2str(tau,3),'}'), 'Color','r')
            xmin = min(edgeCen);            
            xmax = max(edgeCen);  
            
        end

    end
    
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('S (events)')
    ylabel('P(S)')

    prob = N;

end