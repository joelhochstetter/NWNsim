function [alpha, dal, xmin, xmax, p, pcrit, ks, bins, prob, MLcompare] = plotAvalancheLifetime(lifeAv, fitP)
%{
    Plots the avalanche size distribution
    Inputs:
    sizeAv: Avalanche sizes
     fitP: A struct containing parameters to fit
        fitP.lc:    Lower cut-off of IEI
        fitP.uc:    Upper cut-off of IEI
        fitP.logBin: plot data with logarithmic binning 


    Option to fit if we provide cut-offs

%}
    
    MLcompare = struct();
    alpha = 0.0;
    dal = 0.0;

    if nargin == 1
        fitPL = 0;
    else
        fitPL = 1;
    end
    

    if fitPL
        %add defaults for cut-offs for PL
        if ~isfield(fitP, 'lc')
            fitP.lc = 1;
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
            fitP.logBin = false;
         end                
        
         if ~isfield(fitP, 'minBinEvents')
            fitP.minBinEvents = 3;
         end           
    end

    xmin = 1.0;
    xmax = Inf;
    p    = 0.0;
    pcrit = 0.0;
    ks = 0.0;
    
    if fitP.logBin
%         nbins = ceil(2*iqr(lifeAv)/(numel(lifeAv)^(1/3))); %calculated by Freeman Diaconis rule
        [bins, N, edges] = LogBin(lifeAv);
    else 
        [N,edges] = histcounts(lifeAv, 'Normalization', 'probability','BinMethod','integers');
        bins = (edges(1:end-1) + edges(2:end))/2;
    end
    
    %Exclude bins which have too few events
%     lifeAv = sort(lifeAv);
%     upperCut = min(edges(N < fitP.minBinEvents));
%     lifeAv(lifeAv >= upperCut) = [];    
    
    %% Extract region of distribution that is strictly decreasing
%     [~, firstMin] = findpeaks(-N);
%     firstMin = firstMin(1);
%     if firstMin < numel(N)
%         fitP.uc = edges(firstMin + 1);
%     end

    %%

    loglog(bins, N, 'bx')
    hold on;

    if fitPL
        
        if fitP.useML
            if numel(unique(lifeAv)) > 2            
%                 fitP.uc = 60;
                MLcompare = mlFit(lifeAv(lifeAv <= fitP.uc), fitP.fitTrun);
                alpha   = MLcompare.PL.tau;
                xmin = MLcompare.PL.xmin;
                xmax = MLcompare.PL.xmax;
                if isfield(MLcompare.PL, 'dtau')
                    dal    = MLcompare.PL.dtau; 
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
                y = A*x.^(-alpha);
                loglog(x, y, 'r--');
                text(x(2), y(2)/3, strcat('T^{-', num2str(alpha,3),'}'), 'Color','r')
            end
            
        else
            if numel(unique(lifeAv)) > 2
                %only include bins within include range to fit
                fitEdges = edges((edges >= fitP.lc) & (edges <= fitP.uc));
                cutFront = numel(edges(edges < fitP.lc));
                cutEnd   = numel(edges(edges > fitP.uc));
                edgeCen  = (fitEdges(1:end-1)  + fitEdges(2:end))/2;
                fitN     = N(1 + cutFront : end - cutEnd);         

                %fit power law
                [alpha, dal] = fitPowerLawLinearLogLog(edgeCen, fitN);         

                x = min(edgeCen):min(edgeCen)/100:max(edgeCen);
                A = N(find(edges >= min(fitEdges), 1));
                y = A*(x/min(x)).^(-alpha);
                loglog(x, y, 'r--');
                text(x(1), y(1)/3, strcat('T^{-', num2str(alpha,3),'}'), 'Color','r')
                xmin = min(edgeCen);            
                xmax = max(edgeCen);            
            end
        end
    end

    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('T (bins)')
    ylabel('P(T)')
    
    prob = N;
    
end