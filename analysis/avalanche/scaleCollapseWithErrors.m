function [invsignutau, errEst] = scaleCollapseWithErrors(binnedEvents, Tmax, minTime, minFreq, nBootstraps, sampleFraction, ncpu)
%{
    Input: 
  binnedEvents: events per bin after avlanches are temporally binned:
        Tmax: maximum avalanche life-time considered
     minTime: 
     minFreq: 
       ToPlot: plots output

    Ouput:
        Determines scale-collapse parameter gamma

    Uses Marshall (2016) method to fitting scaling function

    Written by Joel Hochstetter


%}

    if nargin < 7
        ncpu = 0;
    end
    
    %%
    if nBootstraps <= 0
        [~, size_t, time_t, ~, ~, lifeAv] = avalancheShape(binnedEvents, sampleFraction, -1);
        lives = sort(unique(lifeAv))';
        lifeFreq = sum(lifeAv == lives);
        [gamma, gamma_vals, RMS_errors, ~, ~, ~] = scaleCollapse(lives, lifeFreq, time_t, size_t, Tmax, minTime, minFreq, false);
        invsignutau = 1 + gamma;
        errEst = range(gamma_vals(RMS_errors < min(RMS_errors)/0.90))/2;
        gammaVals = gamma;
    else  
        gammaVals = zeros(nBootstraps, 1);
        if ncpu > 0 
            if gcp('nocreate') == 0
                parpool(ncpu)
            end
            parfor s = 1:nBootstraps
                [~, size_t, time_t, ~, sample, lifeAv] = avalancheShape(binnedEvents, sampleFraction, s);
                lives = sort(unique(lifeAv))';
                lifeFreq = sum(lifeAv == lives);
                [gammaVals(s), ~, ~, ~, ~, ~] = scaleCollapse(lives, lifeFreq, time_t, size_t, Tmax, minTime, minFreq, false);
            end
        else
            for s = 1:nBootstraps
                [~, size_t, time_t, ~, sample, lifeAv] = avalancheShape(binnedEvents, sampleFraction, s);
                lives = sort(unique(lifeAv))';
                lifeFreq = sum(lifeAv == lives);
                [gammaVals(s), ~, ~, ~, ~, ~] = scaleCollapse(lives, lifeFreq, time_t, size_t, Tmax, minTime, minFreq, false);
            end            
        end
        invsignutau = 1 + mean(gammaVals);
        errEst = std(gammaVals);        
    end
    
    
    
    %% Plot distribution of estimates
%     figure; 
%     histogram(gammaVals + 1, 'Normalization', 'probability');
%     hold on;
%     vline(invsignutau, 'k')
%     vline(invsignutau - errEst, 'r')
%     vline(invsignutau + errEst, 'r')
%     xlabel('1/\sigma\nu\tau')
%     ylabel('P(1/\sigma\nu\tau)')
    
    
    
end