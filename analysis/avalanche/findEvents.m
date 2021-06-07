function events =  findEvents(G, eventDetect)
    defEvDet.method = 'hybrid';
    defEvDet.window = 1;
    defEvDet.thresh = 5.0e-8; %threshold of form dG >= thr or dG./G >= thr
    defEvDet.k_std  = 0.0; %threshold of form dG >= k*std(dG) or dG./G = k*std(dG./G)
    defEvDet.k_mean = 0.0; %threshold of form dG >= k*mean(dG) or dG./G = k*mean(dG./G)
    defEvDet.relThresh = 0.01;
    defEvDet.noiseFloor = 5e-10;     
    
    fields = fieldnames(defEvDet);
    for i = 1:numel(fields)
        if isfield(eventDetect, fields{i}) == 0
            eventDetect.(fields{i}) = defEvDet.(fields{i});
        end
    end
    dG = abs(gradient(G));
    dG(isnan(dG)) = 0;
    events = zeros(size(dG));
    dG(dG < eventDetect.noiseFloor) = 0.0;
    
    if numel(dG) <= 5
        return
    end
    
    switch eventDetect.method 
        case 'threshold'
            dG = abs(dG);
            events = (dG >= eventDetect.thresh) & (dG >= eventDetect.k_std*std(dG)) & (dG >= eventDetect.k_mean*mean(dG));
        case 'ratioThreshold'
            dGG = abs(dG./G);
            events = (dGG >= eventDetect.relThresh) & (dGG >= eventDetect.k_std*std(dGG)) & (dGG >= eventDetect.k_mean*mean(dGG));
        case 'hybrid'
            dG = abs(dG);
            dGG = abs(dG./G);
            events = (dGG >= eventDetect.relThresh) | (dG >= eventDetect.thresh);
        case 'stationaryPt'
            %find peaks, check bigger than threshold ratio
            %threshold ratio
            dGG = abs(dG./G);
            [pks, locs] = findpeaks(dGG);
            events(locs(pks > eventDetect.relThresh)) = true;             
            %threshold
            dG = abs(dG);
            [pks, locs] = findpeaks(dG);
            events(locs(pks > eventDetect.thresh)) = true;
        case 'thresholdPeak'
            % find crossing of threshold. Take peak from that sequence as event time.
            events = events + thresholdCrossingPeaks(dG, eventDetect.thresh);
            events = events + thresholdCrossingPeaks(-dG, eventDetect.thresh);
            events = events + thresholdCrossingPeaks(dG./G, eventDetect.relThresh);
            events = events + thresholdCrossingPeaks(-dG./G, eventDetect.relThresh);
            events = events > 0;
            
        case 'kirchoff' %use kirchoff laws to estimate an event given network conductance.
            %to do
            events = zeros(size(dG));
    end

    


end