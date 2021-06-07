function [gamma, gamma_vals, RMS_errors, new_t, size_t_ave, size_t_var] = scaleCollapse(lives, lifeFreq, time_t, size_t, Tmax, minTime, minFreq, ToPlot)
%{
    Input: 
        lifeAv:   lifetime of avalanches
       freq: frequency of avalanche
        re_tm: rescaled time
       size_t: size as a function of time
        Tmax: maximum avalanche life-time considered
       ToPlot: plots output

    Ouput:
        Determines scale-collapse parameter gamma

    Uses Marshall method

%}

    %% set default values
    if nargin < 5
        Tmax = inf;
    end

    if nargin < 6
        minTime = 10;
    end
    
    if nargin < 7
      minFreq  = 50;
    end
    
    if nargin < 8
        ToPlot = true;
    end
    
    %% Snip zeros off data
    if time_t{1}(1) == 0
        snip = 1; %if contains zeroes we snip
    else
        snip = 0;
    end
    
    for i = 1:numel(time_t)
        time_t{i} = time_t{i}(1+snip:end-snip);
        size_t{i}  = size_t{i}(1+snip:end-snip);        
    end
    
    
    %% Reshape to row vectors
    lives = reshape(lives, [numel(lives), 1]);
    lifeFreq = reshape(lifeFreq, [numel(lifeFreq), 1]);
    
    %% Calculate re-scaled time
    assert(numel(lives) == numel(time_t));
    re_tm = cell(size(time_t));
    for i = 1:numel(time_t)
        re_tm{i} = time_t{i}/(lives(i));
    end
    
    
    
    %% Prune data according to appropriate datasets
   
    selected = (lives >= minTime) & (lifeFreq >= minFreq) & (lives <= Tmax);
    lives     = lives(selected);
    re_tm = re_tm(selected);
    size_t    = size_t(selected);
    time_t = time_t(selected);
    
    
    
    %% Interpolate data
    min_t = 1/minTime;
    dt = 0.001;
    new_t = min_t:dt:1;
    size_t_int = cell(size(size_t));
    
    for i = 1:numel(size_t)
        size_t_int{i} = interp1(re_tm{i}, size_t{i}, new_t);
    end


    %% Loop through gamma values and times
    gamma_vals = 0:0.01:2.0;
    size_t_var =  zeros(numel(gamma_vals, numel(new_t)));
    size_t_ave =  zeros(numel(gamma_vals, numel(new_t)));
    RMS_errors = zeros(size(gamma_vals));
    
    
    for g = 1:numel(gamma_vals)
        for t = 1:numel(new_t)
            m_sz       = 0; %average size
            m_sq_sz = 0; %average size squared  
            maxVal = 0;
            minVal  = inf;
            for d = 1:numel(lives)
                m_sz       = m_sz +size_t_int{d}(t)/lives(d)^(gamma_vals(g));
                m_sq_sz = m_sq_sz + (size_t_int{d}(t)/lives(d)^(gamma_vals(g)))^2;
%                 if size_t_int{d}(t) < minVal
%                     minVal = size_t_int{d}(t);%/lifeAv(d)^(gamma_vals(g));
%                 end
%                 if size_t_int{d}(t) > maxVal
%                     maxVal = size_t_int{d}(t);%/lifeAv(d)^(gamma_vals(g));
%                 end                
            end
            m_sz         = m_sz/numel(lives);
            m_sq_sz   = m_sq_sz/numel(lives);
            size_t_var(g, t) = m_sq_sz - m_sz^2;
            size_t_ave(g,t) = m_sz;
            RMS_errors(g) = RMS_errors(g) + size_t_var(g, t);
        end
        maxVal = max(size_t_ave(g,:));
        minVal  = min(size_t_ave(g,:));
        RMS_errors(g) = RMS_errors(g)/(maxVal - minVal)^2;%*lifeAv(d)^(2*gamma_vals(g));
    end
    
    %% Plot fitting results
%     if ToPlot
%         figure;
%         plot(gamma_vals, RMS_errors);
%         xlabel('\gamma')
%         ylabel('RMS error')
%     end
    
    %% Get gamma
    [~, gm] = min(RMS_errors); %index for best gamma 

    gamma = gamma_vals(gm);
    
    
    %% Plot scaling function
    if ToPlot
        subplot(1,2,1);
    %     boundedline(new_t, size_t_ave(gm,:), sqrt(size_t_var(gm,:)), 'transparency', 0.1)%, 'k--',)    
        for i = 1:numel(lives)
            hold on;        
            plot(time_t{i}, size_t{i})
        end
    %     plot(new_t, size_t_ave(gm,:), 'k', 'Linewidth', 2.5)
        xlabel('t')
        ylabel('s(t,T)')
    %     title(strcat2({'1/\sigma\nu z = ', gamma + 1}))
        axis square;
        box on;

        subplot(1,2,2);
        boundedline(new_t, size_t_ave(gm,:), sqrt(size_t_var(gm,:)), 'transparency', 0.1)%, 'k--',)    
        for i = 1:numel(lives)
            hold on;        
            plot(re_tm{i}, size_t{i}/lives(i)^gamma)
        end
        plot(new_t, size_t_ave(gm,:), 'k', 'Linewidth', 2.5)
        xlabel('t/T')
        ylabel('s(t,T) T^{1 - 1/\sigma\nu z}')
        title(strcat2({'1/\sigma\nu z = ', num2str(gamma + 1,2)}))
        axis square;
        box on;    
    end
    
end