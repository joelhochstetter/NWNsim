function [re_tm,re_sz, scale, acoeff, gamma] = avalancheShapeCollapse(dur, size_t, time_t)
%{
  Computes shape collapse function

    Input:
        dur (Nx1 array) - each unique avalanche duration
        size_t (Nx1 cell) - average avalanche size as a function of time
                              for each unique duration. Each element is a
                              time-series vector
        time_t (Nx1 cell) - stores the time vectors for the duration of
                              each event in the analysis 


    Output:
        re_sz (Nx1 cell) - average avalanche size as a function of time
                              for each unique duration. Each element is a
                              time-series vector. Rescaled for time.
        
        re_tm (Nx1 cell) - stores the time vectors for the duration of
                              each event in the analysis. Rescaled for
                              time.

        scale (Nx1 cell) - stores the scale function
        
%}

    N = numel(dur);
    
    re_tm = cell(N,1);
    re_sz = cell(N,1);
    scale = cell(N,1);
    
    %% rescale durations  to follow convention of Mehta 2002
    for i = 1:N
        dur(i) = dur(i);
    end
    
    
    %% calculate rescaled time
    for i = 1:N
        re_tm{i} = time_t{i}/(dur(i) + 1);
    end
    
    %% rescaled size function
    for i = 1:N
        re_sz{i} = size_t{i};
    end
    
   
    %% calculate scale function
    %Fits shape function to find exponent
    
    [acoeff, gamma] = fitScaleFunction(dur, re_tm, re_sz);    
    
    for i = 1:N
        scale{i} = size_t{i}*dur(i)^(-gamma);
    end

end