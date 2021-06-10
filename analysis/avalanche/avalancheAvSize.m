function [mSize, mLife] = avalancheAvSize(sizeAv, lifeAv)
%{
    Input:
        sizeAv (Ax1 array) - number of events in given avalanche
        lifeAv (Ax1 array) - number of bins avalanche goes for 

    Avalanche defined such that events happen in subsequent time bins and
    no events happen in preceding and after time-bin

    A avalanches are recorded

    N unique lifetimes


    Output:
        mLife (Nx1 array) - life-times
        mSize (Nx1 array) - mean number of events for a given life-time

    Written by Joel Hochstetter       
%}
    
    mLife = unique(lifeAv, 'sorted');
    mSize = zeros(size(mLife));
    
    for i = 1:numel(mLife)  
        mSize(i) = mean(sizeAv(find(lifeAv == mLife(i))));
    end



end