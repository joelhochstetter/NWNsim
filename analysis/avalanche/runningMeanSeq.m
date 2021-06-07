function mx =  runningMeanSeq(x, times, N, dt)
%{
    Calculates the running mean of a given variable:
            x: some variable (such as event, avalanche size, branch)
    times: times when the value x is measured
           N: N is total number of time-points
          dt: Discretisation used for running bin
%}


    seq                      = binSequence(x, times, N);
    [means, counts] = seqMeanCounts(seq);
     mx = runningWeightedMean(means, counts, dt);

end