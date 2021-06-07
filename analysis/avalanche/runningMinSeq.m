function mx =  runningMinSeq(x, times, N, dt)
%{
    Calculates the running min of a given variable:
            x: some variable (such as event, avalanche size, branch)
    times: times when the value x is measured
           N: N is total number of time-points
          dt: Discretisation used for running bin
%}


    seq      = binSequence(x, times, N);
    mins    = seqMinCounts(seq);
    mx       = runningMin(mins, dt);

end