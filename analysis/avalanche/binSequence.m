function seqs = binSequence(x, times, N)
%{
    Inputs:
            x: some variable (such as event, avalanche size, branch)
    times: times when the value x is measured
           N: N is total number of time-points


    Outputs:
        seqs: cell of all elements of x at each time-point

    Written by Joel Hochstetter

%}

    assert(numel(times) == numel(x));

    seqs = cell(N,1);
    for i = 1:numel(x)
        seqs{times(i)} = [seqs{times(i)}, x(i)];
    end


end