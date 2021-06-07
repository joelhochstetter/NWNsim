function [weights, mse, rnmse, y] = NLT(target, nwV, theseNWs)
%{
    Performs non-linear transformation task:
        inputs: 
            nwV (TxN array): nanowire voltages
            theseNWs: array a subset of all nanowires (If theseNWs = 0 then
            all are used)
            target: Tx1 array
            
    Written by Joel Hochstetter
%}
    
    if theseNWs == 0
       theseNWs = 1:size(nwV, 2); 
    end
    Xmat = nwV(:, theseNWs);
    weights = regress(target, Xmat);
    y = Xmat*weights;
    
    [mse, rnmse] = calcMSE(y, target);

end
