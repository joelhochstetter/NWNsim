function [mse, rnmse] = calcMSE(y, T)
%{
    Calculates mean square error (mse) and root mean square error (rnmse)
    given vector y and target vector T

%}

    mse   = mean((y-T).^2);
    rnmse = sqrt(mse/mean(T.^2));
end
