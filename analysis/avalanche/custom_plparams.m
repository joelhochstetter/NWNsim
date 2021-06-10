function [fullResults] = custom_plparams(x, xmin)
%{
    This function is equivalent to pl_params.m from Marshall 2016.
    This is modified for a given xmin loops over the possible xmax
        values and performs power-law fit

    For attribution see pl_params.m
    Modified by Joel Hochstetter
%}

%% Parse command line for parameters

nSamples = 500;
pCrit = 1.0;
likelihood = 1e-2;

%% Process data

% ensure data is oriented vertically
nX = numel(x);
x = reshape(x, nX, 1);

% find unique x values
unqX = unique(x);

%% Sort data in decreasing order by normalized Euclidean distance

% get all support pairs
% support = nchoosek(unqX,2);
support = [xmin*ones(sum(unqX > xmin + 2),1) , unqX(unqX > xmin + 2)];

% get the amount of data covered by each support pair (nData)
nSupport = size(support, 1);


%% Store fullResults
fullResults = struct('xmin', zeros(nSupport,1), ...
    'xmax', zeros(nSupport,1), 'p', zeros(nSupport,1), ...
    'tau', zeros(nSupport,1));


%% Initiate greedy search for optimal support
xmins = zeros(nSupport,1);
xmaxs = zeros(nSupport,1);
ps        = zeros(nSupport,1);
taus     = zeros(nSupport,1);


parfor iSupport = 1:nSupport
    
    
    % print the current support pair to inform user
    iSupport
    
    xmin = support(iSupport,1); 
    xmax = support(iSupport,2);
    
    % MLE of truncated distribution
    [tau, ~, ~, ~] = plmle(x, 'xmin', xmin, 'xmax', xmax);
    
    % p-value for MLE
    p = pvcalc(x, tau, 'xmin', xmin, 'xmax', xmax, 'samples', nSamples,...
        'threshold', pCrit, 'likelihood', likelihood);
    
    xmins(iSupport)  = xmin;
    xmaxs(iSupport) = xmax;
    ps(iSupport)       = p;
    taus(iSupport)    = tau;
    

end

 fullResults.xmin = xmins;
 fullResults.xmax = xmaxs;
 fullResults.p        = ps;
 fullResults.tau     = taus;



end