%{
    This produces avalanche statistics for experiments
        to plot Figures 4d-f

    You can skip to line 51 to use pre-processed analysis.

%}

%% Load time-series data
load('experiments/experimental/experimental_network_2_timeseries_fixed_voltage.mat');


%% Pre-process data for analysis
%select datasets to use for analysis
G = netC([1,2,4]);

% get time-vectors, and cut when current is too low
times = cell(size(G)); %time-vectors
V  = voltage*ones(numel(G)); %specify voltage for each

Icut = 1e-8;

for i = 1:numel(G)
    times{i} = dt*[1:numel(G{i})];
    I = G{i}*V(i);
    G{i}(I < Icut) = Icut*V(i);
end


%% Process avalanches from data
Gt = 5e-8; %threshold on delta G to be an event
r    = 0.03; %threshold on deltaG/G to be an event

fitML = false; %fit using maximum likelihood method
saveFolder = 'expAvalanches';
binSize = -1; %use average inter-event-interval as bin-size

%we use threshold Peak event detection method. For other methods see findEvents.m
eventDetect = struct('method', 'thresholdPeak', 'thresh', Gt, 'relThresh', r); 

% conditions allow us to cut from first to last event. See applyConditions.m
conditions =  struct('type','eventInterval', 'thresh', Gt, 'ratio', r); 

combinedCritAnalysis(G, V, times, strcat2({saveFolder, '/bs', binSize, '/'}), eventDetect, fitML, binSize, conditions)



%% SKIP HERE FOR QUICK ANALYSIS 
%% Import files of processed avalanches
% Exp = load('expAvalanches/bs-1/critResults.mat'); %for run as above
Exp = load('experiments/experimental/mainCritResults.mat'); %pre-processed

Exp = Exp.results;


%% Plot Figure 4
nb = 20; %number of bins
figure;
subplot(2,3,4);
sizeAv = Exp.avalanche.sizeAv;
xmin = Exp.avalanche.sizeFit.lc;
xmax = Exp.avalanche.sizeFit.uc;
tau    = Exp.avalanche.sizeFit.tau;
dtau  = Exp.avalanche.sizeFit.dTau;
% [N,edges] = histcounts(sizeAv, 'Normalization', 'probability'); %linear bins
% bins = (edges(1:end-1) + edges(2:end))/2;
[bins, N, edges] = LogBin(sizeAv, nb);

x = xmin:0.01:xmax;
A = N(find(edges <= xmin, 1));
y = A*(x/xmin).^(-tau);
loglog(bins, N, 'r.')
hold on;
loglog(x, y, 'k-');
xlabel('S')
ylabel('P(S)')
text(100, 1e-1, strcat('S^{-', num2str(tau,2),'}'), 'Color','k')
title('Experiment')
xlim([1,1000])
xticks([1,10,100,1000])
ylim([1e-4, 1])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'd','fontweight','bold','fontsize',12)

subplot(2,3,5);
lifeAv = Exp.avalanche.lifeAv;
xmin = Exp.avalanche.timeFit.lc;
xmax = Exp.avalanche.timeFit.uc;
alpha    = Exp.avalanche.timeFit.alpha;
% [N,edges] = histcounts(lifeAv, 'Normalization', 'probability'); %linear bins
% bins = (edges(1:end-1) + edges(2:end))/2;

[bins, N, edges] = LogBin(lifeAv, nb);

x = xmin:0.01:xmax;
A = N(find(edges <= xmin, 1))/2;
y = A*(x/xmin).^(-alpha);
loglog(bins, N, 'r.')
hold on;
loglog(x, y, 'k-');
xlabel('T')
ylabel('P(T)')
text(20, 1e-1, strcat('T^{-', num2str(alpha,2),'}'), 'Color','k')
xlim([1,100])
xticks([1,10,100])
ylim([1e-4, 1])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'e','fontweight','bold','fontsize',12)

subplot(2,3,6);
[mSize, mLife] = avalancheAvSize(sizeAv, lifeAv);
gamma_m_1 = Exp.avalanche.gamma.x2;
loglog(mLife, mSize, 'r.')
A = mSize(find(mLife > xmin, 1));
y = A*(x/xmin).^(gamma_m_1);
hold on;
loglog(x, y, 'k-');
xlabel('T')
ylabel('\langle S \rangle (T)')
text(2, 100, strcat('T^{', num2str(gamma_m_1,2),'}'), 'Color','k')
xlim([1,100])
xticks([1,10,100])
ylim([1,1000])
yticks([1,10,100,1000])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'f','fontweight','bold','fontsize',12)


%% Experimental exponents
disp(strcat2({'\tau = ', Exp.avalanche.sizeFit.tau, '+/-', Exp.avalanche.sizeFit.dTau}));
disp(strcat2({'\alpha = ', Exp.avalanche.timeFit.alpha, '+/-', Exp.avalanche.timeFit.dAlpha}));
dx1 = ((Exp.avalanche.timeFit.dAlpha/(Exp.avalanche.timeFit.alpha - 1)) + (Exp.avalanche.sizeFit.dTau/(Exp.avalanche.sizeFit.tau - 1))) *Exp.avalanche.gamma.x1;
disp(strcat2({'S,T: 1/\sigma \nu z = ', Exp.avalanche.gamma.x1, '+/-', dx1}));
disp(strcat2({'<S>(T): 1/\sigma \nu z = ', Exp.avalanche.gamma.x2, '+/-', Exp.avalanche.gamma.dx2}));
disp(strcat2({'Suc = ', Exp.avalanche.sizeFit.uc, ', Tuc = ', Exp.avalanche.timeFit.uc}))
