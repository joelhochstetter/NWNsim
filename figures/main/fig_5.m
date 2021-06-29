%{
    Pre-processed file is provided in NWNsim/experiments/DC/ 

    To use set avFolder = 'NWNsim/experiments/DC/' and run code
        block on line 47
%}


%% Generate relevant networks
% L = 100x100, density = 0.1 nw/(um)^2
%From bash: 
% python multi_generate_networks.py --Lx 100 --seedMax 1000 --density 0.10 --folder Density0.10ChangeSize

netFolder = 'Density0.10ChangeSize';

%ensure all networks  have source-drain paths
GetConnectedNWNs(netFolder)


%% Run avalanches (as in Fig. 4)
saveFolder = 'sims/Density0.10ChangeSize';
NSims = 1000;
L = 100; %set network size (side-length of square)
density = 0.1;
Vstar   = [0.7, 1.0, 1.05, 1.8];
T = 30;

%loop over seed 
parfor vv = 1:numel(Vstar)
    for seed = 1:NSims
        v = Vstar(vv);
        DC_vary_seed_by_ensemble(seed, netFolder, L, saveFolder, v, T)
    end
end


%% Process avalanches
binSize = 160; %use Average inter event interval at Vstar = 1
baseFolder = '.'; %'sims';
avFolder = 'simAvalanches';
fitML = false;

simAvAnalysis(baseFolder, avFolder, Vstar, L, density, binSize, NSims, T, fitML);



%% Import avalanches
d = 0.10;
L = 100;
bs = 160; %this is the <IEI> at Vstar = 1 for this network. Change 

% ensure bins are the same for each distribution
dS = 3;
edgesS =[0:dS:5000] + 0.5*mod(dS,2);
dT = 1;
edgesT =[0:dT:5000] + 0.5*mod(dT,2);

V0pt7 =  load(strcat(avFolder, '/density', num2str(d, '%.2f'), '/Vstar0.7/Lx', num2str(L), '/bs', num2str(bs), '/critResults.mat'));
V0pt7 = V0pt7.critResults;
Sprob{1} = V0pt7.avalanche.sizeFit.prob;
Sbins{1} = V0pt7.avalanche.sizeFit.bins;
Tprob{1} = V0pt7.avalanche.timeFit.prob;
Tbins{1} = V0pt7.avalanche.timeFit.bins;
[Sprob{1}, e1] = histcounts(V0pt7.avalanche.sizeAv, edgesS, 'Normalization', 'probability');
Sbins{1} = (e1(1:end-1) + e1(2:end))/2;
[Tprob{1}, e1] = histcounts(V0pt7.avalanche.lifeAv, edgesT, 'Normalization', 'probability');
Tbins{1} = (e1(1:end-1) + e1(2:end))/2;

V1 =  load(strcat(avFolder, '/density',num2str(d, '%.2f'), '/Vstar1/Lx', num2str(L), '/bs', num2str(bs), '/critResults.mat'));
V1 = V1.critResults;
Sprob{2} = V1.avalanche.sizeFit.prob;
Sbins{2} = V1.avalanche.sizeFit.bins;
Tprob{2} = V1.avalanche.timeFit.prob;
Tbins{2} = V1.avalanche.timeFit.bins;
[Sprob{2}, e1] = histcounts(V1.avalanche.sizeAv, edgesS, 'Normalization', 'probability');
Sbins{2} = (e1(1:end-1) + e1(2:end))/2;
[Tprob{2}, e1] = histcounts(V1.avalanche.lifeAv, edgesT, 'Normalization', 'probability');
Tbins{2} = (e1(1:end-1) + e1(2:end))/2;

V1pt05 =  load(strcat(avFolder, '/density',num2str(d, '%.2f'), '/Vstar1.05/Lx', num2str(L), '/bs', num2str(bs), '/critResults.mat'));
V1pt05 = V1pt05.critResults;
Sprob{3} = V1pt05.avalanche.sizeFit.prob;
Sbins{3} = V1pt05.avalanche.sizeFit.bins;
Tprob{3} = V1pt05.avalanche.timeFit.prob;
Tbins{3} = V1pt05.avalanche.timeFit.bins;
[Sprob{3}, e1] = histcounts(V1pt05.avalanche.sizeAv, edgesS, 'Normalization', 'probability');
Sbins{3} = (e1(1:end-1) + e1(2:end))/2;
[Tprob{3}, e1] = histcounts(V1pt05.avalanche.lifeAv, edgesT, 'Normalization', 'probability');
Tbins{3} = (e1(1:end-1) + e1(2:end))/2;

V1pt8 =  load(strcat(avFolder, '/density',num2str(d, '%.2f'), '/Vstar1.8/Lx', num2str(L), '/bs', num2str(bs), '/critResults.mat'));
V1pt8 = V1pt8.critResults;
Sprob{4} = V1pt8.avalanche.sizeFit.prob;
Sbins{4} = V1pt8.avalanche.sizeFit.bins;
Tprob{4} = V1pt8.avalanche.timeFit.prob;
Tbins{4} = V1pt8.avalanche.timeFit.bins;
[Sprob{4}, e1] = histcounts(V1pt8.avalanche.sizeAv, edgesS, 'Normalization', 'probability');
Sbins{4} = (e1(1:end-1) + e1(2:end))/2;
[Tprob{4}, e1] = histcounts(V1pt8.avalanche.lifeAv, edgesT, 'Normalization', 'probability');
Tbins{4} = (e1(1:end-1) + e1(2:end))/2;



%% Plot figure 5
ms = 8;
figure;
set(gcf, 'color', 'w');
figSize = [0 0 17 8.5];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); 
set(gcf, 'Units', 'centimeters', 'Position',figSize)
subplot(1,2,1)
loglog(Sbins{1}, Sprob{1}, '.', 'Markersize', ms, 'color', 'k')
hold on;
loglog(Sbins{2}, Sprob{2}, '.', 'Markersize', ms, 'color', 'r')
loglog(Sbins{3}, Sprob{3}, '.', 'Markersize', ms, 'color', 'c')
loglog(Sbins{4}, Sprob{4}, '.', 'Markersize', ms, 'color', 'b')
xlabel('S')
ylabel('P(S)')

xmin = V1.avalanche.sizeFit.lc;
xmax = V1.avalanche.sizeFit.uc;
tau    = V1.avalanche.sizeFit.tau;
x = xmin:0.01:xmax;
A = Sprob{2}(find(Sbins{2} <= xmin + 1, 1))*2;
y = A*(x/xmin).^(-tau);
loglog(x, y, 'r-');
text(20, 0.1, strcat('S^{-', num2str(tau,2),'}'), 'Color','r')

xlim([1,2000])
xticks([1,10,100,1000])
ylim([1e-5, 1])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
xrange = xlim;
yrange = ylim;
shift = -0.07;
axis square;
leg =  legend('0.7', '1.0',  '1.05', '1.8', 'location', 'northeast');
title(leg, 'V^*', 'fontweight', 'normal')
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'a','fontweight','bold','fontsize',12)


subplot(1,2,2);
loglog(Tbins{1}, Tprob{1}, '.', 'Markersize', ms, 'color', 'k')
hold on;
loglog(Tbins{2}, Tprob{2}, '.', 'Markersize', ms, 'color', 'r')
loglog(Tbins{3}, Tprob{4}, '.', 'Markersize', ms, 'color', 'c')
loglog(Tbins{4}, Tprob{4}, '.', 'Markersize', ms, 'color', 'b')

xmin = V1.avalanche.timeFit.lc;
xmax = V1.avalanche.timeFit.uc;
alpha    = V1.avalanche.timeFit.alpha;
x = xmin:0.01:xmax;
A = Tprob{2}(find(Tbins{2} <= xmin, 1))/2;
y = A*(x/xmin).^(-alpha);
loglog(x, y, 'r-');
text(5, 2e-1, strcat('T^{-', num2str(alpha,2),'}'), 'Color','r')

xlabel('T')
ylabel('P(T)')
xlim([1,100])
ylim([1e-5, 1])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
axis square;
leg =  legend('0.7', '1.0',  '1.05', '1.8', 'location', 'southwest');
title(leg, 'V^*', 'fontweight', 'normal')
xrange = xlim;
yrange = ylim;
shift = -0.07;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'b','fontweight','bold','fontsize',12)
