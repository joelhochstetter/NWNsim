%% Generate networks
% Run python commands from bash to generate networks
% Note we require extra number of simulations as some have no s-d paths
%
% python multi_generate_networks.py --Lx 50 --LxMax 200 --numSims 4 --seedMax 3100 --density 0.06 --folder nets/Density0.06ChangeSize
% python multi_generate_networks.py --Lx 50 --LxMax 200 --numSims 4 --seedMax 1050 --density 0.10 --folder nets/Density0.10ChangeSize
% python multi_generate_networks.py --Lx 50 --LxMax 200 --numSims 4 --seedMax 1010 --density 0.14 --folder nets/Density0.14ChangeSize


%% Run DC simulations for avalanches
% Vstar = 1.8 simulations are run for SFig 8 as well

baseFolder = '.';%'sims';
netFolders = {'Density0.06ChangeSize', 'Density0.10ChangeSize', 'Density0.14ChangeSize'};
NSims = 1000;
Lvals = [50, 100, 150, 200]; %set network size (side-length of square)
Vstar     = [1, 1.8];
T = 30;


%%
%loop over seed 
for seed = 1:1000
    for d = 1:3
        for L = Lvals
            for v = Vstar
                saveFolder = strcat(baseFolder, '/', netFolders{d});
                netFolder = strcat('nets/', netFolders{d});
                DC_vary_seed_by_ensemble(seed, netFolder, L, saveFolder, v)
            end
        end
    end
end


%% Process avalanches from simulations
binSize = -1; %use Average inter event interval
density = [0.06, 0.10, 0.14];
fitML = false;
avFolder = 'simAvalanches';
simAvAnalysis(baseFolder, avFolder, Vstar, Lvals, density, binSize, NSims, T, fitML);


%% network size dependence: import sims
V =1.0; %voltage
d = 0.10; %density
mxS = 1000; %maximum avalanche size for plotting
mxT = 200; %maximum avalanche life-time for plotting
nb = 25;   %number of bins

Small =  load(strcat(avFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx50/bs-1/critResults.mat'));
Small = Small.critResults;
Sprob{1} = Small.avalanche.sizeFit.prob;
Sbins{1} = Small.avalanche.sizeFit.bins;
Tprob{1} = Small.avalanche.timeFit.prob;
Tbins{1} = Small.avalanche.timeFit.bins;
mSze{1} = Small.avalanche.avSizeFit.mSize;
mLfe{1} = Small.avalanche.avSizeFit.mLife;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb, mxS);
[Tbins{1}, Tprob{1}] = LogBin(Small.avalanche.lifeAv, nb, mxT);


Mid =  load(strcat(avFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx100/bs-1/critResults.mat'));
Mid = Mid.critResults;
Sprob{2} = Mid.avalanche.sizeFit.prob;
Sbins{2} = Mid.avalanche.sizeFit.bins;
Tprob{2} = Mid.avalanche.timeFit.prob;
Tbins{2} = Mid.avalanche.timeFit.bins;
mSze{2} = Mid.avalanche.avSizeFit.mSize;
mLfe{2} = Mid.avalanche.avSizeFit.mLife;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb, mxS);
[Tbins{2}, Tprob{2}] = LogBin(Mid.avalanche.lifeAv, nb, mxT);


Big =  load(strcat(avFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx150/bs-1/critResults.mat'));
Big = Big.critResults;
Sprob{3} = Big.avalanche.sizeFit.prob;
Sbins{3} = Big.avalanche.sizeFit.bins;
Tprob{3} = Big.avalanche.timeFit.prob;
Tbins{3} = Big.avalanche.timeFit.bins;
mSze{3} = Big.avalanche.avSizeFit.mSize;
mLfe{3} = Big.avalanche.avSizeFit.mLife;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb, mxS);
[Tbins{3}, Tprob{3}] = LogBin(Big.avalanche.lifeAv, nb, mxT);


Fat =  load(strcat(avFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx200/bs-1/critResults.mat'));
Fat = Fat.critResults;
Sprob{4} = Fat.avalanche.sizeFit.prob;
Sbins{4} = Fat.avalanche.sizeFit.bins;
Tprob{4} = Fat.avalanche.timeFit.prob;
Tbins{4} = Fat.avalanche.timeFit.bins;
mSze{4} = Fat.avalanche.avSizeFit.mSize;
mLfe{4} = Fat.avalanche.avSizeFit.mLife;
[Sbins{4}, Sprob{4}] = LogBin(Fat.avalanche.sizeAv, nb, mxS);
[Tbins{4}, Tprob{4}] = LogBin(Fat.avalanche.lifeAv, nb, mxT);


%% Plot SF7(a-c)
ms = 8;
shift = -0.12;

figure;
set(gcf, 'color', 'w');
figSize = [0 0 26 15];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); 
set(gcf, 'Units', 'centimeters', 'Position',figSize)
subplot(2,3,1)
loglog(Sbins{1}, Sprob{1}, '.', 'Markersize', ms, 'Color', [0.9290, 0.6940, 0.1250])
hold on;
loglog(Sbins{2}, Sprob{2}, 'r.', 'Markersize', ms)
loglog(Sbins{3}, Sprob{3}, 'b.', 'Markersize', ms)
loglog(Sbins{4}, Sprob{4}, 'k.', 'Markersize', ms)
xlabel('S')
ylabel('P(S)')

xmin = Fat.avalanche.sizeFit.lc;
xmax = Fat.avalanche.sizeFit.uc;
tau    =Fat.avalanche.sizeFit.tau;
x = xmin:0.01:xmax;
A = 1.1*Sprob{2}(round(xmin));
y = A*(x/xmin).^(-tau);
loglog(x, y, 'r-');
text(10, 0.2, strcat('S^{-2.0}'), 'Color','r')

xlim([1,1000])
xticks([1,10,100,1000])
ylim([1e-6, 1])
yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
xrange = xlim;
yrange = ylim;
axis square;
leg = legend('50', '100', '150', '200', 'location', 'southwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')

text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'a','fontweight','bold','fontsize',12)

subplot(2,3,2);
loglog(Tbins{1}, Tprob{1}, '.', 'Markersize', ms, 'Color', [0.9290, 0.6940, 0.1250])
hold on;
loglog(Tbins{2}, Tprob{2}, 'r.', 'Markersize', ms)
loglog(Tbins{3}, Tprob{3}, 'b.', 'Markersize', ms)
loglog(Tbins{4}, Tprob{4}, 'k.', 'Markersize', ms)

xmin = Fat.avalanche.sizeFit.lc;
xmax = Fat.avalanche.sizeFit.uc;
tau    = Fat.avalanche.sizeFit.tau;
x = xmin:0.01:xmax;
A = 1.4*Sprob{2}(4);
y = A*(x/xmin).^(-tau);
loglog(x, y, 'r-');
text(5, 2e-1, 'T^{-2.3}', 'Color','r')

xlabel('T')
ylabel('P(T)')
xlim([1,100])
xticks([1,10,100,1000])
ylim([1e-5, 1])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
axis square;
leg = legend('50', '100', '150', '200', 'location', 'southwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'b','fontweight','bold','fontsize',12)

subplot(2,3,3);
loglog(mLfe{1}, 0.9*mSze{1}, '.', 'Markersize', ms, 'Color', [0.9290, 0.6940, 0.1250])
hold on;
loglog(mLfe{2}, mSze{2}, 'r.', 'Markersize', ms)
loglog(mLfe{3}, mSze{3}, 'k.', 'Markersize', ms)
loglog(mLfe{4}, mSze{4}, 'k.', 'Markersize', ms)
axis square;
xlim([1,100])
xticks([1,10,100])

xmin = Fat.avalanche.timeFit.lc;
xmax = Fat.avalanche.timeFit.uc;
gamma    = Fat.avalanche.avSizeFit.gamma_m_1;
x = xmin:0.01:xmax;
A = 0.70*mSze{4}(find(mLfe{4} <= xmin + 1, 1))*3;
y = A*(x/xmin).^(gamma);
loglog(x, y, 'r-');
text(5, 5, strcat('T^{', num2str(gamma,2),'}'), 'Color','r')
xlabel('T')
ylabel('<S>(T)')
leg = legend('50', '100', '150', '200', 'location', 'northwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'c','fontweight','bold','fontsize',12)

%% density dependence: import sims
V = 1;  %voltage
L = 200; %network size
bs = -1; %bin size
mxS = 1000;  %maximum avalanche size for plotting
mxT = 200;%maximum avalanche life-time for plotting
nb = 25; %number of bins

Small =  load(strcat('Av/density0.06/Vstar', num2str(V), '/Lx',  num2str(L), '/bs', num2str(bs), '/critResults.mat'));
Small = Small.critResults;
Sprob{1} = Small.avalanche.sizeFit.prob;
Sbins{1} = Small.avalanche.sizeFit.bins;
Tprob{1} = Small.avalanche.timeFit.prob;
Tbins{1} = Small.avalanche.timeFit.bins;
mSze{1} = Small.avalanche.avSizeFit.mSize;
mLfe{1} = Small.avalanche.avSizeFit.mLife;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb, mxS);
[Tbins{1}, Tprob{1}] = LogBin(Small.avalanche.lifeAv, nb, mxS);

Mid =  load(strcat('Av/density0.10/Vstar', num2str(V), '/Lx',  num2str(L), '/bs', num2str(bs), '/critResults.mat'));
Mid = Mid.critResults;
Sprob{2} = Mid.avalanche.sizeFit.prob;
Sbins{2} = Mid.avalanche.sizeFit.bins;
Tprob{2} = Mid.avalanche.timeFit.prob;
Tbins{2} = Mid.avalanche.timeFit.bins;
mSze{2} = Mid.avalanche.avSizeFit.mSize;
mLfe{2} = Mid.avalanche.avSizeFit.mLife;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb, mxS);
[Tbins{2}, Tprob{2}] = LogBin(Mid.avalanche.lifeAv, nb, mxS);

Big =  load(strcat('Av/density0.14/Vstar', num2str(V), '/Lx',  num2str(L), '/bs', num2str(bs), '/critResults.mat'));
Big = Big.critResults;
Sprob{3} = Big.avalanche.sizeFit.prob;
Sbins{3} = Big.avalanche.sizeFit.bins;
Tprob{3} = Big.avalanche.timeFit.prob;
Tbins{3} = Big.avalanche.timeFit.bins;
mSze{3} = Big.avalanche.avSizeFit.mSize;
mLfe{3} = Big.avalanche.avSizeFit.mLife;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb, mxS);
[Tbins{3}, Tprob{3}] = LogBin(Big.avalanche.lifeAv, nb, mxS);


%%
subplot(2,3,4)
hold off;
loglog(Sbins{1}, Sprob{1}, '.', 'Markersize', ms)
hold on;
loglog(Sbins{2}, Sprob{2}, 'r.', 'Markersize', ms)
loglog(Sbins{3}, Sprob{3}, 'k.', 'Markersize', ms)
xlabel('S')
ylabel('P(S)')

xmin = Mid.avalanche.sizeFit.lc;
xmax =Mid.avalanche.sizeFit.uc;
tau    =Mid.avalanche.sizeFit.tau;
x = xmin:0.01:xmax;
A = 1.1*Sprob{2}(3);
y = A*(x/xmin).^(-tau);
loglog(x, y, 'r-');
text(10, 0.2, strcat('S^{-2.0}'), 'Color','r')

xlim([1,1000])
xticks([1,10,100,1000])
ylim([1e-5, 1])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
xrange = xlim;
yrange = ylim;
axis square;
leg = legend('0.06', '0.10', '0.14', 'location', 'northeast');
title(leg, 'nw/\mum^2', 'fontweight', 'normal')

text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'd','fontweight','bold','fontsize',12)

subplot(2,3,5);
hold off;
loglog(Tbins{1}, Tprob{1}, '.', 'Markersize', ms)
hold on;
loglog(Tbins{2}, Tprob{2}, 'r.', 'Markersize', ms)
loglog(Tbins{3}, Tprob{3}, 'k.', 'Markersize', ms)

xmin = 3;
xmax = 25;
tau    = 2.3;
x = xmin:0.01:xmax;
A = 0.9*Sprob{2}(3);
y = A*(x/xmin).^(-tau);
loglog(x, y, 'r-');
text(5, 2e-1, 'T^{-2.3}', 'Color','r')

xlabel('T')
ylabel('P(T)')
xlim([1,100])
xticks([1,10,100])
ylim([1e-5, 1])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
axis square;
leg = legend('0.06', '0.10', '0.14', 'location', 'southwest');
title(leg, 'nw/\mum^2', 'fontweight', 'normal')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'e','fontweight','bold','fontsize',12)

subplot(2,3,6);
hold off;
loglog(mLfe{1}, 0.9*mSze{1}, '.', 'Markersize', ms)
hold on;
loglog(mLfe{2}, mSze{2}, 'r.', 'Markersize', ms)
loglog(mLfe{3}, mSze{3}, 'k.', 'Markersize', ms)
axis square;
xlim([1,100])
xticks([1,10,100])

xmin = Mid.avalanche.timeFit.lc;
xmax = Mid.avalanche.timeFit.uc;
gamma    = Mid.avalanche.avSizeFit.gamma_m_1;
x = xmin:0.01:xmax;
A = 0.70*mSze{2}(find(mLfe{2} <= xmin + 1, 1))*3;
y = A*(x/xmin).^(gamma);
loglog(x, y, 'r-');
text(5, 50, strcat('T^{', num2str(gamma,2),'}'), 'Color','r')
xlabel('T')
ylabel('<S>(T)')
leg = legend('0.06', '0.10', '0.14', 'location', 'southeast');
title(leg, 'nw/\mum^2', 'fontweight', 'normal')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'f','fontweight','bold','fontsize',12)
