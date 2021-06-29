%% Set-up (Follow SFig7.m)
% Generate networks (multi_generate_networks.py)
% Run DC simulations for avalanches (DC_vary_seed_by_ensemble.m)
% Run avalanche analysis on data (simAvAnalysis.m)



%% Scaling collapse for changing L
L = [50, 100, 150, 200]; %network sizes
ms =2; %marker size for plotting

figure;
figSize = [0 0 30 20];
set(gcf, 'Units', 'centimeters', 'Position',figSize);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); 



%%
V =1; %voltage
d = 0.06; %density
nb = 15; %number of bins for plotting
bs = -1.0; %bin-size = <IEI> for time discretisation
AvFolder = 'simAvalanches';

Small =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx50/bs', num2str(bs), '/critResults.mat'));
Small = Small.critResults;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb);

Mid =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx100/bs', num2str(bs), '/critResults.mat'));
Mid = Mid.critResults;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb);

Big =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx150/bs', num2str(bs), '/critResults.mat'));
Big = Big.critResults;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb);

Fat =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx200/bs', num2str(bs), '/critResults.mat'));
Fat = Fat.critResults;
[Sbins{4}, Sprob{4}] = LogBin(Fat.avalanche.sizeAv, nb);

a = 1.6;
tau    = 1.8; 
A = 1.4;

subplot(2,3,1);
hold off;
loglog(Sbins{2}/L(2)^a, 1*A*Sprob{2}.*Sbins{2}.^tau, 'r-', 'LineWidth', ms); hold on;
loglog(Sbins{3}/L(3)^a, 0.95*A*Sprob{3}.*Sbins{3}.^tau, 'b-', 'LineWidth', ms); hold on;
loglog(Sbins{4}/L(4)^a, 1.*A*Sprob{4}.*Sbins{4}.^tau, 'k-', 'LineWidth', ms); hold on;
loglog([1e-4,1], [1,1], 'k:')
xlim([1e-4,1])
xlabel('S/L^{1.6}')
ylabel({'V^* = 1.0', '', 'P(S) S^{1.8}'})
axis square;
leg = legend('100', '150', '200', 'location', 'southwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
title('d = 0.06 nw(\mu m)^{-2}')
xrange = xlim;
yrange = ylim;
shift = -0.07;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'a','fontweight','bold','fontsize',12)




%%
V =1;
d = 0.10;
nb = 20;
bs = -1.0;


Small =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx50/bs', num2str(bs), '/critResults.mat'));
Small = Small.critResults;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb);

Mid =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx100/bs', num2str(bs), '/critResults.mat'));
Mid = Mid.critResults;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb);

Big =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx150/bs', num2str(bs), '/critResults.mat'));
Big = Big.critResults;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb);

Fat =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx200/bs', num2str(bs), '/critResults.mat'));
Fat = Fat.critResults;
[Sbins{4}, Sprob{4}] = LogBin(Fat.avalanche.sizeAv, nb);


a = 1.4;
tau    = 2.0;
A = 1;

subplot(2,3,2);
hold off;
loglog(Sbins{2}/L(2)^a, 0.9*A*Sprob{2}.*Sbins{2}.^tau, 'r-', 'LineWidth', ms); hold on;
loglog(Sbins{3}/L(3)^a, A*Sprob{3}.*Sbins{3}.^tau, 'b-', 'LineWidth', ms); hold on;
loglog(Sbins{4}/L(4)^a, A*Sprob{4}.*Sbins{4}.^tau, 'k-', 'LineWidth', ms); hold on;
loglog([1e-4,1], [1,1], 'k:')
xlim([1e-4,1])
xlabel('S/L^{1.4}')
ylabel('P(S) S^{2.0}')
axis square;
leg = legend('100', '150', '200', 'location', 'southwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
title('d = 0.10 nw(\mu m)^{-2}')

xrange = xlim;
yrange = ylim;
shift = -0.07;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'b','fontweight','bold','fontsize',12)


%%
V =1;
d = 0.14;
nb = 25;
bs = -1.0;

a = 1.2;
tau    = 2.0; 
A = 1;

Small =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx50/bs', num2str(bs), '/critResults.mat'));
Small = Small.critResults;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb);

Mid =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx100/bs', num2str(bs), '/critResults.mat'));
Mid = Mid.critResults;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb);

Big =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx150/bs', num2str(bs), '/critResults.mat'));
Big = Big.critResults;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb);

Fat =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx200/bs', num2str(bs), '/critResults.mat'));
Fat = Fat.critResults;
[Sbins{4}, Sprob{4}] = LogBin(Fat.avalanche.sizeAv, nb);

subplot(2,3,3);
hold off;
loglog(Sbins{2}/L(2)^a, 0.9*A*Sprob{2}.*Sbins{2}.^tau, 'r-', 'LineWidth', ms); hold on;
loglog(Sbins{3}/L(3)^a, A*Sprob{3}.*Sbins{3}.^tau, 'b-', 'LineWidth', ms); hold on;
loglog(Sbins{4}/L(4)^a, 1.1*A*Sprob{4}.*Sbins{4}.^tau, 'k-', 'LineWidth', ms); hold on;
loglog([1e-3,5], [1,1], 'k:')
xlim([1e-3, 5])
xticks([1e-3, 1e-2, 1e-1, 1])
xlabel('S/L^{1.2}')
ylabel('P(S) S^{2.0}')
axis square;
leg = legend('100', '150', '200', 'location', 'southwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
title('d = 0.14 nw(\mu m)^{-2}')

xrange = xlim;
yrange = ylim;
shift = -0.07;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'c','fontweight','bold','fontsize',12)


%%
V =1.8;
d = 0.06;
nb = 25;
bs = -1.0;

Small =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx50/bs', num2str(bs), '/critResults.mat'));
Small = Small.critResults;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb);

Mid =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx100/bs', num2str(bs), '/critResults.mat'));
Mid = Mid.critResults;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb);

Big =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx150/bs', num2str(bs), '/critResults.mat'));
Big = Big.critResults;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb);


a = 1.7;
tau    = 2.4; 
A = 0.5;

subplot(2,3,4);
hold off;
loglog(Sbins{1}/L(1)^a, 1.05*A*Sprob{1}.*Sbins{1}.^tau, '-', 'LineWidth', ms, 'Color', [0.9290, 0.6940, 0.1250]); hold on;
loglog(Sbins{2}/L(2)^a, A*Sprob{2}.*Sbins{2}.^tau, 'r-', 'LineWidth', ms); hold on;
loglog(Sbins{3}/L(3)^a, A*Sprob{3}.*Sbins{3}.^tau, 'b-', 'LineWidth', ms); hold on;
loglog([1e-4,1], [1,1], 'k:')
xlim([1e-4, 1])
xlabel('S/L^{1.7}')
ylabel({'V^* = 1.8', '', 'P(S) S^{2.4}'})
axis square;
leg = legend('50', '100', '150', 'location', 'northwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
title('d = 0.06 nw(\mu m)^{-2}')

xrange = xlim;
yrange = ylim;
shift = -0.07;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'd','fontweight','bold','fontsize',12)

%%
V =1.8;
d = 0.10;
nb = 25;
bs = -1.0;

Small =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx50/bs', num2str(bs), '/critResults.mat'));
Small = Small.critResults;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb);

Mid =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx100/bs', num2str(bs), '/critResults.mat'));
Mid = Mid.critResults;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb);

Big =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx150/bs', num2str(bs), '/critResults.mat'));
Big = Big.critResults;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb);

a = 2.1;
tau    = 2.1;
A = 1.1;

subplot(2,3,5);
hold off;
h(3) = loglog(Sbins{3}/L(3)^a, A*Sprob{3}.*Sbins{3}.^tau, 'b-', 'LineWidth', ms); hold on;
h(2) = loglog(Sbins{2}/L(2)^a, A*Sprob{2}.*Sbins{2}.^tau, 'r-', 'LineWidth', ms); hold on;
h(1) = loglog(Sbins{1}/L(1)^a, A*Sprob{1}.*Sbins{1}.^tau, '-', 'LineWidth', ms, 'Color', [0.9290, 0.6940, 0.1250]); hold on;

loglog([1e-4,1], [1,1], 'k:')
xlim([1e-4,1])
xlabel('S/L^{2.1}')
ylabel('P(S) S^{2.1}')
axis square;
leg = legend([h(1), h(2), h(3)],{'50', '100', '150'}, 'location', 'northwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
title('d = 0.10 nw(\mu m)^{-2}')

xrange = xlim;
yrange = ylim;
shift = -0.07;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'e','fontweight','bold','fontsize',12)


%%
V =1.8;
d = 0.14;
nb = 25;
bs = -1.0;

Small =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx50/bs', num2str(bs), '/critResults.mat'));
Small = Small.critResults;
[Sbins{1}, Sprob{1}] = LogBin(Small.avalanche.sizeAv, nb);

Mid =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx100/bs', num2str(bs), '/critResults.mat'));
Mid = Mid.critResults;
[Sbins{2}, Sprob{2}] = LogBin(Mid.avalanche.sizeAv, nb);

Big =  load(strcat(AvFolder, '/density', num2str(d, '%.2f'), '/Vstar', num2str(V), '/Lx150/bs', num2str(bs), '/critResults.mat'));
Big = Big.critResults;
[Sbins{3}, Sprob{3}] = LogBin(Big.avalanche.sizeAv, nb);


a = 2.2;
tau    = 2.1;
A = 1.12;

subplot(2,3,6);
hold off;
h(3) = loglog(Sbins{3}/L(3)^a, A*Sprob{3}.*Sbins{3}.^tau, 'b-', 'LineWidth', ms); hold on;
h(2) = loglog(Sbins{2}/L(2)^a, A*Sprob{2}.*Sbins{2}.^tau, 'r-', 'LineWidth', ms); hold on;
h(1) = loglog(Sbins{1}/L(1)^a, A*Sprob{1}.*Sbins{1}.^tau, '-', 'LineWidth', ms, 'Color', [0.9290, 0.6940, 0.1250]); hold on;

loglog([1e-4,1], [1,1], 'k:')
xlim([1e-4,1])
xlabel('S/L^{2.2}')
ylabel('P(S) S^{2.1}')
axis square;
leg = legend([h(1), h(2), h(3)],{'50', '100', '150'}, 'location', 'northwest');
title(leg, 'L (\mum)', 'fontweight', 'normal')
title('d = 0.14 nw(\mu m)^{-2}')

xrange = xlim;
yrange = ylim;
shift = -0.07;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'f','fontweight','bold','fontsize',12)