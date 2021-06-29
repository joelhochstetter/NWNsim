%% Run simulations as in Fig 5



%% Loop over bin-sizes and produce criticality analysis
baseFolder = '.'; %'sims';
density = 0.10;
L = 100; %network size
binSize = [80:40:320]'; %in units of ms
NSims = 1000;
fitML = false;
avFolder = 'simAvalanches';
Vstar =[0.7, 1.0, 1.8];

simAvAnalysis(baseFolder, avFolder, Vstar, L, density, binSize, NSims, 30, fitML);



%% Load criticality analysis results
Vstars = [0.7, 1.0, 1.8];

Szbins   = cell(numel(binSize), numel(Vstars));
Szprob  = cell(numel(binSize), numel(Vstars));
Tmbins  = cell(numel(binSize), numel(Vstars));
Tmprob  = cell(numel(binSize), numel(Vstars));
ASlife      = cell(numel(binSize), numel(Vstars));
ASsize    = cell(numel(binSize), numel(Vstars));
Nbs = numel(binSize);

for v = 1:numel(Vstars)
    for j = 1:Nbs
        Vstar = Vstars(v);
        bs = binSize(j);

        saveName = strcat2({avFolder, '/density', num2str(density, '%.2f'), '/Vstar', Vstar, '/Lx', L, '/bs',  num2str(bs), '/critResults.mat'});
        critResults = load(saveName);
        critResults = critResults.critResults;
        Szbins{j,v} = critResults.avalanche.sizeFit.bins;
        Szprob{j,v} = critResults.avalanche.sizeFit.prob;
        Tmbins{j,v} = critResults.avalanche.timeFit.bins;
        Tmprob{j,v} = critResults.avalanche.timeFit.prob;
        ASlife{j,v} = critResults.avalanche.avSizeFit.mLife;
        ASsize{j,v} = critResults.avalanche.avSizeFit.mSize;        
    end
end



%% Time-step dependence of avalanches
sdx1 = 4;
sdx2 = 2;
shift = -0.12;

figure;
figSize = [0 0 24 21];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm
set(gcf, 'Units', 'centimeters', 'Position',figSize)

p = parula(numel(binSize));
subplot(3,3,1)
for i = 1:numel(binSize)
    loglog(Szbins{i,1}, Szprob{i,1}, 'color', p(i,:))
    hold on;
end
xlabel('S')
ylabel('P(S)')
title('V^* = 0.7')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'a','fontweight','bold','fontsize',12)

subplot(3,3,2)
for i = 1:numel(binSize)
    loglog(Tmbins{i,1}, Tmprob{i,1}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('P(T)')

xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'b','fontweight','bold','fontsize',12)

subplot(3,3,3)
for i = 1:numel(binSize)
    loglog(ASlife{i,1}, ASsize{i,1}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('<S>(T)')
leg = legend(num2str(round(binSize)));
title(leg, '\Delta t (ms)')
xrange = xlim;
yrange = ylim;
box on;

text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'c','fontweight','bold','fontsize',12)

subplot(3,3,4)
for i = 1:numel(binSize)
    loglog(Szbins{i,2}, Szprob{i,2}, 'color', p(i,:))
    hold on;
end
xlabel('S')
ylabel('P(S)')
title('V^* = 1.0')
xlim([1,1000])
xticks([1,10,100,1000])
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'd','fontweight','bold','fontsize',12)

subplot(3,3,5)
for i = 1:numel(binSize)
    loglog(Tmbins{i,2}, Tmprob{i,2}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('P(T)')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'e','fontweight','bold','fontsize',12)

subplot(3,3,6)
for i = 1:numel(binSize)
    loglog(ASlife{i,2}, ASsize{i,2}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('<S>(T)')
box on;
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'f','fontweight','bold','fontsize',12)

subplot(3,3,7)
for i = 1:numel(binSize)
    loglog(Szbins{i,3}, Szprob{i,3}, 'color', p(i,:))
    hold on;
end
xlabel('S')
ylabel('P(S)')
title('V^* = 1.8')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'g','fontweight','bold','fontsize',12)

subplot(3,3,8)
for i = 1:numel(binSize)
    loglog(Tmbins{i,3}, Tmprob{i,3}, 'color', p(i,:))
    hold on;
end
xlim([1,100])
xlabel('T')
ylabel('P(T)')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'g','fontweight','bold','fontsize',12)

subplot(3,3,9)
for i = 1:numel(binSize)
    loglog(ASlife{i,3}, ASsize{i,3}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('<S>(T)')
xlim([1,100])
xrange = xlim;
yrange = ylim;
box on;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'i','fontweight','bold','fontsize',12)