%% Run simulations as in Fig 4.



%% Loop over bin-sizes and produce criticality analysis
binSize = [80:40:320]; %use Average inter event interval
simAvAnalysis(saveFolder, 'simAvalanches', 1.0, L, density, binSize, NSims, 30);



%% Load criticality analysis results
for j = 1:Nbs
    bs = binSize(j);

    saveFolder = strcat(baseFolder, '/SimAvalanches/bs', num2str(bs),'/');
    mkdir(saveFolder);


    for i = 1:numel(vals)
        critResults{i,j} = load(strcat2({baseFolder, '/', subtype, num2str(vals(i), fmt), '/bs', bs, '/critResults.mat'}));
        critResults{i,j} = critResults{i,j}.critResults;
    end
end



%% Plot SFig 16
compAvalanche( '/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar0.7/', 100, 'Lx', 'Lx', 80:40:320)
compAvalanche( '/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar1/', 100, 'Lx', 'Lx', 80:40:320)
compAvalanche( '/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar1.8/', 100, 'Lx', 'Lx', 80:40:320)



%% Time-step dependence of avalanches
sdx1 = 4;
sdx2 = 2;
shift = -0.12;
copyfile '/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar0.7/AvCompare/AvComp.mat' '/import/silo2/joelh/Criticality/Avalanche/PaperAvalanches/AvCompVstar0.7.mat'
copyfile '/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar1/AvCompare/AvComp.mat' '/import/silo2/joelh/Criticality/Avalanche/PaperAvalanches/AvCompVstar1.mat'
copyfile '/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar1.8/AvCompare/AvComp.mat' '/import/silo2/joelh/Criticality/Avalanche/PaperAvalanches/AvCompVstar1.8.mat'

avFolder = '/import/silo2/joelh/Criticality/Avalanche/PaperAvalanches/';
figure;
figSize = [0 0 24 21];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm
set(gcf, 'Units', 'centimeters', 'Position',figSize)

load(strcat(avFolder, 'AvCompVstar0.7.mat'))
p = parula(numel(Szbins));
subplot(3,3,1)
for i = 1:numel(Szbins)
    loglog(Szbins{i}, Szprob{i}, 'color', p(i,:))
    hold on;
end
xlabel('S')
ylabel('P(S)')
title('V^* = 0.7')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'a','fontweight','bold','fontsize',12)

subplot(3,3,2)
for i = 1:numel(Tmbins)
    loglog(Tmbins{i}, Tmprob{i}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('P(T)')

xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'b','fontweight','bold','fontsize',12)

subplot(3,3,3)
for i = 1:numel(ASlife)
    loglog(ASlife{i}, ASsize{i}, 'color', p(i,:))
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

load(strcat(avFolder, 'AvCompVstar1.mat'))
p = parula(numel(Szbins));
subplot(3,3,4)
for i = 1:numel(Szbins)
    loglog(Szbins{i}, Szprob{i}, 'color', p(i,:))
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
for i = 1:numel(Tmbins)
    loglog(Tmbins{i}, Tmprob{i}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('P(T)')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'e','fontweight','bold','fontsize',12)

subplot(3,3,6)
for i = 1:numel(ASlife)
    loglog(ASlife{i}, ASsize{i}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('<S>(T)')
box on;
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'f','fontweight','bold','fontsize',12)

load(strcat(avFolder, 'AvCompVstar1.8.mat'))
p = parula(numel(Szbins));
subplot(3,3,7)
for i = 1:numel(Szbins)
    loglog(Szbins{i}, Szprob{i}, 'color', p(i,:))
    hold on;
end
xlabel('S')
ylabel('P(S)')
title('V^* = 1.8')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'g','fontweight','bold','fontsize',12)

subplot(3,3,8)
for i = 1:numel(Tmbins)
    loglog(Tmbins{i}, Tmprob{i}, 'color', p(i,:))
    hold on;
end
xlim([1,100])
xlabel('T')
ylabel('P(T)')
xrange = xlim;
yrange = ylim;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'g','fontweight','bold','fontsize',12)

subplot(3,3,9)
for i = 1:numel(ASlife)
    loglog(ASlife{i}, ASsize{i}, 'color', p(i,:))
    hold on;
end
xlabel('T')
ylabel('<S>(T)')
xlim([1,100])
xrange = xlim;
yrange = ylim;
box on;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'i','fontweight','bold','fontsize',12)