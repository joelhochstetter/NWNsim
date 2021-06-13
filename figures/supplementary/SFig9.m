%% Run simulations / experiments and avalanche results (as in Fig. 4-5)
%{
    experimental processed files is provided for quick processing

%} 

%% Import avalanche results
Sim = load('/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar1/Lx100/bs120/critResults.mat');
results = Sim.critResults;
%included is experimental critResults.mat file
Exp = load('bs-2/Avalanche_Joint_thresholdPeak_Gt5e-08_rt0.01_nf0/critResults.mat');


%% Simulation: Calculate avalanche shape collapse exponent with errors
baseFolder = '/import/silo2/joelh/Criticality/Avalanche/FixDensity/Av/density0.10/Vstar1/Lx100';
Sim = load(strcat(baseFolder, '/bs120/critResults.mat'));
events = load(strcat(baseFolder, '/events.mat'));
events = events.events;
results = Sim.critResults;
Tmax = results.avalanche.timeFit.uc;
joinPeriod = 30000;
binSize = round(results.IEI.meanIEI);
[binned, ~, ~] = binWithJoins(events, binSize, joinPeriod);
[SimInvsignutau, SimErrEst] = scaleCollapseWithErrors(binned, Tmax, 5, 50, 1000, 1.0, 15);



%% Experiment: Calculate avalanche shape collapse exponent with errors
results = Exp.results;
events = results.events.eventTrain;
lifeAv = results.avalanche.lifeAv;
Tmax = results.avalanche.timeFit.uc;
joinPeriod = -1;
binSize = round(results.IEI.meanIEI);
[binned, ~, ~] = binWithJoins(events, binSize, joinPeriod);
[ExpInvsignutau, ExpErrEst] = scaleCollapseWithErrors(binned, Tmax, 5, 50, 1000, 1.0, 15);



%% Plot supp Fig 9a-b
ds = results.avalanche.shape.dur;
Lifes = ds';
lifeAv = results.avalanche.lifeAv;
Freq = sum(lifeAv == Lifes);
results.avalanche.timeFit.uc

figure;
scaleCollapse(Lifes, Freq, results.avalanche.shape.time_t, results.avalanche.shape.size_t, results.avalanche.timeFit.uc, 5, 50);

title(strcat('1/\sigma\nu z = ',  SimInvsignutau, ' \pm ', SimErrEst))

subplot(1,2,1);
xrange = xlim;
yrange = ylim;
shift = -0.10;
text(((1-shift)*(xrange(1)) + shift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'a','fontweight','bold','fontsize',12)
title('Simulation')


subplot(1,2,2);
xrange = xlim;
yrange = ylim;
shift = -0.10;
text(((1-shift)*(xrange(1)) + shift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'b','fontweight','bold','fontsize',12)
set(gcf, 'visible', 'on');
print(gcf,strcat('/import/silo2/joelh/Criticality/Avalanche/PaperAvalanches/SimShapeCollapse.png'), '-dpng', '-r300', '-painters')


%% Plot supp Fig 9c-d
results = Exp.results;
ds = results.avalanche.shape.dur;
Lifes = ds';
lifeAv = results.avalanche.lifeAv;
Freq = sum(lifeAv == Lifes);
results.avalanche.timeFit.uc
figure;
scaleCollapse(Lifes, Freq, results.avalanche.shape.time_t, results.avalanche.shape.size_t, 318, 5, 50);
title(strcat('1/\sigma\nu z = ',  ExpInvsignutau, ' \pm ', ExpErrEst))

subplot(1,2,1);
xrange = xlim;
yrange = ylim;
shift = -0.10;
text(((1-shift)*(xrange(1)) + shift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'c','fontweight','bold','fontsize',12)
title('Experiment')


subplot(1,2,2);
xrange = xlim;
yrange = ylim;
shift = -0.10;
text(((1-shift)*(xrange(1)) + shift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'd','fontweight','bold','fontsize',12)



