%% Run simulations / experiments and avalanche results (as in Fig. 4-5)
%{
    experimental processed files is provided for quick processing
        in 'NWNsim/experiments/combinedCritResults.mat'

    simulation processed files for quick processing in 
        NWNsim/DC/L100
%} 

%% Import avalanche results
Sim = load('/import/silo2/joelh/Criticality/Avalanche/FixDensity/AvNew1/density0.10/Vstar1/Lx100/bs120/critResults.mat');
results = Sim.critResults;
%included is experimental critResults.mat file
Exp = load('experiments/experimental/combinedCritResults.mat');


%% Simulation: Calculate avalanche shape collapse exponent with errors
baseFolder = '/import/silo2/joelh/Criticality/Avalanche/FixDensity/AvNew/density0.10/Vstar1/Lx100';
Sim = load(strcat(baseFolder, '/bs120/critResults.mat'));
results = Sim.critResults;
events = results.events.eventTrain;
Tmax = results.avalanche.timeFit.uc;
joinPeriod = 30000;
binSize = round(results.IEI.meanIEI);
[binned, ~, ~] = binWithJoins(events, binSize, joinPeriod);
[SimInvsignutau, SimErrEst] = scaleCollapseWithErrors(binned, Tmax, 5, 50, 1000, 1.0, 15);



%% Plot supp Fig 9a-b
ds = results.avalanche.shape.dur;
Lifes = ds';
lifeAv = results.avalanche.lifeAv;
Freq = sum(lifeAv == Lifes);

figure;
scaleCollapse(Lifes, Freq, results.avalanche.shape.time_t, results.avalanche.shape.size_t, results.avalanche.timeFit.uc, 5, 50);

title(strcat('1/\sigma\nu z = ',  num2str(SimInvsignutau, '%.2f'), ' \pm ', num2str(SimErrEst, '%.2f')))

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


%% Experiment: Calculate avalanche shape collapse exponent with errors
results = Exp.results;
events = results.events.eventTrain;
lifeAv = results.avalanche.lifeAv;
Tmax = results.avalanche.timeFit.uc;
joinPeriod = -1;
binSize = round(results.IEI.meanIEI);
[binned, ~, ~] = binWithJoins(events, binSize, joinPeriod);
[ExpInvsignutau, ExpErrEst] = scaleCollapseWithErrors(binned, Tmax, 5, 50, 1000, 1.0, 15);


%% Plot supp Fig 9c-d
results = Exp.results;
ds = results.avalanche.shape.dur;
Lifes = ds';
lifeAv = results.avalanche.lifeAv;
Freq = sum(lifeAv == Lifes);

figure;
scaleCollapse(Lifes, Freq, results.avalanche.shape.time_t, results.avalanche.shape.size_t, 318, 5, 50);
title(strcat('1/\sigma\nu z = ',  num2str(ExpInvsignutau, '%.2f'), ' \pm ', num2str(ExpErrEst, '%.2f')));

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



