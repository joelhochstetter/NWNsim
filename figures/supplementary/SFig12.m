%% Run NLT simulations (as in Fig 8)


%% Change to directory where NLT simulations were saved
NLTfolder  = '.';
cd(NLTfolder);

%%
N = 225;
load( 'NLTresults.mat');
sinRNMSE  = NLTres.sinRNMSE;
dblRNMSE  = NLTres.dblRNMSE;
phsRNMSE = NLTres.phsRNMSE;
squRNMSE = NLTres.squRNMSE;
sinResult  = NLTres.sinResult;
dblResult  = NLTres.dblResult;
phsResult  = NLTres.phsResult;
squResult  = NLTres.squResult;    


%%
%select frequencies and amplitudes
thisAmp = 1;
thisFreq = [0.25, 0.85, 1];
ToPlot = sort(intersect(find(NLTres.Amps == thisAmp), find(sum(NLTres.Freqs == thisFreq, 2) > 0)))'; 


%% Plotting target signal
dt = 1e-3;
T  = 80;
numT = 2;
j = 8;

lw = 2.5;

timeVector = dt:dt:T;
figure;
figSize = [0 0 20 15];
set(gcf, 'Units', 'centimeters', 'Position',figSize);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm

targetStim = struct('Amplitude', 1, 'Frequency', 1, 'Phase', 0.0, 'BiasType', 'AC');
targetSOpt = struct('T', numT, 'dt', 1e-3);
sine = getStimulus(targetStim, targetSOpt);

targetStim = struct('AmplitudeOn', 1, 'AmplitudeOff', -1, 'OffTime', 1/2, 'Phase', 0.0, 'BiasType', 'Square', 'Duty', 50.0);
targetSOpt = struct('T', numT, 'dt', 1e-3);
square = getStimulus(targetStim, targetSOpt);

targetStim = struct('Amplitude', 1, 'Frequency', 2, 'Phase', 0.0, 'BiasType', 'ACsaw');
targetSOpt = struct('T', numT, 'dt', 1e-3);
dbl = getStimulus(targetStim, targetSOpt);

targetStim = struct('Amplitude', 1, 'Frequency', 1, 'Phase', pi/2, 'BiasType', 'ACsaw');
targetSOpt = struct('T', numT, 'dt', 1e-3);
phs = getStimulus(targetStim, targetSOpt);

subplot(2,2,1);
plot(dt:dt:numT, sine.Signal, 'k--', 'LineWidth', lw)
hold on;
plot(timeVector'.*(thisFreq.*ones(80000,1)), sinResult(ToPlot, :)')
xlim([0,numT])
xrange = xlim;
yrange = ylim;
shift = -0.09;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'a','fontweight','bold','fontsize',12)
xlabel('t/T')
ylabel('Signal')
title('Sinusoidal')
subplot(2,2,2);
plot(dt:dt:numT, square.Signal, 'k--', 'LineWidth', lw)
hold on;
plot(timeVector'.*(thisFreq.*ones(80000,1)), squResult(ToPlot, :)')
xlim([0,numT])
xrange = xlim;
yrange = ylim;
shift = -0.09;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'b','fontweight','bold','fontsize',12)
xlabel('t/T')
ylabel('Signal')
title('Square')
subplot(2,2,3);
plot(dt:dt:numT, phs.Signal, 'k--', 'LineWidth', lw)
hold on;
plot(timeVector'.*(thisFreq.*ones(80000,1)), phsResult(ToPlot, :)')
xlim([0,numT])
xrange = xlim;
yrange = ylim;
shift = -0.09;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'c','fontweight','bold','fontsize',12)
xlabel('t/T')
ylabel('Signal')
title('\pi/2 phase shifted')
subplot(2,2,4);
plot(dt:dt:numT, dbl.Signal, 'k--', 'LineWidth', lw)
hold on;
plot(timeVector'.*(thisFreq.*ones(80000,1)), dblResult(ToPlot, :)')
xlim([0,numT])
xrange = xlim;
yrange = ylim;
shift = -0.09;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'd','fontweight','bold','fontsize',12)
xlabel('t/T')
ylabel('Signal')
title('Double frequency')