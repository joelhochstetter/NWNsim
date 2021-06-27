%% Run simulations (as in Fig 5)


%% Specify parameters
%using bin-size = 160ms as this corresponds to average inter-event interval
    % at V^* = 1
binSize = 160; 

    
    
%% Sub-critical simulation: Load corresponding conductance time-series and events file 
baseFolder = 'Av/density0.10/Vstar0.7/Lx100/';
G = load(strcat(baseFolder, '/netC.mat'));
G1 = G.netC;
events = load(strcat(baseFolder, '/events.mat'));
events = events.events;
joinperiod = 30000; %each simulation was 30s = 30000 time-steps
[binned1, binTimes1, binJoinPeriod1] = binWithJoins(events, binSize, joinperiod);



%% Critical simulation: Load corresponding conductance time-series and events file
baseFolder = 'Av/density0.10/Vstar1/Lx100/';
G = load(strcat(baseFolder, '/netC.mat'));
G2 = G.netC;
events = load(strcat(baseFolder, '/events.mat'));
events = events.events;
joinperiod = 30000; %each simulation was 30s = 30000 time-steps
[binned2, binTimes2, binJoinPeriod2] = binWithJoins(events, binSize, joinperiod);




%% Super-critical simulation: Load corresponding conductance time-series and events file
baseFolder = 'Av/density0.10/Vstar1.8/Lx100/';
G = load(strcat(baseFolder, '/netC.mat'));
G3 = G.netC;
events = load(strcat(baseFolder, '/events.mat'));
events = events.events;
joinperiod = 30000;
[binned3, binTimes3, binJoinPeriod3] = binWithJoins(events, binSize, joinperiod);



%% Experiments: Load corresponding conductance time-series and events 
%Load critResults.mat file for experiments as this saves G and events
Exp = load('ExpFolder/critResults.mat');
Exp = Exp.results;
Gx   = Exp.net.G;

critResultsx = Exp;
events = Exp.events.eventTrain;
[binnedx, binTimesx, binJoinPeriodx] = binWithJoins(events, Exp.avalanche.binSize, -1);


%%
j = 150;
h = 0.9;


figure;

%sub-critical
subplot(2,2,1);
binTimes      = binTimes1;
binned             = binned1;
G                               = G1;
time = [(1:1:numel(G))*1e-3]';
yl = ylim;
zx = [4500, 4515, 4515, 4500];
zy = [yl(1), yl(1), yl(2), yl(2)];
hold on;
plot(binTimes*1e-3, binned)
hold on;
ylabel('events/bin')
yyaxis right;
plot(time, G)
ylabel('G (S)')
xlabel('t (s)')
xlim([j*30 + 0,(j + 1)*30])
box on;
title('V^* < 1')
xrange = xlim;
yrange = ylim;
shift = -0.07;
text(((1-shift)*(xrange(1)) + shift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'a','fontweight','bold','fontsize',12)


%critical
subplot(2,2,2);
binTimes      = binTimes2;
binned             = binned2;
G                               = G2;
time = [(1:1:numel(G))*1e-3]';
yl = ylim;
zx = [4500, 4515, 4515, 4500];
zy = [yl(1), yl(1), yl(2), yl(2)];
hold on;
plot(binTimes*1e-3, binned)
hold on;
ylabel('events/bin')
yyaxis right;
plot(time, G)
ylabel('G (S)')
xlabel('t (s)')
xlim([j*30 + 0,(j + 1)*30])
ylim([1.2e-8, 3.21e-8])
box on;
title('V^* = 1')
xrange = xlim;
yrange = ylim;
xshift = -0.07;
yshift = -0.07;
text(((1-xshift)*(xrange(1)) + xshift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'b','fontweight','bold','fontsize',12)


%super-critical
subplot(2,2,3);
binTimes      = binTimes3;
binned             = binned3;
G                               = G3;
time = [(1:1:numel(G))*1e-3]';

yl = ylim;
zx = [4500, 4515, 4515, 4500];
zy = [yl(1), yl(1), yl(2), yl(2)];
hold on;
plot(binTimes*1e-3, binned)
hold on;
ylabel('events/bin')
yyaxis right;
plot(time, G)
ylabel('G (S)')
xlabel('t (s)')
xlim([j*30 + 0,(j + 1)*30])
box on;
title('V^* > 1')
xrange = xlim;
yrange = ylim;
shift = -0.07;
text(((1-shift)*(xrange(1)) + shift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'c','fontweight','bold','fontsize',12)


%Experiment
subplot(2,2,4);
binTimes      = binTimesx;
binned             = binnedx;
G                               = Gx;
time = [(1:1:numel(G))*1e-3]';

yl = ylim;
zx = [4500, 4515, 4515, 4500];
zy = [yl(1), yl(1), yl(2), yl(2)];
hold on;
plot(binTimes*1e-3, binned)
hold on;
ylabel('events/bin')
yyaxis right;
plot(time, G)
ylabel('G (S)')
xlabel('t (s)')
xlim([j*30 + 0,(j + 1)*30])
box on;
title('Experiment')
xrange = xlim;
yrange = ylim;
shift = -0.07;
text(((1-shift)*(xrange(1)) + shift*(xrange(2))) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'd','fontweight','bold','fontsize',12)

