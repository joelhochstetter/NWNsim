%% Generate networks (as in Fig. 4)
netSize = 150;
netFolder = 'nets/Density0.10ChangeSize';


%% Get network names
% for a given electrode pairing some networks at low sizes will have no source-drain paths 
% in this case (eg at L = 50, or density = 0.06) run 'GetConnectedNWNs'
% and this will give a mapping between seed and seed of connected nwn

if exist(strcat2({netFolder, '/conn_lx_', netSize, '.mat'}), 'file')
    load(strcat2({netFolder, '/conn_lx_', netSize, '.mat'}), 'conSeeds')
    actualSeed = conSeeds(seedIdx);
    nets = dir(strcat(netFolder, '/*_seed_', num2str(actualSeed,'%03.f'), '*lx_', num2str(netSize), '*.mat'))';                
else
    nets = dir(strcat(netFolder, '/*_seed_', num2str(seedIdx - 1,'%03.f'), '*lx_', num2str(netSize), '*.mat'))';
end

connFile = nets(1).name; %this file must exist. but if everything before is done connectly it will
addpath(netFolder);


%% Run AC attractor
dt = 5e-4; %set time-step for simulations
T  = 3000; %Set simulation time: Must run an integer number of periods for all frequencies
Amps = [5]; %Specify stimulus amplitudes
Freqs = [0.25, 0.5, 1]; %Specify stimulus frequencies
attractorFolder = 'S15attractors'; %set save folder for attractor
runACAttractor(Amps, Freqs, dt, T, attractorFolder, connFile)


%% Calculate Av exponents




%% Run Avalanche simulations





%% Run Lyapunov simulations




%% import AC avalanche
mxS = 50000;
nb = 60;

load('AvChangeFreq/f0.25/bs-1/critResults.mat');
[bins{1}, prob{1}] = LogBin(critResults.avalanche.sizeAv, nb, mxS);

load('AvChangeFreq/f0.5/bs-1/critResults.mat');
[bins{2}, prob{2}] = LogBin(critResults.avalanche.sizeAv, nb, mxS);


load('AvChangeFreq/f1/bs-1/critResults.mat');
[bins{3}, prob{3}] = LogBin(critResults.avalanche.sizeAv, nb, mxS);


Sim = load('simAvalanches/density0.10/Vstar1/Lx100/bs-1/critResults.mat');
Sim = Sim.critResults;
[bins{4}, prob{4}] = LogBin(Sim.avalanche.sizeAv, nb, mxS);


%%
figure;
figSize = [0 0 15 15];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm
set(gcf, 'Units', 'centimeters', 'Position',figSize);

loglog(bins{1}, prob{1}, 'g.', 'Markersize', 20)
hold on;
loglog(bins{2}, prob{2}, 'm.', 'Markersize', 20)
loglog(bins{3}, prob{3}, 'k.', 'Markersize', 20)
loglog(bins{4}, prob{4}, 'r.', 'Markersize', 20)

x = 4:0.01:82;
y = 1*(x/1).^(-2);
loglog(x, y, 'r', 'Linewidth', 2);

xlabel('S', 'Fontsize', 15)
ylabel('P(S)', 'Fontsize', 15)
xlim([1,5000])
ylim([1e-6, 1])
axis square;
lg = legend('\lambda < 0', '\lambda \approx 0', '\lambda > 0', 'DC V^* = 1');
lg.FontSize = 12;
