%% Set-up simulation parameters
params = struct();

% Set Simulation Options
params.SimOpt.useWorkspace    = false;
params.SimOpt.saveSim         = true;
params.SimOpt.T               = 35.0; 
params.SimOpt.dt              = 1e-2;
params.SimOpt.useParallel     = false;
params.SimOpt.hdfSave         = true;
params.importSwitch = false;
params.SimOpt.stopIfDupName = true;

%Set Stimulus
params.Stim.BiasType     = 'DC'; % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'ACsaw'
params.Stim.Amplitude    = 0.22; 


%Set Components
params.Comp.onConductance  = 0.7*7.77e-5;
params.Comp.offConductance  = 0.7*1e-9;
params.Comp.setVoltage     = 1e-2;
params.Comp.resetVoltage   = 5e-3;
params.Comp.criticalFlux   =  0.01*3.1;
params.Comp.maxFlux        = 0.015*3.1;
params.Comp.penalty        =    1;
params.Comp.boost          =    2;

params.Conn.filename = 'asn_nw_00992_nj_03008_seed_2821_avl_10.00_disp_01.00_lx_100.00_ly_100.00.mat';
params.SimOpt.nameComment = '_big';
params.SimOpt.RectElectrodes = true;


%% Run simulations
params.Comp.ComponentType  = 'atomicSwitch';
multiRun(params);
params.Comp.ComponentType  = 'tunnelSwitchL';
multiRun(params);


%% Import simulations
params.Comp.ComponentType  = 'atomicSwitch';
a = multiImport(params);

params.Comp.ComponentType  = 'tunnelSwitchL';
t = multiImport(params);


%% Load experimental dataset
load('NWNsim/experiments/experimental/experimental_network_1_timeseries_change_voltage.mat')
G = netC{1};
V = voltage(1);
tm = dt*[1:numel(G)];


%% Plot Supp Fig 20
tend = 35;
tvec = a{1}.Stim.TimeAxis;

figure('color','w', 'units', 'centimeters', 'OuterPosition', [5 5 15 15]);
h = [];
h(2) = semilogy(tvec, a{1}.netC);
hold on;
h(1) = semilogy(tm, runningMean(G, 20));
h(3) = semilogy(tvec, t{1}.netC);
xlabel('t (s)')
ylabel('G (S)')
myleg = legend(h, {'Exp', 'Bin', 'Tun'}, 'location', 'southeast');
set(findall(gca, 'Type', 'Line'),'LineWidth',2.0);
xlim([0,tend])
hold off



%% Plot Supp Fig 20 inset
figure('color','w', 'units', 'centimeters', 'OuterPosition', [5 5 8 8]);
h = [];
h(2) = semilogy(tvec, a{1}.netC);
hold on;
h(1) = semilogy(tm, runningMean(G,20));
h(3) = semilogy(tvec, t{1}.netC);
xlabel('t (s)')
ylabel('G (S)')
set(findall(gca, 'Type', 'Line'),'LineWidth',2.0);
set(gca, 'yscale', 'linear')
xlim([0,tend])
xticks([0:5:25])
axis square
hold off

