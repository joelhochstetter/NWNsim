%% Set-up parameters to run simulation
%To see full list of possible parameters go to runSim.m
params = struct();

% Set Simulation Options
params.SimOpt.saveSim         = true;

params.SimOpt.T                =  7;
params.SimOpt.dt               = 1e-3;

%Set Stimulus
params.Stim.BiasType     = 'DC'; % for stimulus options see getStimulus.m
params.Stim.Amplitude    =  0.1; %Vstar = 1.1

%set connectivity file
params.Conn.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';

%Set Components
params.Comp.ComponentType  = 'tunnelSwitchL';
params.Comp.onConductance  = 7.77e-5;
params.Comp.offConductance  = 7.77e-8;
params.Comp.setVoltage     = 1e-2;
params.Comp.resetVoltage  = 5e-3;
params.Comp.criticalFlux   =  0.01;
params.Comp.maxFlux        = 0.015;
params.Comp.penalty        =    1;
params.Comp.boost          =   1;


%% Run simulation
multiRun(params);


%% Import saved simulation file
t = multiImport(params);


%% Get connectivity 
sim = t{1};
contacts = sim.ContactNodes;
Connectivity          = struct('filename', sim.ConnectFile);
Connectivity          = getConnectivity(Connectivity);

adjMat       = Connectivity.weights;
sp  = kShortestPath(adjMat, contacts(1), contacts(2), 10); % Get shortest path nodes 
spE = getPathEdges(sp{1}, Connectivity.EdgeList); % Get shortest path eges 

dt        = sim.dt;
timeVec   = dt:dt:sim.T;



%% Plot Fig 3a (Junction Conductance)
figure;
tend = 7.0;
cmap = parula(10);

G0 = sim.Comp.onG;

%Ordered by activitation time
ordAct = [1,2,9,3,4,5,8,7,6];
spA = spE(ordAct);

hold on;
h = [];

for i = 9:-1:1
    h(i) = semilogy(timeVec,sim.swC(:,spA(i))/G0, '-', 'Color', cmap(10-ordAct(i), :));
end

xlabel('t (s)', 'FontWeight', 'bold')
ylabel('G_{jn} (G_0)', 'FontWeight', 'bold')
set(gca, 'YScale', 'log')
ylim([0.8e-3,1.5])
set(gca,'XLim',[0,tend])
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
yyaxis right
h(10) = semilogy(timeVec, 9*sim.netC/G0,'--','Linewidth',2., 'Color', 'r');
ylabel('G_{nw}^*', 'Color','r', 'FontWeight', 'bold')
set(gca,'YColor','red');
leg = [string(1:9), 'net'];
ordLeg = [1,2,4,5,6,9,8,7,3, 10];
myleg = legend(h(ordLeg), leg, 'location','northwest');
title(myleg, 'Junction #', 'FontWeight', 'bold')
set(gca, 'YScale', 'log')
ylim([0.8e-3,1.5])
axis square;
box on;


%% Plot Fig 3b (Voltage distribution)
timeVector = timeVec;

figure('color','w', 'units', 'centimeters', 'OuterPosition', [5 5 18 18]);
contourf(timeVector(1:1:end),1:9, abs(sim.swV(1:1:end,spE))'/1e-2)
ylabel('Distance from source','fontweight','bold','fontsize',10);
xlabel('t (s)','fontsize',10);
colormap(parula)
hcb = colorbar;
caxis([0,2.5])
title(hcb,'V_{jn}^*', 'FontSize', 10)
hcb.FontSize = 10;
hcb.FontWeight = 'bold';
xlim([0,5.5])
set(gca,'YDir','reverse')
grid on;
axis square


%% Plot Figs 3c-e Creates snapshot from data - enter snapshots
i = 1;

edgs = [];
nds  = [];
idx = [];
idx(1) = find(timeVec >= 2, 1);
idx(2) = find(timeVec >= 4, 1);
idx(3) = find(timeVec >= 6, 1);
pad = 300;

for i = idx
    edgs  = [];
    whatToPlot            = struct('Dissipation',  false, 'VDrop',  false, 'GraphRep', true, 'Voltages', true, 'Nanowires', true, 'Lambda', false, 'Currents',false, 'Lyapunov', false, 'LogConductance', true);
    axesLimits            = struct('dVCbar',[0; max(max(abs(sim.swV(1,:))))], 'CurrentArrowScaling',8e-3, 'ConCbar', [min(min(sim.swC)), max(max(sim.swC))]);
    snapshot = generateSnapshotFromData(sim.swV(i,:)', sim.swLam(i,:)', sim.swC(i,:)',  sim.Comp.critFlux, sim.Stim.Signal(i), sim.netC(i), i*sim.dt);
    [~, p] = snapshotToFigurePaper(snapshot, sim.ContactNodes, Connectivity, whatToPlot, axesLimits, nds, edgs);
    if i == idx(3) %label shortest path junction of Fig 3e
        labeledge(p, spE, 1:9)
        p.EdgeFontSize = 20;
        p.EdgeFontWeight = 'Bold';
        p.EdgeLabelColor = 'w';
    end
    set(gcf, 'visible','on')
    xlim([-pad,3000+pad])
    ylim([-pad,3000+pad])

    parula1 = parula;
    parula1(end,:) = [1,1,1];
    colormap(parula);
    caxis([-3, -1.7])   
    colorbar('hide');
    axis square;
    text(0, 3150, strcat('t = ', num2str(timeVec(i), 2), 's'), 'Color', 'w', 'FontSize', 23)
    myfig = gcf;
    myfig.InvertHardcopy = 'off'; 
    strcat('snapshot_t',  num2str(timeVec(i), 1), 's.png')
end


%% Get colorbar for Figs 3c-e
figure;
ax = axes;
cbar = colorbar(ax, 'Location', 'west');
ax.Visible = 'off';
colormap(parula);
caxis([-3, -1.7])  
title(cbar, 'G_{jn} (G_0)', 'FontSize', 10, 'FontWeight', 'bold');
cbar.FontSize = 10;
cbar.Ticks = [-3, -2.7, -2.3, -2, -1.7];
cbar.TickLabels = {'1\times10^{-3}', '2\times 10^{-3}', '5\times 10^{-3}','1\times10^{-2}', '2\times10^{-2}'};
