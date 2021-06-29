%{
    This produces avalanche statistics for simulations
        to plot Figures 4a-c

    Pre-processed file is provided in NWNsim/experiments/DC (see line 43)
%}

%% Generate networks
% Follow instructions in generation/README
% L = 150x150, density = 0.1 nw/(um)^2
%From bash: 
% python multi_generate_networks.py --Lx 150 --seedMax 1000 --density 0.10 --folder nets/Density0.10ChangeSize
netFolder = 'Density0.10ChangeSize'; %'nets/Density0.10ChangeSize';

%ensure all networks  have source-drain paths
GetConnectedNWNs(netFolder)


%% Run simulations for simulated avalanches
saveFolder = 'sims/Density0.10ChangeSize';
NSims = 1000;
L = 150; %set network size (side-length of square)
density = 0.1;
Vstar   = 1;
T = 30;
EventThreshold = 1e-3; %threshold on deltaG/G

%loop over seed 
parfor seed = 1:NSims
     DC_vary_seed_by_ensemble(seed, netFolder, L, saveFolder, Vstar, T)
end


%% Process avalanches from simulations
baseFolder = 'sims';
binSize = -1; %use Average inter event interval
fitML = true; %use maximum likelihood fitting procedure

simAvAnalysis(baseFolder, 'simAvalanches', 1.0, L, density, binSize, NSims, T, fitML);


%% Import files of processed avalanches
Sim = load('experiments/DC/simAvalanches/density0.10/Vstar1/Lx150/bs-1/critResults.mat'); 
Sim = Sim.critResults;


%% Plot Figure 4
figure;

figSize = [0 0 20 10];
set(gcf, 'Units', 'centimeters', 'Position',figSize);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm

subplot(2,3,1);
sizeAv = Sim.avalanche.sizeAv;
xmin = Sim.avalanche.sizeFit.lc;
xmax = Sim.avalanche.sizeFit.uc;
tau    = Sim.avalanche.sizeFit.tau;

[N,edges] = histcounts(sizeAv, 'Normalization', 'probability');
x = xmin:0.01:xmax;
A = N(find(edges <= xmin, 1));
y = A*(x/xmin).^(-tau);

loglog((edges(1:end-1) + edges(2:end))/2, N, 'r.')
hold on;
loglog(x, y, 'k-');
xlabel('S')
ylabel('P(S)')
text(100, 1e-1, strcat('S^{-', num2str(tau,2),'}'), 'Color','k')
title('Simulation')
xlim([1,1000])
xticks([1,10,100,1000])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'a','fontweight','bold','fontsize',12)


subplot(2,3,2);
lifeAv = Sim.avalanche.lifeAv;
xmin = Sim.avalanche.timeFit.lc - 1;
xmax = Sim.avalanche.timeFit.uc + 1;
alpha    = Sim.avalanche.timeFit.alpha;

[N,edges] = histcounts(lifeAv, 'Normalization', 'probability');
x = xmin:0.01:xmax;
A = N(find(edges <= xmin, 1));
y = A*(x/xmin).^(-alpha);
loglog((edges(1:end-1) + edges(2:end))/2, N, 'r.')
hold on;
loglog(x, y, 'k-');
xlabel('T')
ylabel('P(T)')
text(10, 1e-1, strcat('T^{-', num2str(alpha,2),'}'), 'Color','k')
xlim([1,100])
xticks([1,10,100])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'b','fontweight','bold','fontsize',12)


subplot(2,3,3);
[mSize, mLife] = avalancheAvSize(sizeAv, lifeAv);
gamma_m_1 = Sim.avalanche.gamma.x2;
loglog(mLife, mSize, 'r.')
A = mSize(find(mLife <= xmin, 1))*3;
y = A*(x/xmin).^(gamma_m_1);
hold on;
loglog(x, y, 'k-');
xlabel('T')
ylabel('\langle S \rangle (T)')
text(2, 100, strcat('T^{', num2str(gamma_m_1,2),'}'), 'Color','k')
xlim([1,100])
xticks([1,10,100])
ylim([1,1000])
yticks([1,10,100,1000])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text(10^((1-shift)*log10(xrange(1)) + shift*log10(xrange(2))) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'c','fontweight','bold','fontsize',12)


%% Simulation exponents
disp(strcat2({'\tau = ', Sim.avalanche.sizeFit.tau, '+/-', Sim.avalanche.sizeFit.dTau}));
disp(strcat2({'\alpha = ', Sim.avalanche.timeFit.alpha, '+/-', Sim.avalanche.timeFit.dAlpha}));
dx1 = sqrt((Sim.avalanche.timeFit.dAlpha/(Sim.avalanche.timeFit.alpha - 1))^2 + (Sim.avalanche.sizeFit.dTau/(Sim.avalanche.sizeFit.tau - 1))^2) *Sim.avalanche.gamma.x1;
disp(strcat2({'S,T: 1/\sigma \nu z = ', Sim.avalanche.gamma.x1, '+/-', dx1}));
disp(strcat2({'<S>(T): 1/\sigma \nu z = ', Sim.avalanche.gamma.x2, '+/-', Sim.avalanche.gamma.dx2}));
disp(strcat2({'Suc = ', Sim.avalanche.sizeFit.uc, ', Tuc = ', Sim.avalanche.timeFit.uc}))

