%% Load relevant attractors (as in Fig 7):
%specify frequencies and amplitudes
Freq = [0.1, 0.5, 0.85];
Amp = 1.25*[1,1,1];

attractorFolder = '.'; %set save folder for attractor

%specify relevant filenames
attractorsToUse = cell(3,1);
for i = 1:3
    attractorsToUse{i} = strcat2({'t2_T3000_ACsaw', Amp(i) 'V_f', Freq(i), 'Hz_s0.01_r0.01_c0.01_m0.015_b10_p1.mat'});
end
%importSwitch required to get data to plot. Must have been run with
%params.SimOpt.saveSwitches         = true

params = struct('importSwitch', true, 'SimOpt', struct('saveFolder', attractorFolder));
params.importByName = attractorsToUse;
sims = multiImport(params);


%% get sd shortest path
contacts = sims{1}.ContactNodes; %source-drain electrode pairing
ConnectFile = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
load(ConnectFile, 'adj_matrix');
sp  = kShortestPath(double(adj_matrix), contacts(1), contacts(2), 1);
Connectivity          = struct('filename', ConnectFile);
Connectivity          = getConnectivity(Connectivity);
spE = getPathEdges(sp{1}, Connectivity.EdgeList);



%% Get junction conductances and voltages
spCon = cell(3,1);
spVol   = cell(3,1);
netV     = cell(3,1);
netC     = cell(3,1);

for j = 1:3
    spCon{j} =  sims{j}.swC(:,spE);
    spVol{j}   =  sims{j}.swV(:,spE); 
    netV{j}     =  sims{j}.Stim.Signal;
    netC{j}     =  sims{j}.netC;
end


%% Plot supplementary Fig 10
Vlim = [-1.5, 1.5];
Glim = [0.5e-4, 1];
timeVec =(1:numel(spCon{1}(:,1)))*sims{1}.dt;

figure;
set(gcf, 'color', 'w');
figSize = [0 0 20 15];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); 
set(gcf, 'Units', 'centimeters', 'Position',figSize);

subplot(3,3,1);
semilogy(timeVec*Freq(1),(spCon{1})/7.77e-5)
xlim([0,1]);
ylim([1e-4,1.5])
yticks([1e-4, 1e-2, 1])
xlabel('t/T','fontweight','normal','fontsize',10)
ylabel('G_{jn} (G_0)','fontweight','normal','fontsize',10)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'a','fontweight','bold','fontsize',12)
yyaxis right;
plot(timeVec*Freq(1), sims{1}. Stim.Signal, 'k-')
ax = gca;
ax.YAxis(2).Color = 'black';
title('\lambda < 0')

subplot(3,3,2);
semilogy(timeVec*Freq(2),spCon{2}/7.77e-5)
xlim([0,1]);
ylim([1e-4,1.5])
yticks([1e-4, 1e-2, 1])
xlabel('t/T','fontweight','normal','fontsize',10)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'b','fontweight','bold','fontsize',12)
yyaxis right;
plot(timeVec*Freq(2), sims{2}. Stim.Signal, 'k-')
ax = gca;
ax.YAxis(2).Color = 'black';
title('\lambda \approx 0')

subplot(3,3,3);
semilogy(timeVec*Freq(3),(spCon{3})/7.77e-5)
set(gca, 'YScale', 'log')
xlim([0,1]);
ylim([1e-4,1.5])
yticks([1e-4, 1e-2, 1])
xlabel('t/T','fontweight','normal','fontsize',10)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'c','fontweight','bold','fontsize',12)
yyaxis right;
plot(timeVec*Freq(3), sims{3}. Stim.Signal, 'k-')
ax = gca;
ax.YAxis(2).Color = 'black';
ylabel('V (V)', 'fontweight', 'normal','fontsize',10)
title('\lambda > 0')


subplot(3,3,4);
plot(timeVec*Freq(1),abs(spVol{1}))
xlim([0,1]);
xlabel('t/T','fontweight','normal','fontsize',10)
ylabel('|V_{jn}| (V)','fontweight','normal','fontsize',10)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'd','fontweight','bold','fontsize',12)
yyaxis right;
plot(timeVec*Freq(3), sims{3}. Stim.Signal, 'k-')
ax = gca;
ax.YAxis(2).Color = 'black';

subplot(3,3,5);
plot(timeVec*Freq(2),abs(spVol{2}))
xlim([0,1]);
xlabel('t/T','fontweight','normal','fontsize',10)
ylim([0,0.6])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'e','fontweight','bold','fontsize',12)
yyaxis right;
plot(timeVec*Freq(3), sims{3}. Stim.Signal, 'k-')
ax = gca;
ax.YAxis(2).Color = 'black';

subplot(3,3,6);
plot(timeVec*Freq(3),abs(spVol{3}))
xlim([0,1]);
xlabel('t/T','fontweight','normal','fontsize',10)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,(shift*(yrange(1)) + (1-shift)*(yrange(2))), 'f','fontweight','bold','fontsize',12)
yyaxis right;
plot(timeVec*Freq(3), sims{3}. Stim.Signal, 'k-')
ax = gca;
ax.YAxis(2).Color = 'black';
ylabel('V (V)', 'fontweight', 'normal','fontsize',10)

subplot(3,3,7);
semilogy(netV{1}, netC{1}/7.77e-5, 'LineWidth', lw)
xlabel('V (V)','fontweight','normal','fontsize',10)
ylabel('G_{nw} (G_0)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim(Glim)
yticks([1e-4, 1e-2, 1])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'g','fontweight','bold','fontsize',12)

subplot(3,3,8);
semilogy(netV{2}, netC{2}/7.77e-5, 'LineWidth', lw)
xlabel('V (V)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim(Glim)
yticks([1e-4, 1e-2, 1])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'h','fontweight','bold','fontsize',12)


subplot(3,3,9);
semilogy(netV{3}, netC{3}/7.77e-5, 'LineWidth', lw)
xlabel('V (V)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim(Glim)
yticks([1e-4, 1e-2, 1])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) ,10^(shift*log10(yrange(1)) + (1-shift)*log10(yrange(2))), 'i','fontweight','bold','fontsize',12)
