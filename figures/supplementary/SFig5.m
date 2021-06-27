%% Run simulations
%Connectivity file is provided in 'neuromorphicNWN/generation'
connFile = '2020-04-02-121953_asn_nw_00540_nj_01079_seed_001_avl_07.00_disp_01.60_gns_05.00_cdisp_700.00.mat';
plengths = 5:15:65; %specify source-drain path lengths
saveFolder = 'n540j1079/'; %name folder to same simulations
Amps = (0.005:0.0005:0.03); %these are re-scaled by path length
rescalePLength = true; %multiplies amplitudes by source-drain path length
T = 1000; % simulation time

for p = plengths
    saveFolder1 = strcat(saveFolder, 'sd', num2str(p));
    DC_Vsweep(saveFolder1, Amps, T, connFile,  0, '', -1, rescalePLength, false, 0, p, false, false, 1)
end
    


%%
plengths = 5:15:65;
nSims = numel(Amps);
Vlist      = zeros(nSims, numel(plengths)); 
Gend    = zeros(nSims, numel(plengths));


for j = 1:numel(plengths)
    p = plengths(j);
    params = struct();
    params.importAll = true;
    params.SimOpt.useParallel = true;
    params.SimOpt.saveFolder = strcat(saveFolder, 'sd', num2str(p));
    sims = multiImport(params);
    
    timeVec = sims{1}.Stim.TimeAxis;
    for i = 1:numel(sims)
        Vlist(i,j)   = sims{i}.Stim.Amplitude;
        Gend(i,j) = sims{i}.netC(end);
    end
end


%% Plot SFig 5
figure;
set(gcf, 'color', 'w');
plot(Vlist/0.01./plengths, Gend/7.77e-5,  'o', 'MarkerSize', 6, 'LineWidth', 1.0)
vline(1, 'k:')
xlabel('V^*','fontweight','normal','fontsize',10)
ylabel('G_{nw}/G_0 (t = \infty)','fontweight','normal','fontsize',10)
leg = legend(string(plengths), 'location', 'northwest');
title(leg, 'n', 'fontweight', 'normal');
axis square;


%% Plot SFig 5 inset
figure;
set(gcf, 'color', 'w');
%set figure size
figSize = [0 0 4 4];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm
set(gcf, 'Units', 'centimeters', 'Position',figSize);
plot(Vlist/0.01./plengths, Gend/7.77e-5.*plengths,  'o', 'MarkerSize', 6, 'LineWidth', 1.0)
xlim([0.5, 1.5])
hline(1, 'k:')
vline(1, 'k:')
xlabel('V^*','fontweight','normal','fontsize',10)
ylabel('G_{nw}^* (t = \infty)','fontweight','normal','fontsize',10)
axis square;