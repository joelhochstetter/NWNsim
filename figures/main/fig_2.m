%% Run simulations varying V 
%% Vsweep with Initial state of all junctions 0
T = 10; %set simulation time
uninitFolder = 'uninit'; %folder for uninitialised network
%make sure these simulations are only sims saved in uninitFolder
%extra datapoints specified near transition
Amps =  unique(sort([0.010:0.005:0.080, 0.09003:0.00001:0.09004, 0.09006:0.00001:0.09009, 0.085:0.0005:0.1000, 0.101:0.001:0.149, 0.15:0.01:0.5, 0.55:0.05:3.0, 0.101:0.005:0.225]));
DC_Vsweep(uninitFolder, Amps, T)
 

%% Run pre-activated DC simulation at V = 1.8V for 10 sec
saveFolder = '.'; %folder to save pre-initialisation of network
DC_Vsweep(saveFolder, 1.8, 10);


%% Run Vsweep initial state from pre-activated network
preinitFolder = 'preinit'; %folder for pre-initialised network
%make sure these simulations are only sims saved in preinitFolder
preFileName = 'tl_T10_DC1.8V_s0.01_r0.005_c0.01_m0.015_b10_p1.mat'; %name of sim file we pre-initalised to
sort(unique([0.01:0.005:0.27, 0.041:0.001:0.044, 0.045:0.001:0.090]));
DC_Vsweep(preinitFolder, Amps, T, -1, preFileName, saveFolder)


%% Import simulations
%must ensure have not saved extra simulations besides relevant in folders
simsU = multiImport(struct('importAll', true, 'SimOpt', struct('saveFolder', uninitFolder)));

simsP = multiImport(struct('importAll', true, 'SimOpt', struct('saveFolder', preinitFolder)));


%% Process data (extract voltages and Gend)
%uninitialised
Vlist    = zeros(size(simsU)); 
Gend = zeros(size(simsU));

for i = 1:numel(simsU)
    Vlist(i) = simsU{i}.Stim.Amplitude;
    Gend(i) = simsU{i}.netC(end);    
end

%pre-initialised
Vlist1    = zeros(size(simsP)); 
Gend1 = zeros(size(simsP));
for i = 1:numel(simsP)
    Vlist1(i)   = simsP{i}.Stim.Amplitude;
    Gend1(i) = simsP{i}.netC(end);    
end


%% Plot Fig 2a
timeVec = 1e-3:1e-3:T;

idx = zeros(9,1);
idx(1) = find(Vlist >= 0.09*0.7,1);
idx(2) = find(Vlist >= 0.0899,1);
idx(3) = find(Vlist >= 0.09*1.1,1);
idx(4) = find(Vlist >= 0.09*1.2,1);
idx(5) = find(Vlist >= 0.09*1.5,1);
idx(6) = find(Vlist >= 0.09*2.0,1);
idx(7) = find(Vlist >= 0.09*3.0,1);
idx(8) = find(Vlist >= 0.09*4.0,1);
idx(9) = find(Vlist >= 0.09*5.0,1);
leg = {};
cmap = parula(10);
figure;
for i =1:numel(idx) 
    semilogy(timeVec, simsU{idx(i)}.netC/7.77e-5, '-', 'color', cmap(i,:))
    leg{i} = num2str(Vlist(idx(i))/0.09, '%.1f');
    hold on;
end

xlim([0,7])
ylim([4e-4, 1])
xlabel('t (s)', 'fontweight','normal','fontsize',10)
ylabel('G_{nw} (G_0)', 'fontweight','normal','fontsize',10)
leg = legend(leg, 'location', 'east');
pos = get(leg,'position');
pos(1) = 0.57;
pos(2) = 0.29;
set(leg,'position', pos)
title(leg, 'V^*', 'fontweight', 'normal')
set(findall(gca, 'Type', 'Line'),'LineWidth',2.0);
axis square;


%% Plot Fig 2b
%   See Fig 3(c-e) (in 'fig_3.m) but use last time-point to make snapshot


%% Plot Fig 2c
figure;
plot(Vlist/0.01/9, Gend/7.77e-5*9, 'ko', 'MarkerSize', 6, 'LineWidth', 1.0)
hold on;
plot(Vlist1/0.01/9, Gend1/7.77e-5*9, 'r^', 'MarkerSize', 6, 'LineWidth', 1.0)
xlim([0,3.0])
ylim([0,4.7])
hline(1, 'k:')
vline(1, 'k:')
vline(0.5, 'k:')
text(0.9, 2, 'n V_{set}', 'Rotation', 90, 'fontsize', 10);
text(0.4, 2.001, 'n V_{reset}', 'Rotation', 90, 'Color', 'r', 'fontsize', 10);
xlabel('V^*','fontweight','normal','fontsize',10)
ylabel('G_{nw}^* (t = \infty)','fontweight','normal','fontsize',10)
leg = legend('Inactive', 'Preactivated', 'location', 'northwest');
axis square;


%% Inset for Fig2c
figure;
set(gcf, 'color', 'w');
figSize = [0 0 10 5];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm
set(gcf, 'Units', 'centimeters', 'Position',figSize);

semilogy(Vlist/9/0.01, Gend/7.77e-5*9, 'ko', 'MarkerSize', 6, 'LineWidth', 1.0)
hline(1, 'k--')
vline(1, 'k--')
% xlabel('t (s)','fontweight','bold','fontsize',10)
% ylabel('G (S)','fontweight','bold','fontsize',10)
xlim([0.7,1.3])
ylim([0.9e-3, 11])
xticks([0.7:0.3:1.3])
yticks([1e-3, 1e-1, 10])
axis square;

