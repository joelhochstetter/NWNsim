%% Run pre-activated DC simulation at V = 1.8V for 10 sec (as in Fig 2)
%{
to reproduce figure exactly must change offConductance to 1e-8.
To do this set 'params.Comp.offConductance   = 1e-8;' 
    on line 204 of DC_Vsweep.m
Must set back to 7.77e-8 for the rest of the simulations
%}

preActFolder ='.'; %folder to save pre-initialisation of network
DC_Vsweep(preActFolder, 1.8, 10, -1, -1, -1, -1, false, false, 0, -1, true, false);
presim = multiImport(struct('SimOpt', struct('saveFolder', preActFolder), 'importByName', 'tl_T10_DC1.8V_s0.01_r0.005_c0.01_m0.015_b10_p1.mat', 'importStateOnly', true));


%% Specify conductances of initial states
numInitStates = 15;
Gmin  = presim{1}.netC(1); %minimum conductance of initial state to use in simulations
Gmax = presim{1}.netC(end); %maximum conductance of initial state to use in simulations
initCons = 10.^linspace(log10(Gmin), log10(Gmax), numInitStates); %make log(netC) of init states equally spaced


%% Run simulations for changing initial state
saveFolder = 'VaryVChangeInitState';
Amps = 0.01*9*[0.5:0.025:3.0];
T = 1; % 1000;

%loop over initial conductance states
for initCon = initCons 
    DC_Vsweep(saveFolder, Amps, T, -1, 'tl_T10_DC1.8V_s0.01_r0.005_c0.01_m0.015_b10_p1.mat' , preActFolder, initCon);
end


%% Import simulations
sims = multiImport(struct('importAll', true, 'SimOpt', struct('saveFolder', saveFolder, 'useParallel', true)));


%% Process simulations
Vlist    = zeros(size(sims)); %voltage
Gend = zeros(size(sims)); %steady state conductance
Ginit   = zeros(size(sims)); %initial conductance

for i = 1:numel(sims)
    Vlist(i) = sims{i}.Stim.Amplitude;
    Gend(i) = sims{i}.netC(end);
    [~,I] = min(abs(sims{i}.netC(1) - initCon));
    fname = split(sims{i}.filename, 'S.mat');
    fname = split(fname(1), 'initCon');
    Ginit(i)   = str2double(char(fname(2)));
end

initCon = sort(unique(Ginit));
Vvals    = sort(unique(Vlist));
Vvals = [0.0025:0.00125:0.0425, Vvals];
initCon_2d  = zeros(numel(Vvals), numel(initCon)); 
Volt_2d       = zeros(numel(Vvals), numel(initCon));
Gend_2d    = zeros(numel(Vvals), numel(initCon));

for i = 1:numel(sims)
    idxV = find(Vlist(i) == Vvals);
    idxG = find(Ginit(i) == initCon);  
    initCon_2d(idxV, idxG) = initCon(idxG);    
    Gend_2d(idxV, idxG)   =Gend(i);   
    Volt_2d(idxV, idxG)      = Vvals(idxV);    
end

Volt_2d = Vvals'*ones(1, numel(initCon));



%%
col =  jet(numel(initCon));
legEntry = {};
figure;
set(gcf, 'color', 'w');
hold on;
for i = 1:numel(initCon)
    plot(Volt_2d(:,i)/0.01/9, Gend_2d(:,i)/7.77e-5*9, '.', 'MarkerSize', 17, 'LineWidth', 1.0, 'Color', col(i,:))
    legEntry{i} = num2str(initCon(i)/7.77e-5*9, '%.1e');
end
xlim([0,3.5])
ylim([0,4.7])
hline(1, 'k:')
vline(1, 'k:')
vline(0.5, 'k:')
xlabel('V^*','fontweight','normal','fontsize',10)
ylabel('G_{nw}^* (t = \infty)','fontweight','normal','fontsize',10)
leg = legend(legEntry, 'location', 'southeast');
title(leg, 'G_{init}^*', 'fontweight', 'normal')
box on;