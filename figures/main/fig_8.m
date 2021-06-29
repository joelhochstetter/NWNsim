%% Load relevant Lyapunov exponents as in (figure 6)
numA = 225; %number of attractor files
E = 261; %number of perturbed junctions that Lyapunov simulations 
li = zeros(numA, E); %junction Lyapunov exponents for each simulation
ml = zeros(numA, 1);
Amp = zeros(numA, 1); %amplitude of simulation stimulus
Freq = zeros(numA, 1); %frequency of simulation stimulus

lyFolder = 'lyapunov'; %folder to save Lyapunov exponent simulations
files = dir(strcat(lyFolder, '/tl*')); %get Lyapunov exponents simulations folders

%need to ensure folders only contains these files otherwise will not work

for i = 1:numA
    if exist(strcat(files(i).folder, '/', files(i).name,'/LyCalc.mat'), 'file')
        f = load(strcat(files(i).folder, '/', files(i).name,'/LyCalc.mat'));
        li(i,:) = f.li; 
        ml(i) = f.ml;     
        Amp(i) = f.params.Stim.Amplitude;
        Freq(i) =f.params.Stim.Frequency;        
    end
end

%% Run stimulus for initial period to relax to attractor (as in figure 6)
attractorFolder = 'attractors';  %name of folder containing networks relaxed to attractor


%% Run non-linear transformation task
saveFolder = 'NLT'; %name of folder to save simulation
NLTres = attractorToNLT(attractorFolder, saveFolder);
save(strcat(attractorFolder, '/', saveFolder, '/NLTresults.mat'), 'NLTres')


%% Import non-linear transformation task
load(strcat(attractorFolder, '/', saveFolder, '/NLTresults.mat'), 'NLTres')

dblAcc_2d = zeros(numel(Amps), numel(Freqs));
phsAcc_2d = zeros(numel(Amps), numel(Freqs));
squAcc_2d  = zeros(numel(Amps), numel(Freqs));
sinAcc_2d  = zeros(numel(Amps), numel(Freqs));

sinAcc = 1- NLTres.sinRNMSE;
squAcc = 1- NLTres.squRNMSE;
dblAcc = 1- NLTres.dblRNMSE;
phsAcc = 1- NLTres.phsRNMSE;


for i = 1:numel(NLTres.Amps)
    idxA = find(NLTres.Amps(i) == Amps);
    idxF = find(NLTres.Freqs(i) == Freqs);    
    dblAcc_2d(idxA, idxF) = 1- NLTres.dblRNMSE(i);
    phsAcc_2d(idxA, idxF) = 1- NLTres.phsRNMSE(i);
    squAcc_2d(idxA, idxF) = 1- NLTres.squRNMSE(i);
    sinAcc_2d(idxA, idxF) = 1- NLTres.sinRNMSE(i);
end

%convert to 2d array
ml1_2d = zeros(numel(Amps), numel(Freqs));
for i = 1:numel(files)
    idxA = find(Amp(i) == Amps);
    idxF = find(Freq(i) == Freqs);    
    ml1_2d(idxA, idxF) = ml(i);
end


%% Plot figure
idx = r >= 10;

figure;
set(gcf, 'color', 'w');
hold on;
lxx = -40:0.1:40;
byy = 0.88*ones(size(lxx));
gyy = 0.50*ones(size(lxx));
kyy = 0.0*ones(size(lxx));
plot([0,0], [0,1], 'k:', 'Linewidth', 1.5, 'HandleVisibility', 'off')

plot(ml(idx), sinAcc(idx), 'b.', 'MarkerSize', 10);
plot(ml(idx), squAcc(idx), 'g.', 'MarkerSize', 10);
plot(ml(idx), phsAcc(idx), 'm.', 'MarkerSize', 10);
plot(ml(idx), dblAcc(idx),  'k.', 'MarkerSize', 10);
xlabel('\lambda (s^{-1})','fontweight','normal','fontsize',10)
ylabel('Accuracy','fontweight','normal','fontsize',10)
lg = legend('sin', 'square',  '\pi/2-shift', '2f', 'location', 'west');
set(lg, 'FontSize', 7);
xlim([-40,40])
ylim([0,1])
xrange = xlim;
yrange = ylim;
pos = get(lg,'Position');
set(lg, 'Position', [pos(1) - 0.005, pos(2) - 0.02, pos(3) - 0.01, pos(4) - 0.01]) 
shift = 0.06;
yticks([0:0.1:1.0]);
ytl = string(yticklabels);
ytl(2:2:end) = ' ';
yticklabels(ytl);
xticks([-40:5:40]);
xtl = string(xticklabels);
xtl(2:4:end) = ' ';
xtl(3:4:end) = ' ';
xtl(4:4:end) = ' ';
xticklabels(xtl);
box on;