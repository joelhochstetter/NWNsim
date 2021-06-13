%% Run network under AC stimulus for long time to allow convergence to attractor
%set-up parameters to run
dt = 5e-4; %set time-step for simulations
T  = 3000; %Set simulation time: Must run an integer number of periods for all frequencies
Amps = [0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3]; %Specify stimulus amplitudes
Freqs = [0.025, 0.05, 0.1, 0.15, 0.25, 0.35, 0.5, 0.65, 0.75, 0.85, 1.0, 1.25, 1.5, 1.75, 2.0]; %Specify stimulus frequencies
attractorFolder = 'attractors'; %set save folder for attractor
runACAttractor(Amps, Freqs, dt, T, attractorFolder) % Run simulations for attractors



%% Calculate Lyapunov exponents of attractor
%set parameters for simulations
eps   = 5e-4; %size of infintesimal perturbation 
dt    = 5e-4; %size of simulation time-step
T     = 200;

files = dir(strcat(attractorFolder, '/*.mat')); %get attractor simulation files
numA = numel(files); %number of attractor files
lyFolder = 'lyapunov'; %folder to save Lyapunov exponent simulations

for i = numA
    calcLyapunov(attractorFolder, files(i).name, lyFolder, eps, dt, T)
end


%% Load relevant Lyapunov exponents
E = 261; %number of junctions simulations are run for
li = zeros(numA, E); %junction Lyapunov exponents for each simulation
ml = zeros(numA, 1);
Amp = zeros(numA, 1); %amplitude of simulation stimulus
Freq = zeros(numA, 1); %frequency of simulation stimulus
files = dir(strcat(lyFolder, '/tl*')); %get Lyapunov exponents simulations folders
%need to ensure files only contains these folders otherwise will not work

for i = 1:numA
    if exist(strcat(files(i).folder, '/', files(i).name,'/LyCalc.mat'), 'file')
        f = load(strcat(files(i).folder, '/', files(i).name,'/LyCalc.mat'));
        li(i,:) = f.li; 
        ml(i) = f.ml;     
        Amp(i) = f.params.Stim.Amplitude;
        Freq(i) =f.params.Stim.Frequency;        
    end
end


%% Converts maximal Lyapunov exponent to 2D array for plotting
ml_2d = zeros(numel(Amps), numel(Freqs));
for i = 1:numel(files)
    idxA = find(Amp(i) == Amps);
    idxF = find(Freq(i) == Freqs);    
    ml_2d(idxA, idxF) = ml(i);  
end


%% Plot Fig 6a
crange = [-15, 15]; %set range for colour-bar
figure;
subplot(1,2,1);
imagesc([0.02, 2], [0.1,3], ml_2d);
ylabel('A (V)','fontsize',10)
xlabel('f (Hz)','fontsize',10)
set(gca,'YDir','normal')
inferno1 = inferno;
colormap(inferno1);
caxis(crange);
cb = colorbar;
cb.Title.String = '\lambda (s^{-1})';
cb.Title.FontSize = 9;
xrange = xlim;
yrange = ylim;
xticks([0:0.5:2]);
shift = -0.07;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'a','fontweight','bold','fontsize',12)


%% Calculate r (average ratio of on and off conductances)
r = zeros(numA, 1);
files = dir(strcat(lyFolder, '/tl*')); %get Lyapunov exponents simulations folders
%need to ensure files only contains these folders otherwise will not work

for i = 1:numel(files)
        f1 =    dir(strcat(files(i).folder, '/', files(i).name, '/*unperturbed.mat'));
        if numel(f1) == 0
            disp('WARNING: no file found')
            continue
        end
        sim = multiImport(struct('importSwitch', false, 'SimOpt', struct('saveFolder', strcat(files(i).folder, '/', files(i).name)), 'importByName', f1(1).name));
        r(i)      = aveMaxMin(sim{1}.netC, sim{1}.Stim.Frequency, sim{1}.dt);
end


%% Plot Fig 6b
subplot(1,2,2);
semilogy(ml, r, '.', 'MarkerSize', 10);
ylabel('r','fontsize',10) %'fontweight','bold'
xlabel('\lambda (s^{-1})','fontsize',10) %'fontweight','bold'
xlim([-20,20])
xrange = xlim;
yrange = ylim;
xticks([-20:10:20])
shift = -0.07;
text((1-shift)*xrange(1) + (shift)*xrange(2), 10^((shift)*log10(yrange(1)) + (1-shift)*log10(yrange(2))) , 'b','fontweight','bold','fontsize',12)


