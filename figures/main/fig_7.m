%% Run attractor simulations as in Fig 6
% runACattractor(Amps, Freqs, dt, T, attractorFolder) % Run simulations for attractors


%% Load relevant attractors
attractorFolder = '.'; %set save folder for attractor
attractorName1 = 'tl_T3000_ACsaw1.25V_f0.1Hz_s0.01_r0.01_c0.01_m0.015_b10_p1.mat';
attractorName2 = 'tl_T3000_ACsaw1.25V_f0.5Hz_s0.01_r0.01_c0.01_m0.015_b10_p1.mat';
attractorName3 = 'tl_T3000_ACsaw1.25V_f0.85Hz_s0.01_r0.01_c0.01_m0.015_b10_p1.mat';

params = struct('SimOpt', struct('saveFolder', attractorFolder));
params.importByName = {attractorName1, attractorName2, attractorName3};
sims = multiImport(params);


%% Extract last 20 periods from attractor
TtoPlot = 20; %specify number of periods to plot
netV = cell(3,1); %voltage time-series
netI = cell(3,1);  %current time-series
netC = cell(3,1);  %conductance time-series

for i = 1:3
    netV{i} = sims{i}.Stim.Signal;
    netI{i}  = sims{i}.netI;
    nstepT = round(1/(sims{i}.dt*sims{i}.Stim.Frequency));  %number of time-steps in a period
    netV{i} = netV{i}(end - TtoPlot*nstepT+1:end); %extract number of periods for voltage
    netI{i}  = netI{i}(  end - TtoPlot*nstepT+1:end); %extract number of periods for current
    netC{i} = netI{i}./netV{i};
end
    

%% Plot figure 6

Vlim = [-1.5, 1.5];
Glim = [0, 5e-5]/7.77e-5;

figure;
lw = 1.5; %line-width

subplot(2,3,1);
plot(netV{1}, netI{1}*1e6, 'LineWidth', lw)
title('\lambda < 0')
xlabel('V (V)','fontweight','normal','fontsize',10)
ylabel('I (\mu A)','fontweight','normal','fontsize',10)
ylim([-50,50])
xlim(Vlim)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'a','fontweight','bold','fontsize',12)
hold on;
x = 0.15+ [0 0 0.15];
y = [-3 3 0];
patch(x,y,'red')
x = -0.15+ [0 0 -0.15];
y = [-3 3 0];
patch(x,y,'red')
[x, y] = arrow([0.35,15], [-1.5, -2], 0.18, 6);
patch(x,y,'red')
[x, y] = arrow([-0.35,-15], [1.5, 2], 0.18, 6);
patch(x,y,'red')


subplot(2,3,2);
plot(netV{2}, netI{2}*1e6, 'LineWidth', lw)
title('\lambda \approx 0')
xlabel('V (V)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim([-50,50])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'b','fontweight','bold','fontsize',12)
hold on;
x = 0.5+ [0 0 0.15];
y = [-3 3 0];
patch(x,y,'red')
x = -0.5+ [0 0 -0.15];
y = [-3 3 0];
patch(x,y,'red')
[x, y] = arrow([0.6,24], [-1.5, -2], 0.18, 6);
patch(x,y,'red')
[x, y] = arrow([-0.6,-24], [1.5, 2], 0.18, 6);
patch(x,y,'red')



subplot(2,3,3);
plot(netV{3}, netI{3}*1e6, 'LineWidth', lw)
title('\lambda > 0')
xlabel('V (V)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim([-50,50])
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'c','fontweight','bold','fontsize',12)
hold on;
[x, y] = arrow([0.7,0], [1, 0], 0.18, 6);
patch(x,y,'red')
[x, y] = arrow([-0.7,0], [-1, 0], 0.18, 6);
patch(x,y,'red')
patch(x,y,'red')
[x, y] = arrow([0.6,20], [-2, -1.5], 0.18, 6);
patch(x,y,'red')
[x, y] = arrow([-0.6,-20], [2, 1.5], 0.18, 6);
patch(x,y,'red')

subplot(2,3,4);
plot(netV{1}, netC{1}/7.77e-5, 'LineWidth', lw)
xlabel('V (V)','fontweight','normal','fontsize',10)
ylabel('G_{nw} (G_0)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim(Glim)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'd','fontweight','bold','fontsize',12)
hold on;
[x, y] = arrow([0.15,0], [1, 0], 0.18, 0.05);
patch(x,y,'red')
[x, y] = arrow([-0.15,0], [-1, 0], 0.18, 0.05);
patch(x,y,'red')
patch(x,y,'red')
[x, y] = arrow([0.6,0.52], [-1, 0], 0.18, 0.05);
patch(x,y,'red')
[x, y] = arrow([-0.6,0.52], [1, 0], 0.18, 0.05);
patch(x,y,'red')

subplot(2,3,5);
plot(netV{2}, netC{2}/7.77e-5, 'LineWidth', lw)
xlabel('V (V)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim(Glim)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'e','fontweight','bold','fontsize',12)
hold on;
[x, y] = arrow([0.6,0], [1, 0], 0.18, 0.05);
patch(x,y,'red')
[x, y] = arrow([-0.6,0], [-1, 0], 0.18, 0.05);
patch(x,y,'red')
patch(x,y,'red')
[x, y] = arrow([0.7,0.51], [-1, 0], 0.18, 0.05);
patch(x,y,'red')
[x, y] = arrow([-0.7,0.51], [1, 0], 0.18, 0.05);
patch(x,y,'red')


subplot(2,3,6);
plot(netV{3}, netC{3}/7.77e-5, 'LineWidth', lw)
xlabel('V (V)','fontweight','normal','fontsize',10)
xlim(Vlim)
ylim(Glim)
xrange = xlim;
yrange = ylim;
shift = -0.17;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'f','fontweight','bold','fontsize',12)
hold on;
[x, y] = arrow([0.8,0], [1, 0], 0.18, 0.05);
patch(x,y,'red')
[x, y] = arrow([-0.8,0], [-1, 0], 0.18, 0.05);
patch(x,y,'red')
patch(x,y,'red')
[x, y] = arrow([0.5,0.44], [-1, 0], 0.18, 0.05);
patch(x,y,'red')
[x, y] = arrow([-0.5,0.45], [1, 0], 0.18, 0.05);
patch(x,y,'red')