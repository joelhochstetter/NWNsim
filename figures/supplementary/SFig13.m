%% Change to directory where NLT simulations were saved
NLTfolder  = '.';
cd(NLTfolder);


%% Select amplitudes and frequencies of curves to plot
thisAmp = 1;
thisFreq = [0.25, 0.85, 1];



%% Plot Supp Fig 13
numT = 1.5;
trange = 1:80000;
nds = 1:10;


figure;
figSize = [0 0 20 6];
set(gcf, 'Units', 'centimeters', 'Position',figSize);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize);

load(strcat2({'tl_T80_ACsaw', thisAmp, 'V_f', thisFreq(1), 'Hz_s0.01_r0.01_c0.01_m0.015_b10_p1.mat'}))
f = thisFreq(1);
subplot(1, 3, 1);
plot(sim.Stim.TimeAxis(trange)*f, sim.nwV(trange,nds)./max(abs(sim.nwV(trange,nds))))
xlabel('t/T')
ylabel('V_{nd} rescaled')
xlim([0, numT])
xrange = xlim;
yrange = ylim;
shift = -0.09;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'a','fontweight','bold','fontsize',12)

title('\lambda < 0')

load(strcat2({'tl_T80_ACsaw', thisAmp, 'V_f', thisFreq(2), 'Hz_s0.01_r0.01_c0.01_m0.015_b10_p1.mat'}))
f = thisFreq(2);
subplot(1, 3, 2);
plot(sim.Stim.TimeAxis(trange)*f, sim.nwV(trange,nds)./max(abs(sim.nwV(trange,nds))))
xlabel('t/T')
xlim([0, numT])
xrange = xlim;
yrange = ylim;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'b','fontweight','bold','fontsize',12)
title('\lambda \approx 0')

load(strcat2({'tl_T80_ACsaw', thisAmp, 'V_f', thisFreq(3), 'Hz_s0.01_r0.01_c0.01_m0.015_b10_p1.mat'}))
f = thisFreq(3);
subplot(1, 3, 3);
plot(sim.Stim.TimeAxis(trange)*f, sim.nwV(trange,nds)./max(abs(sim.nwV(trange,nds))))
xlabel('t/T')
xlim([0, numT])
xrange = xlim;
yrange = ylim;
text((1-shift)*xrange(1) + shift*xrange(2) , shift*yrange(1) + (1-shift)*yrange(2), 'c','fontweight','bold','fontsize',12)
title('\lambda > 0')