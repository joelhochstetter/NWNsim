%% Run simulations and criticality analysis (as in Fig 4):


%% Load criticality results file
Sim = load('Av/density0.10/Vstar1/Lx150/bs-1/critResults.mat');
Sim = Sim.critResults;


%% Obtain p-value and exponent for each xmax for P(S), P(T)
SLs = 7; %lower-cutoff for size distribution
TLs = 7; %lower-cutoff for lfe-time

%extracts xmin, xmax, p, tau for P(S) for each xmax
Sim.avalanche.sizeFit.MLcompare.PL(1).fullResults = custom_plparams(Sim.avalanche.sizeAv, SLs);

sSxmin = Sim.avalanche.sizeFit.MLcompare.PL(1).fullResults.xmin;
sSxmax = Sim.avalanche.sizeFit.MLcompare.PL(1).fullResults.xmax;
sSp        = Sim.avalanche.sizeFit.MLcompare.PL(1).fullResults.p;
sStau     = Sim.avalanche.sizeFit.MLcompare.PL(1).fullResults.tau;

%keep only values for the fixed xmax, the others it did not run for 
sSxmax = sSxmax(sSxmin == SLs);
sStau = sStau(sSxmin == SLs);
sSp = sSp(sSxmin == SLs);
sSxmin = sSxmin(sSxmin == SLs);


%extracts xmin, xmax, p, alpha for P(T) for each xmax
 Sim.avalanche.timeFit.MLcompare.PL(1).fullResults = custom_plparams(Sim.avalanche.lifeAv, TLs);

sTxmin = Sim.avalanche.timeFit.MLcompare.PL(1).fullResults.xmin;
sTxmax = Sim.avalanche.timeFit.MLcompare.PL(1).fullResults.xmax;
sTp        = Sim.avalanche.timeFit.MLcompare.PL(1).fullResults.p;
sTtau     = Sim.avalanche.timeFit.MLcompare.PL(1).fullResults.tau;

%keep only values for the fixed xmax, the others it did not run for 
sTxmax = sTxmax(sTxmin == TLs);
sTtau = sTtau(sTxmin == TLs);
sTp = sTp(sTxmin == TLs);
sTxmin = sTxmin(sTxmin == TLs);




%%
figure;
semilogx(sSxmax, sStau)
hold on;
semilogx(sTxmax, sTtau)
ylim([1.8,2.4])
yticks([1.8:0.2:2.4])
ylabel('\tau, \alpha')

yyaxis right;
semilogx(sSxmax, sSp)
hold on;
semilogx(sTxmax, sTp, 'k')
legend('\tau', '\alpha', 'p_\tau', 'p_\alpha')
xlim([10,500])
ylabel('p-value')
xlabel('x_{max}')
