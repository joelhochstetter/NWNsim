%% Set-up single junction parameters
Vset = 1e-2; %Vset
Vres = 5e-3; %Vreset
b    = 2;
Gon  = 7.77e-5; %on conductance
Goff = 7.77e-8; %off conductance
cLam = 0.01; %lambda crit
mLam = 0.015;  %lambda max
Amp = 0.025; %signal amplitude
Vsec = 0.005; %signal ramp rate (V/s)
freq = Vsec/Amp; %signal frequency
T = 1/freq; %period
dt = 1e-3;
lam1 = 7.5e-3; %initial filament state
Icc = 1e-6;  %compliance current


%% Evolve jnction over time
t = 0:dt:T;
tvec = t';
AC = Amp*sawtooth(2*pi*freq*(tvec-0.75/freq), 0.5); %AC triangular signal
conT = zeros(size(tvec));
lam = zeros(size(tvec));
curr = zeros(size(tvec));
Signal = AC;
lam(1) = lam1;

for i = 1: numel(tvec)
    d = (cLam-abs(lam(i)))*5/cLam;
    if d < 0.0
        d = 0.0;
    end
    
    conT(i) = tunnelSwitchL(d, 0.81, 0.17, Goff, Gon);
    
    V = Signal(i);
    
    if abs(V)*conT(i) > Icc
        V = Icc/conT(i)*sign(Signal(i));
    end
    
    curr(i) = conT(i)*V;    
    
    if i < numel(tvec)
        lam(i+1) = lam(i) + (abs(V)-Vset)*dt*(abs(V) > Vset)*sign(V) + b*(abs(V)-Vres)*dt*(abs(V) < Vres)*sign(lam(i));
        if abs(lam(i+1)) >= mLam 
            lam(i+1) = sign(lam(i+1))*mLam;
        end
    end
end



%% Plot SFig 19b
figure;
figSize = [0 0 13.4 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', figSize); %x_width=10cm y_width=15cm
set(gcf, 'Units', 'centimeters', 'Position',figSize)
plot(Signal*1000,(curr/1e-6),'LineWidth',2);
hold on;
ylim([-1, 1.5])
xlim([-30, 30])
xlabel 'V (mV)'
ylabel 'I (\mu A)'