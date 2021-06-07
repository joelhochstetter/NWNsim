function Stimulus = getStimulus(Stimulus, SimulationOptions, includet0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a structure with the details of an external voltage signal
% applied to the network.

%
% ARGUMENTS: 
% Stimulus - Structure containing the details of the required stimulus. It
%           must contain a field .BiasType, whose value is a string
%           specifying the type of stimulus.:
%           - 'DC' - Constant external voltage. 

%                      .Amplitude
%           - 'AC' - Sinusoidal external voltage.
%                    Required fields:
%                      .Frequency
%                      .Amplitude
%           - 'DCandWait' - A DC signal followed by a much smaller DC
%                           signal (which is meant only for probing, rather 
%                           then for stimulating).
%                           Required fields:
%                             .T 
%                             .dt
%                             .OffTime % Time at which tthe stimulus changes to AmplitudeOff
%                             .AmplitudeOn
%                             .AmplitudeOff
%           - 'Ramp' - A ramping voltage signal.
%                      Required fields:
%                        .T
%                        .dt
%                        .AmplitudeMin
%                        .AmplitudeMax
%           - 'Custom' - An arbitrary voltage signal.
%                        Required fields:
%                          .T
%                          .dt
%                          .TimeAxis
%                          .Signal
% SimulationOptions -  a struct with additional information about the simulation
%                    Required fields:
%                      .T                           (signal duration)
%                      .dt                          (duration of time-step)
%           It is assumed that the units of all input fields are sec, Hz
%           and Volt. Thus, the output fields are in sec and Volt.
%
%
%   includet0 -  if we include  t = 0, and 0 otherwise (default 0)
%
%



% OUTPUT:
% Stimulus - Structure with the details of the time axis and with the 
%            external voltage signal. Fields:
%              .BiasType
%              .T
%              .dt
%              .TimeAxis
%              .Signal
%              .Frequency [Hz] (for AC and Sawtooth)
%
% REQUIRES:
% none
%
% USAGE:
%{
    Options.T          = 1e+1; % (sec)
    Options.dt         = 1e-3; % (sec)
    Stimulus.BiasType  = 'DC';       
    Stimulus.Amplitude = 3;    % (Volt)

    Stimulus = getStimulus(Stimulus, Options);
%}
%
%
% Authors:
% Ido Marcus
% Paula Sanz-Leon
% Joel Hochstetter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if nargin < 3
        includet0 = false;
    end

    % Signal duration and timestep:
    Stimulus.T  = SimulationOptions.T;
    Stimulus.dt = SimulationOptions.dt;
    Stimulus.TimeAxis = (Stimulus.dt:Stimulus.dt:Stimulus.T)';
    
    if includet0
        Stimulus.TimeAxis = [0; Stimulus.TimeAxis];
    end
    
    % = (Stimulus.dt : Stimulus.dt : Stimulus.T)';
    % The first time point we can actually records happens at dt. 
    % Time=0 correspond for the values of the intial conditions
    
    % External voltage signal:
    switch Stimulus.BiasType
        case 'DC'
            Stimulus.Signal = Stimulus.Amplitude*ones(size(Stimulus.TimeAxis));% + normrnd(0,0.2,size(Stimulus.TimeAxis));
            Stimulus.stimName   = strcat(Stimulus.BiasType,num2str( Stimulus.Amplitude,3),'V');
        case 'AC'
            Stimulus.Signal = Stimulus.Amplitude*sin(Stimulus.Phase + 2*pi*Stimulus.Frequency*Stimulus.TimeAxis);
            Stimulus.stimName = strcat(Stimulus.BiasType, num2str( Stimulus.Amplitude,3),'V_f',num2str(Stimulus.Frequency,3),'Hz');
        case 'ACsaw'
            Stimulus.stimName = strcat(Stimulus.BiasType, num2str( Stimulus.Amplitude,3),'V_f',num2str(Stimulus.Frequency,3),'Hz');
            Stimulus.Signal = Stimulus.Amplitude*sawtooth(Stimulus.Phase + 2*pi*Stimulus.Frequency*(Stimulus.TimeAxis-0.75/Stimulus.Frequency) , 0.5);
        case 'DCsaw'
            Stimulus.stimName = strcat(Stimulus.BiasType, num2str( Stimulus.Amplitude,3),'V_f',num2str(Stimulus.Frequency,3),'Hz');
            Stimulus.Signal = abs(Stimulus.Amplitude*sawtooth(2*pi*Stimulus.Frequency*(Stimulus.TimeAxis-0.75/Stimulus.Frequency) , 0.5)) + 1e-8;
        case 'DCandWait'
            % Standard boxcar:
            %Stimulus.Signal = Stimulus.AmplitudeOn*ones(size(Stimulus.TimeAxis));
            % rectangular pulse train:
            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*square(1*pi*Stimulus.TimeAxis/Stimulus.OffTime));
            % sawtooth (for I-V plots):
            %Stimulus.Signal = Stimulus.AmplitudeOn*sawtooth(2*pi*1.0*(Stimulus.TimeAxis-0.75/10.0) , 0.5);
%            % Walsh functions:
%            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*sign(cos(8*pi*Stimulus.TimeAxis).*sin(2*pi*Stimulus.TimeAxis))); % Walsh9
%            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*sign(cos(2*pi*Stimulus.TimeAxis).*cos(4*pi*Stimulus.TimeAxis)));
%            % Walsh6
%            Stimulus.Signal =
%            max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*sign(cos(2*pi*Stimulus.TimeAxis).*sin(2*pi*Stimulus.TimeAxis)));
%            spike:
%            Stimulus.Signal = Stimulus.AmplitudeOff*ones(size(Stimulus.TimeAxis));
%            Stimulus.Signal(Stimulus.TimeAxis >= 5) = max( Stimulus.AmplitudeOff, exp( -1.5*((Stimulus.TimeAxis(Stimulus.TimeAxis >= 5)/5)- 1) )*Stimulus.AmplitudeOn );
            %
            Stimulus.Signal(Stimulus.TimeAxis > Stimulus.OffTime) = Stimulus.AmplitudeOff;
            Stimulus.stimName   = strcat(Stimulus.BiasType,num2str(Stimulus.AmplitudeOn,3),'V_off',num2str(Stimulus.OffTime,3),'s', ...
               'offV', num2str(Stimulus.AmplitudeOff,3), 'V');
           
        case 'Square'
            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*square(1*pi*Stimulus.TimeAxis/Stimulus.OffTime + Stimulus.Phase, Stimulus.Duty));           
            Stimulus.stimName   = strcat(Stimulus.BiasType,num2str(Stimulus.AmplitudeOn,3),'V_off',num2str(Stimulus.OffTime,3),'s', ...
               'offV', num2str(Stimulus.AmplitudeOff,3), 'V');
           
           
        case 'Ramp'
            Stimulus.Signal = linspace(Stimulus.AmplitudeMin, Stimulus.AmplitudeMax, length(Stimulus.TimeAxis))';
            Stimulus.stimName = strcat(Stimulus.BiasType, 'max', num2str( Stimulus.AmplitudeMax,3),'V_min',num2str(Stimulus.AmplitudeMin,3),'V');              
            
        case 'Custom'
            Stimulus.Signal = Stimulus.Signal;
            Stimulus.stimName = 'custom';
            
        case 'Triangle'
            t  = SimulationOptions.TimeVector - SimulationOptions.dt - SimulationOptions.T/2; 
            w  = SimulationOptions.T;
            scaling = abs(Stimulus.AmplitudeMax - Stimulus.AmplitudeMin);
            Stimulus.Signal = scaling * tripuls(t, w) + Stimulus.AmplitudeMin;
            
            
        case 'AlonPulse'
            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*square(2*pi*Stimulus.TimeAxis/Stimulus.Period));
            Stimulus.Signal(Stimulus.TimeAxis >= Stimulus.NumPulse1*Stimulus.Period) = Stimulus.AmplitudeOff;
            SecondSetStart =  Stimulus.NumPulse1*Stimulus.Period + Stimulus.LongWait; %second set of pulses starts here
            Stim1          = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*square(2*pi*(Stimulus.TimeAxis-SecondSetStart)/Stimulus.Period));
            SecondSetEnd   = SecondSetStart + Stimulus.NumPulse2*Stimulus.Period;
            Stimulus.Signal(Stimulus.TimeAxis >= SecondSetStart) = Stim1(Stimulus.TimeAxis >= SecondSetStart);
            Stimulus.Signal(Stimulus.TimeAxis >= SecondSetEnd)   = Stimulus.AmplitudeOff;
            
            
        case 'SinglePulse'
           Stimulus.Signal = Stimulus.AmplitudeOff*ones(size(Stimulus.TimeAxis));
           Stimulus.Signal(Stimulus.TimeAxis > Stimulus.OnTime) = Stimulus.AmplitudeOn;
           Stimulus.Signal(Stimulus.TimeAxis > Stimulus.OffTime) = Stimulus.AmplitudeOff;
           
           %
            
    end
end
