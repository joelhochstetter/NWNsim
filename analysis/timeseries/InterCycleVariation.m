function ICV = InterCycleVariation(netC, dt, T, f, G0) 
%{
netC: network conductance
dt:   time-step
T: Length of simulation
f: Frequency
%}
    
    if nargin < 5
        G0 = 7.77e-5;
    end
    
	tstepT = round(1/dt/f);
	numT = floor(T*f);
	netCMat = zeros(numT, tstepT);
	
	for i = 1:numT
		netCMat(i,:) = netC((1:tstepT) + round((i - 1)*1/dt/f));	
	end
	
	ICV = (nanmean(nanmean((netCMat/G0).^2, 1) - nanmean((netCMat/G0),1).^2))/dt/f;
	ICV(ICV < 0) = 0;
    ICV = sqrt(ICV);
    
end