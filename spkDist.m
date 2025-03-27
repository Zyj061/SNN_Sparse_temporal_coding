function dRossum=spkDist(tSpike1,tSpike2,tmax,tau1)
% vR distance of two filtered spike trains with a kernel

tau2=tau1/4;
dt = tau1/100;

tk = 0:dt:5*tau1;   % kernel window
t = 0:dt:tmax;

spiketrain1 = zeros(size(t));
spiketrain2 = zeros(size(t));
spiketrain1(fix(tSpike1/dt)) = 1;
spiketrain2(fix(tSpike2/dt)) = 1;

Kernel = exp(-tk/tau1)-exp(-tk/tau2);
Kernel = Kernel/max(Kernel);  % normalized
cSpike1 = conv(Kernel,spiketrain1);
cSpike2 = conv(Kernel,spiketrain2);

dRossum =  dt/tau1 * sum((cSpike1-cSpike2).^2) ;
end


