function [DistPtns,actFiring]=Testing_ReSuMe(ptnSet,Wts,targetTime)
% Testing the learning result through the ReSuMe method
% Testing ptnSet= (nAfferents,nPtns) need to convert to cell(1,p);  Wts=(nAfferents,1)
% DistPtns=(nPtn,1)

nWts = 1;

[nAfferents,nPtn] = size(ptnSet);
% actFiring = [];
actFiring = cell(1,nPtn);
count = 0;
% dt = 0.05;

%---- Neuron parameters ---%
ptnTime=25;      % Pattern time 25 ms
dt = 0.1;   % (ms)

tau_m=10;
Rm = 10e6;
Vr = -60;        % membrane potential at rest -60mV
I_ns = 0;           % noise
thr = -55;       % threshold -55mV
V_init = -60;    % t=0 initial membrane potential
V_res=-65;       % reset potential
t_ref = 1;       % refactory time 5ms     @@@@@@@@@@@@@@@@@@@@@@@

tau_d = 0.75;
%---- end Neuron parameters ---%


targetPtn = targetTime;
% load('TargetTime','targetTime');


DistPtns = zeros(nPtn, 1);
for iw = 1:nWts
    
    for iPtn = 1:nPtn
        
        spkPtn=[ptnSet(:,iPtn) (1:nAfferents)'];
        
        firings = [];      % neuron's actual firing times
        WtsFlag = zeros(1,nWts);
        %----- reset states --------
        Vm = V_init;
        %----- end reset states ----
        %    Dw = zeros(size(Wts)); % incremental learning
        
        for t= dt : dt : (ptnTime-dt)
            if WtsFlag(iw) == 1
                    break;
                end
            I_syn = 0;
            %--- caculate current states ---
            if isempty(firings)||t-firings(end)>t_ref
                notHeld = 1;   % held membrane potential in the refactory time.
            else
                notHeld = 0;
            end
            
            if notHeld == 1
                    Tsyn = find(spkPtn(:,1)<=t+0.1*dt);
                    if isempty(Tsyn)
                        continue;
                    else
                        I_syn = I_syn + Wts{iw}(Tsyn)'*(exp(-(t-spkPtn(Tsyn,1))./tau_d));       % nA
                    end                
                Vm = Vm + dt*(Vr-Vm+Rm*(I_syn+I_ns))/tau_m;
            else
                Vm = V_res;
            end
            
            if Vm>=thr && notHeld
                firings = [firings t];
                Vm = V_res;
                WtsFlag(iw) = 1;
            end
            
        end
        actFiring{iPtn} = firings;
        %         DistPtns(iPtn) = spkDist(targetPtn,firings,ptnTime,10);
        
        %         actFiring{iPtn} = firings;
    end % iPtn
end % iw

for iPtn = 1:nPtn
    DistPtns(iPtn) = spkDist(targetPtn,actFiring{iPtn},ptnTime,10);
end

end