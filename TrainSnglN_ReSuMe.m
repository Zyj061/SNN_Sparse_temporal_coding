function Wts = TrainSnglN_ReSuMe(InitWt,ptnSet,ptnLabel,NeuroCls,targetTime)
% ReSuMe training rule
% Training of Single neuron with the ptnSet for one Epoch; Return the
% weights of one epoch training
% the new precise learning rule (Based on Exp3) Trial learning
% InitWt -- 2d matrix of (nAfferent,1);
% ptnSet -- input pattern set (nAfferent,nPtns): need to convert to cell(1,p)
% ptnLabel -- (nPtns,1) values: 0,1,2...9;
% NeuroCls -- single scale for neuron class label: 0,1,2....9 meaning the
% neuron is supposed to learn that class.

% Author: Yajing Zheng
% Time: 2016

nPtn = length(ptnLabel);   % total number of patterns
TemtargetPtn = targetTime;
count = 0;
PtnNum = sum(ptnLabel == NeuroCls);

nAfferents=size(ptnSet,1);      % Number of inputs
eta = 0.1;             % learning rate  0.1
ptnTime=25;      % Pattern time 25 ms
dt = 0.1;   % (ms)

tau_m=10;
Rm = 10e6;
Vr = -60;        % membrane potential at rest -60mV
I_ns = 0;           % noise
thr = -55;       % threshold -55mV
V_init = -60;    % t=0 initial membrane potential
V_res = -65;       % reset potential
t_ref = 1;       % refactory time 5ms     @@@@@@@@@@@@@@@@@@@@@@@

tau_d = 0.75;
a_d = 0.5;
% A_di = 20;
A_di = 20;
tau_di = 1.25;
% time resolution is 1e-5;
%---- end Neuron parameters ---%

% AvgDistEpoch = -1;
nWts = 1;     % associate with 3 output neuron

for iw = 1:nWts
    Wts{iw} = InitWt{iw}(:,NeuroCls+1);   %0.5+0.2*randn(nAfferents,1);   % initial weights****************
end


%        Isyn=0;    % total input synaptic current for postneuron;  Isyn=sum(Wts*PSC)

%        PSC = K1-K2;              % post synaptic current of each afferent neuron
Ins = 0;            % background noise


%     for iepoch=1:maxEpoch

%PtnFlag =[];  pattern flag: 0 for negative pattern; 1 for positive
WtsFlag = zeros(1,nWts);
ActFire = cell(nPtn,1);

spkOutput = [];

for iW = 1:nWts
    
    for iPtn = 1:nPtn
       
        %     CPtnId = TrSeqOrd(iPtn);
        %     spkPtn=[ptnSet(:,CPtnId) (1:nAfferents)'];
        spkPtn=[ptnSet(:,iPtn) (1:nAfferents)'];
        
        %     if ptnLabel(CPtnId)==NeuroCls
        if ptnLabel(iPtn) == NeuroCls
            targetPtn = TemtargetPtn;
            %         targetPtn = cell2mat(targetTime(NeuroCls+1));
            PtnFlag = 1; % positive pattern (Flag)
        else
            targetPtn = ones(1,nWts) * (ptnTime+10); % this means teach neuron never fire
            PtnFlag = 0;
        end
        
        %----- reset states --------
        WtsFlag(iW) = 0;
        firings = [];      % neuron's actual firing times
        Vm = V_init;
        %----- end reset states ----
        
        for t=dt:dt:(ptnTime-dt)
            
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
                        I_syn = I_syn + Wts{iW}(Tsyn)'*(exp(-(t-spkPtn(Tsyn,1))./tau_d));       % nA
                    end                    
                
                Vm = Vm + dt*(Vr-Vm+Rm*(I_syn+I_ns))/tau_m;
            else
                Vm = V_res;
            end
            
            if ~isempty(find(abs(targetPtn-t)<0.1*dt, 1))               % target firing time t
                Sd = 1;             % fire
            else
                Sd = 0;
            end
            
            if Vm>=thr && notHeld
                firings = [firings t];
                So = 1;
                Vm = V_res;
            else
                So = 0;
            end
            
            if Sd || So
                for p = 1:nAfferents
                    Tsyn = spkPtn(p,1)<=t+0.1*dt;
                    if Tsyn==0 || Sd-So == 0
                        adi = 0;
                    else
                        adi = A_di* exp(-(t-spkPtn(p,1))/tau_di);
                    end
                    Wts{iW}(p)=Wts{iW}(p) + eta*(Sd-So)*(a_d+adi);
                end
            end            
        end
        
        
        if ~isempty(firings)  && PtnFlag == 1
            firings;
            
        else if PtnFlag == 1
                count = count + 1;
            end
        end
        
    end % iPtn
    WtsFlag(iW) = 1;
end % iW

rate = count / PtnNum * 100;
fprintf('not Firing Rate on Object %d is %.1f%%\n',NeuroCls+1,rate);

end


