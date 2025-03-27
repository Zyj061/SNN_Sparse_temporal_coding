function [DistPtns,actFiring]=Testing(ptnSet,Wts,targetTime)
% Testing the learning result throught the PSD rule
% Testing ptnSet= (nAfferents,nPtns) need to convert to cell(1,p);  Wts=(nAfferents,1)
% DistPtns=(nPtn,1)

nWts = 1;
ptnTime = 25;      % Pattern time (ms)
dt = 0.1;          % step time - resolution (ms)

[nAfferents,nPtn] = size(ptnSet);
% actFiring = [];
actFiring = cell(1,nPtn);


% dt = 0.05;

%---- Neuron parameters ---%
tau_m = 10;   % (ms)
Rm = 1;       % (Mohm)
Vr = 0;       % rest potential (also as the reset potential) (mV)
Vthr = 10;     % threshold (mV)
%  Vthr = 25;
t_ref = 0.75;    % refectory time (ms)
%  t_ref = 0.5;
tau1=tau_m;   % PSC      ************************************$$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!
tau2=tau1/4;
Sc1 = exp(-dt/tau1);  % scale factor of step state
Sc2 = exp(-dt/tau2);
V0 = 1/max(exp(-(0:dt:5*tau1)/tau1)-exp(-(0:dt:5*tau1)/tau2)); % normalization factor
%---- end Neuron parameters ---%

% nCls = 4;  % number of categories
% kTgtSpks = nCls; % number of spikes in the target train
% targetPtn = [40 80 120 160];        %[30 60 90];     % target firing (ms) -- teacher
targetPtn = targetTime;
% load('TargetTime','targetTime');

Ins = 0;            % background noise


DistPtns = zeros(nPtn, 1);
for iw = 1:nWts
    K1 = zeros(size(Wts{iw}));
    K2 = zeros(size(Wts{iw}));
    
    for iPtn = 1:nPtn
        
        spkPtn=[ptnSet(:,iPtn) (1:nAfferents)'];
        
        firings = [];      % neuron's actual firing times
        WtsFlag = zeros(1,nWts);
        
        %----- reset states --------
        Vm = 0.7*Vthr;
        K1 = 0* K1;
        K2 = 0* K2;
        PSC = K1-K2;
        %----- end reset states ----
        %    Dw = zeros(size(Wts)); % incremental learning
        
        for t=dt:dt:(ptnTime-dt)
            
            %--- caculate current states ---
            Isyn = Wts{iw}'*PSC;
            Vm = Vm + (dt/tau_m)*(Vr-Vm + Rm*(Ins+Isyn));
            SpkTimesIdx = abs(spkPtn(:,1)-t) < 0.1 * dt;
            SpkAfrnt = spkPtn(SpkTimesIdx,2);  % current firing afferents
            if ~isempty(SpkAfrnt)
                K1(SpkAfrnt) = K1(SpkAfrnt) + V0;
                K2(SpkAfrnt) = K2(SpkAfrnt) + V0;
            end
            PSC = K1-K2;   % also used for updating weights
            %--- end of caculate current states ---
            
            noRefFlag = isempty(firings)||t-firings(end)>t_ref; % flag of not in refractory period
            if ~noRefFlag % in refractory period
                Vm = Vr;        % keep
            end
            
            if Vm>Vthr     % actual firing & instanteous update LTD
                Vm = Vr;
                if WtsFlag(iw) == 1
                    break;
                end
                
                len = length(actFiring{iPtn});
                if len == 0
                    actFiring{iPtn} = t;
                    WtsFlag(iw) = 1;
                else if t > actFiring{iPtn}(len)
                        actFiring{iPtn}(len + 1) = t;
                        WtsFlag(iw) = 1;
                    end
                end
                
                firings = [firings t];
            end
            
            K1 = Sc1*K1;  % used for next time step
            K2 = Sc2*K2;
        end
        
        %         DistPtns(iPtn) = spkDist(targetPtn,firings,ptnTime,10);
        
        %         actFiring{iPtn} = firings;
    end % iPtn
end % iw

for iPtn = 1:nPtn
    DistPtns(iPtn) = spkDist(targetPtn,actFiring{iPtn},ptnTime,10);
end

end