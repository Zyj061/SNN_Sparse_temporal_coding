function Wts = TrainSnglN_SPAN(InitWt,ptnSet,ptnLabel,NeuroCls,targetTime)
% SPAN training rule
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

nCls = 10; % corresponding to 10 classes
nPtn = length(ptnLabel);   % total number of patterns
TemtargetPtn = targetTime;

nAfferents=size(ptnSet,1);      % Number of inputs
ptnTime = 25;      % Pattern time (ms)
dt = 0.1;          % step time - resolution (ms)
count = 0;
PtnNum = sum(ptnLabel == NeuroCls);
%---- Neuron parameters ---%
tau_m = 10;   % (ms)
Rm = 1;       % (Mohm)
Vr = 0;       % rest potential (also as the reset potential) (mV)
Vthr = 10;     % threshold (mV)
% Vthr = 25;
t_ref =   0.75;    % refectory time (ms)
% t_ref = 3;
wmax = 10;   % synaptic weight range (nA)
wmin= -wmax;
tau_s = 0.75;   % (ms)
V0 = 1;
Sc1 = exp(-dt/tau_s);
lamda = nCls / nPtn;  % learning rate  ************************************$$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!
%---- end Neuron parameters ---%

% kTgtSpks = 4; % number of spikes in the target train
% targetPtn = round((1:kTgtSpks)/(kTgtSpks+1)*ptnTime);         %[30 60 90];     % target firing (ms) -- teacher


% maxEpoch=200;
ErrorCrit =0.06;  % converge criteria (spike distance)
% ErrorCrit = 1;

% AvgDistEpoch = -1;
nWts = 1;     % associate with 3 output neuron

for iw = 1:nWts
    Wts{iw} = InitWt{iw}(:,NeuroCls+1);   %0.5+0.2*randn(nAfferents,1);   % initial weights****************
end


%        Isyn=0;    % total input synaptic current for postneuron;  Isyn=sum(Wts*PSC)

%        PSC = K1-K2;              % post synaptic current of each afferent neuron
Ins = 0;            % background noise


%     for iepoch=1:maxEpoch

% TrSeqOrd = randperm(nPtn);
DistPtns = zeros(1, nPtn);   % used for stopping criteria
beginVm=zeros(1, round(ptnTime/dt)); % monitor of the membrane potentials

%PtnFlag =[];  pattern flag: 0 for negative pattern; 1 for positive
WtsFlag = zeros(1,nWts);
ActFire = cell(nPtn,1);

spkOutput = [];

for iW = 1:nWts
    
    K1 = zeros(size(Wts{iW}));
    
    for iPtn = 1:nPtn
        Dw = zeros(size(Wts{iW}));   % trial learning
        
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
        
        WtsFlag(iW) = 0;
        firings = [];      % neuron's actual firing times
        
        %----- reset states --------
        Vm = 0.7*Vthr;
        K1 = 0 * K1;
        PSC = K1;
        %----- end reset states ----
        %    Dw = zeros(size(Wts)); % incremental learning
        
        for t=dt:dt:(ptnTime-dt)
            %                     Dw = zeros(size(Wts)); % instanteous learning
            Isyn = Wts{iW}'*PSC;
            Vm = Vm + (dt/tau_m)*(Rm*(Ins+Isyn) - Vm);
            %--- caculate current states ---
           
            SpkTimesIdx = abs(spkPtn(:,1)-t) < 0.1 * dt;
            
            SpkAfrnt = spkPtn(SpkTimesIdx,2);  % current firing afferents
            if ~isempty(SpkAfrnt)
                K1(SpkAfrnt) = K1(SpkAfrnt) + V0;
            end
            PSC = K1;   % also used for updating weights            
            
            
            beginVm(round(t/dt))=Vm;
            noZeroIdx = find(K1);
            t_i = spkPtn(noZeroIdx,1);
            %--- end of caculate current states ---
            
            noRefFlag = isempty(firings)||t-firings(end)>t_ref; % flag of not in refractory period
            if ~noRefFlag % in refractory period
                Vm = Vr;        % keep
            end
            
            if Vm>Vthr     % actual firing & instanteous update LTD
                beginVm(round(t/dt))=1.2*Vm;
                Vm = Vr;
                firings = [firings t];
                % LTD
                %tmpDw = zeros(size(Wts{iW}));
                Dw1 = (abs(t_i-t)+tau_s) .* exp(-abs(t_i-t)/tau_s);
                tmpDw = lamda*1e2/4*Dw1;
                Dw(noZeroIdx) = Dw(noZeroIdx) - tmpDw;
                
                %                 if WtsFlag(iW) == 1
                %                     continue;
                %                 end
                
                
                len = length(ActFire{iPtn});
                if len == 0
                    ActFire{iPtn} = t;
                    %                     WtsFlag(iW) = 1;
                else if t > ActFire{iPtn}(len)
                        ActFire{iPtn}(len+1) = t;
                        %                         WtsFlag(iW) = 1;
                    end
                end
                
            end
            
            %         if ~isempty(find(abs(targetPtn-t)< 0.6, 1))   % desired spike time & LTP
            if ~isempty(find(abs(targetPtn(iW)-t) < 0.1 * dt, 1))   % desired spike time & LTP
                t_d = targetPtn(iW);
                Dw2 = (abs(t_i-t_d)+tau_s) .* exp(-abs(t_i-t_d)/tau_s);
                tmpDw = lamda * 1e2 / 4 * Dw2;
                Dw(noZeroIdx) = Dw(noZeroIdx) + tmpDw;
%                 WtsFlag(iW) = 1;  %当在目标处有发放才将其flag位立起
            end
            
            addTerm = zeros(size(K1));
            addTerm(noZeroIdx) = exp(1-(abs(t_i-t)+dt)/tau_s)*tau_s^(-1)*dt;
            K1 = K1* Sc1 + addTerm;  % used for next time step
            
        end
        
        if PtnFlag==1 % positive pattern
            DistPtns(iPtn) = spkDist(targetPtn,ActFire{iPtn},ptnTime,10);
            if DistPtns(iPtn)<=ErrorCrit  % if correct: don't change Wts
                disp(['Training of Object ' int2str(NeuroCls + 1) ' converged']);
                Dw = zeros(size(Wts{iW}));
            end
        end
        
        if (PtnFlag==0) && isempty(firings)
            %         if (PtnFlag==0) && isempty(ActFire{iPtn})
            Dw = zeros(size(Wts{iW}));
            %         else
            %             WtsFlag(iW) = 1; % signing output neuron iW have been fire
        end
        
        Wts{iW} = Wts{iW} + Dw;
        %     Wts = max(-wmax, min(wmax,Wts) );
        Wts{iW} = max(wmin, min(wmax,Wts{iW}) );
        
        if ~isempty(firings)  && PtnFlag == 1
            firings;           
            
        else if PtnFlag == 1
                %             disp(['Object ' int2str(NeuroCls+1) ' not Fire']);
%                 beginVmFig = figure('position',[10 100 800 100]);
%             figure(beginVmFig);
%             plot([0 dt:dt:ptnTime], [0.7*Vthr beginVm],'LineWidth',1.5);
%             hold on; plot([0 ptnTime],[Vthr Vthr],'--r', 'LineWidth',1);
%             xlim([0 ptnTime]); ylim([0 1.2*Vm]);
%             ylabel('Vm'); xlabel('Time (ms)');
                count = count + 1;
            end
        end
        
    end % iPtn
    WtsFlag(iW) = 1;
end % iW

rate = count / PtnNum * 100;
fprintf('not Firing Rate on Object %d is %.1f%%\n',NeuroCls+1,rate);

%             fprintf('Epoch:%d avgE: %6.3f\n',...
%            iepoch, mean(DistPtns));
%         AvgDistEpoch=mean(DistPtns);

%    Wts = Wts + Dw;  % incremental learning
%             Dist(iepoch) = spkDist(targetPtn,firings,ptnTime,tau1);


%     end % iepoch

%         % Testing
%         disp('Testing...');
%         for iJ=1:length(testJit)
%             TestDist=Testing(ptnSet,Wts,testJit(iJ));
%             CorrectRate(iRun,iJ)=sum(TestDist<=Crt)/length(TestDist);
%         end
%
%         save('Result','AvgDistRunEpoch','CorrectRate');

end


