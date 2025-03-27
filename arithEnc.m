function spiketrain = arithEnc(firingTime,localFiringSpike,firingSpike,picScale,biasIdx)
% generating the spike train based the information gotten from the features
% extraction procedure throught the arithmetic method

% author: Yajing Zheng
% time : 2016.08


nModel = 50;
nSpkTrn = size(firingTime,1);  % num of pictures
nAfrnts = size(firingTime,2) + size(firingSpike,2);  % num of stored neuron of each Object
% nAfrnts = size(firingTime,2);

Twnd = 25;
deltaT = Twnd / nAfrnts;   % 对应每个storedNeuron的时间段
countIdx = 0;
spiketrain = [];
for spk = 1 : nSpkTrn
    
    
    if sum(firingTime(spk,:)) ==0 || sum(firingTime(spk,:)~=0) == 1
        continue;
    else
        countIdx = countIdx + 1;
    end
    
    %     for aff = 1 : nAfrnts
    for m = 1:nModel
        %         if firingTime(spk,aff) == 0
        %         if firingTime(spk,aff) == 0
        if firingTime(spk,m) == 0
%                         spiketrain(countIdx,2*m-1) = Twnd + 10;
                        spiketrain(countIdx,2*m-1) = 2 * (biasIdx(m) - 1) * deltaT;
%             spiketrain(countIdx,aff) = (aff - 1) * deltaT;
        else
            
            if firingTime(spk,m) / localFiringSpike(spk,m) >= deltaT
                %                             spiketrain(spk,2*m-1) = (m -1) * deltaT + (deltaT - 0.01);
                spiketrain(countIdx,2*m-1) = 2*(biasIdx(m) -1) * deltaT + (deltaT - 0.001);
%                   spiketrain(countIdx,2*m-1) = 2*(m - 1) * deltaT + deltaT / 2;
                %                               spiketrain(spk,2*m-1) = 2*(m-1)*deltaT;
            else
                spiketrain(countIdx,2*m-1) = 2 * (biasIdx(m) - 1) * deltaT + firingTime(spk,m) / localFiringSpike(spk,m);
            end
            %             if firingTime(spk,aff) >= deltaT
            %                 firingTime(spk,aff) = deltaT - 0.001;
            %             end
            %             spiketrain(countIdx,aff) = (aff - 1)*deltaT + firingTime(spk,aff);
            %                  spiketrain(spk,2*m-1) =2*(m-1)*deltaT;
        end
        
        %----------------------------firingSpike--------------------------%
        if firingSpike(spk,m) == Inf
%                             spiketrain(countIdx,2*m) = Twnd + 10;
            %               spiketrain(spk,2*m) = ReferTOff(2*m);
            spiketrain(countIdx,2*m) = (2*biasIdx(m) - 1) * deltaT;
        else
            
            if firingSpike(spk,m) / picScale(spk) <= (1.0 - deltaT)
                firingSpike(countIdx,m) = picScale(spk) * (1.01 - deltaT);
            end
            %               spiketrain(spk,2*m) =(2*m-1) * deltaT;
            spiketrain(countIdx,2*m) = (2*biasIdx(m) - 1) * deltaT + (1.0 - firingSpike(spk,m) / picScale(spk));
            %                 spiketrain(spk,2*m) = Twnd - ((m-1) * deltaT + (1.0 - firingSpike(spk,m) / picScale(spk)));
        end
    end
    
end
end
% spiketrain(firingSpike == Inf) = Twnd + 10;

