%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Zheng Yajing
%
% Object recognition learning main function 
% Files:
% 1. arithEnc.m
% 2. TrainSnglN.m
% 3. Testing.m
% 4. ResultAnalysis.m
%
% Tasks:
% 1. Encoding the features, which were extracted from the STDP-based HMAX unsupervised
%    learning rule
% 2. 3d object recognition through the PSD rule
% 
% PARAM.learningRule: 1 (PSD rule)
%                     2 (SPAN rule)
%                     3 (ReSuMe rule)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nCls = 10;
nNeurons = nCls;
maxEpoch = 20; % for each run, perform 'maxEpoch' epochs

nAfferents = PARAM.nNeuron * 2;

% Relative Confidence
TrnAccuracyArray_RC = zeros(1,10);
TstAccuracyArray_RC = zeros(1,10);

% Absolute Confidence
TrnAccuracyArray_AC = zeros(1,10);
TstAccuracyArray_AC = zeros(1,10);

TrnRc = [];
TstRc = [];
TrnAc = [];
TstAc = [];
% Relative Confidience and Absolute Confidence for each category
% respectively
for i = 1 : nCls
    TrnRc{i} = zeros(1,10);
    TstRc{i} = zeros(1,10);
    TrnAc{i} = zeros(1,10);
    TstAc{i} = zeros(1,10);
end

%10-fold cross validation
nRun = 10;
k_fold = 10;

% Each corresponding to each category
PtnEachLen_10 = floor(PtnEachLen ./ k_fold);
PtnEachClsIdx = zeros(1,k_fold);
PtnEachClsIdx(1) = 1;
PtnSubsetLen = sum(PtnEachLen_10);
for ifold = 2:k_fold
    PtnEachClsIdx(ifold) = sum(PtnEachLen(1:ifold-1))+1;
end

PtnSubSet = [];   % initial subset which is used to 10-fold cross validation
PtnSubLabels = [];  % labels corresponding to each subset

for iSubSet = 1:10
    idx = 1;
    PtnSubSet{iSubSet} = [];
    PtnSubLabels{iSubSet} = zeros(PtnSubsetLen,1);
    for ifold = 1:k_fold
        beginIdx = PtnEachClsIdx(ifold) + (iSubSet-1) * PtnEachLen_10(ifold);
        endIdx = PtnEachClsIdx(ifold) + iSubSet * PtnEachLen_10(ifold) - 1;
        PtnSubSet{iSubSet} = [PtnSubSet{iSubSet} PtnSet(:,beginIdx:endIdx)];
        PtnSubLabels{iSubSet}(idx:idx+PtnEachLen_10(ifold)-1) = ifold - 1;
        idx = idx + PtnEachLen_10(ifold);
    end
end


% begin Training
for iRun = 1:nRun
    InitWts = [];
    for i = 1:PARAM.nWts
        InitWts{i}= 0.5+0.2*randn(nAfferents,nNeurons); % random initialization w
        %       InitWts{i} = unifrnd(0,25) * ones(nAfferents,nNeurons);
    end
    Wts = InitWts;
    
    % PtnTrSet = [PtnTrSet phaseEnc(cell2mat(COMMON.tst_neg_firingTime),cell2mat(COMMON.tst_neg_picScale))'];
    PtnTstSet = PtnSubSet{iRun};
    tst_labels = PtnSubLabels{iRun};
    PtnTrSet = [];
    trn_labels = [];
    
    
    for ifold = 1:10
        if ifold == iRun
            continue;
        else
            PtnTrSet = [PtnTrSet PtnSubSet{ifold}];
            trn_labels = [trn_labels; PtnSubLabels{ifold}];
            
        end
    end
    
    nSample = length(PtnTrSet);
    mini_batch = fix(nSample / nCls);
    % mini_batch = 50;
    Twnd = 25;
    targetTime = COMMON.targetTime;
    tic
    for iepoch = 1:maxEpoch
        
        currentOrder = randperm(nSample);
        currentPtnSet = PtnTrSet(:,currentOrder);
        currentTrnLabels = trn_labels(currentOrder,1);
        idx = 1;
        
        
        for iNeuron = 1:nNeurons             % Training Neurons
            if iNeuron == 1
                idx = 1;
            else
                idx = idx +  mini_batch;
                %                       idx = idx + nIns;
            end
            fprintf('epoch:%d Neuron:%d Training...\n',iepoch,iNeuron);
            %        Wts(:,iNeuron) = TrainSnglN(Wts(:,iNeuron),PtnTrSet,trn_labels,iNeuron-1,targetTime{iNeuron});
            %           Wts(:,iNeuron) = TrainSnglN(Wts(:,iNeuron), PtnTrSet(:,idx: idx + nIns - 1), trn_labels(idx:idx + nIns - 1,1), iNeuron-1,targetTime{iNeuron});
            tmpWts = [];
            switch PARAM.learningRule
                case 1
                    % PSD training
                    tmpWts = TrainSnglN(Wts, currentPtnSet(:,idx:idx + mini_batch -1), currentTrnLabels(idx:idx+mini_batch - 1,1), iNeuron-1,targetTime{iNeuron});
                case 2
                    % SPAN training
                    tmpWts = TrainSnglN(Wts, currentPtnSet(:,idx:idx + mini_batch -1), currentTrnLabels(idx:idx+mini_batch - 1,1), iNeuron-1,targetTime{iNeuron});
            for iw = 1:PARAM.nWts
                Wts{iw}(:,iNeuron) = tmpWts{iw};
            end
        end
    end %end iepoch
    % display computational time
    computationTime = toc / 60; % in minute
    if computationTime > 120
        disp(['Computation time = ' num2str(computationTime/60) ' hours']);
    else
        disp(['Computation time = ' num2str(computationTime) ' min']);
    end
    
    
    % Testing
    disp('Testing...');
    
    DistTrnPtns = zeros(length(trn_labels),nCls);
    DistTePtns = zeros(length(tst_labels),nCls);
    
    TeFiring = [];
    TrnFiring = [];
    
    for iNeuron=1:nNeurons % Train pattern set
        tmpWts = [];
        for iw = 1:PARAM.nWts
            tmpWts{iw} = Wts{iw}(:,iNeuron);
        end
        
        % PSD Testing
        [DistTePtns(:,iNeuron),TeFiring{iNeuron}] = Testing(PtnTstSet,tmpWts,targetTime{iNeuron});
        [DistTrnPtns(:,iNeuron),TrnFiring{iNeuron}] = Testing(PtnTrSet,tmpWts,targetTime{iNeuron});
        
    end
    
    TestedLabels = ResultAnalysis(DistTePtns);
    TrainLabels = ResultAnalysis(DistTrnPtns);
    
    TrnlabelsHead = [];
    TrnlabelsMonitor = [];
    TrnlabelsCar = [];
    DistTrnPtns_head = [];
    DistTrnPtns_monitor = [];
    DistTrnPtns_car = [];
    
    for i = 1:nCls
        DistTrnPtns_10{i} = [];
        Trnlabels_10{i} = [];
        DistTstPtns_10{i} = [];
        TestedLabels_10{i} = [];
    end   
    
    for i = 1:nCls
        if i == 1
            beginIdx = 1;
        else
            beginIdx = sum(PtnEachLen_10(1:i-1)) + 1;
        end
        endIdx = beginIdx + PtnEachLen_10(i) - 1;
        TestedLabels_10{i} = TestedLabels(beginIdx:endIdx);
        DistTstPtns_10{i} = DistTePtns(beginIdx:endIdx, i);
        tst_labels10{i} = ones(PtnEachLen_10(i),1) * (i-1);
    end
        
    count = 0;
    for ifold = 1:10
        if ifold == iRun
            continue;
        else
            
            bias = PtnSubsetLen * count;
            for i = 1:nCls
                if i == 1
                    beginIdx = 1;
                else
                    beginIdx = sum(PtnEachLen_10(1:i-1)) + 1 + bias;
                end
                endIdx = beginIdx + PtnEachLen_10(i) - 1;
                Trnlabels_10{i} =  [Trnlabels_10{i}; TrainLabels(beginIdx:endIdx)];
                DistTrnPtns_10{i} = [DistTrnPtns_10{i}; DistTrnPtns(beginIdx:endIdx,i)];
            end
            count = count + 1;
        end
    end
    
    TestAccrcy= sum(TestedLabels==tst_labels)/length(tst_labels) * 100;
    TrainAccrcy = sum(TrainLabels == trn_labels) / length(trn_labels) * 100;
    
    TestAccrcy_10 = zeros(1,nCls);
    TrainAccrcy_10 = zeros(1,nCls);
    
    for i = 1:nCls
        TestAccrcy_10(i) = sum(TestedLabels_10{i} == i-1) / length(TestedLabels_10{i}) * 100;
        TrainAccrcy_10(i) = sum(Trnlabels_10{i} == i-1) / length(Trnlabels_10{i}) * 100;
        TrnRc{i}(iRun) = TrainAccrcy_10(i);
        TstRc{i}(iRun) = TestAccrcy_10(i);
    end
   
    TrnAccry_AC = 0;
    TstAccry_AC = 0;
    
    idx = 1;
    for iSubset = 1 : nRun - 1
        for iCls = 1:nCls
            len = PtnEachLen_10(iCls);
            TrnAccry_AC = TrnAccry_AC + sum(DistTrnPtns(idx:idx+len-1,iCls)<0.05);
            idx = idx+len;
        end
    end
    idx = 1;
    for iCls = 1:nCls
        len = PtnEachLen_10(iCls);
        TstAccry_AC = TstAccry_AC + sum(DistTePtns(idx:idx+len-1,iCls)<0.05);
        idx = idx + len;
    end
    % TestAccrcy= sum(DistTePtns <= 0.06)/length(tst_labels);
    
    TrnAccuracyArray_RC(iRun) = TrainAccrcy;
    TstAccuracyArray_RC(iRun) = TestAccrcy;    
    
    
    TrnAccuracyArray_AC(iRun) = TrnAccry_AC / size(PtnTrSet,2) * 100;
    TstAccuracyArray_AC(iRun) = TstAccry_AC / size(PtnTstSet,2) * 100;
    
    for i = 1:nCls
        TrnAc{i}(iRun) = sum(DistTrnPtns_10{i} < 0.05) / length(Trnlabels_10{i}) * 100;
        TstAc{i}(iRun) = sum(DistTstPtns_10{i} < 0.05) / length(TestedLabels_10{i}) * 100;
    end
   
    disp(['iRun: ' int2str(iRun)]);
    fprintf('Accuracy on Testing Dataset: %.2f%% \n',TestAccrcy);
    fprintf('Accuracy on Training Dataset: %.2f%% \n',TrainAccrcy);
    
    disp('end');
    
end

disp('10-fold cross validation done...');
fprintf('Average Relative Accuracy on Training Dataset: %.2f%% \n', mean(TrnAccuracyArray_RC));
fprintf('Average Relative Accuracy on Testing Dataset: %.2f%% \n', mean(TstAccuracyArray_RC));

fprintf('Average Absolute Accuracy on Training Dataset: %.2f%% \n', mean(TrnAccuracyArray_AC));
fprintf('Average Ansolute Accuracy on Testing Dataset: %.2f%% \n', mean(TstAccuracyArray_AC));

for i = 1:nCls
    disp(COMMON.subdir(i+2).name);
    fprintf('Average Relative Accuracy on Training Dataset: %.2f%% \n', mean(TrnRc{i}));
    fprintf('Average Relative Accuracy on Testing Dataset: %.2f%% \n', mean(TstRc{i}));
    
    fprintf('Average Absolute Accuracy on Training Dataset: %.2f%% \n', mean(TrnAc{i}));
    fprintf('Average Ansolute Accuracy on Testing Dataset: %.2f%% \n\n', mean(TstAc{i}));
    
end
