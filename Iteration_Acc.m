%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Zheng Yajing
% time: 2016.09
% 
% based on the Absolute Cofidence and Relative Cofidence criteria,
% evaluating the performance of learning rule of (1,3,5,10,15,20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iteration = [1,3,5,10,15,20];
len = length(iteration);
times = 10;
tmpTrnAccry_rc = zeros(1,times);
tmpTstAccry_rc = zeros(1,times);
tmpTrnAccry_ac =zeros(1,times);
tmpTstAccry_ac = zeros(1,times);

% Relative Confidence
meanTrnAccry_rc = zeros(1,6);
stdTrnAccry_rc = zeros(1,6);
meanTstAccry_rc = zeros(1,6);
stdTstAccry_rc = zeros(1,6);

% Absolute Confidence
meanTrnAccry_ac = zeros(1,6);
stdTrnAccry_ac = zeros(1,6);
meanTstAccry_ac = zeros(1,6);
stdTstAccry_ac = zeros(1,6);

for index = 1:6
    for iTime = 1:times
        wholeAP = 0;
        maxEpoch = iteration(index);
        OCRmain;
%         idx = 1;
%         for i = 1:10
%             len = PtnTrnEachLen(i);
%             wholeAP = wholeAP + sum(DistTrnPtns(idx:idx+len-1,i)<0.05);
%             idx = len + idx;
%         end
%         wholeTstAP = 0;idx=1;
%         for i = 1:10
%             len = PtnTstEachLen(i);
%             wholeTstAP = wholeTstAP + sum(DistTePtns(idx:idx+len-1,i)<0.05);
%             idx = len + idx;
%         end
        tmpTrnAccry_rc(iTime) = mean(TrnAccuracyArray_RC);
        tmpTstAccry_rc(iTime) =  mean(TstAccuracyArray_RC);
        tmpTrnAccry_ac(iTime) = mean(TrnAccuracyArray_AC);
        tmpTstAccry_ac(iTime) = mean(TstAccuracyArray_AC);
    end
  
    % relative confidience
    meanTrnAccry_rc(index) = mean(tmpTrnAccry_rc);
    stdTrnAccry_rc(index) = std(tmpTrnAccry_rc);
    meanTstAccry_rc(index) = mean(tmpTstAccry_rc);
    stdTstAccry_rc(index) = std(tmpTstAccry_rc);  
    
    % absolute confidience
    meanTrnAccry_ac(index) = mean(tmpTrnAccry_ac);
    stdTrnAccry_ac(index) = std(tmpTrnAccry_ac);
    meanTstAccry_ac(index) = mean(tmpTstAccry_ac);
    stdTstAccry_ac(index) = std(tmpTstAccry_ac); 
end
save(['Accry_iteration' 'RC_AC'],'meanTrnAccry_rc','stdTrnAccry_rc','meanTstAccry_rc','stdTstAccry_rc','meanTrnAccry_ac','stdTrnAccry_ac','meanTstAccry_ac','stdTstAccry_ac');

%------------------------draw classification accuracy-------------
% x = iteration;
% figure
% errorbar(x,meanAccry,stdTrn,'r');
% hold on
% plot1=plot(x,meanAccry,'-k^',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',5);
% hold on;
%
% errorbar(x,meanTstAccry,stdTstAccry,'r');
% hold on
% plot2=plot(x,meanTstAccry,'-b^',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',5);
% xlabel('Epoch');ylabel('Classification Accuracy(%)');ylim([90 100]);
% legend([plot1,plot2],'Train','Test')