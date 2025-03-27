function TestedLabels = ResultAnalysis(DistPtns)
% return tested Labels
% DistPtns -- (nPtns, nCls)

[nPtns, ~] = size(DistPtns);
TestedLabels = zeros(nPtns,1);

for iP = 1: nPtns
    CurDist = DistPtns(iP,:);
    [~,I] = min(CurDist);
    TestedLabels(iP)=mod(I,10) - 1;
    if TestedLabels(iP) == -1
        TestedLabels(iP) = 9;
    end
end

end