function Indices = crossvalindTrials(alabels,atrials,K)

% Remap trials to [1:N]
UniqueTrials = unique(atrials);
ctrials = zeros(length(atrials),1);
for ntr=1:length(UniqueTrials)
    ctrials(find(atrials==UniqueTrials(ntr)))=ntr;
end
atrials = ctrials;

TrLbl=[];
for tr=1:max(atrials)
    tmptr = alabels(find(atrials==tr));
    TrLbl(tr) = tmptr(1);
end
cvind = crossvalind('Kfold',TrLbl,K);

Indices = alabels;
for tr=1:max(atrials)
    Indices(find(atrials==tr)) = cvind(tr);
end