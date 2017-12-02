function AvgAcc = cvk(afeats, alabels, atrials, K, NFeat)

% Accuracy with cross validation (balanced groups, random shuffling)
cvind = crossvalindTrials(alabels,atrials,K);

for i=K:-1:1
    tstfoldind = find(cvind==i);
    trfoldind = setdiff([1:size(afeats,1)],find(cvind==i))';
   
    % Do feature selection on training set only
    for fr=1:size(afeats,2)
        for ch=1:size(afeats,3)
            cfeat = squeeze(afeats(trfoldind,fr,ch));
            clabels = alabels(trfoldind);
            FSmatOutliers(fr,ch) = fisherScore(cfeat,clabels); 
        end
    end

    FSmatOutliers(isnan(FSmatOutliers))=0;
    
    [SortedFS,SortInd] = sort(FSmatOutliers(:),'descend');
    [i1,i2] = ind2sub(size(FSmatOutliers),SortInd);
    i1 = i1(1:NFeat);
    i2 = i2(1:NFeat);
    
    troutInd = tstfoldind;
    trdata = afeats;
    trdata(troutInd,:,:) = [];
    strdata = [];
    strlabels = alabels;
    strlabels(troutInd) = [];
    for s=1:size(trdata,1)
        for j=1:NFeat
            strdata(s,j) = trdata(s,i1(j),i2(j));
        end
    end
    
    tstoutInd = trfoldind;
    ststdata = [];
    tstdata = afeats;
    tstdata(tstoutInd,:,:)=[];
    ststlabels = alabels;
    ststlabels(tstoutInd)=[];
    for s=1:size(tstdata,1)
        for j=1:NFeat
            ststdata(s,j) = tstdata(s,i1(j),i2(j));
        end
    end
    
    [~,~,~, err_tr(i)] = myclassify(strdata,strdata,strlabels,strlabels);
    
    AccCVtr(i) = 100*(1-err_tr(i));
    [~,~,AccCVtst(i)] = myclassify(ststdata,strdata,strlabels,ststlabels);
end
AvgAcc = mean(AccCVtst);
