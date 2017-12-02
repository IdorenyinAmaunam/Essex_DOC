function [ConfMatClass ConfMatAll Accuracy Error] = myclassify(tstdata,trdata,trlbl,tstlbl)

try
    tmpclass = classify(tstdata,trdata,trlbl);
catch
    try
        tmpclass = classify(tstdata,trdata,trlbl,'diaglinear');
    catch
        disp(['Cannot even classify with Naive Bayes!!!']);
        ConfMatClass = NaN;
        ConfMatAll = NaN;
        Accuracy = NaN;
        Error = NaN;
        return;
    end
end
[ConfMatClass ConfMatAll Accuracy Error] = confusion_matrix(tstlbl, tmpclass, length(unique(trlbl)));