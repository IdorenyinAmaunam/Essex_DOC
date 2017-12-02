function [RunResults] = analyzeArtifactPercentage(FilePath, RunResults)

try
    [data, header] = sload(FilePath);
    
    % Interpolate NaNs
    data = proc_naninterpolation(data);
    
    % Some times (rarely) there are still some NaNs remaining in some
    % channels Set those to 0
    data(isnan(data))=0;
    
catch
    disp(['Problem loading file, skipping: ' FilePath]);
    RunResults.fine = 0;
    return;
end


% Remove bullshit triggers
IndRem = find(header.EVENT.TYP > 30000);
header.EVENT.TYP(IndRem) = [];
header.EVENT.POS(IndRem) = [];
header.EVENT.DUR(IndRem) = [];

% Convert 14 to 897 since it seems to be the "hit" in the Sham protocol
header.EVENT.TYP(header.EVENT.TYP==14)=897;

if(sum(ismember(unique(header.EVENT.TYP),[1 769 770 783 781 786 897 898])) ~= length(unique(header.EVENT.TYP)))
    % More weird triggers detected...
    RunResults.fine = 0;
    return
end



% Trial extraction
pos = header.EVENT.POS;
dur = header.EVENT.DUR;
cf = find(header.EVENT.TYP==781);
fix = find(header.EVENT.TYP==786);

if (length(cf) < 15)
    % Too few trials
    RunResults.fine = 0;
    disp('Too few trials, skipping...');
    return;
end

% Find durations of trials
trdur = [];
for tr=1:length(cf)
    try
        if(  (header.EVENT.TYP(cf(tr)+1)==897) || (header.EVENT.TYP(cf(tr)+1)==898) )
            trdur(tr) = pos(cf(tr)+1)-pos(cf(tr));
        elseif(header.EVENT.TYP(cf(tr)+1)== 1)
            % It is impossible to know exactly the end of a timeoue trial since
            % the timeout was not logged (stupid), but it was normally 7
            % seconds
            trdur(tr) = 7*512;
        else
            disp(['Something is wrong with the triggers!']);
            RunResults.fine = 0;
            return;
        end
    catch
        % Apparently that was the last trial and was a timeout
        trdur(tr) = 7*512;
    end
end

% A trial begins with the next window after the one where the trailing (right)
% edge is adjacent to the 781 trigger
% Note that the way I use trials to extract ftrials later on, the two edges
% are not treated in the same way. The first (left) refers to the window's
% leading (left) edge, while the second (right) refers to the window's
% trailing (right) edge.

trials = [pos(cf)-512+32 pos(cf)+trdur'];


%% Filter all channelels below 20 Hz to remove any EMG noise
fs = 512;
cutoff = 20;
[b,a] = butter(4,cutoff/(fs/2));
for ch=1:16
    filtdata(:,ch) = filtfilt(b,a,data(:,ch));
    M(ch) = nanmean(filtdata(:,ch));
    S(ch) = nanstd(filtdata(:,ch));
end

% Remove any EMG channels (some subjects have those)
data = data(:,1:16);

zdata = (data-repmat(M,size(data,1),1))./repmat(S,size(data,1),1); 
zfiltdata = (filtdata-repmat(M,size(filtdata,1),1))./repmat(S,size(filtdata,1),1);

RunResults.artifact.rawtrial = zeros(1,size(data,1));
RunResults.artifact.notrial = ones(1,size(data,1));
for tr=1:size(trials,1)
    RunResults.artifact.rawtrial(trials(tr,1):trials(tr,2))=1;
    RunResults.artifact.notrial(trials(tr,1):trials(tr,2))=0;
end

RunResults.artifact.zth = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0];

for th=1:length(RunResults.artifact.zth)
    RunResults.artifact.clean(th) = 100*sum(sum((abs(zdata) > RunResults.artifact.zth(th))') & RunResults.artifact.rawtrial)/sum(RunResults.artifact.rawtrial);
    RunResults.artifact.noisy(th) = 100*sum(sum((abs(zdata) > RunResults.artifact.zth(th))') & RunResults.artifact.notrial)/sum(RunResults.artifact.notrial);
end
