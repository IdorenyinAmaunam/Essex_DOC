function [RunResults] = analyzeOnlineStroke(FilePath, MATpath, lap, freqs)

RunResults.fine = 1;
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


% Remove overall DC
data = data-repmat(mean(data),size(data,1),1);

% Laplacian spatial filtering
data = data(:,1:16);
data = laplacianSP(data,lap);

% Trial extraction
pos = header.EVENT.POS;
dur = header.EVENT.DUR;
cf = find(header.EVENT.TYP==781);
cue = cf-1;
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

trials = [pos(cf)-512+32 pos(cf)+trdur' pos(fix) pos(fix)+dur(fix)];

% Useful params for PSD extraction with the fast algorithm
psdshift = 512*0.5*0.5;
winshift = 512*0.0625;

if((mod(psdshift,winshift) ~=0) && (mod(winshift,psdshift) ~=0))
    disp(['[eegc3_smr_simloop_fast] The fast PSD method cannot be applied with the current settings!']);
    disp(['[eegc3_smr_simloop_fast] The internal welch window shift must be a multiple of the overall feature window shift (or vice versa)!']);
    return;
end

% Create arguments for spectrogram
spec_win = 512*0.5;
% Careful here: The overlapping depends on whether the winshift or the
% psdshift is smaller. Some calculated internal windows will be redundant,
% but the speed is much faster anyway

if(psdshift <= winshift)
    spec_ovl = spec_win - psdshift;
else
    spec_ovl = spec_win - winshift;
end

% Calculate all the internal PSD windows
for ch=1:16
    %disp(['[eegc3_smr_simloop_fast] Internal PSDs on electrode ' num2str(ch)]);
    [~,f,t,p(:,:,ch)] = spectrogram(data(:,ch), spec_win, spec_ovl, [], 512);
end

% Keep only desired frequencies
p = p(find(ismember(f,freqs)),:,:);

% Setup moving average filter parameters
FiltA = 1;
if(winshift >= psdshift)
    % Case where internal windows are shifted according to psdshift
    MAsize = (512*1.00)/psdshift - 1;   
    FiltB = (1/MAsize)*ones(1,MAsize);
    MAstep = winshift/psdshift;
else
    % Case where internal windows are shifted according to winshift
    FiltStep = psdshift/winshift;
    MAsize = (512*1.00)/winshift - FiltStep;   
    FiltB = zeros(1,MAsize);
    FiltB(1:FiltStep:end-1) = 1;
    FiltB = FiltB/sum(FiltB);
    MAstep = 1;
end

StartInd = find(FiltB~=0);
StartInd = StartInd(end);

afeats = filter(FiltB,FiltA,p,[],2);
afeats = permute(afeats, [2 1 3]);

% Get rid of initial filter byproducts
afeats = afeats(StartInd:end,:,:);

% In case of psdshift, there will be redundant windows. Remove them
if(MAstep > 1)
   afeats = afeats(1:MAstep:end,:,:);
end

% Take the log as final feature values
bfeats = afeats;
afeats = log(afeats);

% Remap trials to PSD space
ftrials(:,1) = floor(trials(:,1)/32)+1;
ftrials(:,2) = floor(trials(:,2)/32)-15;
ftrials(:,3) = floor(trials(:,3)/32)+1;
ftrials(:,4) = floor(trials(:,4)/32)-15;



% Crop to trials
fdata = [];
flbl = [];
ftrlbl = [];
albl = -ones(size(afeats,1),1);
atrlbl = -ones(size(afeats,1),1);
for i=1:size(ftrials,1)
    albl(ftrials(i,1):ftrials(i,2)) = 1;
    atrlbl(ftrials(i,1):ftrials(i,2)) = i;
    albl(ftrials(i,3):ftrials(i,4)) = 0;
    atrlbl(ftrials(i,3):ftrials(i,4)) = i;
    
end
RunResults.afeats = afeats;
RunResults.bfeats = bfeats;
RunResults.albl = albl;
RunResults.atrlbl = atrlbl;



