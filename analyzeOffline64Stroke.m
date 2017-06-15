function [RunResults] = analyzeOffline64Stroke(FilePath, lap64, freqs)

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

if(~isempty(strfind(FilePath,'rhrst')))
    Class = 1;
else
    Class = 2;
end

% Find re-referencing signal
ref = (data(:,63)+data(:,64))/2;
data = data(:,1:62) - repmat(ref,[1 62]);

% Remove non-relevant channels
% Remove overall DC
data = data-repmat(mean(data),size(data,1),1);


% Laplacian spatial filtering
data = laplacianSP(data,lap64);

% Trial extraction
pos = header.EVENT.POS;
dur = header.EVENT.DUR;
cf = find(header.EVENT.TYP==781);
cue = cf-1;
fix = find(header.EVENT.TYP==786);

if (length(cf) < 23)
    % Too few trials, there should be 15 + 8 (rest) = 23
    RunResults.fine = 0;
    disp('Too few trials, skipping...');
    return;
end


% Prefer 781 as start, since for this data I have 4 seconds standard
trials = [pos(cf) pos(cf)+dur(cf) pos(fix) pos(fix)+dur(fix)];
trialsLbl = header.EVENT.TYP(cue);

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
for ch=1:62
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
albl = -ones(size(afeats,1),1);
atrlbl = -ones(size(afeats,1),1);
for i=1:size(ftrials,1)
    if(trialsLbl(i)==783)
        albl(ftrials(i,1):ftrials(i,2)) = 2; % 2 corresponds to "rest"
    else
        albl(ftrials(i,1):ftrials(i,2)) = 1; % 1 corresponds to the "active" class, either left or right hand
    end
    atrlbl(ftrials(i,1):ftrials(i,2)) = i;
    albl(ftrials(i,3):ftrials(i,4)) = 0;
    atrlbl(ftrials(i,3):ftrials(i,4)) = i;
end
RunResults.afeats = afeats;
RunResults.bfeats = bfeats;
RunResults.albl = albl;
RunResults.atrlbl = atrlbl;

% % ERSP
for tr=1:size(trials,1)
    % ERSP (ERSP-TB-Z, equations (7)-(10) in Single-Trial Normalization for Event-Related Spectral Decomposition Reduces Sensitivity to Noisy Trials)
    base_ersp = bfeats(intersect(find(atrlbl==tr),find(albl==0)),:,:);
    act_ersp = bfeats(intersect(find(atrlbl==tr),union(find(albl==1),find(albl==2))),:,:);
    m_base_ersp = squeeze(mean(base_ersp,1));
    s_base_ersp = squeeze(std(base_ersp,1));
    act_ersp_allsamples = (act_ersp - permute(repmat(m_base_ersp,1,1,size(act_ersp,1)),[3 1 2]))./permute(repmat(s_base_ersp,1,1,size(act_ersp,1)),[3 1 2]);
    act_ersp_allsamples(act_ersp_allsamples>5)=NaN;
    act_ersp_tr(tr,:,:) = squeeze(nanmean(act_ersp_allsamples));
end
erspAct = squeeze(nanmean(act_ersp_tr(find(trialsLbl~=783),:,:),1))';
erspRest = squeeze(nanmean(act_ersp_tr(find(trialsLbl==783),:,:),1))';

RunResults.erspAct = erspAct;
RunResults.erspRest = erspRest;

FS = zeros(size(afeats,2),size(afeats,3));
for f=1:size(afeats,2)
    for ch=1:size(afeats,3)
        FS(f,ch) = eegc3_fs2(squeeze(afeats(albl==0,f,ch)),squeeze(afeats(albl==1,f,ch)));
    end
end
RunResults.FS = FS;

% Other important metrics
RunResults.NTrials = length(find(header.EVENT.TYP==781));

[usl sortedFSInd] = sort(FS(:),'descend');
sortedFSInd = sortedFSInd(1:5);
[Bfr Bch] = ind2sub([23,62],sortedFSInd);
Bdata1 = [];
Bdata0 = [];
Bdata2 = [];
for i=1:length(Bch)
    Bdata1 = [Bdata1  squeeze(afeats(albl==1,Bfr(i),Bch(i)))];
    Bdata0 = [Bdata0  squeeze(afeats(albl==0,Bfr(i),Bch(i)))];
    Bdata2 = [Bdata2  squeeze(afeats(albl==2,Bfr(i),Bch(i)))];
end
Bdata01 = [Bdata0 ; Bdata1];
Blbl01 = [zeros(size(Bdata0,1),1) ; ones(size(Bdata1,1),1)];
Bdata02 = [Bdata0 ; Bdata2];
Blbl02 = [zeros(size(Bdata0,1),1) ; ones(size(Bdata2,1),1)];
Bdata12 = [Bdata1 ; Bdata2];
Blbl12 = [zeros(size(Bdata1,1),1) ; ones(size(Bdata2,1),1)];
M1 = mean(Bdata1);
C1 = cov(Bdata1);
M2 = mean(Bdata2);
C2 = cov(Bdata2);
M0 = mean(Bdata0);
C0 = cov(Bdata0);
RunResults.SI01 = KLNormMulti(M1,C1,M0,C0);
RunResults.SI02 = KLNormMulti(M2,C2,M0,C0);
RunResults.SI12 = KLNormMulti(M1,C1,M2,C2);

% Single-sample accuracies
% Simulated single-sample accuracy with best features
try
    Bclass01 = classify(Bdata01,Bdata01,Blbl01);
    [~, ~, RunResults.SimulatedAccBest01] = eegc3_confusion_matrix(Blbl01+1,Bclass01+1);
catch
    Bclass01 = classify(Bdata01,Bdata01,Blbl01,'diaglinear');
    [~, ~, RunResults.SimulatedAccBest01] = eegc3_confusion_matrix(Blbl01+1,Bclass01+1);
end

try
    Bclass02 = classify(Bdata02,Bdata02,Blbl02);
    [~, ~, RunResults.SimulatedAccBest02] = eegc3_confusion_matrix(Blbl02+1,Bclass02+1);
catch
    Bclass02 = classify(Bdata02,Bdata02,Blbl02,'diaglinear');
    [~, ~, RunResults.SimulatedAccBest02] = eegc3_confusion_matrix(Blbl02+1,Bclass02+1);
end

try
    Bclass12 = classify(Bdata12,Bdata12,Blbl12);
    [~, ~, RunResults.SimulatedAccBest12] = eegc3_confusion_matrix(Blbl12+1,Bclass12+1);
catch
    Bclass12 = classify(Bdata12,Bdata12,Blbl12,'diaglinear');
    [~, ~, RunResults.SimulatedAccBest12] = eegc3_confusion_matrix(Blbl12+1,Bclass12+1);
end