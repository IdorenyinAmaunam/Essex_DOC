function [RunResults] = analyzeOnlineStroke(FilePath, MATpath, lap, freqs)

RunResults.fine = 1;
try
    [data, header] = sload(FilePath);
catch
    disp(['Problem loading file, skipping: ' FilePath]);
    RunResults.fine = 0;
    return;
end

try
    analysis = load(MATpath);
    analysis = analysis.analysis;
    if(~isempty(strfind(MATpath,'rhrst')))
        Class = 1;
    else
        Class = 2;
    end
catch
    disp(['Problem loading classifier mat file: ' MATpath ', skipping: ' FilePath]);
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


%trials = [pos(cf) pos(cf)+dur(cf) pos(fix) pos(fix)+dur(fix)];
% Prefer cue as start, otherwise some trials are super short
trials = [pos(cue) pos(cue)+dur(cue)+dur(cf) pos(fix) pos(fix)+dur(fix)];

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

% Find features used
UsedFeat = [];
fInd = 0;

for ch=1:length(analysis.tools.features.channels)
    for fr=1:length(analysis.tools.features.bands{analysis.tools.features.channels(ch)})
        fInd = fInd + 1;

        UsedFeat = [UsedFeat ; analysis.tools.features.channels(ch) ...
            analysis.tools.features.bands{analysis.tools.features.channels(ch)}(fr)];

    end
end

% Crop all data to selected features
sfeats = zeros(size(afeats,1),size(UsedFeat,1));
for s=1:size(afeats,1)
    for f=1:size(UsedFeat,1)
        sfeats(s,f) = afeats(s,(UsedFeat(f,2)-4)/2+1,UsedFeat(f,1));
    end
end

% Crop to trials
fdata = [];
flbl = [];
ftrlbl = [];
albl = -ones(size(afeats,1),1);
atrlbl = -ones(size(afeats,1),1);
for i=1:size(ftrials,1)
    fdata = [fdata; sfeats(ftrials(i,1):ftrials(i,2),:)];
    flbl = [flbl; ones(ftrials(i,2)-ftrials(i,1)+1,1)];
    albl(ftrials(i,1):ftrials(i,2)) = 1;
    atrlbl(ftrials(i,1):ftrials(i,2)) = i;
    fdata = [fdata; sfeats(ftrials(i,3):ftrials(i,4),:)];
    flbl = [flbl; zeros(ftrials(i,4)-ftrials(i,3)+1,1)];
    albl(ftrials(i,3):ftrials(i,4)) = 0;
    atrlbl(ftrials(i,3):ftrials(i,4)) = i;
    ftrlbl = [ftrlbl; i*ones(ftrials(i,2)-ftrials(i,1)+1+ftrials(i,4)-ftrials(i,3)+1,1)];
end
RunResults.afeats = afeats;
RunResults.bfeats = bfeats;
RunResults.albl = albl;
RunResults.atrlbl = atrlbl;
RunResults.sfeats = sfeats;
RunResults.fdata = fdata;
RunResults.flbl = flbl;
RunResults.ftrlbl = ftrlbl;

% Single-sample accuracies

% Online "real" single-sample accuracy (more like detection rate)
% Take care of some bad classifiers (of jy18) where 769 is first while
% normally the taskset is [783 769]
if( (Class==2) && (analysis.settings.task.classes_old(1)==769) )
    analysis.tools.net.gau.M = flip(analysis.tools.net.gau.M);
    analysis.tools.net.gau.C = flip(analysis.tools.net.gau.C);
end
prob = [];
actdata = fdata(flbl==1,:);
for i=1:size(actdata,1)
    [usl prob(i,:)] = gauClassifier(analysis.tools.net.gau.M, analysis.tools.net.gau.C,...
        actdata(i,:));
end
[maxV maxI] = max(prob');
RunResults.OnDetectionRate = 100*sum(maxI==Class)/size(actdata,1);

% Simulated single-sample accuracy with same features
try
    fclass = classify(fdata,fdata,flbl);
    [~, ~, RunResults.SimulatedAccSelected] = eegc3_confusion_matrix(flbl+1,fclass+1);
catch
    % Sometimes selected features are too correlated for LDA, switch to Naive Bayes...
    fclass = classify(fdata,fdata,flbl,'diaglinear');
    [~, ~, RunResults.SimulatedAccSelected] = eegc3_confusion_matrix(flbl+1,fclass+1);
end

% % ERD/ERS
mbase = zeros(23,16);
bmbase = zeros(23,16);
mact = zeros(23,16);
bmact = zeros(23,16);
for tr=1:size(trials,1)
    baseline = squeeze(mean(afeats(intersect(find(atrlbl==tr),find(albl==0)),:,:),1));
    bbaseline = squeeze(mean(bfeats(intersect(find(atrlbl==tr),find(albl==0)),:,:),1));
    mbase = mbase + baseline;
    bmbase = bmbase + bbaseline;
    activity = squeeze(mean(afeats(intersect(find(atrlbl==tr),find(albl==1)),:,:),1));
    bactivity = squeeze(mean(bfeats(intersect(find(atrlbl==tr),find(albl==1)),:,:),1));
    mact = mact+activity;
    bmact = bmact+bactivity;
    erds(tr,:,:) = (activity-baseline)./baseline;
    berds(tr,:,:) = (bactivity-bbaseline)./bbaseline;
    erdsdiff(tr,:,:) = (activity-baseline);
    berdsdiff(tr,:,:) = (bactivity-bbaseline);
    % ERSP (ERSP-TB-Z, equations (7)-(10) in Single-Trial Normalization for Event-Related Spectral Decomposition Reduces Sensitivity to Noisy Trials)
    base_ersp = bfeats(intersect(find(atrlbl==tr),find(albl==0)),:,:);
    act_ersp = bfeats(intersect(find(atrlbl==tr),find(albl==1)),:,:);
    m_base_ersp = squeeze(mean(base_ersp,1));
    s_base_ersp = squeeze(std(base_ersp,1));
    act_ersp_allsamples = (act_ersp - permute(repmat(m_base_ersp,1,1,size(act_ersp,1)),[3 1 2]))./permute(repmat(s_base_ersp,1,1,size(act_ersp,1)),[3 1 2]);
    act_ersp_allsamples(act_ersp_allsamples>5)=NaN;
    act_ersp_tr(tr,:,:) = squeeze(nanmean(act_ersp_allsamples));
end
mbase = mbase/size(trials,1);
bmbase = bmbase/size(trials,1);
mact = mact/size(trials,1);
bmact = bmact/size(trials,1);
erdsavg = (mact-mbase)./mbase;
berdsavg = (bmact-bmbase)./bmbase;
ersp = squeeze(mean(act_ersp_tr,1))';

RunResults.logERDS = squeeze(mean(erds,1));
RunResults.ERDS = squeeze(mean(berds,1));
RunResults.logERDSd = squeeze(mean(erdsdiff,1));
RunResults.ERDSd = squeeze(mean(berdsdiff,1));
RunResults.logERDSa = erdsavg;
RunResults.ERDSa = berdsavg;
RunResults.ersp = ersp;

FS = zeros(size(afeats,2),size(afeats,3));
for f=1:size(afeats,2)
    for ch=1:size(afeats,3)
        FS(f,ch) = eegc3_fs2(squeeze(afeats(albl==0,f,ch)),squeeze(afeats(albl==1,f,ch)));
    end
end
RunResults.FS = FS;

% Other important metrics
RunResults.NTrials = length(find(header.EVENT.TYP==781));
RunResults.NHit = length(find(header.EVENT.TYP==897));
RunResults.TrAccA = 100*(RunResults.NHit/RunResults.NTrials);
RunResults.DurAll = header.EVENT.DUR(find(header.EVENT.TYP==781))/512;
RunResults.AvgDurAll = mean(RunResults.DurAll);
tmpTYP = [header.EVENT.TYP; 1];
RunResults.IndHit = find(tmpTYP(find(tmpTYP==781)+1)==897);
RunResults.AvgDurHit = mean(RunResults.DurAll(RunResults.IndHit));
RunResults.UsedFeat = UsedFeat;

[usl sortedFSInd] = sort(FS(:),'descend');
sortedFSInd = sortedFSInd(1:5);
[Bfr Bch] = ind2sub([23,16],sortedFSInd);
Bdata1 = [];
Bdata0 = [];
for i=1:length(Bch)
    Bdata1 = [Bdata1  squeeze(afeats(albl==1,Bfr(i),Bch(i)))];
    Bdata0 = [Bdata0  squeeze(afeats(albl==0,Bfr(i),Bch(i)))];
end
Bdata = [Bdata0 ; Bdata1];
Blbl = [zeros(size(Bdata0,1),1) ; ones(size(Bdata1,1),1)];
M1 = mean(Bdata1);
C1 = cov(Bdata1);
M0 = mean(Bdata0);
C0 = cov(Bdata0);
RunResults.SI = KLNormMulti(M1,C1,M0,C0);

% Simulated single-sample accuracy with best features
try
    Bclass = classify(Bdata,Bdata,Blbl);
    [~, ~, RunResults.SimulatedAccBest] = eegc3_confusion_matrix(Blbl+1,Bclass+1);
catch
    Bclass = classify(Bdata,Bdata,Blbl,'diaglinear');
    [~, ~, RunResults.SimulatedAccBest] = eegc3_confusion_matrix(Blbl+1,Bclass+1);    
end