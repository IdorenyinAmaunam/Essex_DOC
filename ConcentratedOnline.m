addpath(genpath('/home/ido/New DOC Literature/DOCpipeline/Code/biosig/'));
FigureSavePath = '~/New DOC Literature/DOCpipeline/Code/code/SavedData/';
Path = '/home/ido/New DOC Literature/DOCpipeline/Code/code/EEG';

SubID  = {'qq62','ckg8','ds86','fh47','jy18','ma93','qv39','rj31','wu60','ya00',...
    'odr2','ji34','ao48','pk72','lm90', 'rai1','mo17','cy97','qs36','vxd9',...
    'eo60','ia41','b2sc','och9','mdl9','se34','jv20','mnv4'};
Group  = {'BCI','BCI','Sham','Sham','BCI','BCI','Sham','BCI','BCI','BCI','BCI','Sham',...
    'BCI','BCI','Sham','BCI','BCI','Sham','Sham','Sham','Sham','Sham','Sham',...
    'BCI','BCI','BCI','Sham','Sham'};
Side = {'L','L','L','L','L','L','L','L','L','R','L','R','R','L','L','R','L','L',...
    'R','L','R','L','R','R','L','L','L','L'};

Electrodes16 = {'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'};

% Colors for BCI and Sham
Color = {[255 40 40]/255,[102 184 184]/255};

nBCI=0;
nSham=0;
for s=1:length(Group)
    if(strcmp(Group{s},'BCI'))
        nBCI=nBCI+1;
        SubPaperID{s} = ['BCI' num2str(nBCI,'%02d')];
    else
        nSham=nSham+1;
        SubPaperID{s} = ['sham' num2str(nSham,'%02d')];
    end
end

Electrodes16 = {'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'};

% Clinical outcomes
FMAUE{1} = [NaN 8 0 0 0 0  0 0 23  6 12 30  0  0 4 21 10  0 25  4 NaN NaN 10 8 18 NaN NaN 4];
FMAUE{2} = [NaN 9 0 8 0 1 10 0 47 32 40 43  3 35 9 38 10 34 40 34 NaN NaN 17 7 18 NaN NaN 4];
FMAdiff = FMAUE{2}-FMAUE{1};

ESS{1} = [NaN 48 42 39 48 57 52 57 74 57 72 70 61 57 42 NaN NaN NaN NaN 4 NaN NaN 45 55 62 NaN NaN 57 ];
ESS{2} = [NaN 48 54 72 64 67 82 61 84 81 84 80 73 93 56 NaN NaN NaN NaN 71 NaN NaN 56 62 74 NaN NaN 58];
ESSdiff = ESS{2}-ESS{1};

BarthelIndex{1} = [NaN 25 55 5  75 50 10 35 70 50 50 70 15  60  20 NaN NaN NaN NaN 15 NaN NaN 10 45 50 NaN NaN 45];
BarthelIndex{2} = [NaN 40 60 65 85 60 60 65 85 55 80 70 100 100 40 NaN NaN NaN NaN 50 NaN NaN 5 55 55 NaN NaN 65];
BarthelIndexdiff = BarthelIndex{2}-BarthelIndex{1};

% Print info for clinical outcome statistical testing
printCO(FMAUE{1}, FMAUE{2}, Group, 'FMA-UE');
printCO(ESS{1}, ESS{2}, Group, 'ESS');
printCO(BarthelIndex{1}, BarthelIndex{2}, Group, 'Barthel Index');

for sub=1:length(SubID)
    
    if(strcmp(Side{sub},'L'))
        affectedtask = 'lhrst';
        probindex = 2;
    else
        affectedtask = 'rhrst';
        probindex = 1;
    end
       
    % Initialize per-group variables here
    PG.HitSum = 0;

    if(exist([Path '/' SubID{sub} '/' SubID{sub} '_Res.mat'],'file'))
        load([Path '/' SubID{sub} '/' SubID{sub} '_Res.mat']);
    else
        
        % Arrange sessions in temporal order
        RunDir = dir([Path '/' SubID{sub}]);
        RunDir = RunDir(3:end);
        
        % Keep only online 16 runs
        KeepInd = [];
        for r=1:length({RunDir.name})
            if( (~isempty(strfind(RunDir(r).name,'online'))) && (~isempty(strfind(RunDir(r).name,affectedtask))))
                KeepInd = [KeepInd; r];
            end
        end
        RunDir = RunDir(KeepInd);        
        FileNames = {RunDir.name};
        
        Date = {};
        for f=1:length(FileNames)
            Dots = strfind(FileNames{f},'.');
            Date{f} = FileNames{f}(Dots(1)+1:Dots(2)-1);
        end
        SesDates = unique(Date);
        SesDatesNum = datenum(SesDates,'yyyymmdd');
        [SortedSesDatesNum SortedSesDatesNumInd] = sort(SesDatesNum,'ascend');
        SesDates = SesDates(SortedSesDatesNumInd);
        Nses = length(SesDates);
        RunSesInd = [];
        for f=1:length(FileNames)
            RunSesInd(f) = find(strcmp(SesDates,Date{f}));
        end
        
        %% Initialize ALL variables
        AllHit = [];
        AllHitRate = [];
        AllTrialDuration = [];
        AllHitTrialDuration = [];
        AllDetectionRate = [];
        AllSimAccSel = [];
        AllSimAccBest = [];
        AllSI = [];
        AllFS = {};
        AlllogERDSd = {};
        AlllogERDSdpr = {};
        AllERSP = {};
        
        AllERSPLastHit = zeros(23,16,0);
        AllERSPLastMiss = zeros(23,16,0);
        
        AllPSDLastHit = [];
        AllPSDLastMiss = [];
        AllPSDRest = [];

        AllPSDLastHitClassif = [];
        AllPSDLastMissClassif = [];
            
        AllPSDHitClassif = [];
        AllPSDMissClassif = [];
        
        AllPSDLastClassifHitMiss = [];
        AllFESHitMiss = [];
        
        AllFSSel = {};
        
        AllArtifactPrctClean = []; 
        AllArtifactPrctNoisy = []; 
        
        AllTrIntegProbs = nan(0,400);
        
        AllSesInd = [];
        
        % Analyze per session
        for ses=1:length(SesDates)
            ThisSesRunind = find(RunSesInd==ses);

            % Initialize SES variables
            SesHit = [];
            SesHitRate = [];
            SesTrialDuration = [];
            SesHitTrialDuration = [];
            SesDetectionRate = [];
            SesSimAccSel = [];
            SesSimAccBest = [];
            SesSI = [];
            SesFS = {};
            SeslogERDSd = {};
            SeslogERDSdpr = {};
            SesERSP = {};

            for run=1:length(ThisSesRunind)
                load([Path '/' SubID{sub} '/' FileNames{ThisSesRunind(run)}]);
                disp(['Processing Subject:' SubID{sub} ', Session: ' num2str(ses) ', Run: ' num2str(run)]);
                AllSesInd = [AllSesInd ses];
                
                SesHit = [SesHit RunResults.NHit];
                AllHit = [AllHit RunResults.NHit];
                SesHitRate = [SesHitRate RunResults.TrAccA];
                AllHitRate = [AllHitRate RunResults.TrAccA];            
                SesTrialDuration = [SesTrialDuration ; RunResults.DurAll];
                SesHitTrialDuration = [SesHitTrialDuration ; RunResults.DurAll(RunResults.IndHit)];
                AllTrialDuration = [AllTrialDuration ; RunResults.DurAll];
                AllHitTrialDuration = [AllHitTrialDuration ; RunResults.DurAll(RunResults.IndHit)];
                AllDetectionRate = [AllDetectionRate ; RunResults.OnDetectionRate];
                SesDetectionRate = [SesDetectionRate ; RunResults.OnDetectionRate];
                SesSimAccSel = [SesSimAccSel ; RunResults.SimulatedAccSelected];
                SesSimAccBest = [SesSimAccBest ; RunResults.SimulatedAccBest];
                AllSimAccSel = [AllSimAccSel ; RunResults.SimulatedAccSelected];
                AllSimAccBest = [AllSimAccBest ; RunResults.SimulatedAccBest];
                AllSI = [AllSI; RunResults.SI];
                SesSI = [SesSI; RunResults.SI];
                SesFS = [SesFS ; RunResults.FS'];
                AllFS = [AllFS ; RunResults.FS'];
                SeslogERDSd = [SeslogERDSd ; RunResults.logERDSd'];
                AlllogERDSd = [AlllogERDSd ; RunResults.logERDSd'];
                act = squeeze(mean(RunResults.afeats(RunResults.albl==1,:,:),1));
                base = squeeze(mean(RunResults.afeats(RunResults.albl==0,:,:),1));
                SeslogERDSdpr = [SeslogERDSdpr ; (act-base)'];
                AlllogERDSdpr = [AlllogERDSdpr ; (act-base)'];
                SesERSP = [SesERSP ; RunResults.ersp];
                AllERSP = [AllERSP ; RunResults.ersp];
                
                % TPR and FPR of last decision
%                IndRealHit = find( RunResults.TrialOutcome > 30000 | RunResults.TrialOutcome==14 | RunResults.TrialOutcome==897);
%                IndRealMiss = setdiff([1:max(RunResults.atrlbl)],IndRealHit);
%                AllPSDLastHitClassif = [AllPSDLastHitClassif; RunResults.TrialLastAct(IndRealHit)'];
%                AllPSDLastMissClassif = [AllPSDLastMissClassif; RunResults.TrialLastAct(IndRealMiss)'];

               IndRealHit = find( RunResults.TrialOutcome > 30000 | RunResults.TrialOutcome==14 | RunResults.TrialOutcome==897);
               IndRealMiss = setdiff([1:max(RunResults.atrlbl)],IndRealHit);
               AllPSDLastHitClassif = [AllPSDLastHitClassif; RunResults.TrialLastAct(IndRealHit)'];
               AllPSDLastMissClassif = [AllPSDLastMissClassif; RunResults.TrialLastAct(IndRealMiss)'];

                % Latest test for Nat Comm revision
                FESHitMiss = zeros(1,length(RunResults.TrialOutcome));
                FESHitMiss(IndRealHit)=1;
                PSDHitMiss = RunResults.TrialLastAct;
                AllPSDLastClassifHitMiss = [AllPSDLastClassifHitMiss PSDHitMiss];
                AllFESHitMiss = [AllFESHitMiss FESHitMiss];               
               
               
                tmpERSP = [];
                for atr=1:max(RunResults.atrlbl)
                    tmp = squeeze(RunResults.afeats(intersect(find(RunResults.atrlbl==1),find(RunResults.albl==1)),:,:));
                    tmpM = squeeze(mean(tmp,1));
                    tmpS = squeeze(std(tmp,1));
                    tmpERSP(:,:,atr) = (squeeze(RunResults.LastPSDAll(atr,:,:)) - tmpM)./tmpS;
                end
                
                AllERSPLastHit = cat(3,AllERSPLastHit,tmpERSP(:,:,IndRealHit));
                AllERSPLastMiss = cat(3,AllERSPLastMiss,tmpERSP(:,:,IndRealMiss));
                
                RunResults.LastPSDAll = permute(RunResults.LastPSDAll,[2 3 1]);

                AllPSDLastHit = cat(3,AllPSDLastHit,squeeze(RunResults.LastPSDAll(:,:,IndRealHit)));
                AllPSDLastMiss = cat(3,AllPSDLastMiss,squeeze(RunResults.LastPSDAll(:,:,IndRealMiss)));
                AllPSDRest = cat(1,AllPSDRest,squeeze(RunResults.afeats(intersect(find(RunResults.atrlbl==atr),find(RunResults.albl==0)),:,:)));

                %% TPR and FPR on all, not only last sec
                for ss=1:size(RunResults.actdata,1)
                    [~, prob(ss,:)] = gauClassifier(RunResults.analysis.tools.net.gau.M,...
                        RunResults.analysis.tools.net.gau.C,...
                        RunResults.actdata(ss,:));
                end
                [maxV maxI] = max(prob');
                Class = find(RunResults.analysis.settings.task.classes_old~=783);
                ActAcc = zeros(1,length(maxI));
                ActAcc(maxI==Class)=1;
                RunTrialInd = RunResults.atrlbl(RunResults.albl==1);
                TPR = ActAcc(ismember(RunTrialInd,IndRealHit));
                FPR = ActAcc(ismember(RunTrialInd,IndRealMiss));
                AllPSDHitClassif = [AllPSDHitClassif; TPR'];
                AllPSDMissClassif = [AllPSDMissClassif; FPR'];
                
                % Fisher scores of selected features
                FSsel = [];
                for f=1:size(RunResults.UsedFeat,1)
                    FSsel = [FSsel RunResults.FS((RunResults.UsedFeat(f,2)-4)/2+1,RunResults.UsedFeat(f,1))];
                end
                AllFSSel{end+1} = FSsel;
                
                %% Percentage of artifact
                %AllArtifactPrctClean = [AllArtifactPrctClean ; RunResults.artifact.clean];
                %AllArtifactPrctNoisy = [AllArtifactPrctNoisy ; RunResults.artifact.noisy];
                
                
                % Simulated hit rates
                for tr=1:max(RunResults.atrlbl)
                    trcprobs = RunResults.probs(RunResults.atrlbl(find(RunResults.albl==1))==tr,probindex);
                    AllTrIntegProbs(end+1,1) = 0.5;
                    for p=1:length(trcprobs)
                        AllTrIntegProbs(end,p+1) = 0.96*AllTrIntegProbs(end,p)+0.04*trcprobs(p);
                    end
                end
                
            end

            %% Per Session averages
            %% Hit
            SubResults.HitSes(ses) = nanmean(SesHit);
            SubResults.HitRateSes(ses) = nanmean(SesHitRate);
            %% Trial durations
            SubResults.TrialDurationSes(ses) = nanmean(SesTrialDuration);
            SubResults.HitTrialDurationSes(ses) = nanmean(SesHitTrialDuration);
            %% Accuracies
            SubResults.DetectionRateSes(ses) = nanmean(SesDetectionRate);
            SubResults.SimAccSelSes(ses) = nanmean(SesSimAccSel);
            SubResults.SimAccBestSes(ses) = nanmean(SesSimAccBest);
            % SI
            SubResults.SISes(ses) = nanmean(SesSI);
            %FS
            Nrun = length(SesFS);
            SubResults.FSSes{ses} = squeeze(nanmean(reshape(cell2mat(SesFS),[Nrun 16 23]),1));
            % ERDSd
            SubResults.logERDSdSes{ses} = squeeze(nanmean(reshape(cell2mat(SeslogERDSd),[Nrun 16 23]),1));
            % ERDSdpr (first average across trials, then take the difference)
            SubResults.logERDSdprSes{ses} = squeeze(nanmean(reshape(cell2mat(SeslogERDSdpr),[Nrun 16 23]),1));
            %ERSP
            Nrun = length(SesERSP);
            SubResults.ERSPSes{ses} = squeeze(nanmean(reshape(cell2mat(SesERSP),[Nrun 16 23]),1));
        end

        SubResults.RunSesInd = RunSesInd; 
        %% Hit
        SubResults.HitAll = AllHit;
        SubResults.HitSum = sum(AllHit);
        SubResults.HitRateAll = AllHitRate;
        %% Trial durations
        SubResults.TrialDurationAll = AllTrialDuration;
        SubResults.HitTrialDurationAll = AllHitTrialDuration;   
        %% Accuracies
        SubResults.DetectionRateAll = AllDetectionRate;
        SubResults.SimAccSelAll = AllSimAccSel;
        SubResults.SimAccBestAll = AllSimAccBest;
        % SI
        SubResults.SIAll = AllSI;
        %FS
        Nrun = length(AllFS);
        SubResults.FSAll = squeeze(nanmean(reshape(cell2mat(AllFS),[Nrun 16 23]),1));
        % ERDSd
        SubResults.logERDSdAll = squeeze(nanmean(reshape(cell2mat(AlllogERDSd),[Nrun 16 23]),1));
        % ERDSdpr
        SubResults.logERDSdprAll = squeeze(nanmean(reshape(cell2mat(AlllogERDSdpr),[Nrun 16 23]),1));
        % ERSP
        Nrun = length(AllERSP);
        SubResults.ERSPAll = squeeze(nanmean(reshape(cell2mat(AllERSP),[Nrun 16 23]),1));        
        %% ERSP last second hit
        %SubResults.ERSPAllLastHit = AllERSPLastHit;

        SubResults.LastTPR = 100*sum(AllPSDLastHitClassif)/length(AllPSDLastHitClassif);
        SubResults.LastFPR = 100*sum(AllPSDLastMissClassif)/length(AllPSDLastMissClassif);
        
        SubResults.TPR = 100*sum(AllPSDHitClassif)/length(AllPSDHitClassif);
        SubResults.FPR = 100*sum(AllPSDMissClassif)/length(AllPSDMissClassif);        
        
        SubResults.ERSPLastHit = squeeze(mean(AllERSPLastHit,3));
        SubResults.ERSPLastMiss = squeeze(mean(AllERSPLastMiss,3));
        
        AllPSDLastHit = permute(AllPSDLastHit,[3 1 2]);
        AllPSDLastMiss = permute(AllPSDLastMiss,[3 1 2]);
        tmpFSLblHit = [ones(size(AllPSDLastHit,1),1) ; 2*ones(size(AllPSDRest,1),1)];
        tmpFSLblMiss = [ones(size(AllPSDLastMiss,1),1) ; 2*ones(size(AllPSDRest,1),1)];

        SubResults.FSHit = fisherScore(cat(1,AllPSDLastHit,AllPSDRest),tmpFSLblHit);
        SubResults.FSMiss = fisherScore(cat(1,AllPSDLastMiss,AllPSDRest),tmpFSLblMiss);
        
        SubResults.AllFSSel = AllFSSel;
        SubResults.AllPSDLastClassifHitMiss = AllPSDLastClassifHitMiss;
        SubResults.AllFESHitMiss = AllFESHitMiss;
              
        SubResults.AllArtifactPrctClean = AllArtifactPrctClean;
        SubResults.AllArtifactPrctNoisy = AllArtifactPrctNoisy;
        
        SubResults.AllTrIntegProbs = AllTrIntegProbs;
        
        % Create Feature Selection map
        SubResults.AllSelFeat = zeros(16,23);
        for f=1:size(RunResults.UsedFeat,1)
            SubResults.AllSelFeat(RunResults.UsedFeat(f,1),(RunResults.UsedFeat(f,2)-4)/2+1) = 1; 
        end
        % Mirrored
        if(strcmp(Side{sub},'L'))
            doflip = 1;
        else
            doflip = 0;
        end 
        SubResults.AllSelFeatMirror = flipSide16Matrix(SubResults.AllSelFeat,doflip);
        
        SubResults.AllSesInd = AllSesInd;
        
        % Save subject results
        if(~exist([Path '/' SubID{sub} '/' SubID{sub} '_Res.mat'],'file'))
            save([Path '/' SubID{sub} '/' SubID{sub} '_Res.mat'],'SubResults');
        end
    end
    
    % Per-subject variables
    PS.SesNum(sub) = length(SubResults.HitSes);
    PS.HitSum(sub) = SubResults.HitSum;
    PS.HitPerRun(sub) = nanmean(SubResults.HitAll);
    PS.TrialDur(sub) = nanmean(SubResults.TrialDurationAll);
    PS.HitDur(sub) = nanmean(SubResults.HitTrialDurationAll);
    tmpHit = [];
    tmpRunsPerSes = [];
    for s=1:max(SubResults.RunSesInd)
        tmpHit(s) = sum(SubResults.HitAll(SubResults.RunSesInd==s));
        tmpRunsPerSes(s) = sum(SubResults.RunSesInd==s);
    end
    PS.RunsPerSes(sub) = nanmean(tmpRunsPerSes);
    PS.HitPerSession(sub) = nanmean(tmpHit);
    PS.HitRate(sub) = nanmean(SubResults.HitRateAll);
    PS.DetectionRate(sub) = nanmean(SubResults.DetectionRateAll);
    PS.SimAccSel(sub) = nanmean(SubResults.SimAccSelAll);
    PS.SimAccBest(sub) = nanmean(SubResults.SimAccBestAll);
    PS.SI(sub) = nanmean(SubResults.SIAll);
    PS.MaxFS(sub) = max(SubResults.FSAll(:));
    
    PS.DetectionRateSes{sub} = SubResults.DetectionRateSes;
    PS.SimAccSelSes{sub} = SubResults.SimAccSelSes;
    PS.SimAccBestSes{sub} = SubResults.SimAccBestSes;
    PS.SISes{sub} = SubResults.SISes;
    PS.FSSes{sub} = SubResults.FSSes;
    PS.logERDSdSes{sub} = SubResults.logERDSdSes;
    PS.logERDSdprSes{sub} = SubResults.logERDSdprSes;
    PS.ERSPSes{sub} = SubResults.ERSPSes;
    PS.FSavg{sub} = SubResults.FSAll';
    PS.ERSPavg{sub} = SubResults.ERSPAll';
    %PS.ERSPAllLastHit{sub} = SubResults.ERSPAllLastHit;
    %PS.LastPSDFS{sub} = SubResults.AllPSDBaseFS; 
    PS.LastTPR(sub) = SubResults.LastTPR;
    PS.LastFPR(sub) = SubResults.LastFPR;
    
    PS.TPR(sub) = SubResults.TPR;
    PS.FPR(sub) = SubResults.FPR;    
    
    PS.ERSPLastHit{sub} = SubResults.ERSPLastHit;
    PS.ERSPLastMiss{sub} = SubResults.ERSPLastMiss;
    PS.FSLastHit{sub} = SubResults.FSHit;
    PS.FSLastMiss{sub} = SubResults.FSMiss;
    
    PS.FSSel{sub} = SubResults.AllFSSel;
    PS.DetectionRateAll{sub} = SubResults.DetectionRateAll;
    
    
    PS.AllPSDLastClassifHitMiss{sub} = SubResults.AllPSDLastClassifHitMiss;
    PS.AllFESHitMiss{sub} = SubResults.AllFESHitMiss;
%     [PS.LastPSDConfMat(sub,:,:) PS.LastPSDConfMatAll(sub,:,:) PS.LastPSDAcc(sub)] = ...
%         confusion_matrix(SubResults.AllPSDLastClassifHitMiss, SubResults.AllFESHitMiss,2);
%     
    %% This is the matrix I finally want to use, but I need to make 1=positive 2=negative
    % Also, Ground Truth I consider the BCI classification and prediction
    % the FES Hit/Miss
    [PS.LastPSDConfMat(sub,:,:) PS.LastPSDConfMatAll(sub,:,:) PS.LastPSDAcc(sub)] = ...
        confusion_matrix(1-SubResults.AllPSDLastClassifHitMiss, 1-SubResults.AllFESHitMiss, 2);
    
    
    PS.ArtPrctC(sub,:) = mean(SubResults.AllArtifactPrctClean);
    PS.ArtPrctN(sub,:) = mean(SubResults.AllArtifactPrctNoisy);
    
    PS.AllTrIntegProbs{sub} = SubResults.AllTrIntegProbs;
    
    PS.SelFeatMap(sub,:,:) = SubResults.AllSelFeat;
    PS.SelFeatMirrorMap(sub,:,:) = SubResults.AllSelFeatMirror;
    
    PS.AllSesInd{sub} = SubResults.AllSesInd;
    clear SubResults;
end

IndBCI = find(strcmp(Group,'BCI'));
IndSham = find(strcmp(Group,'Sham'));
% Present results
disp(['Number of sessions ==> BCI: ' num2str(mean(PS.SesNum(IndBCI))) ' +- ' num2str(std2(PS.SesNum(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.SesNum(IndSham))) ' +- ' num2str(std2(PS.SesNum(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.SesNum(IndBCI),PS.SesNum(IndSham)))]);

disp(['Number of runs per session ==> BCI: ' num2str(mean(PS.RunsPerSes(IndBCI))) ' +- ' num2str(std2(PS.RunsPerSes(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.RunsPerSes(IndSham))) ' +- ' num2str(std2(PS.RunsPerSes(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.RunsPerSes(IndBCI),PS.RunsPerSes(IndSham)))]);

disp(['Total hits (FES stimulations) ==> BCI: ' num2str(mean(PS.HitSum(IndBCI))) ' +- ' num2str(std2(PS.HitSum(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.HitSum(IndSham))) ' +- ' num2str(std2(PS.HitSum(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.HitSum(IndBCI),PS.HitSum(IndSham)))]);

disp(['Hits per run ==> BCI: ' num2str(mean(PS.HitPerRun(IndBCI))) ' +- ' num2str(std2(PS.HitPerRun(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.HitPerRun(IndSham))) ' +- ' num2str(std2(PS.HitPerRun(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.HitPerRun(IndBCI),PS.HitPerRun(IndSham)))]);

disp(['Hits per session ==> BCI: ' num2str(mean(PS.HitPerSession(IndBCI))) ' +- ' num2str(std2(PS.HitPerSession(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.HitPerSession(IndSham))) ' +- ' num2str(std2(PS.HitPerSession(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.HitPerSession(IndBCI),PS.HitPerSession(IndSham)))]);

disp(['Hit rate ==> BCI: ' num2str(mean(PS.HitRate(IndBCI))) ' +- ' num2str(std2(PS.HitRate(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.HitRate(IndSham))) ' +- ' num2str(std2(PS.HitRate(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.HitRate(IndBCI),PS.HitRate(IndSham)))]);

disp(['Trial duration ==> BCI: ' num2str(mean(PS.TrialDur(IndBCI))) ' +- ' num2str(std2(PS.TrialDur(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.TrialDur(IndSham))) ' +- ' num2str(std2(PS.TrialDur(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.TrialDur(IndBCI),PS.TrialDur(IndSham)))]);

disp(['Trial duration (Hit Only) ==> BCI: ' num2str(mean(PS.HitDur(IndBCI))) ' +- ' num2str(std2(PS.HitDur(IndBCI))) ...
    '  Sham: ' num2str(mean(PS.HitDur(IndSham))) ' +- ' num2str(std2(PS.HitDur(IndSham))) ...
    '  p-val: ' num2str(ranksum(PS.HitDur(IndBCI),PS.HitDur(IndSham)))]);


%% Figures for revision of chronic stroke paper in Nature Communications


%% Find threshold from integrated probabilities
for sub=1:length(SubID)
    FoundTh{sub} = [];
    for run=1:floor(size(PS.AllTrIntegProbs{sub},1)/15)% There are two subjects with one trial less...
        RunFoundTh = [];
        for tr=1:15
            if(PS.AllFESHitMiss{sub}((run-1)*15+tr)==1) % Exclude miss trials
                EndProbInd = find(PS.AllTrIntegProbs{sub}((run-1)*15+tr,:)~=0);
                EndProbBefore = PS.AllTrIntegProbs{sub}((run-1)*15+tr,EndProbInd(end-1));
                RunFoundTh(end+1) = round(100*EndProbBefore)/100;
            end
        end
        if(~isempty(RunFoundTh))
            FoundTh{sub}(end+1) = mean(RunFoundTh);
        else
            FoundTh{sub}(end+1) = NaN;
        end
    end
end

for s=1:length(IndBCI)
    MaxRun(s) = length(FoundTh{s});
end
MaxRun = max(MaxRun);

FTh = nan(length(IndBCI),MaxRun);
for s=1:length(IndBCI)
    Fth(s,1:length(FoundTh{s})) = FoundTh{s};
end
Fth(Fth==0)=NaN;
Fth = Fth(:,1:60) % Stop at 60 runs becase above that there are too few patients

figure(132);
shadedErrorBar(1:60,nanmean(Fth),nanstd(Fth),'b');
hold on; plot(nanmean(Fth),'-b','LineWidth',3);
P = polyfit([1:60],nanmean(Fth),1);
f=@(x)(P(1).*x+P(2));
plot([1:60],f([1:60]),'-r','LineWidth',3);
hold off;
axis([1 60 0.5 1.0]);
xlabel('Run index (chronological)','FontSize',20);
ylabel('Average Confidence Threshold value','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);


figure(134);
DetectRateSes = nan(length(IndBCI),12);
for s=1:length(IndBCI)
    DetectRateSes(s,1:length(PS.DetectionRateSes{s})) = PS.DetectionRateSes{s};
end
DetectRateSes = DetectRateSes(:,1:10);
shadedErrorBar(1:10,nanmean(DetectRateSes),nanstd(DetectRateSes),'b');
hold on; plot(nanmean(DetectRateSes),'-b','LineWidth',3);
P = polyfit([1:10],nanmean(DetectRateSes),1);
f=@(x)(P(1).*x+P(2));
plot([1:10],f([1:10]),'-r','LineWidth',3);
hold off;
axis([1 10 20 120]);
xlabel('Session index (chronological)','FontSize',20);
ylabel('Detection Rate (%)','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);

figure(133);
DetectRate = nan(length(IndBCI),60);
for s=1:length(IndBCI)
    if(length(PS.DetectionRateAll{s}) > 60)
        DetectRate(s,1:60) = PS.DetectionRateAll{s}(1:60);
    else
        DetectRate(s,1:length(PS.DetectionRateAll{s})) = PS.DetectionRateAll{s};
    end
end
shadedErrorBar(1:60,nanmean(DetectRate),nanstd(DetectRate),'b');
hold on; plot(nanmean(DetectRate),'-b','LineWidth',3);
P = polyfit([1:60],nanmean(DetectRate),1);
f=@(x)(P(1).*x+P(2));
plot([1:60],f([1:60]),'-r','LineWidth',3);
hold off;
axis([1 60 20 100]);
xlabel('Run index (chronological)','FontSize',20);
ylabel('Detection Rate (%)','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);

%% For Supplementary Figure S3, recalculate Hit Rate with conservative threshold
% because before I was showing the detection rate averaged over
% runs/sessions
% Simulated hit rates per threshold
GAHitRateSim055 = nan(length(SubID),12);
GAHitRateSim055Run = nan(length(SubID),100);
GAThSes = nan(length(IndBCI),12);
Decth = 0.55;
for s=1:length(SubID)
    for ses=1:max(PS.AllSesInd{s})
        MaxIProb = [];
        IndRun = find(PS.AllSesInd{s}==ses);
        TrInd = (IndRun(1)-1)*15+1:IndRun(end)*15;
        GAHitRateSim055(s,ses) = 100*sum(max(PS.AllTrIntegProbs{s}(TrInd,:)') > Decth)/length(TrInd);
        if(s<=max(IndBCI))
            GAThSes(s,ses) = nanmean(FoundTh{s}(IndRun));
        end
    end
    for run=1:floor(size(PS.AllTrIntegProbs{s},1)/15)
        TrIndRun = (run-1)*15+1:run*15;
        if(TrIndRun(end) < size(PS.AllTrIntegProbs{s},1))
            GAHitRateSim055Run(s,run) = 100*sum(max(PS.AllTrIntegProbs{s}(TrIndRun,:)') > Decth)/length(TrIndRun);
        else
            break;
        end
    end
end


figure(134)
shadedErrorBar(1:10,nanmean(GAHitRateSim055(:,1:10)),nanstd(GAHitRateSim055(:,1:10)),'b');
hold on; plot(nanmean(GAHitRateSim055(:,1:10)),'-b','LineWidth',3);
P = polyfit([1:10],nanmean(squeeze(GAHitRateSim055(:,1:10))),1);
f=@(x)(P(1).*x+P(2));
plot([1:10],f([1:10]),'-r','LineWidth',3);
hold off;
axis([1 10 40 110]);
xlabel('Session index (chronological)','FontSize',20);
ylabel('Simulated Hit Rate (%)','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);


figure(135)
shadedErrorBar(1:60,nanmean(GAHitRateSim055Run(:,1:60)),nanstd(GAHitRateSim055Run(:,1:60)),'b');
hold on; plot(nanmean(GAHitRateSim055Run(:,1:60)),'-b','LineWidth',3);
P = polyfit([1:60],nanmean(squeeze(GAHitRateSim055Run(:,1:60))),1);
f=@(x)(P(1).*x+P(2));
plot([1:60],f([1:60]),'-r','LineWidth',3);
hold off;
axis([1 60 40 120]);
xlabel('Run index (chronological)','FontSize',20);
ylabel('Simulated Hit Rate (%)','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);


figure(136);
shadedErrorBar(1:10,nanmean(GAThSes(:,1:10)),nanstd(GAThSes(:,1:10)),'b');
hold on; plot(nanmean(GAThSes(:,1:10)),'-b','LineWidth',3);
P = polyfit([1:10],nanmean(GAThSes(:,1:10)),1);
f=@(x)(P(1).*x+P(2));
plot([1:10],f([1:10]),'-r','LineWidth',3);
hold off;
axis([1 10 0.55 0.9]);
xlabel('Session index (chronological)','FontSize',20);
ylabel('Average Confidence Threshold value','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);


% Effort to explain exceptions with features selected
figure(11);imagesc(squeeze(mean(PS.SelFeatMirrorMap(IndBCI,:,:))),[0 0.5])
figure(12);imagesc(squeeze(mean(PS.SelFeatMirrorMap(FMAdiff(IndBCI)==0,:,:))),[0 0.5])
figure(13);imagesc(squeeze(mean(PS.SelFeatMirrorMap(FMAdiff(IndBCI)>0,:,:))),[0 0.5])
% 
% % Artifact
% figure(5);
% Zth = [0:0.5:3.5];
% bar([mean(PS.ArtPrctC)' mean(PS.ArtPrctN)']);
% 
% legend({'Motor Attempt', 'Run'},'FontSize',20);
% set(gca,'XTick',[1:length(Zth)]);
% set(gca,'XTickLabel',Zth,'FontSize',20);
% axis([0 length(Zth)+1 0 80]);
% xlabel('z-score threshold','FontSize',20);
% ylabel('Artifact impact (%)','FontSize',20);
% set(gca,'FontSize',20);
% set(gca,'LineWidth',3);

% Grand average of trials
GAIProb = [];
for s=1:length(SubID)
    GAIProb(:,s) = nanmean(PS.AllTrIntegProbs{s});
end
FirstZero = find(sum(GAIProb')==0);
FirstZero = FirstZero(1);
GAIProb = GAIProb(1:FirstZero-1,:);
t=[1/16:1/16:(FirstZero-1)/16];

figure(6);
plot(t,GAIProb(:,IndBCI),'-','LineWidth',3);
hold on;
plot(t,GAIProb(:,IndSham),'--','LineWidth',3);
hold off;
legend(SubID,'Orientation','Vertical','FontSize',10);
%axis([0 FirstZero/16 0 1.0]);
axis([0 7 0 1.0]);
ylabel('Integrated Probability','FontSize',20);
xlabel('Time (s)','FontSize',20);
title('Grand Average of trial probability','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);


% Simulated hit rates per threshold
SimHitRate = [];
Decth = [0.55:0.05:0.95];
for s=1:length(SubID)
    MaxIProb = max(PS.AllTrIntegProbs{s}');
    for th=1:length(Decth)
        SimHitRate(s,th) =  100*sum(MaxIProb > Decth(th))/length(MaxIProb);
        tmphit = MaxIProb > Decth(th);
        [CorrSimHitRealHit(s,th) PvalCorrSimHitRealHit(s,th)] = corr(tmphit',PS.AllFESHitMiss{s}');
    end
end
clc;
for th=1:length(Decth)
     disp(['Threshold ' num2str(Decth(th)) ' : BCI-> ' num2str(mean(SimHitRate(IndBCI,th))) ' +/-' num2str(std2(SimHitRate(IndBCI,th)))...
    ' , Sham-> ' num2str(mean(SimHitRate(IndSham,th))) ' +/-' num2str(std2(SimHitRate(IndSham,th))) ' , p=' num2str(ranksum(SimHitRate(IndSham,th),SimHitRate(IndBCI,th)))])
end
for th=1:length(Decth)
     disp(['Correlation Simulated vs Real, Threshold ' num2str(Decth(th)) ' : BCI-> ' num2str(mean(CorrSimHitRealHit(IndBCI,th))) ' +/-' num2str(std2(CorrSimHitRealHit(IndBCI,th)))...
    ' , Sham-> ' num2str(mean(CorrSimHitRealHit(IndSham,th))) ' +/-' num2str(std2(CorrSimHitRealHit(IndSham,th))) ' , p=' num2str(ranksum(CorrSimHitRealHit(IndSham,th),CorrSimHitRealHit(IndBCI,th)))])
end

close all;
figure(432);
SimHitRate055 = SimHitRate(:,1);
bar(1,mean(PS.DetectionRate(IndBCI)),'FaceColor',Color{1});hold on;
bar(2,mean(PS.DetectionRate(IndSham)),'FaceColor',Color{2});
errorbar(1, mean(PS.DetectionRate(IndBCI)), std2(PS.DetectionRate(IndBCI)),'.k','LineWidth',3);
errorbar(2, mean(PS.DetectionRate(IndSham)), std2(PS.DetectionRate(IndSham)),'.k','LineWidth',3);
bar(4,mean(SimHitRate055(IndBCI)),'FaceColor',Color{1});
bar(5,mean(SimHitRate055(IndSham)),'FaceColor',Color{2});
errorbar(4, mean(SimHitRate055(IndBCI)), std2(SimHitRate055(IndBCI)),'.k','LineWidth',3);
errorbar(5, mean(SimHitRate055(IndSham)), std2(SimHitRate055(IndSham)),'.k','LineWidth',3);
%bar(7,mean(PS.HitRate(IndBCI)),'FaceColor',Color{1});
%errorbar(7, mean(PS.HitRate(IndBCI)), std2(PS.HitRate(IndBCI)),'.k','LineWidth',3);
%bar(8,mean(PS.HitRate(IndSham)),'FaceColor',Color{2});
%errorbar(8, mean(PS.HitRate(IndSham)), std2(PS.HitRate(IndSham)),'.k','LineWidth',3);
[p1 h] = ranksum(PS.DetectionRate(IndBCI),PS.DetectionRate(IndSham));
[p2 h] = ranksum(SimHitRate055(IndBCI),SimHitRate055(IndSham));
axis([0 6 0 120]);
%sigstar({[1,2] , [4,5]},[p1 p2]);
hold off;
legend({'BCI','Sham'},'FontSize',20);
set(gca,'XTick',[1.5 4.5]);
set(gca,'XTickLabel',{'Detection Rate' , 'Simulated Hit Rate'},'FontSize',20);
xlabel('BCI Accuracy Type','FontSize',20);
ylabel('BCI Accuracy (%)','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);


% Correlation of accuracy with discriminncy, in order to claim no
% influence of artifacts (Comment 1.8)
for s=1:length(SubID)
    MeanFSSel = [];
    MaxFSSel = [];
    DR = [];
    for r=1:length(PS.DetectionRateAll{s})
        MeanFSSel(r) = mean(PS.FSSel{s}{r});
        MaxFSSel(r) = max(PS.FSSel{s}{r});
        DR(r) = PS.DetectionRateAll{s}(r);
    end     
    % This is for more detailed (per subject) correlations, still make the 
    % point but need to report too mamy numbers + explain 3 excpetions...
    [MaxR(s) MaxP(s)] = corr(MaxFSSel',DR'); %%
    [MeanR(s) MeanP(s)] = corr(MeanFSSel',DR');
    
    % Better report per population correlations (1 value per subject)
    SMeanFSSel(s) = mean(MeanFSSel);
    SMaxFSSel(s) = mean(MaxFSSel);
    
    % Callculate for Fz
    MeanFSFz(s) = mean(PS.FSavg{s}(:,1));
    MaxFSFz(s) = max(PS.FSavg{s}(:,1));
end

% Final corr to report
%[r p] = corr(SMaxFSSel',PS.DetectionRate');
[r p] = corr(SMeanFSSel',PS.DetectionRate');
[r p] = corr(MeanFSFz',PS.DetectionRate');
%[r p] = corr(MaxFSFz',PS.DetectionRate');


%% Figures for comment 1.13, on last-sec PSD
% % THis is the normalized per column version
% figure(100);
% Metric = PS.LastPSDConfMatInv(:,1,1)+PS.LastPSDConfMatInv(:,2,2)-PS.LastPSDConfMatInv(:,1,2)-PS.LastPSDConfMatInv(:,2,1);
% bar(1,mean(Metric(IndBCI)),'FaceColor',Color{1});hold on;
% bar(2,mean(Metric(IndSham)),'FaceColor',Color{2});
% errorbar([1 2], [mean(Metric(IndBCI)) mean(Metric(IndSham))], [std2(Metric(IndBCI)) std2(Metric(IndSham))],'.k','LineWidth',3);
% [h p] = ttest2(Metric(IndBCI),Metric(IndSham));
% axis([0 3 -10 150]);
% sigstar({[1,2]},[p])
% hold off;
% legend({'BCI','Sham'},'FontSize',20);
% set(gca,'XTick',[1:2]);
% set(gca,'XTickLabel',{'BCI' , 'Sham'},'FontSize',20);
% xlabel('Treatment','FontSize',20);
% ylabel('TPR+TNR-FPR-FNR (%)','FontSize',20);
% set(gca,'FontSize',20);
% set(gca,'LineWidth',3);

%% Recalculate everything according to proper definitions
RemInd = find(isnan(FMAdiff));
FMAdiffNoNan = FMAdiff;
FMAdiffNoNan(RemInd)=[];
SubIDNoNaN = SubID;
SubIDNoNaN(RemInd)=[];
SideNoNaN = Side;
SideNoNaN(RemInd)=[];
GroupNoNaN = Group;
GroupNoNaN(RemInd) = [];
IndBCINoNan = find(strcmp(GroupNoNaN,'BCI'));
IndShamNoNan = find(strcmp(GroupNoNaN,'Sham'));

TPR = PS.LastPSDConfMatAll(:,1,1)./(PS.LastPSDConfMatAll(:,1,1) + PS.LastPSDConfMatAll(:,2,1));
TPR(RemInd)=[];
[r p] = corr(TPR,FMAdiffNoNan');
disp(['TPR, r = ' num2str(r) ' , p = ' num2str(p)]);
TNR = PS.LastPSDConfMatAll(:,2,2)./(PS.LastPSDConfMatAll(:,2,2) + PS.LastPSDConfMatAll(:,1,2));
TNR(RemInd)=[];
[r p] = corr(TNR,FMAdiffNoNan');
disp(['TNR, r = ' num2str(r) ' , p = ' num2str(p)]);
PPV = PS.LastPSDConfMatAll(:,1,1)./(PS.LastPSDConfMatAll(:,1,1) + PS.LastPSDConfMatAll(:,1,2));
PPV(RemInd)=[];
[r p] = corr(PPV,FMAdiffNoNan');
disp(['PPV, r = ' num2str(r) ' , p = ' num2str(p)]);
NPV = PS.LastPSDConfMatAll(:,2,2)./(PS.LastPSDConfMatAll(:,2,1) + PS.LastPSDConfMatAll(:,2,2));
NPV(RemInd)=[];
[r p] = corr(NPV,FMAdiffNoNan');
disp(['NPV, r = ' num2str(r) ' , p = ' num2str(p)]);
ACC = (PS.LastPSDConfMatAll(:,1,1)+PS.LastPSDConfMatAll(:,2,2))./(PS.LastPSDConfMatAll(:,1,1) + PS.LastPSDConfMatAll(:,1,2)+PS.LastPSDConfMatAll(:,2,1) + PS.LastPSDConfMatAll(:,2,2));
ACC(RemInd)=[];
[r p] = corr(ACC,FMAdiffNoNan');
disp(['ACC, r = ' num2str(r) ' , p = ' num2str(p)]);
F1 = 2./(1./TPR + 1./PPV);
[r p] = corr(F1,FMAdiffNoNan');
disp(['F1, r = ' num2str(r) ' , p = ' num2str(p)]);


TP = PS.LastPSDConfMatAll(:,1,1);
TP(RemInd)=[];
[r p] = corr(TP,FMAdiffNoNan');
disp(['TP, r = ' num2str(r) ' , p = ' num2str(p)]);
FP = PS.LastPSDConfMatAll(:,1,2);
FP(RemInd)=[];
[r p] = corr(FP,FMAdiffNoNan');
disp(['FP, r = ' num2str(r) ' , p = ' num2str(p)]);
FN = PS.LastPSDConfMatAll(:,2,1);
FN(RemInd)=[];
[r p] = corr(FN,FMAdiffNoNan');
disp(['FN, r = ' num2str(r) ' , p = ' num2str(p)]);
TN = PS.LastPSDConfMatAll(:,2,2);
TN(RemInd)=[];
[r p] = corr(TN,FMAdiffNoNan');
disp(['TN, r = ' num2str(r) ' , p = ' num2str(p)]);


[r p] = corr(PPV+ACC,FMAdiffNoNan');
disp(['PPV+ACC, r = ' num2str(r) ' , p = ' num2str(p)]);

% 
% close all;
% Metric = 100*(PPV+NPV)';
% 
% figure(86);
% bar(1,mean(Metric(IndBCI)),'FaceColor',Color{1});hold on;
% bar(2,mean(Metric(IndSham)),'FaceColor',Color{2});
% errorbar([1 2], [mean(Metric(IndBCI)) mean(Metric(IndSham))], [std2(Metric(IndBCI)) std2(Metric(IndSham))],'.k','LineWidth',3);
% [h p] = ttest2(Metric(IndBCI),Metric(IndSham));
% axis([0 3 95 200]);
% sigstar({[1,2]},[p]);
% hold off;
% legend({'BCI','Sham'},'FontSize',30);
% set(gca,'XTick',[1:2]);
% set(gca,'XTickLabel',{'BCI' , 'Sham'},'FontSize',30);
% xlabel('Treatment','FontSize',30);
% ylabel('PPV+NPV (%)','FontSize',30);
% set(gca,'FontSize',30);
% set(gca,'LineWidth',3);
% 
% figure(87);
% hold on;
% h = plot(Metric,FMAdiff,'.','Color','k','MarkerSize',40);
% h1 = plot(Metric(IndBCI),FMAdiff(IndBCI)','.','Color',Color{1},'MarkerSize',40);
% h2 = plot(Metric(IndSham),FMAdiff(IndSham)','.','Color',Color{2},'MarkerSize',40);
% hls = lsline;
% set(hls,'LineWidth',4);
% hold off;
% axis([90 180 -3 20]);
% legend([h h1 h2],{'All','BCI','Sham'},'Location','NorthWest','FontSize',40);
% set(gca,'FontSize',40);
% xlabel('PPV+NPV (%)','FontSize',40);
% ylabel('FMA POST-PRE','FontSize',40);
% [r p] = corr(Metric',FMAdiff');
% 
% 
% Metric = 100*ACC';
% figure(88);
% bar(1,mean(Metric(IndBCI)),'FaceColor',Color{1});hold on;
% bar(2,mean(Metric(IndSham)),'FaceColor',Color{2});
% errorbar([1 2], [mean(Metric(IndBCI)) mean(Metric(IndSham))], [std2(Metric(IndBCI)) std2(Metric(IndSham))],'.k','LineWidth',3);
% [h p] = ttest2(Metric(IndBCI),Metric(IndSham));
% axis([0 3 45 100]);
% sigstar({[1,2]},[p]);
% hold off;
% legend({'BCI','Sham'},'FontSize',40);
% set(gca,'XTick',[1:2]);
% set(gca,'XTickLabel',{'BCI' , 'Sham'},'FontSize',40);
% xlabel('Treatment','FontSize',40);
% ylabel('Accuracy (%)','FontSize',40);
% set(gca,'FontSize',40);
% set(gca,'LineWidth',3);
% 
% 
% figure(89);
% hold on;
% h = plot(Metric,FMAdiff,'.','Color','k','MarkerSize',50);
% h1 = plot(Metric(IndBCI),FMAdiff(IndBCI)','.','Color',Color{1},'MarkerSize',50);
% h2 = plot(Metric(IndSham),FMAdiff(IndSham)','.','Color',Color{2},'MarkerSize',50);
% %hls = lsline;
% %set(hls,'LineWidth',4);
% Pbci = polyfit(Metric(IndBCI),FMAdiff(IndBCI),1);
% Psham = polyfit(Metric(IndSham),FMAdiff(IndSham),1);
% Pall = polyfit(Metric,FMAdiff,1);
% fbci=@(x)(Pbci(1).*x+Pbci(2));
% fsham=@(x)(Psham(1).*x+Psham(2));
% fall=@(x)(Pall(1).*x+Pall(2));
% plot(Metric(IndBCI),fbci(Metric(IndBCI)),'-','LineWidth',3,'Color',Color{1});
% plot(Metric(IndSham),fsham(Metric(IndSham)),'-','LineWidth',3,'Color',Color{2});
% plot(Metric,fall(Metric),'-k','LineWidth',3);
% hold off;
% %axis([90 180 -3 20]);
% %legend([h h1 h2],{'All','BCI','Sham'},'Location','NorthWest','FontSize',40);
% set(gca,'FontSize',40);
% xlabel('Accuracy (%)','FontSize',40);
% ylabel('\DeltaFMA (score)','FontSize',40);
% [r p] = corr(Metric',FMAdiff');
% 

% close all;
% figure(1000);
% bar(IndBCI,PS.LastTPR(IndBCI),'FaceColor','b');hold on;
% bar(IndSham,PS.LastTPR(IndSham),'FaceColor','r');
% bar(length(SubID)+1, mean(PS.LastTPR(IndBCI)),'FaceColor','b');
% bar(length(SubID)+2, mean(PS.LastTPR(IndSham)),'FaceColor','r');
% set(gca,'XTick',[1:length(SubID)+2]);
% set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);
% axis([0 length(SubID)+3 0 100]);
% xticklabel_rotate([],45,[]);
% legend({'BCI','Sham'});
% ylabel('True positive rate (%)');
% 
% 
% figure(1001);
% bar(IndBCI,PS.LastFPR(IndBCI),'FaceColor','b');hold on;
% bar(IndSham,PS.LastFPR(IndSham),'FaceColor','r');
% bar(length(SubID)+1, mean(PS.LastFPR(IndBCI)),'FaceColor','b');
% bar(length(SubID)+2, mean(PS.LastFPR(IndSham)),'FaceColor','r');
% set(gca,'XTick',[1:length(SubID)+2]);
% set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);
% axis([0 length(SubID)+3 0 100]);
% xticklabel_rotate([],45,[]);
% legend({'BCI','Sham'});
% ylabel('False positive rate (%)');
% 
% 
% figure(1002);
% plot(PS.LastTPR-PS.LastFPR,FMAdiff','*k',PS.LastTPR(IndBCI)-PS.LastFPR(IndBCI),FMAdiff(IndBCI)','*b',PS.LastTPR(IndSham)-PS.LastFPR(IndSham),FMAdiff(IndSham)','*r');
% lsline;
% axis([-10 80 -5 20]);
% legend({'All','BCI','Sham'});
% xlabel('TPR - FPR (%)');
% ylabel('\DeltaFMA Post-Pre');
% [r p] = corr((PS.LastTPR-PS.LastFPR)',FMAdiff');
% text(0,17,['r = ' num2str(r) ', p = ' num2str(p)]);
% 
% 
% 
% figure(1003);
% plot(PS.LastTPR./(PS.LastTPR + (100 - PS.LastFPR)),FMAdiff','*k');
% lsline;
% %axis([-10 80 -5 20]);
% legend({'All','BCI','Sham'});
% xlabel('TPR/(TPR+FPR) ');
% ylabel('\DeltaFMA Post-Pre');
% [r p] = corr((PS.LastTPR./(PS.LastTPR + (100 - PS.LastFPR)))',FMAdiff');
% text(0,17,['r = ' num2str(r) ', p = ' num2str(p)]);


%% Search for metrics that correlate well wit FMA diff
%findCorrWithFMAdiff(PS.ERSPLastHit,FMAdiff,SubID, Side, IndBCI, IndSham);
findCorrWithFMAdiff(cellfun(@minus,PS.ERSPLastHit,PS.ERSPLastMiss,'UniformOutput',false),FMAdiff,SubID, Side, IndBCI, IndSham);
%findCorrWithFMAdiff(PS.FSLastHit,FMAdiff,SubID, Side, IndBCI, IndSham);
%findCorrWithFMAdiff(cellfun(@minus,PS.FSLastHit,PS.FSLastMiss,'UniformOutput',false),FMAdiff,SubID, Side, IndBCI, IndSham);
%findCorrWithFMAdiff(PS.FSLastMiss,FMAdiff,SubID, Side, IndBCI, IndSham);
%findCorrWithFMAdiff(PS.ERSPLastMiss,FMAdiff,SubID, Side, IndBCI, IndSham);
%findCorrWithFMAdiff(PS.FSavg,FMAdiff,SubID, Side, IndBCI, IndSham);
%findCorrWithFMAdiff(PS.ERSPavg,FMAdiff,SubID, Side, IndBCI, IndSham);


%% Figures

%% Detection Rate plots
%% Supplementary figure for Nat Comm revision, comment 1.13
% Average detection rate per subject

close all;clc;
%% Group colors to look like the existing figures
h=figure(1);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.DetectionRate(sub),'FaceColor',Color{1},'LineWidth',3);
    else
        bar(sub,PS.DetectionRate(sub),'FaceColor',Color{2},'LineWidth',3);
    end
end
hbci = bar(length(SubID)+1, mean(PS.DetectionRate(find(strcmp(Group,'BCI')))),'FaceColor',Color{1},'LineWidth',3);
errorbar(length(SubID)+1, mean(PS.DetectionRate(find(strcmp(Group,'BCI')))),std2(PS.DetectionRate(find(strcmp(Group,'BCI')))),'Color','k','LineWidth',3);
hsham = bar(length(SubID)+2, mean(PS.DetectionRate(find(strcmp(Group,'Sham')))),'FaceColor',Color{2},'LineWidth',3);
errorbar(length(SubID)+2, mean(PS.DetectionRate(find(strcmp(Group,'Sham')))),std2(PS.DetectionRate(find(strcmp(Group,'Sham')))),'Color','k','LineWidth',3);
axis([0.5 length(SubID)+2.5 0 101]);
%line([0.5 length(SubID)+2.5],[50 50],'Color','k','LineWidth',3,'LineStyle','--');
%line([0.5 length(SubID)+2.5],[58 58],'Color','k','LineWidth',3,'LineStyle','--');
line([length(SubID)+0.5 length(SubID)+0.5],[0 101],'Color','k','LineWidth',3,'LineStyle','--');
hold off;
legend([hbci(1) hsham(1)],{'BCI','Sham'},'FontSize',20);
xlabel('Participants');
ylabel('Single-sample online detection rate (%)');
set(gca,'LineWidth',3,'FontSize',20);
set(gca,'XTick',[1:length(SubID)+2]);
%set(gca,'XTickLabel',[SubPaperID 'BCI' 'Sham']);
set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);
% Maximize the figure before running this
xticklabel_rotate([],45,[],'Fontsize',20);


%% Hit rate plot
h=figure(2);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.HitRate(sub),'FaceColor',Color{1},'LineWidth',3);
    else
        bar(sub,PS.HitRate(sub),'FaceColor',Color{2},'LineWidth',3);
    end
end
hbci = bar(length(SubID)+1, mean(PS.HitRate(find(strcmp(Group,'BCI')))),'FaceColor',Color{1},'LineWidth',3);
errorbar(length(SubID)+1, mean(PS.HitRate(find(strcmp(Group,'BCI')))),std2(PS.HitRate(find(strcmp(Group,'BCI')))),'Color','k','LineWidth',3);
hsham = bar(length(SubID)+2, mean(PS.HitRate(find(strcmp(Group,'Sham')))),'FaceColor',Color{2},'LineWidth',3);
errorbar(length(SubID)+2, mean(PS.HitRate(find(strcmp(Group,'Sham')))),std2(PS.HitRate(find(strcmp(Group,'Sham')))),'Color','k','LineWidth',3);
axis([0.5 length(SubID)+2.5 0 101]);
%line([0.5 length(SubID)+2.5],[50 50],'Color','k','LineWidth',3,'LineStyle','--');
%line([0.5 length(SubID)+2.5],[58 58],'Color','k','LineWidth',3,'LineStyle','--');
line([length(SubID)+0.5 length(SubID)+0.5],[40 101],'Color','k','LineWidth',3,'LineStyle','--');
hold off;
legend([hbci(1) hsham(1)],{'BCI','Sham'},'FontSize',20);
xlabel('Participants');
ylabel('Hit rate (%)');
set(gca,'LineWidth',3,'FontSize',20);
set(gca,'XTick',[1:length(SubID)+2]);
%set(gca,'XTickLabel',[SubPaperID 'BCI' 'Sham']);
set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);
% Maximize the figure before running this
xticklabel_rotate([],45,[],'Fontsize',20);

%% Simulated Hit rate plot for decision threshold 0.55
SimHitRate055 = SimHitRate(:,1);
h=figure(234);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,SimHitRate055(sub),'FaceColor',Color{1},'LineWidth',3);
    else
        bar(sub,SimHitRate055(sub),'FaceColor',Color{2},'LineWidth',3);
    end
end
hbci = bar(length(SubID)+1, mean(SimHitRate055(find(strcmp(Group,'BCI')))),'FaceColor',Color{1},'LineWidth',3);
errorbar(length(SubID)+1, mean(SimHitRate055(find(strcmp(Group,'BCI')))),std2(SimHitRate055(find(strcmp(Group,'BCI')))),'Color',Color{1},'LineWidth',3);
hsham = bar(length(SubID)+2, mean(SimHitRate055(find(strcmp(Group,'Sham')))),'FaceColor',Color{2},'LineWidth',3);
errorbar(length(SubID)+2, mean(SimHitRate055(find(strcmp(Group,'Sham')))),std2(SimHitRate055(find(strcmp(Group,'Sham')))),'Color',Color{2},'LineWidth',3);
axis([0.5 length(SubID)+2.5 0 120]);
%line([0.5 length(SubID)+2.5],[50 50],'Color','k','LineWidth',3,'LineStyle','--');
%line([0.5 length(SubID)+2.5],[58 58],'Color','k','LineWidth',3,'LineStyle','--');
line([length(SubID)+0.5 length(SubID)+0.5],[40 120],'Color','k','LineWidth',3,'LineStyle','--');
hold off;
legend([hbci(1) hsham(1)],{'BCI','Sham'},'FontSize',20);
xlabel('Participants');
ylabel('Simulated Hit rate (%)');
set(gca,'LineWidth',3,'FontSize',20);
set(gca,'XTick',[1:length(SubID)+2]);
set(gca,'XTickLabel',[SubPaperID 'BCI' 'Sham']);
% Maximize the figure before running this
xticklabel_rotate([],45,[],'Fontsize',20);



%% Simulated Accuracy plots
% Average simulated accuracy (selected) per subject
h=figure(3);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.SimAccSel(sub),'FaceColor',Color{1},'LineWidth',3);
    else
        bar(sub,PS.SimAccSel(sub),'FaceColor',Color{2},'LineWidth',3);
    end
end
bar(length(SubID)+1, mean(PS.SimAccSel(find(strcmp(Group,'BCI')))),'FaceColor',Color{1},'LineWidth',3);
errorbar(length(SubID)+1, mean(PS.SimAccSel(find(strcmp(Group,'BCI')))),std2(PS.SimAccSel(find(strcmp(Group,'BCI')))),'k','LineWidth',3);
bar(length(SubID)+2, mean(PS.SimAccSel(find(strcmp(Group,'Sham')))),'FaceColor',Color{2},'LineWidth',3);
errorbar(length(SubID)+2, mean(PS.SimAccSel(find(strcmp(Group,'Sham')))),std2(PS.SimAccSel(find(strcmp(Group,'Sham')))),'k','LineWidth',3);
axis([0.5 length(SubID)+2.5 40 101]);
line([0.5 length(SubID)+2.5],[50 50],'Color','k','LineWidth',3,'LineStyle','--');
line([0.5 length(SubID)+2.5],[58 58],'Color','k','LineWidth',3,'LineStyle','--');
line([length(SubID)+0.5 length(SubID)+0.5],[40 101],'Color','k','LineWidth',3,'LineStyle','--');
hold off;
legend({'BCI','Sham'},'FontSize',20);
xlabel('Participants');
ylabel('Simulated detection rate (selected) (%)');
set(gca,'LineWidth',3,'FontSize',20);
set(gca,'XTick',[1:length(SubID)+2]);
set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);
% Maximize the figure before running this
xticklabel_rotate([],45,[],'Fontsize',20);

% Average simulated accuracy (best) per subject
h=figure(4);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.SimAccBest(sub),'FaceColor',Color{1},'LineWidth',3);
    else
        bar(sub,PS.SimAccBest(sub),'FaceColor',Color{2},'LineWidth',3);
    end
end
bar(length(SubID)+1, mean(PS.SimAccBest(find(strcmp(Group,'BCI')))),'FaceColor',Color{1},'LineWidth',3);
errorbar(length(SubID)+1, mean(PS.SimAccBest(find(strcmp(Group,'BCI')))),std2(PS.SimAccBest(find(strcmp(Group,'BCI')))),'k','LineWidth',3);
bar(length(SubID)+2, mean(PS.SimAccBest(find(strcmp(Group,'Sham')))),'FaceColor',Color{2},'LineWidth',3);
errorbar(length(SubID)+2, mean(PS.SimAccBest(find(strcmp(Group,'Sham')))),std2(PS.SimAccBest(find(strcmp(Group,'Sham')))),'k','LineWidth',3);
axis([0.5 length(SubID)+2.5 40 101]);
line([0.5 length(SubID)+2.5],[50 50],'Color','k','LineWidth',3,'LineStyle','--');
line([0.5 length(SubID)+2.5],[58 58],'Color','k','LineWidth',3,'LineStyle','--');
line([length(SubID)+0.5 length(SubID)+0.5],[40 101],'Color','k','LineWidth',3,'LineStyle','--');
hold off;
legend({'BCI','Sham'},'FontSize',20);
xlabel('Participants');
ylabel('Simulated detection rate (best) (%)');
set(gca,'LineWidth',3,'FontSize',20);
set(gca,'XTick',[1:length(SubID)+2]);
set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);
% Maximize the figure before running this
xticklabel_rotate([],45,[],'Fontsize',20);


% 
% %% ERS/ERS topoplots
% % Assume always left-side affected (right hemisphere is the affected)
% for sub=1:length(SubID)
%     AvgERDS{sub} = PS.logERDSdprSes{sub}{1}; 
%     numses = 1;
%     for r=2:length(PS.logERDSdprSes{sub})
%         if(sum(sum(isnan(PS.logERDSdprSes{sub}{r}))) == 0)
%             numses = numses + 1;
%             AvgERDS{sub} = AvgERDS{sub} + PS.logERDSdprSes{sub}{r};
%         end
%     end
%     AvgERDS{sub} = AvgERDS{sub}/numses;
%     if(strcmp(Side{sub},'R'))
%         doflip = 1;
%     else
%         doflip = 0;
%     end
%     TopoMu{sub} = convChans(flipSide16(mean(AvgERDS{sub}(:,3:6),2),doflip));
%     TopoBeta{sub} = convChans(flipSide16(mean(AvgERDS{sub}(:,8:11),2),doflip));
%     MinERDSMu(sub) = min(TopoMu{sub});
%     MaxERDSMu(sub) = max(TopoMu{sub});
%     MinERDSBeta(sub) = min(TopoBeta{sub});
%     MaxERDSBeta(sub) = max(TopoBeta{sub});
%     LimMu(sub) = max(abs(MaxERDSMu(sub)),abs(MinERDSMu(sub)));
%     LimBeta(sub) = max(abs(MaxERDSBeta(sub)),abs(MinERDSBeta(sub)));
% end
% 
% 
% %% ERSP topoplots
% % Assume always left-side affected (right hemisphere is the affected)
% for sub=1:length(SubID)
%     AvgERSP{sub} = PS.ERSPSes{sub}{1}; 
%     numses = 1;
%     for r=2:length(PS.ERSPSes{sub})
%         if(sum(sum(isnan(PS.ERSPSes{sub}{r}))) == 0)
%             numses = numses + 1;
%             AvgERSP{sub} = AvgERSP{sub} + PS.ERSPSes{sub}{r};
%         end
%     end
%     AvgERSP{sub} = AvgERSP{sub}/numses;
%     if(strcmp(Side{sub},'R'))
%         doflip = 1;
%     else
%         doflip = 0;
%     end
%     TopoMuERSP{sub} = convChans(flipSide16(mean(AvgERSP{sub}(:,3:6),2),doflip));
%     TopoBetaERSP{sub} = convChans(flipSide16(mean(AvgERSP{sub}(:,8:11),2),doflip));
%     MinERSPMu(sub) = min(TopoMuERSP{sub});
%     MaxERSPMu(sub) = max(TopoMuERSP{sub});
%     MinERSPBeta(sub) = min(TopoBetaERSP{sub});
%     MaxERSPBeta(sub) = max(TopoBetaERSP{sub});
%     LimERSPMu(sub) = max(abs(MaxERSPMu(sub)),abs(MinERSPMu(sub)));
%     LimERSPBeta(sub) = max(abs(MaxERSPBeta(sub)),abs(MinERSPBeta(sub)));
% end
% 
%% Fisher Score
for sub=1:length(SubID)
    MinFS(sub) = min(PS.FSavg{sub}(:));
    MaxFS(sub) = max(PS.FSavg{sub}(:));
    if(strcmp(Side{sub},'L'))
        doflip = 1;
    else
        doflip = 0;
    end
    %TopoFSMu{sub} = convChans(flipSide16(mean(PS.FSavg{sub}(:, 3:6),2),doflip));
    %TopoFSBeta{sub} = convChans(flipSide16(mean(PS.FSavg{sub}(:,8:11),2),doflip));
    TopoFSMu{sub} = convChans(flipSide16(max(PS.FSavg{sub}(:,3:6)'),doflip));
    TopoFSBeta{sub} = convChans(flipSide16(max(PS.FSavg{sub}(:,8:11)'),doflip));
end
%MinFS = min(MinFS);
%MaxFS = max(MaxFS);

% Load colormaps
cmerds = load('erdersmap.mat');cmerds = cmerds.cm;
cmfs = load('fsmap.mat');cmfs = cmfs.cm;
% Load chanlocs 64
chanlocs = load('chanlocs.mat');chanlocs=chanlocs.chanlocs;

% 
% close all;
% % ERDS
% for sub=1:length(SubID)
%     figure(5);
%     topoplot(TopoMu{sub},chanlocs,'maplimits',[-LimMu(sub) LimMu(sub)]);colormap(cmerds);drawnow;
%     %print([FigureSavePath '/erds/' SubID{sub} '_Mu'],'-r600','-dpng');
%     figure(6);
%     topoplot(TopoBeta{sub},chanlocs,'maplimits',[-LimBeta(sub) LimBeta(sub)]);colormap(cmerds);drawnow;
%     %print([FigureSavePath '/erds/' SubID{sub} '_Beta'],'-r600','-dpng');
% end
% 
% % ERSP
% for sub=1:length(SubID)
%     figure(7);
%     %topoplot(TopoMuERSP{sub},chanlocs,'maplimits',[-LimMu(sub) LimMu(sub)]);colormap(cmerds);drawnow;
%     topoplot(TopoMuERSP{sub},chanlocs);colormap(cmerds);drawnow;
%     %print([FigureSavePath '/erds/' SubID{sub} '_Mu'],'-r600','-dpng');
%     figure(8);
%     %topoplot(TopoBetaERSP{sub},chanlocs,'maplimits',[-LimBeta(sub) LimBeta(sub)]);colormap(cmerds);drawnow;
%     topoplot(TopoBetaERSP{sub},chanlocs);colormap(cmerds);drawnow;
%     %print([FigureSavePath '/erds/' SubID{sub} '_Beta'],'-r600','-dpng');
% end
% 
% FS
close all;
for sub=1:length(IndBCI)
    figure(sub);
    subplot(1,2,1);topoplot(TopoFSMu{sub},chanlocs,'maplimits',[0 max(TopoFSMu{sub})]);%colorbar;
    subplot(1,2,2);topoplot(TopoFSBeta{sub},chanlocs,'maplimits',[0 max(TopoFSBeta{sub})]);%colorbar;
end

% 
% %% Per session DP maps for the QR Appendix
% MaxAbsERSP = 0;
% for sub=1:length(SubID)
%     figure(9);
%     %topoplot(TopoFSMu{sub},chanlocs,'maplimits',[MinFS(sub) MaxFS(sub)]);colormap(cmfs);drawnow;
%     topoplot(TopoFSMu{sub},chanlocs);
%     %print([FigureSavePath '/fs/' SubID{sub} '_Mu'],'-r600','-dpng');
%     figure(10);
%     topoplot(TopoFSBeta{sub},chanlocs,'maplimits',[MinFS(sub) MaxFS(sub)]);colormap(cmfs);drawnow;
%     %print([FigureSavePath '/fs/' SubID{sub} '_Beta'],'-r600','-dpng');
% end


close all;
for sub=1:length(SubID)
    figure(sub);
    Nrows = ceil(PS.SesNum(sub)/5);
    if(strcmp(Side{sub},'R'))
        doflip = 1;
    else
        doflip = 0;
    end
    for ses=1:PS.SesNum(sub)
        subplot(Nrows, 5, ses);imagesc(flipSide16Matrix(PS.ERSPSes{sub}{ses},doflip),[-1.0 1.0]);
        if(ses==1)
            colorbar;
        end
        title(['Session '  num2str(ses)]);
        xlabel('Bands (Hz)');
        ylabel('Channels');
        set(gca,'LineWidth',3,'FontSize',10);
        set(gca,'XTick',[1:4:23]);
        set(gca,'XTickLabel',[4:8:48]);
        set(gca,'YTick',[1:16]);
        set(gca,'YTickLabel',Electrodes16);
    end
end




% 
% close all;
% %% Get some ERSP statistics using the first and last sessions
% % Use only first and last session
% PreERSP = [];
% PostERSP = [];
% for sub=1:length(SubID)
%     if(strcmp(Side{sub},'R'))
%         doflip = 1;
%     else
%         doflip = 0;
%     end
%     PreERSP(sub,:,:) = flipSide16Matrix(PS.ERSPSes{sub}{1},doflip);
%     PostERSP(sub,:,:) = flipSide16Matrix(PS.ERSPSes{sub}{end},doflip);
% end
% 
% % Test for desired conditions
% p = [];
% h = [];
% for ch=1:size(PreERSP,2)
%     for fr=1:size(PreERSP,3)
%         
%         % BCI should be significantly different pre-post
%         [h(1,ch,fr) p(1,ch,fr)] = ttest2(PreERSP(IndBCI,ch,fr),PostERSP(IndBCI,ch,fr));
%         % Sham should NOT be significantly different pre-post
%         [h(2,ch,fr) p(2,ch,fr)] = ttest2(PreERSP(IndSham,ch,fr),PostERSP(IndSham,ch,fr));
%         
%         % BCI vs Sham should NOT be different Pre
%         [h(3,ch,fr) p(3,ch,fr)] = ttest2(PreERSP(IndBCI,ch,fr),PreERSP(IndSham,ch,fr));
%         % BCI vs Sham should be different Post        
%         [h(4,ch,fr) p(4,ch,fr)] = ttest2(PostERSP(IndBCI,ch,fr),PostERSP(IndSham,ch,fr));
%         
%     end
% end
% 
% SelFeatCh=15;
% SelFeatFr=17; % Index, not Hz
% erddiffplot(PreERSP,PostERSP,IndBCI,IndSham,SelFeatCh,SelFeatFr);
% 
% SelFeatCh=10;
% SelFeatFr=8; % Index, not Hz
% erddiffplot(PreERSP,PostERSP,IndBCI,IndSham,SelFeatCh,SelFeatFr);
% 
% % LocReg{1} = [2 3 7 8 12 13];
% % LocReg{2} = [4 9 14];
% % LocReg{3} = [5 6 10 11 15 16];
% % FrReg{1} = [3 4 5 6];
% % FrReg{2} = [7 8 9 10];
% % FrReg{3} = [11 12 13 14];
% % 
% % RegPreERSP = loweResERSP(PreERSP,LocReg,FrReg);
% % RegPostERSP = loweResERSP(PostERSP,LocReg,FrReg);
% 
% 
