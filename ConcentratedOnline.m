FigureSavePath = '/home/sperdikis/Subversion/simis/trunk/WyssAcuteStroke/img/';
Path = '/home/sperdikis/Data/Results/AcuteStrokeResults/';
SubID  = {'ckg8','ds86','fh47','jy18','ma93','qv39','rj31','wu60','ya00','odr2','ji34','ao48','pk72','lm90', 'rai1','mo17'};
Group  = {'BCI','Sham','Sham','BCI','BCI','Sham','BCI','BCI','BCI','BCI','Sham','BCI','BCI','Sham','BCI','BCI'};

Side = {'L','L','L','L','L','L','L','L','R','L','L','R','L','L','R','L'};

Electrodes16 = {'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'};

% Clinical outcomes
FMAUE{1} = [8 0 0 0 0 0  0 23 6  12 30  0 0 4 21 NaN];
FMAUE{2} = [9 0 8 0 1 10 0 47 32 40 NaN 3 35 9 38 NaN];
FMAdiff = FMAUE{2}-FMAUE{1};

ESS{1} = [48 42 39 48 57 52 57 74 57 72 70 61 57 42 NaN NaN];
ESS{2} = [48 54 72 64 67 82 61 84 81 84 80 73 93 56 NaN NaN];
ESSdiff = ESS{2}-ESS{1};

BarthelIndex{1} = [25 55 5  75 50 10 35 70 50 50 70 15  60  20 NaN NaN];
BarthelIndex{2} = [40 60 65 85 60 60 65 85 55 80 70 100 100 40 NaN NaN];
BarthelIndexdiff = BarthelIndex{2}-BarthelIndex{1}; 

% Print info for clinical outcome statistical testing
printCO(FMAUE{1}, FMAUE{2}, Group, 'FMA-UE');
printCO(ESS{1}, ESS{2}, Group, 'ESS');
printCO(BarthelIndex{1}, BarthelIndex{2}, Group, 'Barthel Index');

for sub=1:length(SubID)
    
    if(strcmp(Side{sub},'L'))
        affectedtask = 'lhrst';
    else
        affectedtask = 'rhrst';
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
    PS.FSavg{sub} = SubResults.FSAll;
    PS.ERSPavg{sub} = SubResults.ERSPAll;
    
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


%% Figures
close all;

%% Detection Rate plots
% Average detection rate per subject
h=figure(1);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.DetectionRate(sub),'b','LineWidth',3);
    else
        bar(sub,PS.DetectionRate(sub),'r','LineWidth',3);
    end
end
bar(length(SubID)+1, mean(PS.DetectionRate(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
errorbar(length(SubID)+1, mean(PS.DetectionRate(find(strcmp(Group,'BCI')))),std2(PS.DetectionRate(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
bar(length(SubID)+2, mean(PS.DetectionRate(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
errorbar(length(SubID)+2, mean(PS.DetectionRate(find(strcmp(Group,'Sham')))),std2(PS.DetectionRate(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
axis([0.5 length(SubID)+2.5 40 101]);
line([0.5 length(SubID)+2.5],[50 50],'Color','k','LineWidth',3,'LineStyle','--');
line([0.5 length(SubID)+2.5],[58 58],'Color','k','LineWidth',3,'LineStyle','--');
line([length(SubID)+0.5 length(SubID)+0.5],[40 101],'Color','k','LineWidth',3,'LineStyle','--');
hold off;
legend({'BCI','Sham'},'FontSize',20);
xlabel('Participants');
ylabel('Single-sample online detection rate (%)');
set(gca,'LineWidth',3,'FontSize',20);
set(gca,'XTick',[1:length(SubID)+2]);
set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);



%% Hit rate plot
% Average detection rate per subject
h=figure(2);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.HitRate(sub),'b','LineWidth',3);
    else
        bar(sub,PS.HitRate(sub),'r','LineWidth',3);
    end
end
bar(length(SubID)+1, mean(PS.HitRate(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
errorbar(length(SubID)+1, mean(PS.HitRate(find(strcmp(Group,'BCI')))),std2(PS.HitRate(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
bar(length(SubID)+2, mean(PS.HitRate(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
errorbar(length(SubID)+2, mean(PS.HitRate(find(strcmp(Group,'Sham')))),std2(PS.HitRate(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
axis([0.5 length(SubID)+2.5 0 101]);
line([0.5 length(SubID)+2.5],[50 50],'Color','k','LineWidth',3,'LineStyle','--');
line([length(SubID)+0.5 length(SubID)+0.5],[40 101],'Color','k','LineWidth',3,'LineStyle','--');
hold off;
legend({'BCI','Sham'},'FontSize',20);
xlabel('Participants');
ylabel('Hit rate (%)');
set(gca,'LineWidth',3,'FontSize',20);
set(gca,'XTick',[1:length(SubID)+2]);
set(gca,'XTickLabel',[SubID 'BCI' 'Sham']);


%% Simulated Accuracy plots
% Average simulated accuracy (selected) per subject
h=figure(3);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.SimAccSel(sub),'b','LineWidth',3);
    else
        bar(sub,PS.SimAccSel(sub),'r','LineWidth',3);
    end
end
bar(length(SubID)+1, mean(PS.SimAccSel(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
errorbar(length(SubID)+1, mean(PS.SimAccSel(find(strcmp(Group,'BCI')))),std2(PS.SimAccSel(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
bar(length(SubID)+2, mean(PS.SimAccSel(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
errorbar(length(SubID)+2, mean(PS.SimAccSel(find(strcmp(Group,'Sham')))),std2(PS.SimAccSel(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
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


% Average simulated accuracy (best) per subject
h=figure(4);
hold on;
for sub=1:length(SubID)
    if(strcmp(Group{sub},'BCI'))
        bar(sub,PS.SimAccBest(sub),'b','LineWidth',3);
    else
        bar(sub,PS.SimAccBest(sub),'r','LineWidth',3);
    end
end
bar(length(SubID)+1, mean(PS.SimAccBest(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
errorbar(length(SubID)+1, mean(PS.SimAccBest(find(strcmp(Group,'BCI')))),std2(PS.SimAccBest(find(strcmp(Group,'BCI')))),'b','LineWidth',3);
bar(length(SubID)+2, mean(PS.SimAccBest(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
errorbar(length(SubID)+2, mean(PS.SimAccBest(find(strcmp(Group,'Sham')))),std2(PS.SimAccBest(find(strcmp(Group,'Sham')))),'r','LineWidth',3);
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
% %% Fisher Score
% for sub=1:length(SubID)
%     MinFS(sub) = min(PS.FSavg{sub}(:));
%     MaxFS(sub) = max(PS.FSavg{sub}(:));
%         if(strcmp(Side{sub},'R'))
%         doflip = 1;
%     else
%         doflip = 0;
%     end
%     %TopoFSMu{sub} = convChans(flipSide16(mean(PS.FSavg{sub}(:, 3:6),2),doflip));
%     %TopoFSBeta{sub} = convChans(flipSide16(mean(PS.FSavg{sub}(:,8:11),2),doflip));
%     TopoFSMu{sub} = convChans(flipSide16(max(PS.FSavg{sub}(:,3:6)'),doflip));
%     TopoFSBeta{sub} = convChans(flipSide16(max(PS.FSavg{sub}(:,8:11)'),doflip));
% end
% %MinFS = min(MinFS);
% %MaxFS = max(MaxFS);
% 
% % Load colormaps
% cmerds = load('erdersmap.mat');cmerds = cmerds.cm;
% cmfs = load('fsmap.mat');cmfs = cmfs.cm;
% % Load chanlocs 64
% chanlocs = load('chanlocs.mat');chanlocs=chanlocs.chanlocs;
% 
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
% % FS
% for sub=1:length(SubID)
%     figure(9);
%     %topoplot(TopoFSMu{sub},chanlocs,'maplimits',[MinFS(sub) MaxFS(sub)]);colormap(cmfs);drawnow;
%     topoplot(TopoFSMu{sub},chanlocs);
%     %print([FigureSavePath '/fs/' SubID{sub} '_Mu'],'-r600','-dpng');
%     figure(10);
%     topoplot(TopoFSBeta{sub},chanlocs,'maplimits',[MinFS(sub) MaxFS(sub)]);colormap(cmfs);drawnow;
%     %print([FigureSavePath '/fs/' SubID{sub} '_Beta'],'-r600','-dpng');
% end




% 
% %% Per session DP maps for the QR Appendix
% MaxAbsERSP = 0;
% for sub=1:length(SubID)
%     % Find biggest ERSP overall
%     for ses=1:PS.SesNum(sub)
%         MaxAbsERSP = max(MaxAbsERSP,max(abs(PS.ERSPSes{sub}{ses}(:))));
%     end
% end
% 
% close all;
% for sub=1:length(SubID)
%     figure(sub);
%     Nrows = ceil(PS.SesNum(sub)/5);
%     if(strcmp(Side{sub},'R'))
%         doflip = 1;
%     else
%         doflip = 0;
%     end
%     for ses=1:PS.SesNum(sub)
%         subplot(Nrows, 5, ses);imagesc(flipSide16Matrix(PS.ERSPSes{sub}{ses},doflip),[-1.0 1.0]);
%         if(ses==1)
%             colorbar;
%         end
%         title(['Session '  num2str(ses)]);
%         xlabel('Bands (Hz)');
%         ylabel('Channels');
%         set(gca,'LineWidth',3,'FontSize',10);
%         set(gca,'XTick',[1:4:23]);
%         set(gca,'XTickLabel',[4:8:48]);
%         set(gca,'YTick',[1:16]);
%         set(gca,'YTickLabel',Electrodes16);
%     end
% end
% 




close all;
%% Get some ERSP statistics using the first and last sessions
% Use only first and last session
PreERSP = [];
PostERSP = [];
for sub=1:length(SubID)
    if(strcmp(Side{sub},'R'))
        doflip = 1;
    else
        doflip = 0;
    end
    PreERSP(sub,:,:) = flipSide16Matrix(PS.ERSPSes{sub}{1},doflip);
    PostERSP(sub,:,:) = flipSide16Matrix(PS.ERSPSes{sub}{end},doflip);
end

% Test for desired conditions
p = [];
h = [];
for ch=1:size(PreERSP,2)
    for fr=1:size(PreERSP,3)
        
        % BCI should be significantly different pre-post
        [h(1,ch,fr) p(1,ch,fr)] = ttest2(PreERSP(IndBCI,ch,fr),PostERSP(IndBCI,ch,fr));
        % Sham should NOT be significantly different pre-post
        [h(2,ch,fr) p(2,ch,fr)] = ttest2(PreERSP(IndSham,ch,fr),PostERSP(IndSham,ch,fr));
        
        % BCI vs Sham should NOT be different Pre
        [h(3,ch,fr) p(3,ch,fr)] = ttest2(PreERSP(IndBCI,ch,fr),PreERSP(IndSham,ch,fr));
        % BCI vs Sham should be different Post        
        [h(4,ch,fr) p(4,ch,fr)] = ttest2(PostERSP(IndBCI,ch,fr),PostERSP(IndSham,ch,fr));
        
    end
end

SelFeatCh=15;
SelFeatFr=17; % Index, not Hz
erddiffplot(PreERSP,PostERSP,IndBCI,IndSham,SelFeatCh,SelFeatFr);

SelFeatCh=10;
SelFeatFr=8; % Index, not Hz
erddiffplot(PreERSP,PostERSP,IndBCI,IndSham,SelFeatCh,SelFeatFr);

% LocReg{1} = [2 3 7 8 12 13];
% LocReg{2} = [4 9 14];
% LocReg{3} = [5 6 10 11 15 16];
% FrReg{1} = [3 4 5 6];
% FrReg{2} = [7 8 9 10];
% FrReg{3} = [11 12 13 14];
% 
% RegPreERSP = loweResERSP(PreERSP,LocReg,FrReg);
% RegPostERSP = loweResERSP(PostERSP,LocReg,FrReg);


