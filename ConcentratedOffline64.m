FigureSavePath = '/home/sperdikis/Subversion/simis/trunk/WyssAcuteStroke/img2/';
Path = '/home/sperdikis/AcuteStrokeResults/';
%SubID  = {'ckg8','ds86','fh47','jy18','ma93','qv39','rj31','wu60','ya00','odr2'};
SubID  = {'odr2','ckg8'};
Group  = {'BCI','BCI'};
Side = {'L','L'};
Electrodes64 = {'FP1','FPz','FP2','AF7','AF3','AF4','AF8','F7','F5','F3',...
    'F1','Fz','F2','F4','F6','F8','FT7','FC5','FC3','FC1',...
    'FCz','FC2','FC4','FC6','FT8','T7','C5','C3','C1','Cz',...
    'C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPz','CP2',...
    'CP4','CP6','TP8','P7','P5','P3','P1','Pz','P2','P4',...
    'P6','P8','PO7','PO3','POz','PO4','PO8','O1','Oz','O2',...
    'F9','F10'};

for sub=1:length(SubID)
    
    if(strcmp(Side{sub},'L'))
        affectedtask = 'lhrst';
    else
        affectedtask = 'rhrst';
    end
    
    % Initialize per-group variables here

    if(exist([Path '/' SubID{sub} '/' SubID{sub} '_ResOff62.mat'],'file'))
        load([Path '/' SubID{sub} '/' SubID{sub} '_ResOff62.mat']);
    else
        
        % Arrange sessions in temporal order
        RunDir = dir([Path '/' SubID{sub}]);
        RunDir = RunDir(3:end);
        % Keep only offline 62 runs
        KeepInd = [];
        for r=1:length({RunDir.name})
            if( (~isempty(strfind(RunDir(r).name,'offline'))) && (~isempty(strfind(RunDir(r).name,affectedtask))))
                KeepInd = [KeepInd; r];
            end
        end
        RunDir = RunDir(KeepInd);
        
        %IndDir = [RunDir(:).isdir];
        %FileNames = {RunDir(find(IndDir~=1)).name};
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
        AllSimAccBest01 = [];
        AllSimAccBest02 = [];
        AllSimAccBest12 = [];        
        AllSI01 = [];
        AllSI02 = [];
        AllSI12 = [];
        AllFS = {};
        AllERSPAct = {};
        AllERSPRest = {};
        
        % Analyze per session
        for ses=1:length(SesDates)
            ThisSesRunind = find(RunSesInd==ses);

            % Initialize SES variables
            SesSimAccBest01 = [];
            SesSimAccBest02 = [];
            SesSimAccBest12 = [];
            SesSI01 = [];
            SesSI02 = [];
            SesSI12 = [];
            SesFS = {};
            SesERSPAct = {};
            SesERSPRest = {};
            
            for run=1:length(ThisSesRunind)
                load([Path '/' SubID{sub} '/' FileNames{ThisSesRunind(run)}]);
                disp(['Processing Subject:' SubID{sub} ', Session: ' num2str(ses) ', Run: ' num2str(run)]);
                SesSimAccBest01 = [SesSimAccBest01 ; RunResults.SimulatedAccBest01];
                AllSimAccBest01 = [AllSimAccBest01 ; RunResults.SimulatedAccBest01];
                SesSimAccBest02 = [SesSimAccBest02 ; RunResults.SimulatedAccBest02];
                AllSimAccBest02 = [AllSimAccBest02 ; RunResults.SimulatedAccBest02];
                SesSimAccBest12 = [SesSimAccBest12 ; RunResults.SimulatedAccBest12];
                AllSimAccBest12 = [AllSimAccBest12 ; RunResults.SimulatedAccBest12];                
                AllSI01 = [AllSI01; RunResults.SI01];
                SesSI01 = [SesSI01; RunResults.SI01];
                AllSI02 = [AllSI02; RunResults.SI02];
                SesSI02 = [SesSI02; RunResults.SI02];
                AllSI12 = [AllSI12; RunResults.SI12];
                SesSI12 = [SesSI12; RunResults.SI12];                
                SesFS = [SesFS ; RunResults.FS'];
                AllFS = [AllFS ; RunResults.FS'];
                SesERSPAct = [SesERSPAct ; RunResults.erspAct];
                AllERSPAct = [AllERSPAct ; RunResults.erspAct];
                SesERSPRest = [SesERSPRest ; RunResults.erspRest];
                AllERSPRest = [AllERSPRest ; RunResults.erspRest];                
            end

            %% Per Session averages
            %% Accuracies
            SubResults.SimAccBestSes01(ses) = nanmean(SesSimAccBest01);
            SubResults.SimAccBestSes02(ses) = nanmean(SesSimAccBest02);
            SubResults.SimAccBestSes12(ses) = nanmean(SesSimAccBest12);            
            % SI
            SubResults.SISes01(ses) = nanmean(SesSI01);
            SubResults.SISes02(ses) = nanmean(SesSI02);
            SubResults.SISes12(ses) = nanmean(SesSI12);
            
            %FS
            Nrun = length(SesFS);
            SubResults.FSSes{ses} = squeeze(nanmean(reshape(cell2mat(SesFS),[Nrun 62 23]),1));
            %ERSP
            Nrun = length(SesERSPAct);
            SubResults.ERSPSesAct{ses} = squeeze(nanmean(reshape(cell2mat(SesERSPAct),[Nrun 62 23]),1));
            Nrun = length(SesERSPRest);
            SubResults.ERSPSesRest{ses} = squeeze(nanmean(reshape(cell2mat(SesERSPRest),[Nrun 62 23]),1));            

        end

        SubResults.RunSesInd = RunSesInd; 
        SubResults.SimAccBestAll01 = AllSimAccBest01;
        SubResults.SimAccBestAll02 = AllSimAccBest02;
        SubResults.SimAccBestAll12 = AllSimAccBest12;
        
        % SI
        SubResults.SIAll01 = AllSI01;
        SubResults.SIAll02 = AllSI02;
        SubResults.SIAll12 = AllSI12;        
        
        %FS
        Nrun = length(AllFS);
        SubResults.FSAll = squeeze(nanmean(reshape(cell2mat(AllFS),[Nrun 62 23]),1));
        % ERSP
        Nrun = length(AllERSPAct);
        SubResults.ERSPAllAct = squeeze(nanmean(reshape(cell2mat(AllERSPAct),[Nrun 62 23]),1));        
        Nrun = length(AllERSPRest);
        SubResults.ERSPAllRest = squeeze(nanmean(reshape(cell2mat(AllERSPRest),[Nrun 62 23]),1));                

        % Save subject results
        if(~exist([Path '/' SubID{sub} '/' SubID{sub} '_ResOff62.mat'],'file'))
            save([Path '/' SubID{sub} '/' SubID{sub} '_ResOff62.mat'],'SubResults');
        end
    end
    
    % Per-subject variables
    PS.SimAccBest01(sub) = nanmean(SubResults.SimAccBestAll01);
    PS.SimAccBest02(sub) = nanmean(SubResults.SimAccBestAll02);
    PS.SimAccBest12(sub) = nanmean(SubResults.SimAccBestAll12);
    
    PS.SI01(sub) = nanmean(SubResults.SIAll01);
    PS.SI02(sub) = nanmean(SubResults.SIAll02);
    PS.SI12(sub) = nanmean(SubResults.SIAll12);
    
    PS.MaxFS(sub) = max(SubResults.FSAll(:));
    
    PS.SimAccBestSes01{sub} = SubResults.SimAccBestSes01;
    PS.SimAccBestSes02{sub} = SubResults.SimAccBestSes02;
    PS.SimAccBestSes12{sub} = SubResults.SimAccBestSes12;
    PS.SISes01{sub} = SubResults.SISes01;
    PS.SISes02{sub} = SubResults.SISes02;
    PS.SISes12{sub} = SubResults.SISes12;
    PS.FSSes{sub} = SubResults.FSSes;
    PS.ERSPSesAct{sub} = SubResults.ERSPSesAct;
    PS.ERSPSesRest{sub} = SubResults.ERSPSesRest;
    PS.FSavg{sub} = SubResults.FSAll;
    PS.ERSPavgAct{sub} = SubResults.ERSPAllAct;
    PS.ERSPavgRest{sub} = SubResults.ERSPAllRest;
    
    clear SubResults;
end

%% Figures
close all;

%% Simulated Accuracy plots

%% Per session DP maps
close all;
load chanlocs64.mat
for sub=1:length(SubID)
    if(strcmp(Side{sub},'R'))
        doflip = 1;
    else
        doflip = 0;
    end
    figure(10*sub+1);imagesc(flipSide64Matrix(PS.ERSPSesAct{sub}{1},doflip),[-1.0 1.0]);
    xlabel('Bands (Hz)');ylabel('Channels');set(gca,'LineWidth',3,'FontSize',10);set(gca,'XTick',[1:2:23]);set(gca,'XTickLabel',[4:4:48]);set(gca,'YTick',[1:62]);set(gca,'YTickLabel',Electrodes64);
    figure(10*sub+2);imagesc(flipSide64Matrix(PS.ERSPSesAct{sub}{2},doflip),[-1.0 1.0]);
    xlabel('Bands (Hz)');ylabel('Channels');set(gca,'LineWidth',3,'FontSize',10);set(gca,'XTick',[1:2:23]);set(gca,'XTickLabel',[4:4:48]);set(gca,'YTick',[1:62]);set(gca,'YTickLabel',Electrodes64);
    figure(10*sub+3);imagesc(flipSide64Matrix(PS.ERSPSesRest{sub}{1},doflip),[-1.0 1.0]);
    xlabel('Bands (Hz)');ylabel('Channels');set(gca,'LineWidth',3,'FontSize',10);set(gca,'XTick',[1:2:23]);set(gca,'XTickLabel',[4:4:48]);set(gca,'YTick',[1:62]);set(gca,'YTickLabel',Electrodes64);
    figure(10*sub+4);imagesc(flipSide64Matrix(PS.ERSPSesRest{sub}{2},doflip),[-1.0 1.0]);
    xlabel('Bands (Hz)');ylabel('Channels');set(gca,'LineWidth',3,'FontSize',10);set(gca,'XTick',[1:2:23]);set(gca,'XTickLabel',[4:4:48]);set(gca,'YTick',[1:62]);set(gca,'YTickLabel',Electrodes64);
    
    % Topos
    figure(100*sub+1);topoplot(flipSide64(mean(PS.ERSPSesAct{sub}{1}(:,3:6),2),doflip),chanlocs64,'maplimits',[-1 1]);
    figure(100*sub+2);topoplot(flipSide64(mean(PS.ERSPSesAct{sub}{1}(:,7:11),2),doflip),chanlocs64,'maplimits',[-1 1]);
    figure(100*sub+3);topoplot(flipSide64(mean(PS.ERSPSesAct{sub}{2}(:,3:6),2),doflip),chanlocs64,'maplimits',[-1 1]);
    figure(100*sub+4);topoplot(flipSide64(mean(PS.ERSPSesAct{sub}{2}(:,7:11),2),doflip),chanlocs64,'maplimits',[-1 1]);
    figure(100*sub+5);topoplot(flipSide64(mean(PS.ERSPSesRest{sub}{1}(:,3:6),2),doflip),chanlocs64,'maplimits',[-1 1]);
    figure(100*sub+6);topoplot(flipSide64(mean(PS.ERSPSesRest{sub}{1}(:,7:11),2),doflip),chanlocs64,'maplimits',[-1 1]);
    figure(100*sub+7);topoplot(flipSide64(mean(PS.ERSPSesRest{sub}{2}(:,3:6),2),doflip),chanlocs64,'maplimits',[-1 1]);
    figure(100*sub+8);topoplot(flipSide64(mean(PS.ERSPSesRest{sub}{2}(:,7:11),2),doflip),chanlocs64,'maplimits',[-1 1]);
end