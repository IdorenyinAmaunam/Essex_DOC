% Parametrization
freqs = [4:2:48]; % Desired frequency bands for PSD features


addpath(genpath('~/Git/cnbi-smrtrain/'));
lap = load('laplacian16.mat');
lap = lap.lap;

lap64 = load('laplacian64HiAmp.mat');
lap64 = lap64.lap.HiAmp64Mat;


%Path = '/mnt/cnbiserver/cnbi-commun/data/processed/CNBI_2016_StrokeMagdeburg_PerdikisSerafeim/good/';
%Path = '/mnt/cnbiserver/cnbi-commun/_INBOX/Data/CNBI_2016_StrokeMagdeburg_PerdikisSerafeim/good/';
%Path = '/mnt/cnbiserver/cnbi-commun/_INBOX/Data/CNBI_2016_AcuteStrokeSUVA_PerdikisSerafeim/';
%Path = '/mnt/cnbiserver/cnbi-commun/_INBOX/Data/CNBI_2016_AcuteStrokeLavigny_PerdikisSerafeim/';
Path = '/mnt/cnbiserver/cnbi-commun/_INBOX/Data/CNBI_2017_AcuteStrokeSanCamillo_PerdikisSerafeim/';
%Path = '~/Desktop/tst/';
SavePath = '~/Data/Results/AcuteStrokeResults/';

SubDir = dir(Path);
SubDir = SubDir(3:end);
isd = [SubDir(:).isdir];
SubDir = SubDir(isd);

for subject = 1:length(SubDir)
    
    Sub = SubDir(subject).name;
    SubSes = dir([Path '/' Sub]);
    SubSes = SubSes(3:end);
    isd = [SubSes(:).isdir];
    SubSes = SubSes(isd);

    
    % Save results file
    if(~exist([SavePath '/' Sub],'dir'))
        % Create subject's folder
        mkdir(SavePath,Sub);
        mkdir([SavePath '/' Sub],'excluded');
    end
    
    onses = 0;
    for ses=1:length(SubSes)
        % Check if it is an online session
        SesName = SubSes(ses).name;
        
        if( strcmp(SesName(1:5),[Sub '.']) && (mean(isstrprop(SesName(6:end),'digit'))==1))
            onses = onses + 1;
            
            % Load log file
            LogFile = dir([Path '/' Sub '/' SesName '/*.log']);
            if(~isempty(LogFile))
                LogFile = LogFile.name;
                fid = fopen([Path '/' Sub '/' SesName '/' LogFile]);
                
                run=0;
                while(true)
                    % Read log line
                    Line = fgetl(fid);
                    if(Line == -1)
                        break;
                    end
                    
                    GDFName = Line(strfind(Line, SesName):strfind(Line, 'gdf')+2);
                    MATname = Line(strfind(Line, [Sub '_']):strfind(Line, 'mat')+2);
                    GDFPath = [Path '/' Sub '/' SesName '/' GDFName];
                    MATPath = [Path '/' Sub '/' MATname];
                    
                    if(exist(GDFPath,'file')>0)
                        ff = dir(GDFPath);
                        if( (ff.bytes/(1024^2)) < 1.0 ) % Get rid of too small GDFs, probably failed attempts to start the loop or interrupted runs
                            continue;
                        end
                    else
                        continue;
                    end
                    
                    if( (exist([SavePath Sub '/' GDFName(1:end-4) '.mat'],'file') == 0) && (exist([SavePath Sub '/excluded/' GDFName(1:end-4) '.mat'],'file') == 0))
                        RunResults = analyzeOnlineStroke(GDFPath,MATPath, lap, freqs);
                        if(RunResults.fine == 1)
                            save([SavePath Sub '/' GDFName(1:end-4) '.mat'],'RunResults');
                        else
                            % Save excluded dummy mat file
                            save([SavePath Sub '/excluded/' GDFName(1:end-4) '.mat'],'RunResults');
                            continue;
                        end
                    else
                        if(exist([SavePath Sub '/' GDFName(1:end-4) '.mat'],'file') == 2)
                            load([SavePath Sub '/' GDFName(1:end-4) '.mat']);
                        else
                            % Faulty run saved, skip it
                            continue;
                        end
                    end
                    
                    disp(['Subject: ' Sub ' , Session: ' num2str(onses) ' , Run: ' num2str(run)]);
                    run = run+1;                    
                end
                fclose(fid);
            end
        end
        
    end

    
    offses64 = 0;
    for ses=1:length(SubSes)
        % Check if it is an offline 64ch session
        SesName = SubSes(ses).name;
        SesName
        if( strcmp(SesName(1:5),[Sub '.']) && strcmp(SesName(end-4:end),'_64ch') )
            offses64 = offses64 + 1;
            
            % Load log file
            LogFile = dir([Path '/' Sub '/' SesName '/*.log']);
            if(~isempty(LogFile))
                LogFile = LogFile.name;
                fid = fopen([Path '/' Sub '/' SesName '/' LogFile]);
                
                run=0;
                while(true)
                    % Read log line
                    Line = fgetl(fid);
                    if(Line == -1)
                        break;
                    end
                    
                    GDFName = Line(strfind(Line, Sub):strfind(Line, 'gdf')+2);
                    GDFPath = [Path '/' Sub '/' SesName '/' GDFName];
                    
                    if(exist(GDFPath,'file')>0)
                        ff = dir(GDFPath);
                        if( (ff.bytes/(1024^2)) < 1.0 ) % Get rid of too small GDFs, probably failed attempts to start the loop or interrupted runs
                            continue;
                        end
                    else
                        continue;
                    end
                    
                    if( (exist([SavePath Sub '/' GDFName(1:end-4) '.mat'],'file') == 0) && (exist([SavePath Sub '/excluded/' GDFName(1:end-4) '.mat'],'file') == 0))
                        RunResults = analyzeOffline64Stroke(GDFPath, lap64, freqs);
                        if(RunResults.fine == 1)
                            save([SavePath Sub '/' GDFName(1:end-4) '.mat'],'RunResults');
                        else
                            % Save excluded dummy mat file
                            save([SavePath Sub '/excluded/' GDFName(1:end-4) '.mat'],'RunResults');
                            continue;
                        end
                    else
                        if(exist([SavePath Sub '/' GDFName(1:end-4) '.mat'],'file') == 2)
                            load([SavePath Sub '/' GDFName(1:end-4) '.mat']);
                        else
                            % Faulty run saved, skip it
                            continue;
                        end
                    end
                    
                    run = run+1;
                    disp(['Subject: ' Sub ' , Offline 64 Session: ' num2str(offses64) ' , Run: ' num2str(run)]);
                end
                fclose(fid);
            end
        end
        
    end
    
end