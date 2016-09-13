% Parametrization
freqs = [4:2:48]; % Desired frequency bands for PSD features


addpath(genpath('~/Git/cnbi-smrtrain/'));
lap = load('laplacian16.mat');
lap = lap.lap;

%Path = '/mnt/cnbiserver/cnbi-commun/_INBOX/Data/CNBI_2016_StrokeMagdeburg_PerdikisSerafeim/';
Path = '~/Desktop/tst/';
SavePath = '~/tmp/';

SubDir = dir(Path);
SubDir = SubDir(3:end);
isd = [SubDir(:).isdir];
SubDir = SubDir(isd);

for subject = 1:length(SubDir)
    Acc= {};
    TrAcc = {};
    Sub = SubDir(subject).name;
    SubSes = dir([Path '/' Sub]);
    SubSes = SubSes(3:end);
    isd = [SubSes(:).isdir];
    SubSes = SubSes(isd);
    
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
                    
                    if(exist([SavePath Sub '/' GDFName(1:end-4) '.mat'],'file') == 0)
                        RunResults = analyzeOnlineStroke(GDFPath,MATPath, lap, freqs);
                    else
                        load([SavePath Sub '/' GDFName(1:end-4) '.mat']); 
                    end
                    
                    if(RunResults.fine == 1)
                        % Save results file
                        if(~exist([SavePath '/' Sub],'dir'))
                            % Create subject's playback folder
                            mkdir(SavePath,Sub);
                        end

                        run = run+1;

                        disp(['Subject: ' Sub ' , Session: ' num2str(onses) ' , Run: ' num2str(run)]);

                        save([SavePath Sub '/' GDFName(1:end-4) '.mat'],'RunResults');
                        
                    end
                end
                
                fclose(fid);
                
            end
        end
    end
    % Save concentrated subject results here
    %save([SavePath Sub '/' Sub '_Acc.mat'],'Sum');
end