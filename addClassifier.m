%Path = '/mnt/cnbiserver/cnbi-commun/data/raw/inbox/TOBI_BCI_FES_Stroke/'; % Andrea's part 
Path = '/mnt/cnbiserver/cnbi-commun/data/raw/inbox/CNP_VS_BCI_FES_Stroke/'; % Robert's part
SavePath = '~/Data/Results/ChronicStrokeResults/';

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
        
        if( strcmp(SesName(1:length(Sub)+1),[Sub '_']) && (mean(isstrprop(SesName(length(Sub)+2:end),'digit'))==1))
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
                    
                    GDFName = Line(strfind(Line, [SesName(1:length(Sub)) '.' SesName(length(Sub)+2:end)]):strfind(Line, 'gdf')+2);
                    MATname = Line(strfind(Line, [Sub '_']):strfind(Line, 'mat')+2);
                    if(isempty(MATname))
                        % Sometimes extension is missing
                        MATname = [Line(strfind(Line, [Sub '_']):strfind(Line, [Sub '_'])+14+length(Sub)) '.mat'];
                    end
                    GDFPath = [Path '/' Sub '/' SesName '/' GDFName];
                    MATPath = [Path '/' Sub '/' MATname];
                    
                    if(isempty(GDFName))
                        % It means that this line contained no GDF file
                        % name
                        continue;
                    end
                    
                    if( exist([SavePath Sub '/' GDFName(1:end-4) '.mat'],'file') == 2 )
                        clear RunResults;
                        
                        % Load Run Results
                        load([SavePath Sub '/' GDFName(1:end-4) '.mat']);
                        
                        % Try to load classifier
                        try
                            analysis = load(MATPath);
                            analysis = analysis.analysis;
                            RunResults.analysis = analysis;
                        catch
                            RunResults.analysis = [];
                        end
                        
                        % Save back with updated info
                        save([SavePath Sub '/' GDFName(1:end-4) '.mat'],'RunResults');
                    end
                    
                    disp(['Subject: ' Sub ' , Session: ' num2str(onses) ' , Run: ' num2str(run+1)]);
                    run = run+1;                    
                end
                fclose(fid);
            end
        end
        
    end
    
end