clear;close all;clc; 
% Parametrization of all doc code
% another test comment
%change
freqs = [4:2:48]; % Desired frequency bands for PSD features

addpath(genpath('/home/ido/New DOC Literature/DOCpipeline/Code/biosig/'));
lap = load('laplacian16.mat');
lap = lap.lap;

% lap64 = load('laplacian64HiAmp.mat');
% lap64 = lap64.lap.HiAmp64Mat;
% 


Path = '/home/ido/New DOC Literature/DOCpipeline/Code/code/EEG';

SavePath = '~/New DOC Literature/DOCpipeline/Code/code/SavedData/';

SubDir = dir(Path);
SubDir = SubDir(3:end);
isd = [SubDir(:).isdir];
SubDir = SubDir(isd);
%SubDir = SubDir(1);


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
                
                GDFName = Line(1:strfind(Line, 'gdf')+2);
                GDFPath = [Path '/' Sub '/' SesName '/' GDFName];
              
             
                if(exist(GDFPath,'file')>0)
                    ff = dir(GDFPath);
                    if( (ff.bytes / (1024^2)) < 1.0 ) % Get rid of too small GDFs, probably failed attempts to start the loop or interrupted runs
                        continue;
                    end
                else
                    continue;
                end
                
                if( (exist([SavePath Sub '/savedData/' GDFName(1:end-4) '.mat'],'file') == 0) && (exist([SavePath Sub '/excluded/' GDFName(1:end-4) '.mat'],'file') == 0))
                    RunResults = analyzeOnlineStroke(GDFPath,lap, freqs);
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
    

