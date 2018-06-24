%% TIME EXTRACTOR FOR SEPARATED CONDITIONS 
% Conditions random, pmax, gmax, lmin
% The default for this project is to run 5 blocks
% It creates three separate arrays for conditions, onsets and durations
% This creates an output for the contrast creator for SPM

%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% June 2018
% working on SPM8

%% Clean workspace
clc; clear

%% Base Paths
cd('..')
folder.Root      = pwd;
folder.Processed = fullfile(folder.Root, 'Processed');
folder.Behavior  = fullfile(folder.Root, 'Behavioral');
folder.Time      = fullfile(folder.Root, 'Behavioral', 'Timing');
folder.Scripts   = fullfile(folder.Root, 'Scripts'); 

%% Get all subject paths
folder.ProcessedPaths      = dir(folder.Processed);
folder.ProcessedPaths(1:2) = [];

%% Ask for which subjects to run
[options.Group, ~] = listdlg('ListString',{'Individual Elements','All Subjects'},'Name','No. Subjects to Process?');

%% Set subject paths for input and output images according to subjects selected
if options.Group == 1 % Customized list
    [options.Subjects, ~] = listdlg('ListString',char(folder.ProcessedPaths.name),'Name','Which subjects do you want?');
    group.SubjectsList    = folder.ProcessedPaths(options.Subjects);

elseif options.Group == 2 % All subjects
    group.SubjectsList = folder.ProcessedPaths;
end

for iFolder = 1:size(group.SubjectsList,1)
    group.SubjectsPaths{iFolder,:}  = fullfile(folder.Processed, group.SubjectsList(iFolder).name);
    group.TimePaths{iFolder,:} = fullfile(folder.Processed, group.SubjectsList(iFolder).name, 'timing');
end

if ~exist(folder.Time,'dir')
    mkdir(folder.Time);
end

%% Loop throught the subject list
for iSubj = 1:size(group.SubjectsPaths,1)   
    
    subject.TimeFolder = group.TimePaths{iSubj,:};
    subject.ID         = char(extractAfter(group.SubjectsList(iSubj).name, 's'));
    
    if ~exist(subject.TimeFolder, 'dir')
        mkdir(subject.TimeFolder);
    end
    
    %% Loop through all the functional runs
    for iRun = 1:5
        %% Clear run values to avoid overwritting issues
        clear run
        clear names
        clear onsets
        clear durations
        
        %% Set file names to load and save
        run.BehaviorFile = fullfile(folder.Behavior, ['/Gamble_ET_1_S' subject.ID '_block' num2str(iRun) '.mat']); 
        run.TimeFile     = fullfile(group.TimePaths{iSubj}, ['s' subject.ID '_run' num2str(iRun) '_timing.mat']);
        run.TimeBackup   = fullfile(folder.Time, ['s' subject.ID '_run' num2str(iRun) '_timing.mat']);
        
        %% Load participant responses
        load(run.BehaviorFile, 'stim_choice')
        run.TrialType     = {stim_choice.type}.';
        run.TrialOnset    = [stim_choice.start_time].';
        run.TrialDuration = [stim_choice.rt].';
        run.TrialResponse = [stim_choice.resp_num].';
        
        run.Pmax = [stim_choice.pmax].';
        run.Gmax = [stim_choice.gmax].';
        run.Lmin = [stim_choice.lmin].';
        
        clear stim_choice
        
        %% Create arrays for each condition with trial onset (start time) and duration (RT)
        run.RandLoc  = (contains(run.TrialType, 'rand') & ~isnan(run.TrialResponse));
        run.PmaxLoc  = run.Pmax == run.TrialResponse;
        run.GmaxLoc  = run.Gmax == run.TrialResponse;
        run.LminLoc  = run.Lmin == run.TrialResponse;
        
        run.RandOnset = run.TrialOnset(run.RandLoc);
        run.PmaxOnset = run.TrialOnset(run.PmaxLoc); 
        run.GmaxOnset = run.TrialOnset(run.GmaxLoc); 
        run.LminOnset = run.TrialOnset(run.LminLoc); 
        
        run.RandDuration = run.TrialDuration(run.RandLoc);
        run.PmaxDuration = run.TrialDuration(run.PmaxLoc); 
        run.GmaxDuration = run.TrialDuration(run.GmaxLoc); 
        run.LminDuration = run.TrialDuration(run.LminLoc);
        
        %% Set names for conditions 
        if any(isempty (run.RandOnset), 1)
            run.RandName = {};
        else
            run.RandName ='random';
        end

        if any(isempty (run.PmaxOnset), 1)
            run.PmaxName = {};
        else
            run.PmaxName ='pmax';
        end

        if any(isempty (run.GmaxOnset), 1)
            run.GmaxName = {};
        else
            run.GmaxName ='gmax';
        end

        if any(isempty (run.LminOnset), 1)
            run.LminName = {};
        else
            run.LminName = 'lmin';
        end
        
        %% Create final arrays for the run
        names     = {run.RandName, run.PmaxName, run. GmaxName, run.LminName};
        onsets    = {run.RandOnset, run.PmaxOnset, run. GmaxOnset, run.LminOnset};
        durations = {run.RandDuration, run.PmaxDuration, run. GmaxDuration, run.LminDuration};
        
        names     = names(~cellfun('isempty', names));
        onsets    = onsets(~cellfun('isempty', onsets));
        durations = durations(~cellfun('isempty', durations));
        
        %% Save the arrays
        save(run.TimeFile, 'names', 'onsets', 'durations')
        save(run.TimeBackup, 'names', 'onsets', 'durations')
        
    end

end

%% Return to scripts folder
cd(folder.Scripts)


