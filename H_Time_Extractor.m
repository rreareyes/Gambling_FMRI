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
        run.RandLoc = (contains(run.TrialType, 'rand') & ~isnan(run.TrialResponse));
              
        run.PmaxZeroLoc = contains(run.TrialType, 'trad_zero') & run.Pmax == run.TrialResponse;
        run.GmaxZeroLoc = contains(run.TrialType, 'trad_zero') & run.Gmax == run.TrialResponse;
        run.LminZeroLoc = contains(run.TrialType, 'trad_zero') & run.Lmin == run.TrialResponse;
        
        run.PmaxNegLoc = contains(run.TrialType, 'trad_neg') & run.Pmax == run.TrialResponse;
        run.GmaxNegLoc = contains(run.TrialType, 'trad_neg') & run.Gmax == run.TrialResponse;
        run.LminNegLoc = contains(run.TrialType, 'trad_neg') & run.Lmin == run.TrialResponse;        
        
        run.RandOnset = run.TrialOnset(run.RandLoc);
        
        run.PmaxZeroOnset = run.TrialOnset(run.PmaxZeroLoc); 
        run.GmaxZeroOnset = run.TrialOnset(run.GmaxZeroLoc); 
        run.LminZeroOnset = run.TrialOnset(run.LminZeroLoc); 
        
        run.PmaxNegOnset = run.TrialOnset(run.PmaxNegLoc); 
        run.GmaxNegOnset = run.TrialOnset(run.GmaxNegLoc); 
        run.LminNegOnset = run.TrialOnset(run.LminNegLoc); 
        
        run.RandDuration = run.TrialDuration(run.RandLoc);
        
        run.PmaxZeroDuration = run.TrialDuration(run.PmaxZeroLoc); 
        run.GmaxZeroDuration = run.TrialDuration(run.GmaxZeroLoc); 
        run.LminZeroDuration = run.TrialDuration(run.LminZeroLoc);
        
        run.PmaxNegDuration = run.TrialDuration(run.PmaxNegLoc); 
        run.GmaxNegDuration = run.TrialDuration(run.GmaxNegLoc); 
        run.LminNegDuration = run.TrialDuration(run.LminNegLoc);
        
        %% Set names for conditions 
        if any(isempty (run.RandOnset), 1)
            run.RandName = {};
        else
            run.RandName = 'random';
        end
        
        % Zero Trials
        if any(isempty (run.PmaxZeroOnset), 1)
            run.PmaxZeroName = {};
        else
            run.PmaxZeroName = 'pmax_zero';
        end

        if any(isempty (run.GmaxZeroOnset), 1)
            run.GmaxZeroName = {};
        else
            run.GmaxZeroName ='gmax_zero';
        end

        if any(isempty (run.LminZeroOnset), 1)
            run.LminZeroName = {};
        else
            run.LminZeroName = 'lmin_zero';
        end
        
        % Neg trials
        if any(isempty (run.PmaxNegOnset), 1)
            run.PmaxNegName = {};
        else
            run.PmaxNegName = 'pmax_neg';
        end

        if any(isempty (run.GmaxNegOnset), 1)
            run.GmaxNegName = {};
        else
            run.GmaxNegName ='gmax_neg';
        end

        if any(isempty (run.LminNegOnset), 1)
            run.LminNegName = {};
        else
            run.LminNegName = 'lmin_neg';
        end
        
        %% Create final arrays for the run
        names     = {run.RandName, run.PmaxZeroName, run.GmaxZeroName, run.LminZeroName, run.PmaxNegName, run. GmaxNegName, run.LminNegName};
        onsets    = {run.RandOnset, run.PmaxZeroOnset, run.GmaxZeroOnset, run.LminZeroOnset, run.PmaxNegOnset, run.GmaxNegOnset, run.LminNegOnset};
        durations = {run.RandDuration, run.PmaxZeroDuration, run.GmaxZeroDuration, run.LminZeroDuration, run.PmaxNegDuration, run.GmaxNegDuration, run.LminNegDuration};
        
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


