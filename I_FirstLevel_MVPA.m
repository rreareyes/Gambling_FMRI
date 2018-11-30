%% SPM FIRST LEVEL ANALYSIS
% This script does 1st level design analysis on smoothed nii files
% It considers that the files have the "sm_" prefix

% It uses the timing outputs from the TimeExtractor.m script to get onsets
% names and durations.

% It offers the option of running all subjects at once or specific
% subjects, which can be selected from a list dialog.

% It also allows to model either odd runs (1, 3, 5) or even runs (2, 4)
% Or to create both models one after the other

% Most variables are stored as fields in the following structures:
% folder: contains the base directories names and paths

% group: mostly has list of files, folders and paths from the group of
%        subjects selected to be processed

% subject: holds the information from each individual iteration of the
%          loops through all the subjects in the selected group

% It assumes the following organization from the source files

%   Processed
%       Subject
%           func ----> Functional 4D nii images (sm_sXX_Gamble_run-X.nii)
%           anat ----> Anatomical T1 mprage image (sXX_T1w.nii)

% And it creates the folders
%   1stLevel
%       AllConditions
%           Subject
%               beta_XXXX
%               mask
%               RPV
%               ResMS
%               SPM.mat
%       Responses
%           Subject
%               beta_XXXX
%               mask
%               RPV
%               ResMS
%               SPM.mat

%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% November 2018
% working on SPM8

%% Clean workspace
clc; clear

%% Base Paths
cd('..')
folder.Root       = pwd;
folder.Processed  = fullfile(folder.Root, 'Processed');
folder.Scripts    = fullfile(folder.Root, 'Scripts'); 
folder.Results{1} = fullfile(folder.Root, 'MVPA', 'Odd', '1stLevel');
folder.Results{2} = fullfile(folder.Root, 'MVPA', 'Even', '1stLevel');

%% Get all subject paths
folder.ProcessedPaths      = dir(folder.Processed);
folder.ProcessedPaths(1:2) = [];

%% Define your group
%% Ask for which subjects to run
[options.Group, ~] = listdlg('ListString', {'Individual Elements', 'All Subjects'}, 'Name', 'No. Subjects to Process?');

%% Set subject paths for input and output images according to subjects selected
if options.Group == 1 % Customized list
    [options.Subjects, ~] = listdlg('ListString', char(folder.ProcessedPaths.name), 'Name', 'Which subjects do you want?');
    group.SubjectsList    = folder.ProcessedPaths(options.Subjects);

elseif options.Group == 2 % All subjects
    group.SubjectsList = folder.ProcessedPaths;
end

%% Ask which runs you want to model
[options.Runs, ~] = listdlg('ListString',{'Odd runs', 'Even Runs', 'Even and Odd'}, 'Name','Runs to model?');

for iFolder = 1:size(group.SubjectsList,1)
    group.SubjectsPaths{iFolder, :} = fullfile(folder.Processed, group.SubjectsList(iFolder).name);
    group.ResultPaths{iFolder, 1}   = fullfile(folder.Results{1}, group.SubjectsList(iFolder).name);
    group.ResultPaths{iFolder, 2}   = fullfile(folder.Results{2}, group.SubjectsList(iFolder).name);
end

%% If both odd and even runs models will be created
if options.Runs == 3
   
    options.nModel = 2;
else
    options.nModel = 1;
    
end

%% Create the models requested
for iModel = 1:options.nModel

    if options.Runs == 3
        options.Model = iModel; %first iteration creates odd model and second even model.

    else
        options.Model = options.Runs;

    end    
    
    %% Loop through the subjects to get onsets and perform 1st level analysis
    for iSubject = 1:size(group.SubjectsPaths, 1)
        %% Clear values to avoid overwriting issues
        clear subject matlabbatch

        %% Set folder paths to load files and save results
        subject.ResultPath = group.ResultPaths{iSubject, options.Model};
        subject.SourcePath = group.SubjectsPaths{iSubject, :};

        if ~exist(subject.ResultPath,'dir')
            mkdir(subject.ResultPath);
        
        end

        %% Get all files from functional runs
        subject.FuncFolder = fullfile(subject.SourcePath, 'func');
        subject.FuncList   = dir(fullfile(subject.FuncFolder, 'sm_*.nii'));

        %% Get the multiple motion regressors from the text file inside
        subject.RegressorsFilesList = dir(fullfile(subject.FuncFolder, 'rp_*.txt'));

        %% Get the timing files path from their storage folder to load in each iteration
        subject.TimingPath    = fullfile(subject.SourcePath, 'timing');
        subject.TimeFilesList = dir(fullfile(subject.TimingPath, '*.mat'));

        %% Set the runs to be modeled
        if options.Model == 1 % Odd Runs
            options.RunsList = [1 3 5];

        elseif options.Model == 2 % Even Runs
            options.RunsList = [2 4];

        end
        
        subject.TimeFilesList = subject.TimeFilesList(options.RunsList, :);
        subject.RegressorsFilesList = subject.RegressorsFilesList(options.RunsList, :);
        subject.FuncList = subject.FuncList(options.RunsList, :);
        
        %% Get the absolute path to each file and store them in {Run x 1} cell matrices
        for iRun = 1:size(options.RunsList, 2)
            % Multiple regressors files 
            subject.RegressorsPaths{iRun,:} = fullfile(subject.FuncFolder, subject.RegressorsFilesList(iRun).name);

            % Timing files
            subject.TimePaths{iRun,:} = fullfile(subject.TimingPath, subject.TimeFilesList(iRun).name);

            % Paths for each functional image volume
            for iVolume = 1:200
                subject.FuncVolumesPaths{iRun,:}{iVolume,:} = fullfile(subject.FuncFolder, [subject.FuncList(iRun).name ',' num2str(iVolume)]);
            end

        end

        %% FIRST LEVEL ANALYSIS
        % Initial set up
        spm('defaults','fmri');
        spm_jobman('initcfg');

        % Set base parameters
        matlabbatch{1}.spm.stats.fmri_spec.dir            = {subject.ResultPath}; %Directory where spm folder will be saved
        matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = 1.5;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

        %% PARAMETERS TO MODEL EXPERIMENTAL CONDITIONS (event-related to each decision making event)
        % Loop through all runs to establish the individual parameters for each
        for iRun = 1:size(subject.FuncVolumesPaths,1)
            % Indicate absolute paths to functional images, the mat files with timing, and the multiple motion regression file
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).scans     = subject.FuncVolumesPaths{iRun}; %data files
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond      = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi     = subject.TimePaths(iRun); %timing file
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).regress   = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi_reg = subject.RegressorsPaths(iRun); %multiple regressor file
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).hpf       = 128;
        
        end

        matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';

        %% Model Estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1)                      = cfg_dep;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname                = 'Select SPM.mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname                = 'fMRI model specification: SPM.mat File';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch         = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output           = substruct('.','spmmat');
        matlabbatch{2}.spm.stats.fmri_est.method.Classical               = 1;

        %% Run the batch
        spm_jobman('run',matlabbatch);

    end % End of subject loop
    
end % End of models loop

%% Return to scripts folder
cd(folder.Scripts)
