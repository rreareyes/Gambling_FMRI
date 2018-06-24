%% SPM 4D Creator

% This script creates 4D nii images for each functional run. This allows us
% to fuse all the 3D volumes in a single nii image (for logistic purposes)

% It assumes the following organization from the nii images
%   NII
%       sXX
%           SessionX
%               Images.nii
%           SessionX
%               Images.nii
%           SessionX
%               Images.nii
%           ...

% The resulting images will be saved in the folder "Processed"
% The funcional images will be saved in the folder "func"
% The mprage image is saved in the folder "anat"

% This results in the following organization
%   Processed
%       sXX
%           func   ----> Functional 4D nii images (sXX_Gamble_run-X.nii)
%           anat   ----> Anatomical T1 mprage image (sXX_T1w.nii)

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
folder.Nii       = fullfile(folder.Root,'NII');
folder.Processed = fullfile(folder.Root,'Processed');
folder.Scripts   = fullfile(folder.Root, 'Scripts'); %location of the scripts

%% Get all subject folders
folder.NiiPaths      = dir(folder.Nii);
folder.NiiPaths(1:2) = [];

%% Ask for which subjects to run
[options.Group, ~] = listdlg('ListString',{'Individual Elements','All Subjects'},'Name','No. Subjects to Process?');

%% Set subject paths for input and output images according to subjects selected
if options.Group == 1 % Customized list
    [options.Subjects, ~] = listdlg('ListString',char(folder.NiiPaths.name),'Name','Which subjects do you want?');
    group.NiiList    = folder.NiiPaths(options.Subjects);

elseif options.Group == 2 % All subjects
    group.NiiList = folder.NiiPaths;

end

for iFolder = 1:size(group.NiiList,1)
    group.NiiPaths{iFolder,:}       = fullfile(folder.Nii, group.NiiList(iFolder).name);
    group.ProcessedPaths{iFolder,:} = fullfile(folder.Processed, group.NiiList(iFolder).name);
end

for iSubj = 1:length(group.NiiPaths)     
    clear subject
    currentSubj = ['Working on ' group.NiiList(iSubj).name];
    disp(currentSubj)
    
    %% Set the final directories
    subject.Functional = fullfile(group.ProcessedPaths{iSubj,:}, 'func');
    subject.Anatomical = fullfile(group.ProcessedPaths{iSubj,:}, 'anat');
    
    if ~exist(subject.Functional,'dir')
        mkdir(subject.Functional);
    end
    
    if ~exist(subject.Anatomical,'dir')
        mkdir(subject.Anatomical);
    end
    
    %% Get a list of the relevant folders (functional runs and anatomical)
    subject.FolderList     = dir(group.NiiPaths{iSubj,:});
    subject.RunsLocation   = contains({subject.FolderList.name}, 'run').';
    subject.RunsFolderList = subject.FolderList(subject.RunsLocation);
    subject.RunsImagesList = dir(fullfile(group.NiiPaths{iSubj}, 'epi_bold_run*', '*.nii'));
    subject.RunsImagesPath = fullfile({subject.RunsImagesList.folder},{subject.RunsImagesList.name}).';
    
    subject.AnatLocation = contains({subject.FolderList.name}, 'mprage').';
    subject.AnatList     = subject.FolderList(subject.AnatLocation);
    subject.AnatFile     = dir(fullfile(subject.AnatList.folder, subject.AnatList.name, '*.nii'));
    subject.AnatPath     = fullfile(subject.AnatFile.folder, subject.AnatFile.name);
    
    %% Create 4D images for each run and save them to the func folder
    spm('defaults','FMRI');
    spm_jobman('initcfg');
    
    for iRun = 1:5
        subject.Volumes     = dir(fullfile(subject.RunsFolderList(iRun).folder, subject.RunsFolderList(iRun).name, '*.nii'));
        subject.RunPath     = fullfile({subject.Volumes.folder}, {subject.Volumes.name}).';
        subject.Image4dPath = fullfile(subject.Functional, [group.NiiList(iSubj).name '_Gamble_run-' num2str(iRun) '.nii']);
        
        %% Set images to cast 
        matlabbatch{iRun}.spm.util.cat.vols  = subject.RunPath; %#ok<SAGROW>
        matlabbatch{iRun}.spm.util.cat.name  = subject.Image4dPath; %#ok<SAGROW>
        matlabbatch{iRun}.spm.util.cat.dtype = 4; %#ok<SAGROW>
    end
        
    spm_jobman('run',matlabbatch);
    
    %% Copy the T1 mprage file to the anat folder
    subject.AnatTarget = fullfile(subject.Anatomical, [group.NiiList(iSubj).name '_T1w.nii']);
    copyfile(subject.AnatPath, subject.AnatTarget)

    
end

%% Return to script folder
cd(folder.Scripts)


