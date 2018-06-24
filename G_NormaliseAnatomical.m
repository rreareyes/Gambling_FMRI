%% Anatomical Normaliser

% This script normalises the anatomical nii file
% It assumes that your image is in a folder with mprage on its name
% The output is an image with the "norm_" prefix

%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% August 2017
% working on SPM8

%% Clean workspace
clc; clear

%% Base Paths
cd('..')
folder.Root      = pwd;
folder.Processed = fullfile(folder.Root,'Processed');
folder.Scripts   = fullfile(folder.Root, 'Scripts'); %location of the scripts

%% Get all subject folders
folder.ProcessedPaths      = dir(folder.Processed);
folder.ProcessedPaths(1:2) = [];

%% Ask for which subjects to run
[options.Group,~] = listdlg('ListString',{'Individual Elements','All Subjects'},'Name','N° Subjects to Process?');

%% Set subject paths for input and output images according to subjects selected
if options.Group == 1 % Customized list
    [options.Subjects, ~] = listdlg('ListString',char(folder.ProcessedPaths.name),'Name','Which subjects do you want?');
    group.SubjectsList    = folder.ProcessedPaths(options.Subjects);

elseif options.Group == 2 % All subjects
    group.SubjectsList = folder.ProcessedPaths;

end

for iFolder = 1:size(group.SubjectsList,1)
    group.SubjectsPaths{iFolder,:} = fullfile(folder.Processed, group.SubjectsList(iFolder).name);
end

%% Ask for spm folder location
folder.Spm = uigetdir ('','Please select the SPM folder');

%% Loop the normalization through all selected subjects
for iSubj = 1:length(group.SubjectsPaths)
    %% Clear values to avoid overwritting and show current loop
    clear subject
    currentSubj = ['Working on ' group.SubjectsList(iSubj).name];
    disp(currentSubj)
    
    %% Get anatomic image an data
    subject.AnatFolder = fullfile(group.SubjectsPaths{iSubj}, 'anat');
    subject.ImageList  = dir(fullfile(subject.AnatFolder, 's*.nii'));
    subject.ImagePath  = fullfile(subject.AnatFolder, subject.ImageList.name);
    subject.MprageMat  = dir(fullfile(subject.AnatFolder, 's*_seg_sn.mat')); %reference file
    subject.MatPath    = fullfile(subject.AnatFolder, subject.MprageMat.name);
    
    %% Run SPM Normalization
    spm('defaults','fmri');
    spm_jobman('initcfg');
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.matname      = {subject.MatPath}; 
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample     = {subject.ImagePath}; %the anatomical file we are going to warp
    matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.bb       = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.roptions.vox      = [1 1 1]; %changed voxel size to fit original resolution
    matlabbatch{1}.spm.spatial.normalise.write.roptions.interp   = 1;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap     = [0 0 0];
    matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix   = 'norm_';

    %% Run the batch
    spm_jobman('run',matlabbatch);
  
end

%% Get back to script folder
cd(folder.Scripts)
