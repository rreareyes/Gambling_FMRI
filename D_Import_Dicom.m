%% SPM Dicom importer

% This script allows to import the data from the gamble fmri Dicom files

% You require to create a folder called "DICOM" and store there the folders
% with your subject's scans.

% Also, it requires having the script files in a folder called "Scripts"

% You can modify this if you want, just go to the section "Base Paths" and
% change to your preferred folder names if necessary.

% It uses the original organization of container folders from the scanner:
%   DICOM
%       sXX
%           ProjectCode
%               Session(Run)
%                   Images.DICOM

% It creates a folder per subject with the following organization:
%   NII
%       sXX
%           Session(Run)
%               Images.nii

%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% June 2018
% working on SPM8

%% Clean workspace
clc; clear

%% Base Paths
cd('..')
folder.Root    = pwd;
folder.Dicom   = fullfile(folder.Root, 'DICOM');
folder.Nii     = fullfile(folder.Root, 'NII');
folder.Scripts = fullfile(folder.Root, 'Scripts'); %location of the scripts

%% Get all subject folders
folder.DicomPaths      = dir(folder.Dicom);
folder.DicomPaths(1:2) = [];

%% Ask for which subjects to run
[options.Group, ~] = listdlg('ListString', {'Individual Elements', 'All Subjects'}, 'Name', 'No. Subjects to Process?','ListSize',[280,250]);

%% Set subject paths for input and output images according to subjects selected
if options.Group == 1 % Customized list
    [options.Subjects, ~] = listdlg('ListString', char(folder.DicomPaths.name), 'Name','Which subjects do you want to process?', 'ListSize',[350,250]);
    group.DicomList  = folder.DicomPaths(options.Subjects);

elseif options.Group == 2 % All subjects
    group.DicomList = folder.DicomPaths;

end

for iFolder = 1:size(group.DicomList,1)
    group.DicomPaths{iFolder,:} = fullfile(folder.Dicom, group.DicomList(iFolder).name);
    group.NiiPaths{iFolder,:}   = fullfile(folder.Nii, group.DicomList(iFolder).name);
end

for iSubj = 1:length(group.DicomPaths)
    %% Clear values to avoid overwritting and show current loop
    clear subject
    currentSubj = ['Working on ' group.DicomList(iSubj).name];
    disp(currentSubj)
    
    %% Get inside each participant's folder
    subject.ImageList  = dir([group.DicomPaths{iSubj} filesep '**' filesep '*.IMA']);
    subject.ImagePaths = fullfile({subject.ImageList.folder}, {subject.ImageList.name}).';
    
    if ~exist(group.NiiPaths{iSubj, :}, 'dir')
        mkdir(group.NiiPaths{iSubj, :});
    end
    
    %% Import all IMA files and convert to Nii
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    
    matlabbatch{1}.spm.util.dicom.data             = subject.ImagePaths; 
    matlabbatch{1}.spm.util.dicom.root             = 'series';
    matlabbatch{1}.spm.util.dicom.outdir           = group.NiiPaths(iSubj, :);
    matlabbatch{1}.spm.util.dicom.convopts.format  = 'nii';
    matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;
 
    spm_jobman('run',matlabbatch);

end

%% Return to script folder
cd(folder.Scripts)

