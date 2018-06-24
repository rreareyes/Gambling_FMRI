%% PREPROCESSING

% This script allows to preprocess the imported files as follows:
%-Timing correction ----> Prefix "tr_"
%-Realignment       ----> Prefix "r_"
%-Coregistration
%-Segmentation   
%-Normalization     ----> Prefix "norm_"
%-Smoothing         ----> Prefix "sm_"

% It assumes that all the nii images are organized as follows
%   Processed
%       Subject
%           func ----> Functional 4D nii images (sXX_Gamble_run-X.nii)
%           anat ----> Anatomical T1 mprage image (sXX_T1w.nii)

% The resulting images from each step are saved using the same name from
% the original, and adding the prefix indicated before.

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
folder.Scripts   = fullfile(folder.Root, 'Scripts'); %location of the scripts

%% Get all subject folders
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
    group.SubjectsPaths{iFolder,:} = fullfile(folder.Processed, group.SubjectsList(iFolder).name);
end

%% Ask for spm folder location
folder.Spm = uigetdir ('','Please select the SPM folder');

for iSubj = 1:length(group.SubjectsPaths) 
    %% Clear values to avoid overwritting and show current loop
    clear subject
    currentSubj = ['Working on ' group.SubjectsList(iSubj).name];
    disp(currentSubj)
    
    %% Get paths for runs
    subject.FuncFolder = fullfile(group.SubjectsPaths{iSubj}, 'func');
    subject.FuncList   = dir(fullfile(subject.FuncFolder, '*.nii'));
    
    for iRun = 1:5
        for iVolume = 1:200
           subject.FuncPaths{1 ,iRun}{iVolume,:} = fullfile(subject.FuncFolder, [subject.FuncList(iRun).name ',' num2str(iVolume)]);
        end
    end
    
    subject.AnatFolder = fullfile(group.SubjectsPaths{iSubj}, 'anat');
    subject.AnatList   = dir(fullfile(subject.AnatFolder, '*.nii'));
    subject.AnatPath   = fullfile(subject.AnatFolder, subject.AnatList.name);
    
    %% Run preprocessing
    spm('defaults','FMRI');
    spm_jobman('initcfg');

    %% Time correction    
    matlabbatch{1}.spm.temporal.st.scans    = subject.FuncPaths; %files for each run
    matlabbatch{1}.spm.temporal.st.nslices  = 44;
    matlabbatch{1}.spm.temporal.st.tr       = 1.5;
    matlabbatch{1}.spm.temporal.st.ta       = 1.46590909090909;
    matlabbatch{1}.spm.temporal.st.so       = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44];
    matlabbatch{1}.spm.temporal.st.refslice = 22;
    matlabbatch{1}.spm.temporal.st.prefix   = 'tc_';

    %% Realign
    for iRun = 1:5
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1)                      = cfg_dep;
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).tname                = 'Session';
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).tgt_spec{1}(1).value = 'image';
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).sname                = ['Slice Timing: Slice Timing Corr. Images (Sess ' num2str(iRun) ')'];
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).src_exbranch         = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{2}.spm.spatial.realign.estwrite.data{iRun}(1).src_output           = substruct('()',{iRun}, '.','files');
    end
    
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight  = '';
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which   = [2 1];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp  = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask    = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix  = 'r_';

    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight  = '';
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which   = [2 1];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp  = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask    = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix  = 'r_';   
    
    %% Coregistration
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1)                      = cfg_dep;
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).tname                = 'Reference Image';
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(1).name  = 'filter';
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(2).name  = 'strtype';
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).sname                = 'Realign: Estimate & Reslice: Mean Image';
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).src_exbranch         = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).src_output           = substruct('.','rmean');
    matlabbatch{3}.spm.spatial.coreg.estimate.source                      = {[subject.AnatPath ',1']}; % source on the anatomical folder
    matlabbatch{3}.spm.spatial.coreg.estimate.other                       = {''};
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun           = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep                = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol                = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm               = [7 7];
    
    %% Segmentation
    matlabbatch{4}.spm.spatial.preproc.data(1)                      = cfg_dep;
    matlabbatch{4}.spm.spatial.preproc.data(1).tname                = 'Data';
    matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(1).name  = 'filter';
    matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(2).name  = 'strtype';
    matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.spatial.preproc.data(1).sname                = 'Coregister: Estimate: Coregistered Images';
    matlabbatch{4}.spm.spatial.preproc.data(1).src_exbranch         = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{4}.spm.spatial.preproc.data(1).src_output           = substruct('.','cfiles');
    matlabbatch{4}.spm.spatial.preproc.output.GM                    = [0 0 1];
    matlabbatch{4}.spm.spatial.preproc.output.WM                    = [0 0 1];
    matlabbatch{4}.spm.spatial.preproc.output.CSF                   = [0 0 0];
    matlabbatch{4}.spm.spatial.preproc.output.biascor               = 1;
    matlabbatch{4}.spm.spatial.preproc.output.cleanup               = 0;
    matlabbatch{4}.spm.spatial.preproc.opts.tpm                     = {fullfile(folder.Spm, 'tpm', 'grey.nii');fullfile(folder.Spm, 'tpm', 'white.nii');fullfile(folder.Spm, 'tpm', 'csf.nii')};
    matlabbatch{4}.spm.spatial.preproc.opts.ngaus                   = [2; 2; 2; 4];
    matlabbatch{4}.spm.spatial.preproc.opts.regtype                 = 'mni';
    matlabbatch{4}.spm.spatial.preproc.opts.warpreg                 = 1;
    matlabbatch{4}.spm.spatial.preproc.opts.warpco                  = 25;
    matlabbatch{4}.spm.spatial.preproc.opts.biasreg                 = 0.0001;
    matlabbatch{4}.spm.spatial.preproc.opts.biasfwhm                = 60;
    matlabbatch{4}.spm.spatial.preproc.opts.samp                    = 3;
    matlabbatch{4}.spm.spatial.preproc.opts.msk                     = {''};
    
    %% Normalization
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1)                      = cfg_dep;
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).tname                = 'Parameter File';
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).name  = 'filter';
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).name  = 'strtype';
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).sname                = 'Segment: Norm Params Subj->MNI';
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).src_exbranch         = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).src_output           = substruct('()',{1}, '.','snfile', '()',{':'});

    for iRun = 1:5
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun)                      = cfg_dep;
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).tname                = 'Images to Write';
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).tgt_spec{1}(1).value = 'image';
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).tgt_spec{1}(2).value = 'e';
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).sname                = ['Realign: Estimate & Reslice: Resliced Images (Sess ' num2str(iRun) ')'];
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).src_exbranch         = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(iRun).src_output           = substruct('.','sess', '()',{iRun}, '.','rfiles');
    end

    matlabbatch{5}.spm.spatial.normalise.write.roptions.preserve = 0;
    matlabbatch{5}.spm.spatial.normalise.write.roptions.bb       = [-78 -112 -70; 78 76 85];
    matlabbatch{5}.spm.spatial.normalise.write.roptions.vox      = [3 3 3]; %changed voxel size to fit original resolution
    matlabbatch{5}.spm.spatial.normalise.write.roptions.interp   = 1;
    matlabbatch{5}.spm.spatial.normalise.write.roptions.wrap     = [0 0 0];
    matlabbatch{5}.spm.spatial.normalise.write.roptions.prefix   = 'norm_';

    %% Smoothing
    matlabbatch{6}.spm.spatial.smooth.data(1)                      = cfg_dep;
    matlabbatch{6}.spm.spatial.smooth.data(1).tname                = 'Images to Smooth';
    matlabbatch{6}.spm.spatial.smooth.data(1).tgt_spec{1}(1).name  = 'filter';
    matlabbatch{6}.spm.spatial.smooth.data(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{6}.spm.spatial.smooth.data(1).tgt_spec{1}(2).name  = 'strtype';
    matlabbatch{6}.spm.spatial.smooth.data(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{6}.spm.spatial.smooth.data(1).sname                = 'Normalise: Write: Normalised Images (Subj 1)';
    matlabbatch{6}.spm.spatial.smooth.data(1).src_exbranch         = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{6}.spm.spatial.smooth.data(1).src_output           = substruct('()',{1}, '.','files');
    matlabbatch{6}.spm.spatial.smooth.fwhm                         = [8 8 8];
    matlabbatch{6}.spm.spatial.smooth.dtype                        = 0;
    matlabbatch{6}.spm.spatial.smooth.im                           = 0;
    matlabbatch{6}.spm.spatial.smooth.prefix                       = 'sm_';
    
    %% Run the batch
    spm_jobman('run',matlabbatch);

end

%% Return to script folder
cd(folder.Scripts)
