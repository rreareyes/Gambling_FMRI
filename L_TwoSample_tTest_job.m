%% TWO SAMPLE T TEST
% This script runs a group level analysis on the contrast images generated
% in the first level processing. It estimates a two sample t test comparing
% adolescents and adults, and creates the respective contrasts.

% It gives you the option of running a specific contrast or all at the same
% time.

% It considers that the files have the "con_" prefix and are organized as
% follows:

%   1stLevel
%       AllConditions
%           Subject
%               con_XXXX
%               spmT_ XXXX

% The output from this script is saved in the following way:

%   2ndLevel
%       AultsVSAdo ----------------> Experimental conditions
%           Condition --------> Particular condition
%               con_XXXX
%               spmT_ XXXX
%               SPM.mat

%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% June 2018
% working on SPM8

%% Clean workspace
clc; clear

%% Base paths
cd('..')
folder.Root      = pwd;
folder.Scripts   = fullfile (folder.Root, 'Scripts');
folder.Results   = fullfile (folder.Root, '1stLevel', 'AllConditions');
folder.SecondLvl = fullfile (folder.Root, '2ndLevel');

if ~exist(folder.SecondLvl,'dir')
    mkdir(folder.SecondLvl);
end

%% Names for the options
contrast.Pmax = {'pmax', 'pmax > rand', 'pmax > gmax', 'pmax > lmin', 'pmax > non pmax'};
contrast.Gmax = {'gmax', 'gmax > rand', 'gmax > pmax', 'gmax > lmin', 'gmax > non gmax'};
contrast.Lmin = {'lmin', 'lmin > rand', 'lmin > pmax', 'lmin > gmax', 'lmin > non lmin'};
contrast.Rand = {'rand','rand > all', 'all > rand'};
contrast.Mix  = {'non pmax > pmax', 'non gmax > gmax', 'non lmin > lmin', 'non pmax', 'structured', 'nongmax', 'non_lmin', 'non_pmax > random', 'non_gmax > random'};

contrast.Options = horzcat(contrast.Pmax, contrast.Gmax, contrast.Lmin, contrast.Rand, contrast.Mix).';

%% Names for the contrasts in the file
contrast.PmaxFile = {'pmax', 'pmaxVSrand', 'pmaxVSgmax', 'pmaxVSlmin', 'pmaxVSnon_pmax'};
contrast.GmaxFile = {'gmax', 'gmaxVSrand', 'gmaxVSpmax', 'gmaxVSlmin', 'gmaxVSnon_gmax'};
contrast.LminFile = {'lmin', 'lminVSrand', 'lminVSpmax', 'lminVSgmax', 'lminVSnon_lmin'};
contrast.RandFile = {'rand','randVSall', 'allVSrand'};
contrast.MixFile  = {'non_pmaxVSpmax', 'non_gmaxVSgmax', 'non_lminVSlmin', 'non_pmax', 'structured', 'nongmax', 'non_lmin', 'non_pmaxVSrandom', 'non_gmaxVSrandom'};

contrast.Files = horzcat(contrast.PmaxFile, contrast.GmaxFile, contrast.LminFile, contrast.RandFile, contrast.MixFile).';

%% Ask if you want to create all contrast, an individual contrast for all subjects, or all contrast for a given group
[contrast.Number,~] = listdlg('ListString',{'Single Contrast','All Contrasts'},'Name','Analysis?');

if contrast.Number == 1 %Create single contrast
    % Ask which contrast from the list you want to create
    [contrast.Selection,~] = listdlg('ListString',contrast.Options,'Name','Select the contrast you want to analyze');
    
    contrast.Target = contrast.Files{contrast.Selection}; % Which individual contrast will be generated
    contrast.List   = {contrast.Target}; % Make this contrast the only one in the list (so the for loop iterates just one time
else
    % Create all the contrast (the for loop iterates for n number of contrast in the list)
    contrast.List = contrast.Files;
end

%% Get folders specific for every Group
folder.List        = dir(folder.Results);
folder.List(1:2,:) = []; % first rows of list are meaningless ([], ., ...]
folder.List(contains({folder.List.name},'100') | contains({folder.List.name},'120') | contains({folder.List.name},'130')) = []; % remove subjects with braces

% Locate Adolescents and Adults
folder.ReferenceAdolescent = regexp({folder.List.name}.','s..\d'); % Adolescents have 3 digits as ID: results_sXXX
folder.MapAdult            = cellfun('isempty',folder.ReferenceAdolescent); % location of non Adolescents
folder.MapAdolescent       = folder.MapAdult == 0;

% Create the folder list
folder.Adolescent = folder.List(folder.MapAdolescent,:);
folder.Adult      = folder.List(folder.MapAdult,:);

%% Get the folder paths for both groups
for iSubj = 1:length(folder.Adult)
    folder.AdultPath(iSubj,:) = fullfile(folder.Results,folder.Adult(iSubj).name);
end

for iSubj = 1:length(folder.Adolescent)
    folder.AdolescentPath(iSubj,:) = fullfile(folder.Results,folder.Adolescent(iSubj).name);
end

%% ITERATE OVER ALL THE CON FILES TO CREATE 2ND LEVEL CONTRASTS
for iContrast = 1:length(contrast.List)
    
   %% Assign folder for the results (create if it doesn't exist)
   folder.Contrast = fullfile(folder.SecondLvl, 'AdultsVSAdo', contrast.List{iContrast}); %2ndLevel/ContrastName
        if ~exist(folder.Contrast,'dir')
            mkdir(folder.Contrast);
        end
     
    %% Get the contrast images
    if contrast.Number == 2 
        contrast.Target = contrast.Files{iContrast}; % Creates one contrast from the list in every iteration
    end

    for iFolder = 1:size(folder.AdultPath,1) %% Using just the subjects in the group 
        contrast.AdultImages = dir([folder.AdultPath(iFolder,:) filesep 'con_' contrast.Target '.img']); % Look for all images in the folder list, matching the contrast name (search for ContrastName.img)
        guide = ~isempty(contrast.AdultImages); % Check if that contrast exists for a particular subject

        if  guide == 1 %If it exists, add it to the list
            contrast.AdultPaths{iFolder,:} = fullfile(contrast.AdultImages.folder, contrast.AdultImages.name);
        end
    end
    
    for iFolder = 1:size(folder.AdolescentPath,1) %% Using just the subjects in the group 
        contrast.AdolescentImages = dir([folder.AdolescentPath(iFolder,:) filesep 'con_' contrast.Target '.img']); % Look for all images in the folder list, matching the contrast name (search for ContrastName.img)
        guide = ~isempty(contrast.AdolescentImages); % Check if that contrast exists for a particular subject

        if  guide == 1 %If it exists, add it to the list
            contrast.AdolescentPaths{iFolder,:} = fullfile(contrast.AdolescentImages.folder, contrast.AdolescentImages.name);
        end
    end
    
    %% Remove empty cells from the paths
    contrast.AdultPaths      = contrast.AdultPaths(~cellfun('isempty',contrast.AdultPaths));
    contrast.AdolescentPaths = contrast.AdolescentPaths(~cellfun('isempty',contrast.AdolescentPaths));
    
    %% Save the images used in each 2nd level contrast    
    file.Summary = fullfile(folder.Contrast, 'imgsIncluded');
    imgsIncluded = [{contrast.AdultPaths}, {contrast.AdolescentPaths}];
    save(file.Summary, 'imgsIncluded')
        
    try
        %% Factorial design
        spm('defaults','fmri');
        spm_jobman('initcfg');

        matlabbatch{1}.spm.stats.factorial_design.dir                    = {folder.Contrast};
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1          = contrast.AdultPaths;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2          = contrast.AdolescentPaths;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.dept            = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.variance        = 1; %Assume unequal variance
        matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca           = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova          = 0;
        matlabbatch{1}.spm.stats.factorial_design.cov                    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;
    
        %% Model Estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1)                      = cfg_dep;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname                = 'Select SPM.mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname                = 'Factorial design specification: SPM.mat File';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch         = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output           = substruct('.','spmmat');
        matlabbatch{2}.spm.stats.fmri_est.method.Classical               = 1;

        %% Write the contrasts in the SPM file
        matlabbatch{3}.spm.stats.con.spmmat(1)                      = cfg_dep;
        matlabbatch{3}.spm.stats.con.spmmat(1).tname                = 'Select SPM.mat';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{3}.spm.stats.con.spmmat(1).sname                = 'Model estimation: SPM.mat File';
        matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch         = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{3}.spm.stats.con.spmmat(1).src_output           = substruct('.','spmmat');
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name           = 'Adults > Adolescents';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec         = [1 -1];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name           = 'Adolescents > Adults';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec         = [-1 1];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep        = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;

        %% Run the batch
        spm_jobman('run',matlabbatch);

    catch % handle behavior if an error is found executing Try statements
        %% Show and error message and store it for reference
        contrast.errors{iContrast} = ['No significant voxels found in the contrast: ' contrast.List{iContrast}];
        contrast.message            = errordlg(contrast.errors{iContrast},'Model not estimated');
        
        %% Clear path names to avoid mixing names with previous iterations
        if contrast.Number == 2 
            contrast = rmfield(contrast,'AdultPaths');
            contrast = rmfield(contrast,'AdultImages');

            contrast = rmfield(contrast,'AdolescentPaths');
            contrast = rmfield(contrast,'AdolescentImages');
        end
        
        %% Continue to the nex iteration of the loop
        continue; 
        
    end
    
    %% Clear path names to avoid mixing names with previous iterations
    if contrast.Number == 2 
        contrast = rmfield(contrast,'AdultPaths');
        contrast = rmfield(contrast,'AdultImages');

        contrast = rmfield(contrast,'AdolescentPaths');
        contrast = rmfield(contrast,'AdolescentImages');
    end
end

%% Get back to the starting point
cd(folder.Scripts)
