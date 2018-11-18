%% ONE SAMPLE T TEST
% This script runs a group level analysis on the contrast images generated
% in the first level processing. It estimates a one sample t test at the
% group level and creates the respective contrasts.

% It gives you the option of generating contrast to analyze all
% participants or groups (adolescents or adults) in the different
% experimental conditions, run a specific contrast for all subjects or
% model just the button responses. For the experimental conditions it also
% allows to incorporate covariates in the model (WASI and Age)

% It considers that the files have the "con_" prefix and are organized as
% follows:

%   1stLevel
%       AllConditions
%           Subject
%               con_XXXX
%               spmT_ XXXX
%       Responses
%           Subject
%               con_XXXX
%               spmT_ XXXX

% The output from this script is saved in the following way:

%   2ndLevel
%       AllConditions ----------------> Experimental conditions
%           Type ---------------------> Including or not covariates
%               Group ----------------> Adolescents, Adults or Everyone
%                   Condition --------> Particular condition
%                     con_XXXX
%                     spmT_ XXXX
%                     SPM.mat
%       Responses
%           con_XXXX
%           spmT_ XXXX
%           SPM.mat


%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% June 2018
% working on SPM8

%% Clean workspace
clc; clear

%% Base paths
cd('..')
folder.Root       = pwd;
folder.Scripts    = fullfile (folder.Root, 'Scripts');
folder.Behavior   = fullfile (folder.Root, 'Behavioral', 'Results');
folder.SecondLvl  = fullfile (folder.Root, '2ndLevel');

if ~exist(folder.SecondLvl,'dir')
    mkdir(folder.SecondLvl);
end

%% Ask the type of analysis to perform
[options.ContrastGroup, ~] = listdlg('ListString',{'All Contrasts','Single Contrast','Adult All Contrasts','Adolescent All Contrasts', 'Button Response'},'Name','Select the analysis', 'ListSize',[250,400]);

%% Set specifics for the group analysis using ondly button responses
if options.ContrastGroup == 5
    % Set the source folder for the con images
    folder.Results = fullfile (folder.Root, '1stLevel', 'Responses');
    
    %% Get all the participants folders
    folder.List = dir(folder.Results);
    folder.List(1:3,:) = []; % first rows of list are meaningless ([], ., ...]

    %% Get IDs and remove incomplete subjects 
    folder.ID       = str2double(extractAfter({folder.List.name}, 's')).';
    folder.Complete = sum(folder.ID == [100 120 130], 2) == 0;
    folder.List     = folder.List(folder.Complete,:);

    %% Define subject folders
    for iSubj = 1:size(folder.List,1)
        folder.GroupPath{iSubj,:} = fullfile(folder.Results,folder.List(iSubj).name);
    end

    for iFolder = 1:size(folder.GroupPath,1) %% Number of subjects selected
        contrast.Images = dir([folder.GroupPath{iFolder,:} filesep 'con_0001.img']); % Look for all images in the folder list, matching the contrast name (search for ContrastName.img)
        guide = ~isempty(contrast.Images); % Check if that contrast exists for a particular subject

        if  guide == 1 %If it exists, add it to the list
            contrast.Paths{iFolder,:} = fullfile(contrast.Images.folder, contrast.Images.name); % create the path to the image from a valid subject
        end
    end

    %% Remove empty cells from the paths
    contrast.Paths = contrast.Paths(~cellfun('isempty',contrast.Paths));
    
    contrast.ResultsPath = fullfile(folder.SecondLvl, 'Responses');
    if ~exist(contrast.ResultsPath, 'dir')
        mkdir(contrast.ResultsPath);
    end
    
    %% Factorial design
    spm('defaults','fmri');
    spm_jobman('initcfg');

    matlabbatch{1}.spm.stats.factorial_design.dir = {contrast.ResultsPath};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = contrast.Paths;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    %% Model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %% Write the contrasts in the SPM file
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
    matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
    matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'ButtonPress';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1;

    %% Run the batch
    spm_jobman('run',matlabbatch);
    
%% Set specifics for the group analysis using the experimental conditions
else 
    % Set the source folder for the con images
    folder.Results = fullfile (folder.Root, '1stLevel', 'AllConditions');

    %% Files to load
    file.WASI       = fullfile(folder.Behavior, 'WASI.mat');
    file.Survey     = fullfile(folder.Behavior, 'Demographic.mat');

    %% Load Demographics and WASI
    load(file.WASI, 'Scores');
    load(file.Survey, 'Age');

    %% Names for the options
    
    contrast.Pmax = {'pmax', 'pmax > rand', 'pmax > gmax', 'pmax > lmin', 'pmax > non pmax'};
    contrast.Gmax = {'gmax', 'gmax > rand', 'gmax > pmax', 'gmax > lmin', 'gmax > non gmax'};
    contrast.Lmin = {'lmin', 'lmin > rand', 'lmin > pmax', 'lmin > gmax', 'lmin > non lmin'};
    contrast.Rand = {'rand','rand > all', 'all > rand'};
    contrast.Mix  = {'non pmax > pmax', 'non gmax > gmax', 'non lmin > lmin', 'non pmax', 'structured', 'nongmax', 'non_lmin', 'non_pmax > random', 'non_gmax > random', 'zero', 'negative', 'zero > negative'};

    contrast.Options = horzcat(contrast.Pmax, contrast.Gmax, contrast.Lmin, contrast.Rand, contrast.Mix).';

    %% Names for the contrasts in the file
    contrast.PmaxFile = {'pmax', 'pmaxVSrand', 'pmaxVSgmax', 'pmaxVSlmin', 'pmaxVSnon_pmax'};
    contrast.GmaxFile = {'gmax', 'gmaxVSrand', 'gmaxVSpmax', 'gmaxVSlmin', 'gmaxVSnon_gmax'};
    contrast.LminFile = {'lmin', 'lminVSrand', 'lminVSpmax', 'lminVSgmax', 'lminVSnon_lmin'};
    contrast.RandFile = {'rand','randVSall', 'allVSrand'};
    contrast.MixFile  = {'non_pmaxVSpmax', 'non_gmaxVSgmax', 'non_lminVSlmin', 'non_pmax', 'structured', 'nongmax', 'non_lmin', 'non_pmaxVSrandom', 'non_gmaxVSrandom', 'zero', 'negative', 'zeroVSnegative'};

    contrast.Files = horzcat(contrast.PmaxFile, contrast.GmaxFile, contrast.LminFile, contrast.RandFile, contrast.MixFile).';

    % Ask if you want to include covariates in group analysis
    [options.Covariates, ~] = listdlg('ListString',{'Yes','No'},'Name','Include covariates?', 'ListSize',[250,400]);

    if options.ContrastGroup == 2 %Create single contrast
        % Ask which contrast from the list you want to create
        [options.ContrastSelection, ~] = listdlg('ListString',contrast.Options,'Name','Select which contrast you want to analyze', 'ListSize',[350,400]);
        contrast.Target = contrast.Files{options.ContrastSelection}; % Which individual contrast will be generated
        contrast.List = {contrast.Target}; % Make this contrast the only one in the list (so the for loop iterates just one time)
    else
        % Create all the contrast (the for loop iterates for n number of contrast in the list)
        contrast.List = contrast.Files;
    end
    
    %% Get all the participants folders
    folder.List = dir(folder.Results);
    folder.List(1:3,:) = []; % first rows of list are meaningless ([], ., ...]

    %% Get IDs and remove incomplete subjects 
    folder.ID       = str2double(extractAfter({folder.List.name}, 's')).';
    folder.Complete = sum(folder.ID == [100 120 130], 2) == 0;
    folder.List     = folder.List(folder.Complete,:);
    folder.ID       = folder.ID(folder.Complete,:);

    %% Sort folders according to subject number
    [folder.ID, idx] = sortrows(folder.ID);
    folder.List      = folder.List(idx,:);

    %% Get which individuals completed the task
    data.ValidScores = ismember(Scores(:,1),folder.ID(:,1));
    data.ValidAge    = ismember(Age(:,1),folder.ID(:,1));

    data.WASI = Scores(data.ValidScores,:);
    data.Age = Age(data.ValidAge,:);

    %% Get folders for specific groups
    folder.ReferenceAdolescent = regexp({folder.List.name},'s..\d'); % Adolescents have 3 digits as ID: results_sXXX
    folder.MapAdult = cellfun('isempty',folder.ReferenceAdolescent).'; % location of non Adolescents
    folder.MapAdolescent = folder.MapAdult == 0;
    % Create the folder list
    folder.Adolescent = folder.List(folder.MapAdolescent);
    folder.Adult = folder.List(folder.MapAdult);
    
    %% Define which folder list would be used, depending on selection from the menu
    if options.ContrastGroup == 3 % Adults
        for iSubj = 1:size(folder.Adult,1)
            folder.GroupPath{iSubj,:} = fullfile(folder.Results, folder.Adult(iSubj).name);
        end

    elseif options.ContrastGroup == 4 % Adolescents
        for iSubj = 1:size(folder.Adolescent,1)
            folder.GroupPath{iSubj,:} = fullfile(folder.Results,folder.Adolescent(iSubj).name);
        end

    elseif options.ContrastGroup == 5
        for iSubj = 1:size(folder.Adolescent,1)
            folder.GroupPath{iSubj,:} = fullfile(folder.Results,folder.Adolescent(iSubj).name);
        end

    else
        for iSubj = 1:size(folder.List,1) % Everyone
            folder.GroupPath{iSubj,:} = fullfile(folder.Results,folder.List(iSubj).name);
        end
    end
        
    %% ITERATE OVER ALL THE CON FILES TO CREATE 2ND LEVEL CONTRASTS
    for iContrast = 1:length(contrast.List)

        %% Assign folder for the results (create if it doesn't exist)
        if options.Covariates == 2 

            if options.ContrastGroup == 1 || options.ContrastGroup == 2 %All Subjects no covariates
                folder.Contrast = fullfile(folder.SecondLvl, 'AllConditions', 'NoCovariates', 'Everyone', contrast.List{iContrast}); %2ndLevel/ContrastName

            elseif options.ContrastGroup == 3 %Adults no covariates
                folder.Contrast = fullfile(folder.SecondLvl, 'AllConditions', 'NoCovariates', 'Adults', contrast.List{iContrast}); %2ndLevel/Adults/ContrastName

            elseif options.ContrastGroup == 4 %Adolescents no covariates
                folder.Contrast = fullfile(folder.SecondLvl, 'AllConditions', 'NoCovariates', 'Adolescents', contrast.List{iContrast}); %2ndLevel/Adolescents/ContrastName

            end

        elseif options.Covariates == 1

            if options.ContrastGroup == 1 || options.ContrastGroup == 2  %All Subjects with covariates
                folder.Contrast = fullfile(folder.SecondLvl, 'AllConditions', 'Covariates', 'Everyone', contrast.List{iContrast}); %2ndLevel/ContrastName

            elseif options.ContrastGroup == 3 %Adults with covariates
                folder.Contrast = fullfile(folder.SecondLvl, 'AllConditions', 'Covariates', 'Adults', contrast.List{iContrast}); %2ndLevel/Adults/ContrastName

            elseif options.ContrastGroup == 4 %Adolescents with covariates
                folder.Contrast = fullfile(folder.SecondLvl, 'AllConditions', 'Covariates', 'Adolescents', contrast.List{iContrast}); %2ndLevel/Adolescents/ContrastName
            end

        end

        if ~exist(folder.Contrast,'dir')
            mkdir(folder.Contrast);
        end

        %% Get the contrast images
        if options.ContrastGroup ~= 2 %% To all contrasts for the specific case (All subjects, just Adults or just Adolescents)
            contrast.Target = contrast.Files{iContrast}; % Creates one contrast from the list in every iteration
        end

        for iFolder = 1:size(folder.GroupPath,1) %% Number of subjects selected
            contrast.Images = dir([folder.GroupPath{iFolder,:} filesep 'con_' contrast.Target '.img']); % Look for all images in the folder list, matching the contrast name (search for ContrastName.img)
            guide = ~isempty(contrast.Images); % Check if that contrast exists for a particular subject

            if  guide == 1 %If it exists, add it to the list
                contrast.Paths{iFolder,:} = fullfile(contrast.Images.folder, contrast.Images.name); % create the path to the image from a valid subject
                contrast.ID(iFolder,:)    = str2double(extractAfter({contrast.Images.folder}, [filesep 's'])); % get the ID of subjects with that particular contrast
            end
        end

        %% Remove empty cells from the paths and the ID list
        contrast.Paths = contrast.Paths(~cellfun('isempty',contrast.Paths));
        contrast.ID    = contrast.ID(contrast.ID > 0); % remove spaces filled with 0's in the for loop

        %% Get covariates for the subjects included in this contrast
        contrast.WASI = data.WASI(ismember(data.WASI(:,1), contrast.ID),2); 
        contrast.Age  = data.Age(ismember(data.Age(:,1),contrast.ID),2);

        %% Sort folders according to subject number to match covariates and subject
        [contrast.ID, idxIM] = sortrows(contrast.ID);
        contrast.Paths       = contrast.Paths(idxIM,:);

        %% Save the images used in each 2nd level contrast    
        file.Summary = fullfile(folder.Contrast, 'imgsIncluded');
        imgsIncluded = [num2cell(contrast.ID), num2cell(contrast.Age), num2cell(contrast.WASI), contrast.Paths];
        save(file.Summary, 'imgsIncluded')

        %% Factorial design
        spm('defaults','fmri');
        spm_jobman('initcfg');

        matlabbatch{1}.spm.stats.factorial_design.dir          = {folder.Contrast};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = contrast.Paths;

        if options.Covariates == 1
            matlabbatch{1}.spm.stats.factorial_design.cov(1).c     = contrast.Age;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI  = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC   = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(2).c     = contrast.WASI;
            matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'WASI';
            matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI  = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC   = 1;
        else
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        end

        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

        %% Model estimation
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

        if options.Covariates == 1

            %% Contrast for Covariates
            matlabbatch{3}.spm.stats.con.spmmat(1)                      = cfg_dep;
            matlabbatch{3}.spm.stats.con.spmmat(1).tname                = 'Select SPM.mat';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name  = 'filter';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
            matlabbatch{3}.spm.stats.con.spmmat(1).sname                = 'Model estimation: SPM.mat File';
            matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch         = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
            matlabbatch{3}.spm.stats.con.spmmat(1).src_output           = substruct('.','spmmat');
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name           = contrast.Options{iContrast};
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec         = [1 0 0];
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name           = 'Age';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec         = [0 1 0];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep        = 'none';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.name           = 'WASI';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec         = [0 0 1];
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep        = 'none';
            matlabbatch{3}.spm.stats.con.delete                         = 1;


        else

            %% Contrast without covariates
            matlabbatch{3}.spm.stats.con.spmmat(1)                      = cfg_dep;
            matlabbatch{3}.spm.stats.con.spmmat(1).tname                = 'Select SPM.mat';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name  = 'filter';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
            matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
            matlabbatch{3}.spm.stats.con.spmmat(1).sname                = 'Model estimation: SPM.mat File';
            matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch         = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
            matlabbatch{3}.spm.stats.con.spmmat(1).src_output           = substruct('.','spmmat');
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name           = contrast.Options{iContrast};
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec         = 1;
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
            matlabbatch{3}.spm.stats.con.delete                         = 1;

        end
        %% Run the batch
        spm_jobman('run',matlabbatch);
        
        %% Clear path names to avoid mixing names with previous iterations
        contrast = rmfield(contrast,'Paths');
        contrast = rmfield(contrast,'Images');
        contrast = rmfield(contrast,'ID');
        contrast = rmfield(contrast,'WASI');
        contrast = rmfield(contrast,'Age');
        clear idxIM
               
    end %end of contrast loop
    
end % end of if statement for type of group analysis

%% Get back to the starting point
cd(folder.Scripts)


