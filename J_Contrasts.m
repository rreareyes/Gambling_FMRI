%% CREATE CONTRAST MATRIX BASED ON DATA FROM SPM FILE
% This script writes contrast vectors to the SPM file using the contrast
% manager from the SPM8 package. It also creates the con images for each of
% these comparsions.

% It assumes that the results from the 1st level model estimation are
% stored as follows:

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

% It offers the possibility of creating contrasts for button responses or
% for the different conditions in this paradigm (see contrast.Names). 

% For the experimental conditions it also checks the validity of the
% contrasts (summing either 1 in single conditions or 0 in comparisons).

% Finally, it changes the names of the contrast images, and spmT files
% (from con_0001, con_0002,... to con_pmax, con_pmaxVSrand,...) to make the
% more easily identified and used in group analysis.


%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% June 2018
% working on SPM8

%% Clean workspace
clc; clear

%% Base Paths
cd('..')
folder.Root        = pwd;
folder.Processed   = fullfile(folder.Root, 'Processed');
folder.Scripts     = fullfile(folder.Root, 'Scripts');
folder.ContrBackup = fullfile(folder.Root, 'ContrastVectors');

%% Get all subject paths
folder.ProcessedPaths      = dir(folder.Processed);
folder.ProcessedPaths(1:2) = [];

%% Set subject paths according to selection
%% Type of Model
% Ask what analysis you want to perform
[options.Style, ~] = listdlg('ListString', {'Just Responses', 'Multiple conditions'}, 'Name', 'Type of Analysis');

% Just test for response onsets
if options.Style == 1 
   folder.Results = fullfile(folder.Root, '1stLevel', 'Responses'); 

% Test for responses in all conditions   
elseif options.Style == 2 
   folder.Results = fullfile(folder.Root, '1stLevel', 'AllConditions'); 
end

%% Define your group
%% Ask for which subjects to run
[options.Group, ~] = listdlg('ListString',{'Individual Elements','All Subjects'},'Name','No. Subjects to Process?');

%% Set subject paths for input and output images according to subjects selected
if options.Group == 1 % Customized list
    [options.Subjects, ~] = listdlg('ListString',char(folder.ProcessedPaths.name),'Name','Which subjects do you want?');
    group.SubjectsList = folder.ProcessedPaths(options.Subjects);

elseif options.Group == 2 % All subjects
    group.SubjectsList = folder.ProcessedPaths;
end

for iFolder = 1:size(group.SubjectsList,1)
    group.SubjectsPaths{iFolder,:} = fullfile(folder.Processed, group.SubjectsList(iFolder).name);
    group.ResultPaths{iFolder,:}   = fullfile(folder.Results, group.SubjectsList(iFolder).name);
end

%% Names for the contrasts in the SPM file
contrast.Pmax = {'pmax', 'pmax > rand', 'pmax > gmax', 'pmax > lmin', 'pmax > non pmax'};
contrast.Gmax = {'gmax', 'gmax > rand', 'gmax > pmax', 'gmax > lmin', 'gmax > non gmax'};
contrast.Lmin = {'lmin', 'lmin > rand', 'lmin > pmax', 'lmin > gmax', 'lmin > non lmin'};
contrast.Rand = {'rand','rand > all', 'all > rand'};
contrast.Mix  = {'non pmax > pmax', 'non gmax > gmax', 'non lmin > lmin', 'non pmax', 'non rand'};

contrast.Names = horzcat(contrast.Pmax, contrast.Gmax, contrast.Lmin, contrast.Rand, contrast.Mix).';

%% Names for the contrasts in the SPkeepVariablesM file
file.Pmax = {'pmax', 'pmaxVSrand', 'pmaxVSgmax', 'pmaxVSlmin', 'pmaxVSnon_pmax'};
file.Gmax = {'gmax', 'gmaxVSrand', 'gmaxVSpmax', 'gmaxVSlmin', 'gmaxVSnon_gmax'};
file.Lmin = {'lmin', 'lminVSrand', 'lminVSpmax', 'lminVSgmax', 'lminVSnon_lmin'};
file.Rand = {'rand','randVSall', 'allVSrand'};
file.Mix  = {'non_pmaxVSpmax', 'non_gmaxVSgmax', 'non_lminVSlmin', 'non_pmax', 'non_rand'};

file.Names = horzcat(file.Pmax, file.Gmax, file.Lmin, file.Rand, file.Mix).';

%% Loop through the subjects to get onsets and perform 1st level analysis
for iSubject = 1:size(group.SubjectsPaths,1)
    %% Clear values to avoid overwriting issues
    clear subject matlabbatch errorCheck SPM 

    %% Set paths for SPM.mat source and to save contrast vectors
    subject.ResultFolder   = group.ResultPaths{iSubject,:};
    subject.ContrastFolder = fullfile(subject.ResultFolder, 'contrasts');
    
    if ~exist(folder.ContrBackup,'dir')
        mkdir(folder.ContrBackup);
    end
    
    if ~exist(subject.ContrastFolder,'dir')
        mkdir(subject.ContrastFolder);
    end
    
    
    %% Get the SPM file from result folder
    load(fullfile(subject.ResultFolder, 'SPM.mat'));
    
    %% Get all condition names on the matrix
    subject.Conditions = SPM.xX.name;
    
    %% CONTRAST FOR ALL BUTTON RESPONSES
    if options.Style == 1
        %% Create a contras vector of 1's for each time a button press happened
        subject.ContrastVector = double(contains(subject.Conditions, 'Button'));
        
        %% Generate contrasts in the SPM file
        spm('defaults','fmri');
        spm_jobman('initcfg');
        
        matlabbatch{1}.spm.stats.con.spmmat                  = {fullfile(subject.ResultFolder, 'SPM.mat')};
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name    = 'Button Responses';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec  = subject.ContrastVector;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete                  = 1;
        
        %% Insert constrasts in the SPM file
        spm_jobman('run',matlabbatch);
    
    %% CONTRAST FOR DIFFERENT CONDITIONS
    elseif options.Style == 2 

        %% Locate all the points where a condition appears in the design matrix and create and array of 1's & 0's for each
        subject.Rand = contains(subject.Conditions, 'rand'); 
        subject.Pmax = contains(subject.Conditions, 'pmax'); 
        subject.Gmax = contains(subject.Conditions, 'gmax');
        subject.Lmin = contains(subject.Conditions, 'lmin');

        %% Get the number of runs where a condition was present
        subject.nRand = sum(subject.Rand); 
        subject.nPmax = sum(subject.Pmax); 
        subject.nGmax = sum(subject.Gmax);
        subject.nLmin = sum(subject.Lmin);

        %% Get the actual weight of every condition for the contrast matrix
        %% random
        subject.WeightRand       = subject.Rand ./ subject.nRand;
        subject.WeightRandPair   = subject.Rand ./ (2 * subject.nRand); % When adding 2 conditions, we need to give half of the weight to each
        subject.WeightRandGlobal = subject.Rand ./ (3 * subject.nRand); % When grouping all conditions, we need to give 1/3 of the weight to each

        %% pmax
        subject.WeightPmax       = subject.Pmax ./ subject.nPmax;
        subject.WeightPmaxPair   = subject.Pmax ./ (2 * subject.nPmax); 
        subject.WeightPmaxGlobal = subject.Pmax ./ (3 * subject.nPmax); 

        %% gmax
        subject.WeightGmax       = subject.Gmax ./ subject.nGmax;
        subject.WeightGmaxPair   = subject.Gmax ./ (2 * subject.nGmax); 
        subject.WeightGmaxGlobal = subject.Gmax ./ (3 * subject.nGmax); 

        %% lmin
        subject.WeightLmin       = subject.Lmin ./ subject.nLmin;
        subject.WeightLminPair   = subject.Lmin ./ (2 * subject.nLmin); 
        subject.WeightLminGlobal = subject.Lmin ./ (3 * subject.nLmin); 

        %% mix
        subject.WeightRandPmax     = subject.WeightRandPair + subject.WeightPmaxPair; 
        subject.WeightRandGmax     = subject.WeightRandPair + subject.WeightGmaxPair;
        subject.WeightRandLmin     = subject.WeightRandPair + subject.WeightLminPair;
        subject.WeightPmaxGmax     = subject.WeightPmaxPair + subject.WeightGmaxPair;
        subject.WeightPmaxLmin     = subject.WeightPmaxPair + subject.WeightLminPair;
        subject.WeightGmaxLmin     = subject.WeightGmaxPair + subject.WeightLminPair;
        subject.WeightPmaxGmaxLmin = subject.WeightPmaxGlobal + subject.WeightGmaxGlobal + subject.WeightLminGlobal;

        %% Create the arrays for every contrast
        %% pmax
        subject.ContrPmax        = subject.WeightPmax;
        subject.ContrPmaxRand    = subject.WeightPmax - subject.WeightRand;
        subject.ContrPmaxGmax    = subject.WeightPmax - subject.WeightGmax;
        subject.ContrPmaxLmin    = subject.WeightPmax - subject.WeightLmin;
        subject.ContrPmaxNonpmax = subject.WeightPmax - subject.WeightGmaxLmin;

        %% gmax
        subject.ContrGmax        = subject.WeightGmax;
        subject.ContrGmaxRand    = subject.WeightGmax - subject.WeightRand;
        subject.ContrGmaxPmax    = subject.WeightGmax - subject.WeightPmax;
        subject.ContrGmaxLmin    = subject.WeightGmax - subject.WeightLmin;
        subject.ContrGmaxNongmax = subject.WeightGmax - subject.WeightPmaxLmin;

        %% lmin
        subject.ContrLmin        = subject.WeightLmin;
        subject.ContrLminRand    = subject.WeightLmin - subject.WeightRand;
        subject.ContrLminPmax    = subject.WeightLmin - subject.WeightPmax;
        subject.ContrLminGmax    = subject.WeightLmin - subject.WeightGmax;
        subject.ContrLminNonlmin = subject.WeightLmin - subject.WeightPmaxGmax;

        %% random
        subject.ContrRand        = subject.WeightRand;
        subject.ContrRandNonrand = subject.WeightRand - subject.WeightPmaxGmaxLmin;
        subject.ContrNonrandRand = subject.WeightPmaxGmaxLmin - subject.WeightRand;

        %% mix
        subject.ContrNonpmaxPmax = (subject.WeightGmaxPair + subject.WeightLminPair)- subject.WeightPmax;
        subject.ContrNongmaxGmax = (subject.WeightPmaxPair + subject.WeightLminPair)- subject.WeightGmax;
        subject.ContrNonlminLmin = (subject.WeightPmaxPair + subject.WeightGmaxPair)- subject.WeightLmin;
        subject.ContrNonpmax     = subject.WeightGmaxLmin;
        subject.Nonrand          = subject.WeightPmaxGmaxLmin;

        %% Create contrast matrix
        % Identify the which contrast exist for this subject
        subject.RawContrasts     = vertcat(subject.ContrPmax, subject.ContrPmaxRand, subject.ContrPmaxGmax, subject.ContrPmaxLmin, subject.ContrPmaxNonpmax, subject.ContrGmax, subject.ContrGmaxRand, subject.ContrGmaxPmax, subject.ContrGmaxLmin, subject.ContrGmaxNongmax, subject.ContrLmin, subject.ContrLminRand, subject.ContrLminPmax, subject.ContrLminGmax, subject.ContrLminNonlmin, subject.ContrRand, subject.ContrRandNonrand, subject.ContrNonrandRand, subject.ContrNonpmaxPmax, subject.ContrNongmaxGmax, subject.ContrNonlminLmin, subject.ContrNonpmax, subject.Nonrand);
        subject.NotEmptyContrast = subject.RawContrasts ~= 0 & ~isnan(subject.RawContrasts);
        subject.ValidContrasts   = sum(subject.NotEmptyContrast,2) ~= 0; % mark non-empty contrasts

        % Select vectors, names and file names to save for existing contrast
        subject.ContrastMatrix = subject.RawContrasts(subject.ValidContrasts,:); % contrast value
        subject.ContrastNames  = contrast.Names(subject.ValidContrasts); % contrast names
        subject.ContrastFile   = file.Names(subject.ValidContrasts); % file names

        %% Check integrity of contrasts
        % If the contrast exists (i.e. is different from 0)
        % Comparisons between conditions must sum 0 (rounded to 5 decimals)
        % Single conditions must sum 1
        %% Statements for boxes
        errorCheck.Statements                               = cell(23,1);
        errorCheck.Statements ([1, 6, 11, 16, 22, 23],1)    = {'does not sum 1, check weights'};
        errorCheck.Statements ([2:5, 7:10, 12:15, 17:21],1) = {'does not sum 0, check weights'};
        errorCheck.Messages                                 = horzcat(contrast.Names, errorCheck.Statements);
        errorCheck.Messages                                 = errorCheck.Messages(subject.ValidContrasts,:); 

        %% Create a reference matrix to check the existing contrast values
        errorCheck.ReferenceMatrix                         = zeros(23,1);
        errorCheck.ReferenceMatrix([1, 6, 11, 16, 22, 23]) = 1;
        errorCheck.ReferenceMatrix                         = errorCheck.ReferenceMatrix(subject.ValidContrasts);

        %% Sum values from existing contrasts
        errorCheck.ContrastValue = round(sum(subject.ContrastMatrix,2),5);

        %% Check every contrast and display a message if something is wrong
        for iContrast = 1:size(subject.ContrastNames, 1)
            if errorCheck.ReferenceMatrix (iContrast,1) ~= errorCheck.ContrastValue (iContrast,1)
                msgbox([errorCheck.Messages{iContrast,1} ' ' errorCheck.Messages{iContrast,2}], 'CONTRAST INVALID', 'warn')
                errorCheck.Problem = 1;
            else
                errorCheck.Problem = 0;
            end
        end

        %% If everything is fine, save the contrasts in each subjects folder and in the backup folder
        if errorCheck.Problem == 0
            
            backupMatrix = horzcat(contrast.Names, num2cell(subject.RawContrasts));
            save (fullfile(subject.ContrastFolder, [group.SubjectsList(iSubject).name '_ContrVector.mat']), 'backupMatrix')
            save (fullfile(folder.ContrBackup, [group.SubjectsList(iSubject).name '_ContrVector.mat']), 'backupMatrix')
    
            %% Generate contrasts in the SPM file
            spm('defaults','fmri');
            spm_jobman('initcfg');
            
            for iContrast = 1:size(subject.ContrastMatrix, 1)
                matlabbatch{1}.spm.stats.con.spmmat                          = {fullfile(subject.ResultFolder, 'SPM.mat')}; % SPM file to use
                matlabbatch{1}.spm.stats.con.consess{iContrast}.tcon.name    = char(subject.ContrastNames(iContrast)); % Contrast's name
                matlabbatch{1}.spm.stats.con.consess{iContrast}.tcon.convec  = subject.ContrastMatrix(iContrast,:); % The specified contrast array
                matlabbatch{1}.spm.stats.con.consess{iContrast}.tcon.sessrep = 'none'; % Don't repeat over sessions
                matlabbatch{1}.spm.stats.con.delete                          = 1; %delete pre-existent contrasts; set to 0 if you don't want to do it
                SPM.xCon(iContrast).Vcon.fname                               = char(subject.ContrastFile(iContrast));
            end
            
            %% Insert constrasts in the SPM file
            spm_jobman('run',matlabbatch);

            %% Rename contrast images to have meaningful names
            imgOriginals    = dir(fullfile(subject.ResultFolder, '*con*.img'));
            hdrOriginals    = dir(fullfile(subject.ResultFolder, '*con*.hdr'));
            spmImgOriginals = dir(fullfile(subject.ResultFolder, '*spmT*.img'));
            spmHdrOriginals = dir(fullfile(subject.ResultFolder, '*spmT*.hdr'));

            for iContrast=1:size(subject.ContrastMatrix,1)
                movefile(fullfile(subject.ResultFolder, imgOriginals(iContrast).name), fullfile(subject.ResultFolder, ['con_' char(subject.ContrastFile(iContrast)) '.img']))
                movefile(fullfile(subject.ResultFolder, hdrOriginals(iContrast).name), fullfile(subject.ResultFolder, ['con_' char(subject.ContrastFile(iContrast)) '.hdr']))
                movefile(fullfile(subject.ResultFolder, spmImgOriginals(iContrast).name), fullfile(subject.ResultFolder, ['spmT_' char(subject.ContrastFile(iContrast)) '.img']))
                movefile(fullfile(subject.ResultFolder, spmHdrOriginals(iContrast).name), fullfile(subject.ResultFolder, ['spmT_' char(subject.ContrastFile(iContrast)) '.hdr']))
            end
            
            %% Add the new names to the SPM file
            load(fullfile(subject.ResultFolder, 'SPM.mat'));
            for iContrast = 1:size(subject.ContrastMatrix, 1)
                SPM.xCon(iContrast).Vcon.fname = ['con_' char(subject.ContrastFile(iContrast)) '.img'];
                SPM.xCon(iContrast).Vspm.fname = ['spmT_' char(subject.ContrastFile(iContrast)) '.img'];
            end
            
            save(fullfile(subject.ResultFolder, 'SPM.mat'), 'SPM');

        end % end of if statement on error check


    end  % end of if statement on analysis type

end

%% Return to scripts folder
cd(folder.Scripts)

