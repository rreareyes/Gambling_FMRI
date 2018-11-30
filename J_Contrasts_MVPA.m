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

% It offers the possibility of creating contrasts for even and odd runs for
% the different conditions in this paradigm (see contrast.Names). 

% It also checks the validity of the contrasts (summing either 1 in single 
% conditions or 0 in comparisons).

% Finally, it changes the names of the contrast images, and spmT files
% (from con_0001, con_0002,... to con_pmax, con_pmaxVSrand,...) to make the
% contrasts more easily identified and used in group analysis.

%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% November 2018
% working on SPM8

%% Clean workspace
clc; clear

%% Base Paths
cd('..')
folder.Root        = pwd;
folder.Processed   = fullfile(folder.Root, 'Processed');
folder.Scripts     = fullfile(folder.Root, 'Scripts');
folder.Results{1}  = fullfile(folder.Root, 'MVPA', 'Odd', '1stLevel');
folder.Results{2}  = fullfile(folder.Root, 'MVPA', 'Even', '1stLevel');
folder.ContrBackup{1} = fullfile(folder.Root, 'MVPA', 'Odd', 'ContrastVectors');
folder.ContrBackup{2} = fullfile(folder.Root, 'MVPA', 'Even', 'ContrastVectors');

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

%% Names for the contrasts in the SPM file
contrast.Pmax = {'pmax', 'pmax > rand', 'pmax > gmax', 'pmax > lmin', 'pmax > non pmax'};
contrast.Gmax = {'gmax', 'gmax > rand', 'gmax > pmax', 'gmax > lmin', 'gmax > non gmax'};
contrast.Lmin = {'lmin', 'lmin > rand', 'lmin > pmax', 'lmin > gmax', 'lmin > non lmin'};
contrast.Rand = {'rand','rand > all', 'all > rand'};
contrast.Mix  = {'non pmax > pmax', 'non gmax > gmax', 'non lmin > lmin', 'non pmax', 'structured', 'nongmax', 'non_lmin', 'non_pmax > random', 'non_gmax > random', 'zero', 'negative', 'zero > negative'};

contrast.Names = horzcat(contrast.Pmax, contrast.Gmax, contrast.Lmin, contrast.Rand, contrast.Mix).';

%% Names for the files
file.Pmax = {'pmax', 'pmaxVSrand', 'pmaxVSgmax', 'pmaxVSlmin', 'pmaxVSnon_pmax'};
file.Gmax = {'gmax', 'gmaxVSrand', 'gmaxVSpmax', 'gmaxVSlmin', 'gmaxVSnon_gmax'};
file.Lmin = {'lmin', 'lminVSrand', 'lminVSpmax', 'lminVSgmax', 'lminVSnon_lmin'};
file.Rand = {'rand','randVSall', 'allVSrand'};
file.Mix  = {'non_pmaxVSpmax', 'non_gmaxVSgmax', 'non_lminVSlmin', 'non_pmax', 'structured', 'nongmax', 'non_lmin', 'non_pmaxVSrandom', 'non_gmaxVSrandom', 'zero', 'negative', 'zeroVSnegative'};

file.Names = horzcat(file.Pmax, file.Gmax, file.Lmin, file.Rand, file.Mix).';

%% Create models requested
for iModel = 1:options.nModel
    
    if options.Runs == 3
        options.Model = iModel; %first iteration creates odd model and second even model.

    else
        options.Model = options.Runs;

    end    
    
    %% Loop through the subjects to get onsets and perform 1st level analysis
    for iSubject = 1:size(group.SubjectsPaths,1)
        %% Clear values to avoid overwriting issues
        clear subject matlabbatch errorCheck SPM 

        %% Set paths for SPM.mat source and to save contrast vectors
        subject.ResultFolder   = group.ResultPaths{iSubject, options.Model};
        subject.ContrastFolder = fullfile(subject.ResultFolder, 'contrasts');

        if ~exist(folder.ContrBackup{options.Model},'dir')
            mkdir(folder.ContrBackup{options.Model});
        end

        if ~exist(subject.ContrastFolder,'dir')
            mkdir(subject.ContrastFolder);
        end


        %% Get the SPM file from result folder
        load(fullfile(subject.ResultFolder, 'SPM.mat'));

        %% Get all condition names on the matrix
        subject.Conditions = SPM.xX.name;

        %% Locate all the points where a condition appears in the design matrix and create and array of 1's & 0's for each
        subject.Rand   = contains(subject.Conditions, 'rand');
        subject.Struct = contains(subject.Conditions, {'pmax', 'gmax', 'lmin'});

        subject.PmaxZero = contains(subject.Conditions, 'pmax_zero'); 
        subject.GmaxZero = contains(subject.Conditions, 'gmax_zero');
        subject.LminZero = contains(subject.Conditions, 'lmin_zero');

        subject.PmaxNeg = contains(subject.Conditions, 'pmax_neg'); 
        subject.GmaxNeg = contains(subject.Conditions, 'gmax_neg');
        subject.LminNeg = contains(subject.Conditions, 'lmin_neg');

        subject.Pmax = contains(subject.Conditions, 'pmax'); 
        subject.Gmax = contains(subject.Conditions, 'gmax');
        subject.Lmin = contains(subject.Conditions, 'lmin');

        subject.Zero = contains(subject.Conditions, '_zero'); 
        subject.Neg  = contains(subject.Conditions, '_neg');

        %% Get the number of runs where a condition was present
        subject.nRand   = sum(subject.Rand);
        subject.nStruct = sum(subject.Struct);

        subject.nPmaxZero = sum(subject.PmaxZero); 
        subject.nGmaxZero = sum(subject.GmaxZero);
        subject.nLminZero = sum(subject.LminZero);

        subject.nPmaxNeg = sum(subject.PmaxNeg); 
        subject.nGmaxNeg = sum(subject.GmaxNeg);
        subject.nLminNeg = sum(subject.LminNeg);

        subject.nPmax = sum(subject.Pmax); 
        subject.nGmax = sum(subject.Gmax);
        subject.nLmin = sum(subject.Lmin);

        subject.nZero = sum(subject.Zero);
        subject.nNeg  = sum(subject.Neg);

        %% Get the actual weight of every condition for the contrast matrix
        %% Random Trials
        subject.WeightRand = subject.Rand ./ subject.nRand;

        %% Zero Trials
        % pmax
        subject.WeightPmaxZero = subject.PmaxZero ./ subject.nPmaxZero;

        % gmax
        subject.WeightGmaxZero = subject.GmaxZero ./ subject.nGmaxZero;

        % lmin
        subject.WeightLminZero = subject.LminZero ./ subject.nLminZero;

        % All Zero
        subject.WeightZero = subject.Zero ./ subject.nZero;

        %% Negative Trials
        % pmax
        subject.WeightPmaxNeg = subject.PmaxNeg ./ subject.nPmaxNeg;

        % gmax
        subject.WeightGmaxNeg = subject.GmaxNeg ./ subject.nGmaxNeg;

        % lmin
        subject.WeightLminNeg = subject.LminNeg ./ subject.nLminNeg;

        % All Negative
        subject.WeightNeg = subject.Neg ./ subject.nNeg;

        %% Structured
        subject.WeightPmax = subject.Pmax ./ subject.nPmax;
        subject.WeightGmax = subject.Gmax ./ subject.nGmax;
        subject.WeightLmin = subject.Lmin ./ subject.nLmin;

        subject.WeightStructured = subject.Struct ./ subject.nStruct;

        %% mix
        subject.WeightRandPmax     = (subject.WeightRand + subject.WeightPmax) ./ 2; 
        subject.WeightRandGmax     = (subject.WeightRand + subject.WeightGmax) ./ 2;
        subject.WeightRandLmin     = (subject.WeightRand + subject.WeightLmin) ./ 2;

        subject.WeightNonlmin      = (subject.WeightPmax + subject.WeightGmax) ./ 2;
        subject.WeightNongmax      = (subject.WeightPmax + subject.WeightLmin) ./ 2;
        subject.WeightNonpmax      = (subject.WeightGmax + subject.WeightLmin) ./ 2;

        %% Create the arrays for every contrast
        %% pmax
        subject.ContrPmax        = subject.WeightPmax;
        subject.ContrPmaxRand    = subject.WeightPmax - subject.WeightRand;
        subject.ContrPmaxGmax    = subject.WeightPmax - subject.WeightGmax;
        subject.ContrPmaxLmin    = subject.WeightPmax - subject.WeightLmin;
        subject.ContrPmaxNonpmax = subject.WeightPmax - subject.WeightNonpmax;

        %% zero and negative
        subject.ContrZero    = subject.WeightZero;
        subject.ContrNeg     = subject.WeightNeg;
        subject.ContrZeroNeg = subject.WeightZero - subject.WeightNeg;

        %% gmax
        subject.ContrGmax          = subject.WeightGmax;
        subject.ContrGmaxVSRand    = subject.WeightGmax - subject.WeightRand;
        subject.ContrGmaxVSPmax    = subject.WeightGmax - subject.WeightPmax;
        subject.ContrGmaxVSLmin    = subject.WeightGmax - subject.WeightLmin;
        subject.ContrGmaxVSNongmax = subject.WeightGmax - subject.WeightNongmax;

        %% lmin
        subject.ContrLmin          = subject.WeightLmin;
        subject.ContrLminVSRand    = subject.WeightLmin - subject.WeightRand;
        subject.ContrLminVSPmax    = subject.WeightLmin - subject.WeightPmax;
        subject.ContrLminVSGmax    = subject.WeightLmin - subject.WeightGmax;
        subject.ContrLminVSNonlmin = subject.WeightLmin - subject.WeightNonlmin;

        %% random
        subject.ContrRand             = subject.WeightRand;
        subject.ContrRandVSStructured = subject.WeightRand - subject.WeightStructured;
        subject.ContrStructuredVSRand = subject.WeightStructured - subject.WeightRand;

        %% mix
        subject.ContrNonpmaxPmax = subject.WeightNonpmax - subject.WeightPmax;
        subject.ContrNongmaxGmax = subject.WeightNongmax - subject.WeightGmax;
        subject.ContrNonlminLmin = subject.WeightNonlmin - subject.WeightLmin;
        subject.ContrNonpmax     = subject.WeightNonpmax;
        subject.ContrStructured  = subject.WeightStructured;
        subject.ContrNongmax     = subject.WeightNongmax;
        subject.ContrNonlmin     = subject.WeightNonlmin;
        subject.ContrNonpmaxRand = (subject.WeightNonpmax) - subject.WeightRand;
        subject.ContrNongmaxRand = (subject.WeightNongmax) - subject.WeightRand;

        %% Create contrast matrix
        % Identify which contrasts exist for this subject
        subject.RawContrasts = vertcat(subject.ContrPmax, subject.ContrPmaxRand, subject.ContrPmaxGmax, subject.ContrPmaxLmin,            ...
                                       subject.ContrPmaxNonpmax, subject.ContrGmax, subject.ContrGmaxVSRand, subject.ContrGmaxVSPmax,     ...
                                       subject.ContrGmaxVSLmin, subject.ContrGmaxVSNongmax, subject.ContrLmin, subject.ContrLminVSRand,   ...
                                       subject.ContrLminVSPmax, subject.ContrLminVSGmax, subject.ContrLminVSNonlmin, subject.ContrRand,   ...
                                       subject.ContrRandVSStructured, subject.ContrStructuredVSRand, subject.ContrNonpmaxPmax,            ...
                                       subject.ContrNongmaxGmax, subject.ContrNonlminLmin, subject.ContrNonpmax, subject.ContrStructured, ...
                                       subject.ContrNongmax, subject.ContrNonlmin, subject.ContrNonpmaxRand, subject.ContrNongmaxRand,    ...
                                       subject.ContrZero, subject.ContrNeg, subject.ContrZeroNeg);

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
        errorCheck.Statements                                          = cell(30,1);
        errorCheck.Statements ([1, 6, 11, 16, 22, 23:25, 28:29],1)     = {'does not sum 1, check weights'};
        errorCheck.Statements ([2:5, 7:10, 12:15, 17:21, 26:27, 30],1) = {'does not sum 0, check weights'};
        errorCheck.Messages                                            = horzcat(contrast.Names, errorCheck.Statements);
        errorCheck.Messages                                            = errorCheck.Messages(subject.ValidContrasts,:); 

        %% Create a reference matrix to check the existing contrast values
        errorCheck.ReferenceMatrix                                   = zeros(30, 1);
        errorCheck.ReferenceMatrix([1, 6, 11, 16, 22, 23:25, 28:29]) = 1;
        errorCheck.ReferenceMatrix                                   = errorCheck.ReferenceMatrix(subject.ValidContrasts);

        %% Sum values from existing contrasts
        errorCheck.ContrastValue = round(sum(subject.ContrastMatrix,2), 5);

        %% Check every contrast and display a message if something is wrong
        for iContrast = 1:size(subject.ContrastNames, 1)
            if errorCheck.ReferenceMatrix (iContrast, 1) ~= errorCheck.ContrastValue (iContrast,1)
                msgbox([errorCheck.Messages{iContrast, 1} ' ' errorCheck.Messages{iContrast,2}], 'CONTRAST INVALID', 'warn')
                errorCheck.Problem = 1;
            else
                errorCheck.Problem = 0;
            end
        end

        %% If everything is fine, save the contrasts in each subjects folder and in the backup folder
        if errorCheck.Problem == 0

            backupMatrix = horzcat(contrast.Names, num2cell(subject.RawContrasts));
            save (fullfile(subject.ContrastFolder, [group.SubjectsList(iSubject).name '_ContrVector.mat']), 'backupMatrix')
            save (fullfile(folder.ContrBackup{options.Model}, [group.SubjectsList(iSubject).name '_ContrVector.mat']), 'backupMatrix')

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

            for iContrast=1:size(subject.ContrastMatrix, 1)
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

    end % end of subject loop

end % end of models loop
    
%% Return to scripts folder
cd(folder.Scripts)
