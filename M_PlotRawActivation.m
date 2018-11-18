%% Plot Raw Activation Data from ROI

%% Clean workspace
clc; clear; close all

%% Base paths
cd('..')
folder.Root        = pwd;
folder.Scripts     = fullfile(folder.Root, 'Scripts');
folder.ROI         = fullfile(folder.Root, 'SphericalROI');
folder.ROIclean    = fullfile(folder.ROI, 'Clean');
folder.Activation  = fullfile(folder.ROI, 'RawResults');
folder.GroupDesign = fullfile(folder.Root, '2ndLevel', 'AllConditions', 'NoCovariates');
folder.Adults      = fullfile(folder.GroupDesign, 'Adults');
folder.Adolescents = fullfile(folder.GroupDesign, 'Adolescents');

if ~exist(folder.ROIclean, 'dir')
    mkdir(folder.ROIclean);
end

%% Define contrasts and ROIs names
names.ConFolders = {'pmax', 'non_pmax', 'structured', 'rand', 'zero', 'negative'};
nContrasts       = length(names.ConFolders);

names.ROI      = {'Right AI', 'Right IFG', 'Right DLPFC', 'DMPFC'};
names.Contrast = {'Negative', 'Non Pmax', 'Pmax', 'Random', 'Structured', 'Zero'};

%% Get ROI files and contrasts to extract
roi.List = dir(fullfile(folder.ROI, '*.mat'));

folder.ConAdults      = dir(folder.Adults);
folder.ConAdolescents = dir(folder.Adolescents);

for iContrast = 1:nContrasts
    folder.ConList(iContrast)              = folder.ConAdults(strcmp({folder.ConAdults.name}.', names.ConFolders{iContrast}));
    folder.ConList(iContrast + nContrasts) = folder.ConAdolescents(strcmp({folder.ConAdolescents.name}.', names.ConFolders{iContrast}));
end

%% Prevent problems with NAN's using NN resampling and saving new ROIS
marsbar('on')
for iROI = 1:size(roi.List, 1)
    roi.Path   = fullfile(folder.ROI, roi.List(iROI).name);
    roi.Object = maroi(roi.Path); % make maroi ROI objects
    roi.Clean  = spm_hold(roi.Object, 0); % set NN resampling
    saveroi(roi.Clean, fullfile(folder.ROIclean, roi.List(iROI).name));
end

%% Create absolute paths for SPM files and ROIs
roi.CleanList = dir(fullfile(folder.ROIclean, '*.mat'));

for iROI = 1:size(roi.List, 1)
    roi.CleanPath{iROI, 1} = fullfile(folder.ROIclean, roi.List(iROI).name);
end

for iDesign = 1:size(folder.ConList, 2)
    design.Paths{iDesign, 1} = fullfile(folder.ConList(iDesign).folder, folder.ConList(iDesign).name, 'SPM.mat');
end 

%% Extract activation values from each contrast and ROI
for iDesign = 1:size(design.Paths, 1)
    %roi_files = spm_get(Inf, '*roi.mat', 'Select ROI files');    %yields array of char
    %des_path = spm_get(1, 'SPM.mat', 'Select SPM.mat');
    rois = maroi(roi.CleanPath);
    des = mardo(design.Paths{iDesign});  % make mardo design object
    mY = get_marsy(rois{:}, des, 'mean'); % extract data into marsy data object
    y  = summary_data(mY);
end
 
%% Define colors to use
rgb.LightGreen = [0.443137 0.776471 0.443137];
rgb.DarkGreen  = [0.180392 0.545098 0.341176];

rgb.LightBlue = [0.5765 0.8000 0.9176];
rgb.DarkBlue  = [0.3922 0.5843 0.9294];

rgb.LightGolden = [1.0000 0.9686 0.8000];
rgb.DarkGolden  = [1.0000 0.9020 0.4000];

rgb.LightSalmon = [1.0000 0.8667 0.8000];
rgb.DarkSalmon  = [0.9137 0.5882 0.4784];

rgb.LightPurple = [0.9333 0.8000 1.0000];
rgb.DarkPurple  = [0.8353 0.5020 1.0000];

rgb.LightAqua = [0.8000 1.0000 0.9686];
rgb.DarkAqua  = [0.4000 1.0000 0.9020];

rgb.LightRed = [1.0000 0.8000 0.8000];
rgb.DarkRed  = [1.0000 0.4000 0.4000];


color.Adult      = [rgb.LightBlue; rgb.DarkBlue; rgb.LightSalmon; rgb.DarkSalmon; rgb.LightPurple; rgb.DarkPurple];
color.Adolescent = [rgb.LightGreen; rgb.DarkGreen; rgb.LightGolden; rgb.DarkGolden; rgb.LightRed; rgb.DarkRed];

%% Get files with activation data
file.Names = {'pmax', 'nonpmax', 'struct', 'rand', 'zero', 'neg'};

%% Loop trough files to extract activation from each ROI
for iContrast = 1:size(file.Names, 2)
    %% Load data
    data.ContrastName = file.Names{1, iContrast};
    
    data.AdolescentPath = fullfile(folder.Activation, [data.ContrastName '_Adolescents.mat']);
    data.AdultPath      = fullfile(folder.Activation, [data.ContrastName '_Adults.mat']);
        
    data.Adolescent = load(data.AdolescentPath);
    data.Adult      = load(data.AdultPath);
    
    data.AdolescentVals = data.Adolescent.SPM.marsY.Y;
    data.AdultVals      = data.Adult.SPM.marsY.Y;
        
    %% Get Activation Data
    data.Activation{iContrast, 1} = data.ContrastName;
    data.Error{iContrast, 1} = data.ContrastName;
    
    data.Activation{iContrast, 2} = [mean(data.AdolescentVals), mean(data.AdultVals)]; 
    data.Error{iContrast, 2}      = [std(data.AdolescentVals)./ sqrt(size(data.AdolescentVals, 1)), std(data.AdultVals) ./ sqrt(size(data.AdultVals, 1))];
    
end
%%
names.Plot = {'Pmax', 'Non Pmax', 'Structured', 'Random', 'Zero', 'Negative'};

%% Create plots to compare conditions in Adults and Adolescents in each ROI
for iContrast = [1, 3, 5]
    figure('Name', [names.Plot{iContrast} ' vs ' names.Plot{iContrast + 1}]) %name of the figure
        
    % Plot Pmax Ado
    y1 = data.Activation{iContrast, 2}(1:4);
    e1 = data.Error{iContrast, 2}(1:4);
    x1 = [1:4] - 0.3; % Number of ROIs
    
    bar(x1, y1, 0.20, 'FaceColor', color.Adolescent(iContrast, :),'EdgeColor', color.Adolescent(iContrast, :))
    
    hold on
       
    % Plot NonPmax Ado
    y2 = data.Activation{iContrast + 1, 2}(1:4);
    e2 = data.Error{2,2}(1:4);
    x2 = [1:4] - 0.1;
    
    bar(x2, y2, 0.20, 'FaceColor', color.Adolescent(iContrast + 1, :), 'EdgeColor', color.Adolescent(iContrast + 1, :))
    
    % Plot Pmax Adults
    y3 = data.Activation{iContrast, 2}(5:8);
    e3 = data.Error{iContrast, 2}(5:8);
    x3 = [1:4] + 0.1; % Number of ROIs
    
    bar(x3, y3, 0.20, 'FaceColor', color.Adult(iContrast, :),'EdgeColor', color.Adult(iContrast, :))
    
    hold on
       
    % Plot NonPmax Adults
    y4 = data.Activation{iContrast + 1, 2}(5:8);
    e4 = data.Error{iContrast + 1, 2}(5:8);
    x4 = [1:4] + 0.3;
    
    bar(x4, y4, 0.20, 'FaceColor', color.Adult(iContrast + 1, :), 'EdgeColor', color.Adult(iContrast + 1, :))
    
    legend({['Adolescent ' names.Plot{iContrast}], ['Adolescent ' names.Plot{iContrast + 1}], ['Adult ' names.Plot{iContrast}], ['Adult ' names.Plot{iContrast + 1}]},'FontSize', 10,'Location','northwest')
    
    % Create STD Error Bars
    errorbar(x1, y1, e1, '.', 'Color', 'black')
    errorbar(x2, y2, e2, '.', 'Color', 'black')
    errorbar(x3, y3, e3, '.', 'Color', 'black')
    errorbar(x4, y4, e4, '.', 'Color', 'black')
    
    % Set Labels and axis ticks, and limits
    ax = gca;
    ax.YLabel.String = 'Activation (Raw Value)'; %set label to Y axis
    ax.YLabel.FontSize = 12;
 
    ylim([-0.1, 1])
    xlim([0, 5])
    xticks(1:4)
    xticklabels(names.ROI);
        
    grid on
    grid minor
    
end
%%
cd(folder.Scripts)