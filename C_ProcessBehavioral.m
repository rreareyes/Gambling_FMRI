%% Extract behavioral data Gamble fMRI

%% Authorship
% Created by Eduardo Rea for project Gamble fMRI
% NLP Lab UMass Amherst
% August 2017
% working on SPM8

%% Clean workspace
clear; clc;

%% Base paths
cd('..')
folder.Root     = pwd;
folder.Scripts  = fullfile (folder.Root, 'Scripts');
folder.Behavior = fullfile (folder.Root, 'Behavioral');
folder.Results  = fullfile (folder.Behavior, 'Results');
folder.Graphs   = fullfile(folder.Behavior, 'Results');

if ~exist(folder.Results,'dir')
    mkdir(folder.Results);
end

if ~exist(folder.Graphs,'dir')
    mkdir(folder.Graphs);
end

%% Paths for outputs
graph.Structured = fullfile(folder.Graphs, 'ProportionStructured.fig');
graph.Neg        = fullfile(folder.Graphs, 'ProportionNeg.fig');
graph.Zero       = fullfile(folder.Graphs, 'ProportionZero.fig');
graph.Age	     = fullfile(folder.Graphs, 'Age.fig');
graph.WASI       = fullfile(folder.Graphs, 'WASI.fig');
graph.AgeVSlMin  = fullfile(folder.Graphs, 'AgeVSlMin.fig');
graph.AgeVSgMax  = fullfile(folder.Graphs, 'AgeVSgMax.fig');
graph.AgeVSpMax  = fullfile(folder.Graphs, 'AgeVSpMax.fig');
graph.AgeVSlMin  = fullfile(folder.Graphs, 'AgeVSlMin.fig');

%% Get files to load
file.Behavioral = dir(fullfile(folder.Behavior, '*.mat'));
file.WASI       = fullfile(folder.Graphs, 'WASI.mat');
file.Survey     = fullfile(folder.Graphs, 'Demographic.mat');
file.Summary    = fullfile(folder.Graphs, 'Summary.mat');

for indvFile = 1:length(file.Behavioral)
    file.List{indvFile,:} = fullfile(file.Behavioral(indvFile).folder, file.Behavioral(indvFile).name);
end

%% Set base data for loop
data.nBlocks = 5;
data.nSubjects = length(file.List) ./ data.nBlocks;
data.SubjID = string(unique(extractBetween(file.List,'ET_1_S','_block')));

%% Loop trough all subjects and extract behavioral data
for subject = 1:data.nSubjects
    for block = 1:data.nBlocks
        %% Load files
        file.Name = fullfile(folder.Behavior, ['/Gamble_ET_1_S' char(data.SubjID(subject)) '_block' num2str(block) '.mat']);
        load(file.Name, 'stim_choice')
        data.Raw{subject,1} = data.SubjID{subject};
        data.Raw{subject,block + 1} = stim_choice;
        clear stim_choice
        
        %% Get number of structured trials
        data.StructLoc = contains({data.Raw{subject, block+1}(1:end).type},'trad').';
        data.nStruct(subject,block) = sum(data.StructLoc);
        
        data.TradNegLoc = contains({data.Raw{subject, block+1}(1:end).type},'trad_neg').';
        data.nTradNeg(subject,block) = sum(data.TradNegLoc);
                
        data.TradZeroLoc = contains({data.Raw{subject, block+1}(1:end).type},'trad_zero').';
        data.nTradZero(subject,block) = sum(data.TradZeroLoc);
        
        %% Get participant's selection on structured trials
        data.gMax(subject,1) = str2double(data.SubjID(subject));
        data.pMax(subject,1) = str2double(data.SubjID(subject));
        data.lMin(subject,1) = str2double(data.SubjID(subject));
        
        data.gMax(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(:,:).gmax}) == cell2mat({data.Raw{subject, block+1}(:,:).resp_num}));
        data.pMax(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(:,:).pmax}) == cell2mat({data.Raw{subject, block+1}(:,:).resp_num}));
        data.lMin(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(:,:).lmin}) == cell2mat({data.Raw{subject, block+1}(:,:).resp_num}));
        
        data.gMaxNeg(subject,1) = str2double(data.SubjID(subject));
        data.pMaxNeg(subject,1) = str2double(data.SubjID(subject));
        data.lMinNeg(subject,1) = str2double(data.SubjID(subject));
        
        data.gMaxNeg(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(data.TradNegLoc).gmax}) == cell2mat({data.Raw{subject, block+1}(data.TradNegLoc).resp_num}));
        data.pMaxNeg(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(data.TradNegLoc).pmax}) == cell2mat({data.Raw{subject, block+1}(data.TradNegLoc).resp_num}));
        data.lMinNeg(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(data.TradNegLoc).lmin}) == cell2mat({data.Raw{subject, block+1}(data.TradNegLoc).resp_num}));
        
        data.gMaxZero(subject,1) = str2double(data.SubjID(subject));
        data.pMaxZero(subject,1) = str2double(data.SubjID(subject));
        data.lMinZero(subject,1) = str2double(data.SubjID(subject));
        
        data.gMaxZero(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(data.TradZeroLoc).gmax}) == cell2mat({data.Raw{subject, block+1}(data.TradZeroLoc).resp_num}));
        data.pMaxZero(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(data.TradZeroLoc).pmax}) == cell2mat({data.Raw{subject, block+1}(data.TradZeroLoc).resp_num}));
        data.lMinZero(subject,block+1) = sum(cell2mat({data.Raw{subject, block+1}(data.TradZeroLoc).lmin}) == cell2mat({data.Raw{subject, block+1}(data.TradZeroLoc).resp_num}));
        
        %% Get structured trials with no response
        data.Empty(subject,block)     = sum(isnan(cell2mat({data.Raw{subject, block+1}(data.StructLoc).resp_num})));
        data.EmptyNeg(subject,block)  = sum(isnan(cell2mat({data.Raw{subject, block+1}(data.TradNegLoc).resp_num})));
        data.EmptyZero(subject,block) = sum(isnan(cell2mat({data.Raw{subject, block+1}(data.TradZeroLoc).resp_num})));
                                      
    end
    %% Get the number choices across blocks
    data.gMax(subject,data.nBlocks+2) = sum(data.gMax(subject,2:6),2);
    data.pMax(subject,data.nBlocks+2) = sum(data.pMax(subject,2:6),2);
    data.lMin(subject,data.nBlocks+2) = sum(data.lMin(subject,2:6),2);
    
    data.gMaxNeg(subject,data.nBlocks+2) = sum(data.gMaxNeg(subject,2:6),2);
    data.pMaxNeg(subject,data.nBlocks+2) = sum(data.pMaxNeg(subject,2:6),2);
    data.lMinNeg(subject,data.nBlocks+2) = sum(data.lMinNeg(subject,2:6),2);
    
    data.gMaxZero(subject,data.nBlocks+2) = sum(data.gMaxZero(subject,2:6),2);
    data.pMaxZero(subject,data.nBlocks+2) = sum(data.pMaxZero(subject,2:6),2);
    data.lMinZero(subject,data.nBlocks+2) = sum(data.lMinZero(subject,2:6),2);
    
    
    %% Get number of valid trials across blocks
    data.ValidTrials(subject,1)     = sum(data.nStruct(subject,:),2) - sum(data.Empty(subject,:),2);
    data.ValidTrialsNeg(subject,1)  = sum(data.nTradNeg(subject,:),2) - sum(data.EmptyNeg(subject,:),2);
    data.ValidTrialsZero(subject,1) = sum(data.nTradZero(subject,:),2) - sum(data.EmptyZero(subject,:),2);
    
    %% Calculate proportion of choices across valid trials
    data.gMax(subject,data.nBlocks+3) = data.gMax(subject,data.nBlocks+2) ./ data.ValidTrials(subject);
    data.pMax(subject,data.nBlocks+3) = data.pMax(subject,data.nBlocks+2) ./ data.ValidTrials(subject);
    data.lMin(subject,data.nBlocks+3) = data.lMin(subject,data.nBlocks+2) ./ data.ValidTrials(subject);
    
    data.gMaxNeg(subject,data.nBlocks+3) = data.gMaxNeg(subject,data.nBlocks+2) ./ data.ValidTrialsNeg(subject);
    data.pMaxNeg(subject,data.nBlocks+3) = data.pMaxNeg(subject,data.nBlocks+2) ./ data.ValidTrialsNeg(subject);
    data.lMinNeg(subject,data.nBlocks+3) = data.lMinNeg(subject,data.nBlocks+2) ./ data.ValidTrialsNeg(subject);
    
    data.gMaxZero(subject,data.nBlocks+3) = data.gMaxZero(subject,data.nBlocks+2) ./ data.ValidTrialsZero(subject);
    data.pMaxZero(subject,data.nBlocks+3) = data.pMaxZero(subject,data.nBlocks+2) ./ data.ValidTrialsZero(subject);
    data.lMinZero(subject,data.nBlocks+3) = data.lMinZero(subject,data.nBlocks+2) ./ data.ValidTrialsZero(subject);

end

%% Retrieve proportion of choice
data.Summary     = [data.gMax(:,[1 8]), data.pMax(:,8), data.lMin(:,8),];
data.SummaryNeg  = [data.gMaxNeg(:,[1 8]), data.pMaxNeg(:,8), data.lMinNeg(:,8),];
data.SummaryZero = [data.gMaxZero(:,[1 8]), data.pMaxZero(:,8), data.lMinZero(:,8),];

%% Sort Rows to ease visualization
data.Summary     = sortrows(data.Summary);
data.SummaryNeg  = sortrows(data.SummaryNeg);
data.SummaryZero = sortrows(data.SummaryZero);

%% Save the data for backup
save(file.Summary, '-struct', 'data')

%% Identify individuals from every group
data.Adolescents = data.Summary(:,1) >= 100;
data.Adults = data.Summary(:,1) < 100;

%% Load Demographics and WASI
load(file.WASI, 'Scores');
load(file.Survey, 'Age');

%% Get which individuals completed the task
data.ValidScores = ismember(Scores(:,1),data.Summary(:,1));
data.ValidAge    = ismember(Age(:,1),data.Summary(:,1));

data.WASI = Scores(data.ValidScores,:);
data.Age = Age(data.ValidAge,:);

%% Insert data in the summaries
data.Group(data.Adolescents,1) = {'Adolescents'};
data.Group(data.Adults,1) = {'Adults'};

%% Calculate means and errors for both groups
data.Means(1,:) = mean(data.Summary(data.Adolescents,2:4));
data.Means(2,:) = mean(data.Summary(data.Adults,2:4));

data.Errors(1,:) = std(data.Summary(data.Adolescents,2:4)) ./ sqrt(length(data.Summary(data.Adolescents,:)));
data.Errors(2,:) = std(data.Summary(data.Adults,2:4)) ./ sqrt(length(data.Summary(data.Adults,:)));

data.MeansNeg(1,:) = mean(data.SummaryNeg(data.Adolescents,2:4));
data.MeansNeg(2,:) = mean(data.SummaryNeg(data.Adults,2:4));

data.ErrorsNeg(1,:) = std(data.SummaryNeg(data.Adolescents,2:4)) ./ sqrt(length(data.Summary(data.Adolescents,:)));
data.ErrorsNeg(2,:) = std(data.SummaryNeg(data.Adults,2:4)) ./ sqrt(length(data.Summary(data.Adults,:)));

data.MeansZero(1,:) = mean(data.SummaryZero(data.Adolescents,2:4));
data.MeansZero(2,:) = mean(data.SummaryZero(data.Adults,2:4));

data.ErrorsZero(1,:) = std(data.SummaryZero(data.Adolescents,2:4)) ./ sqrt(length(data.Summary(data.Adolescents,:)));
data.ErrorsZero(2,:) = std(data.SummaryZero(data.Adults,2:4)) ./ sqrt(length(data.Summary(data.Adults,:)));

%% Plot Ages
figure('Name','Age') %name of the figure
boxplot(data.Age(:,2),data.Group)

% Save the plot
saveas (gcf, graph.Age)

%% Plot WASI
figure('Name','WASI') %name of the figure
boxplot(data.WASI(:,2),data.Group)

% Save the plot
saveas (gcf, graph.WASI)

%% Plot Neg trials
errorbar_groups(data.MeansNeg,data.ErrorsNeg, 'bar_width', 0.7, 'errorbar_width', 1,'bar_colors', [0.529412 0.807843 0.921569; 0.27451 0.509804 0.705882],'optional_bar_arguments',{'EdgeColor', 'none'});
title('Proportion of Choice in Neg Trials')
xticklabels({'gmax','pmax','lmin'})
ylabel('Choice Proportion')
legend('Adolescents','Adults')
grid on
grid minor

% Save the plot
saveas (gcf, graph.Neg)

%% Plot Zero trials
errorbar_groups(data.MeansZero,data.ErrorsZero, 'bar_width', 0.7, 'errorbar_width', 1,'bar_colors', [1 0.870588 0.678431; 0.956863 0.643137 0.376471],'optional_bar_arguments',{'EdgeColor', 'none'});
title('Proportion of Choice in Zero Trials')
xticklabels({'gmax','pmax','lmin'})
ylabel('Choice Proportion')
legend('Adolescents','Adults')
grid on
grid minor

% Save the plot
saveas (gcf, graph.Zero)

%% Plot Structured trials
errorbar_groups(data.Means,data.Errors, 'bar_width', 0.7, 'errorbar_width', 1,'bar_colors', [0.443137 0.776471 0.443137; 0.180392 0.545098 0.341176],'optional_bar_arguments',{'EdgeColor', 'none'});
title('Proportion of Choice in Structured Trials')
xticklabels({'gmax','pmax','lmin'})
ylabel('Choice Proportion')
legend('Adolescents','Adults')
grid on
grid minor

% Save the plot
saveas (gcf, graph.Structured)

%% Plot Age vs GMax choices
figure('Name','Age vs GMax') %name of the figure
scatter(data.Age(:,2),data.Summary(:,2), 'Marker', '*','LineWidth',1.5,'MarkerFaceColor', [0.556863 0.219608 0.556863], 'MarkerEdgeColor', [0.556863 0.219608 0.556863])
refLine = lsline;
refLine.LineWidth = 1;
refLine.Color = 'black';
grid on
grid minor

% Save the plot
saveas (gcf, graph.AgeVSgMax)

%% Plot Age vs PMax choices
figure('Name','Age vs PMax') %name of the figure
scatter(data.Age(:,2),data.Summary(:,3), 'Marker', '*','LineWidth',1.5,'MarkerFaceColor', [0.27451 0.509804 0.705882], 'MarkerEdgeColor', [0.27451 0.509804 0.705882])
refLine = lsline;
refLine.LineWidth = 1;
refLine.Color = 'black';
grid on
grid minor

% Save the plot
saveas (gcf, graph.AgeVSpMax)

%% Plot Age vs LMin choices
figure('Name','Age vs LMin') %name of the figure
scatter(data.Age(:,2),data.Summary(:,4), 'Marker', '*','LineWidth',1.5,'MarkerFaceColor', [0.933333 0.509804 0.933333], 'MarkerEdgeColor', [0.933333 0.509804 0.933333])
refLine = lsline;
refLine.LineWidth = 1;
refLine.Color = 'black';
grid on
grid minor

% Save the plot
saveas (gcf, graph.AgeVSlMin)

%% Get back to the starting point
cd(folder.Scripts)