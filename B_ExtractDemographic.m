%% Clear the work space
clear; clc

%% Set data directories
cd('..')
folder.Root    = pwd;
folder.Results = fullfile(folder.Root, 'Behavioral', 'Results'); %where we are saving the output
folder.Scripts = fullfile(folder.Root, 'Scripts'); %location of the scripts

%% Set files to load
file.Survey = fullfile(folder.Results, 'QualtricSurvey.csv');

%% Name for output
file.Demographic = fullfile(folder.Results, 'Demographic.mat');

%% Get the data
survey.Raw = table2cell(readtable(file.Survey));
survey.Age(:,1) = cell2mat(survey.Raw(:,2));
survey.Age(:,2) = cell2mat(survey.Raw(:,5));

%% Save the data
save(file.Demographic, '-struct', 'survey')

%% Go back to where we started
cd(folder.Scripts)