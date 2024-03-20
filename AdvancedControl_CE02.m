%% 2.1 Multiplicative uncertainty

% Define the folder path
folder_path = '/Données EPFL/';

% Get a list of all files in the folder
files = dir(fullfile(folder_path, '*.mat'));

% Loop through each file
for i = 1:length(files)
    % Get the file name
    file_name = files(i).name;
    
    % Load the .mat file
    data = load(fullfile(folder_path, file_name));
end

% Sampling period
Te = 0.002; % seconds
