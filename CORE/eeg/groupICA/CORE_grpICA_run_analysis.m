toolbox_path = 'C:\Data\Matlab\GroupICATv4.0b\icatb';
study_functions_path = 'C:\Data\Matlab\Matlab_files\CORE\eeg\groupICA';
settings_script = 'CORE_icatb_input_settings';

% RUN
close all
addpath(genpath(toolbox_path));
addpath(study_functions_path);
CORE_icatb_batch_file_run(fullfile(study_functions_path,settings_script));