function MNP_addpaths
%% ADD FUNCTIONS/TOOLBOXES TO MATLAB PATH
paths = {
    'C:\Data\Matlab\Matlab_files\MNP',...
    'C:\Data\Matlab\Matlab_files\_generic_eeglab_batch\eeglab_batch_supporting_functions'...
    'C:\Data\Matlab\HGF\HGFv5.0'...
    'C:\Data\Matlab\Matlab_files\_generic_HGF'...
    'C:\Data\Matlab\Matlab_files\_generic_SCIn_data_process'...
    };
subpaths = [1 1 1 1 1]; % add subdirectories too?

for p = 1:length(paths)
    if subpaths(p)
        addpath(genpath(paths{p}));
    else
        addpath(paths{p});
    end
end