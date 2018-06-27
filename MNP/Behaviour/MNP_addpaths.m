function MNP_addpaths
%% ADD FUNCTIONS/TOOLBOXES TO MATLAB PATH
paths = {'C:\Data\Matlab\Matlab_files\MNP','C:\Data\Matlab\Matlab_files\_generic_eeglab_batch\eeglab_batch_supporting_functions'};
subpaths = [1 1]; % add subdirectories too?

for p = 1:length(paths)
    if subpaths(p)
        addpath(genpath(paths{p}));
    else
        addpath(paths{p});
    end
end