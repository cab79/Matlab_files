function CORE_addpaths
%% ADD FUNCTIONS/TOOLBOXES TO MATLAB PATH
paths = {
    0, 'C:\Data\Matlab\Matlab_files';
    1, 'C:\Data\Matlab\Matlab_files\CORE';
    1, 'C:\Data\Matlab\Matlab_files\_generic_eeglab_batch\eeglab_batch_supporting_functions';
    1, 'C:\Data\Matlab\HGF\HGFv5.0';
    1, 'C:\Data\Matlab\Matlab_files\_generic_HGF';
    1, 'C:\Data\Matlab\Matlab_files\_generic_SCIn_data_process';
    1, 'C:\Data\Matlab\VBA-toolbox-master';
    1, 'C:\Data\Matlab\Matlab_files\CORE\EEG\cosmo\modified functions';
    1, 'C:\Data\Matlab\TFCE';
	1, 'C:\Data\Matlab\eeglab14_1_1b';
    1, 'C:\Data\Matlab\PRoNTo_dev-2.0.1\machines\gpml';
    1, 'C:\Data\Matlab\CoSMoMVPA';
    0, 'C:\Data\Matlab\fieldtrip-20180320';
    1, 'C:\Data\Matlab\bayesreg';
    1, 'C:\Data\Matlab\RR_encode';
    0, 'C:\Data\Matlab\spm12';
    };

for p = 1:length(paths)
    if paths{p,1}
        addpath(genpath(paths{p,2}));
    else
        addpath(paths{p,2});
    end
end
try
    rmpath(genpath('C:\Data\Matlab\eeglab14_1_1b\plugins\Biosig3.3.0')); % due to conflict between Fieldtrip and eeglab BIOSIG toolbox
end
try
    ft_defaults; 
end