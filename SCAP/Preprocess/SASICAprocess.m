% Run SASICA and view components
% CLEAR ALL BEFORE FIRST RUN
startf=27;
if exist('EEG','var') && exist('f','var')
    pop_saveset(EEG,'filename',EEG.filename,'filepath',EEG.filepath); 
    save f f
    disp('saving previous EEG and file index');
    clear all
    load f
    f=f+1;
else
    f=startf;
end

%eeglab_path = 'C:\Data\Matlab\eeglab14_1_1b'; % path to eeglab toolbox
%addpath(genpath(eeglab_path)); 
anapath = 'C:\Data\SCAP\EEG\Preprocessed\ICA';
settingspath = 'C:\Data\Matlab\Matlab_files\SCAP\Preprocess';
cd(settingspath)
files = dir(fullfile(anapath,'*_ICA.set'));
fname_ext = '';

file = files(f).name
C = strsplit(file,'_');

EEG = pop_loadset('filename',file,'filepath',anapath);

def=SASICAsettings;
EEG = eeg_SASICA(EEG,def);
