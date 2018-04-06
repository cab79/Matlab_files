% Run SASICA and view components
% CLEAR ALL BEFORE FIRST RUN
startf=1;
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

anapath = 'C:\Data\PET-LEP\Preprocessed\ICA';
files = dir('*_ICA.set');
fname_ext = '';

file = files(f).name;
C = strsplit(file,'_');

EEG = pop_loadset('filename',file,'filepath',anapath);

def=SASICAsettings;
EEG = eeg_SASICA(EEG,def);
