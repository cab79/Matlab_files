% Run SASICA and view components
% CLEAR ALL BEFORE FIRST RUN
startf=35;
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

anapath = 'C:\Data\CORE\Preprocessed_100Hz';
files = dir('*_2_*_1st_ICA.set');
fname_ext = '';
timebins = [-0.2 0.9; % for epoching, TSOT(2)
            -0.2 0.3]; % for epoching, TSOT(4)

file = files(f).name;
C = strsplit(file,'_');

if strcmp(C{2},'2')
    timebin = timebins(1,:);
    %continue
elseif strcmp(C{2},'4')
    timebin = timebins(2,:);
    %continue
end

EEG = pop_loadset('filename',file,'filepath',anapath);

def=SASICAsettings;
EEG = eeg_SASICA(EEG,def);
