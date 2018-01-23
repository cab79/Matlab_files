clear all
rmpath(genpath('C:\Data\Matlab\fieldtrip-20170113'))
addpath('C:\Data\Matlab\spm12')
addpath(genpath('C:\Data\Matlab\spm12\external'));

%% SPECIFY DATA

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
filepath = 'C:\Data\Catastrophising study\Preprocessed'; 
outpath = 'C:\Data\Catastrophising study\SPMdata'; 
outprefix = 'spm12_csd_';

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
fpref = '';
fmid = '';
fsuff = 'cleaned_SPNall.set';

% set data type: 'epoched' or 'continuous'
dattype ='epoched';

% set time window of interest in ms (advise whole epoch)
timewin = [];

% Apply CSD (leave as 0 unless you know what you are doing!)
applyCSD=1; 
% path to CSD montage (only needed if using CSD)
if applyCSD; load('C:\Data\Catastrophising study\Orig\CSDmontage_64.mat'); end

%% LOAD AND CONVERT
files = dir(fullfile(filepath,[fpref '*' fmid  '*' fsuff]));
cd(outpath)

for f = 1:length(files)
    S=struct;
    S.dataset = fullfile(filepath,files(f).name);
    [pth nme ext] = fileparts(files(f).name);
    S.outfile = fullfile(outpath,[outprefix nme]);
    S.mode = dattype;
    S.timewin = timewin;
    
    % list events
    EEG = pop_loadset(files(f).name,filepath);
    
    if applyCSD
        [nchan,nsamp,ntrial] = size(EEG.data);
        EEG.data = reshape(EEG.data,nchan,nsamp*ntrial);
        EEG.data=CSD(EEG.data,G,H); 
        EEG.data = reshape(EEG.data,nchan,nsamp,ntrial);
    end
    
    S.conditionlabels = {EEG.epoch.eventtype};
    
    % check if chanlocs names are all uppercase - if so needs modifying for
    % Fieldtrip functions to recognise the channels (yep, silly I know)
    if ~isempty (cell2mat(strfind({EEG.chanlocs.labels},'CZ')))
        for chan = 1:EEG.nbchan
            chlab = EEG.chanlocs(chan).labels;
            if strfind(chlab,'FP')
                EEG.chanlocs(chan).labels = [chlab(1) lower(chlab(2:end))];
            elseif strfind(chlab,'Z')
                EEG.chanlocs(chan).labels = [chlab(1:end-1) lower(chlab(end))];
            end
        end
        EEG = pop_saveset(EEG,'filename',files(f).name,'filepath',filepath); 
    end
    
    spm_eeg_convert(S);
    clear EEG;
end
