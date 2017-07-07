clear all

%% SPECIFY DATA

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
filepath = 'C:\CORE\Preprocessed_100Hz'; 
outpath = 'C:\CORE\SPMdata'; 
outprefix = 'spm12_';

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
fpref = '';
fmid = '';
%fsuff = '_4_merged_cleaned.set';
fsuff = '_4_cleaned_tm.set';

% set data type: 'epoched' or 'continuous'
dattype ='epoched';
overwrite = 0;

% set time window of interest in ms (advise whole epoch)
timewin = [-200 300];

%% LOAD AND CONVERT
files = dir(fullfile(filepath,[fpref '*' fmid  '*' fsuff]));
cd(outpath)

for f = sort(1:length(files),'descend')
%for 1:length(files)
    S=struct;
    S.dataset = fullfile(filepath,files(f).name);
    [pth nme ext] = fileparts(files(f).name);
    S.outfile = fullfile(outpath,[outprefix nme]);
    S.mode = dattype;
    S.timewin = timewin;
    S.inputformat = 'eeglab_set';
    
    if ~exist([S.outfile '.mat'],'file') || overwrite
        % load data
        EEG = pop_loadset(files(f).name,filepath);

        % check if chanlocs names are all uppercase - if so needs modifying for
        % Fieldtrip functions to recognise the channels (yep, silly I know)
        %S.conditionlabels = {EEG.epoch.eventtype}; % list events
        %if ~isempty (cell2mat(strfind({EEG.chanlocs.labels},'CZ')))
        %    for chan = 1:EEG.nbchan
        %        chlab = EEG.chanlocs(chan).labels;
        %        if strfind(chlab,'FP')
        %            EEG.chanlocs(chan).labels = [chlab(1) lower(chlab(2:end))];
        %        elseif strfind(chlab,'Z')
        %            EEG.chanlocs(chan).labels = [chlab(1:end-1) lower(chlab(end))];
        %        end
        %    end
        %    EEG = pop_saveset(EEG,'filename',files(f).name,'filepath',filepath); 
        %end

        % for data recorded with EGI system and STIM/DIN markers
        [conds, tnums, fnums, bnums] = get_markers(EEG);
        S.conditionlabels = cellfun(@num2str, num2cell(conds), 'UniformOutput', false);
    
        %spm_eeg_convert(S);
        spm_eeg_convert_eeglab_epoched(S,1);  % faster version of spm_eeg_convert, only works with EEGLAB epoched data
    end
    clear EEG;
end
