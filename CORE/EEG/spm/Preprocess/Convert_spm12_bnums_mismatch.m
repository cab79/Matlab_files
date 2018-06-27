clear all

%% conditions

% mismatch conditions
mmconds = [
    1,3;
    2,4;
    5,7;
    6,8;
    9,11; 
    10,12; 
    13,15;
    14,16;
    17,19;
    18,20;
    21,23;
    22,24];

% blocks
bconds = [
    1
    2
    3
];

leftright = {
    [1:4 9:12 17:20]; %left
    [5:8 13:16 21:24]; %right
    };  

%% SPECIFY DATA

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
filepath = 'C:\Data\CORE\Preprocessed_100Hz'; 
outpath = 'C:\Data\CORE\SPMdata'; 
outprefix = 'spm12_blockmismatch_';

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

%for f = sort(1:length(files),'descend')
for f=1:length(files)
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
        
        % create new cond labels
        newconds = NaN(1,length(conds));
        i = 0;
        for b = 1:length(bconds)
            bi = find(ismember(bnums,bconds(b)));
            for lr = 1:length(leftright)
                lri = find(ismember(conds,leftright{lr}));
                for m = 1:size(mmconds,2)
                    mi = find(ismember(conds,mmconds(:,m)'));
                    ii = intersect(intersect(bi,mi),lri);
                    i = i+1;
                    newconds(ii) = i;
                end
            end
        end
        
        S.conditionlabels = cellfun(@num2str, num2cell(newconds), 'UniformOutput', false);
    
        %spm_eeg_convert(S);
        spm_eeg_convert_eeglab_epoched(S,1);  % faster version of spm_eeg_convert, only works with EEGLAB epoched data
    end
    clear EEG;
end
