clear all
% for part 2 data - averaged over the prob factor as data is noisy

dbstop if error
% add toolbox paths
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

%% SPECIFY DATA

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
filepath = 'C:\Data\CORE\EEG\ana\prep\cleaned\part2'; 
outpath = 'C:\Data\CORE\EEG\ana\spm\SPMdata'; 
outprefix = 'spm12_';

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
fpref = '';
fmid = '';
fsuff = '_2_merged_cleaned.set';
%fsuff = '_4_merged_cleaned.set';

% set data type: 'epoched' or 'continuous'
dattype ='epoched';
overwrite = 1;

% set time window of interest in ms (advise whole epoch)
timewin = [-200 900]; % part2
%timewin = [-200 300]; % part4

% conditions within each marker
% mark = {
%     [1 9 17]; % left, Odd, DC1
%     [2 10 18]; % left, Odd, DC3
%     [3 11 19]; % left, stan, DC1
%     [4 12 20]; % left, stan, DC3
%     [5 13 21]; % right, Odd
%     [6 14 22]; % right, Odd
%     [7 15 23]; % right, stan
%     [8 16 24]; % right, Stan
% };
mark = {
    [1 2]; % CP10, left, Odd
    [3 4]; % CP10, left, stan
    [5 6]; % CP10, right, Odd
    [7 8]; % CP10, right, stan
    [9 10]; % CP30, left, Odd
    [11 12]; % CP30, left, stan
    [13 14]; % CP30, right, Odd
    [15 16]; % CP30, right, Stan
    [17 18]; % CP50, left, Odd
    [19 20]; % CP50, left, stan
    [21 22]; % CP50, right, Odd
    [23 24]; % CP50, right, Stan
};

%% LOAD AND CONVERT
fname=[fpref '*' fmid  '*' fsuff];
fname=strrep(fname,'**','*');
files = dir(fullfile(filepath,fname));
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
        cmark = nan(1,length(conds));
        for m = 1:length(mark)
            idx = ismember(conds,mark{m});
            cmark(idx) = m;
        end
        
        S.conditionlabels = cellfun(@num2str, num2cell(cmark), 'UniformOutput', false);
    
        %spm_eeg_convert(S);
        spm_eeg_convert_eeglab_epoched(S,1);  % faster version of spm_eeg_convert, only works with EEGLAB epoched data
    end
    clear EEG;
end
