clear all
dbstop if error
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
%% SPECIFY DATA

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
filepath = 'C:\Data\CORE\EEG\ana\prep\cleaned\part4'; 
outpath = 'C:\Data\CORE\EEG\ana\spm\SPMdata'; 
outprefix = 'spm12_fnums_flip_';
datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
fpref = '';
fmid = '';
%fsuff = '_2_merged_cleaned.set';
fsuff = '_4_merged_cleaned.set';

% set data type: 'epoched' or 'continuous'
dattype ='epoched';
overwrite = 1;

% set time window of interest in ms (advise whole epoch)
timewin = [-200 300];

cond2flip = [5:8];
path.chanlocs = 'C:\Data\CORE\eeg\ana\prep\chanlocs.mat';

% conditions to swap
fcond = {[1:4],[5:8]}; % part4 fnum

% get side affected data from participant file
subhead = 'Subject';
grphead = 'Group';
inchead = 'Include';
CRPSsidehead = 'CRPSlr';
include_codes = [1];
[~,~,pdata] = xlsread(datfile);
grp_col = find(strcmp(pdata(1,:),grphead));
sub_col = find(strcmp(pdata(1,:),subhead));
inc_col = find(strcmp(pdata(1,:),inchead));
side_col = find(strcmp(pdata(1,:),CRPSsidehead));
inc_idx = cellfun(@(x) ismember(x,include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));
subjects = pdata(2:end,sub_col);
CRPSsides = pdata(2:end,side_col);

% add random sides to the healthy subjects (who have no CRPS side)
sidenan = cell2mat(cellfun(@(x) any(x), cellfun(@(x) isnan(x), CRPSsides, 'UniformOutput', 0), 'UniformOutput', 0));
sidepool = CRPSsides(~sidenan);
sn_idx=find(sidenan);
for s=1:length(sn_idx)
    CRPSsides(sn_idx(s)) = sidepool(randi([1 numel(sidepool)]));
end


%% LOAD AND CONVERT
fname=[fpref '*' fmid  '*' fsuff];
fname=strrep(fname,'**','*');
files = dir(fullfile(filepath,fname));
allfiles = {files(:).name};
cd(outpath)

%for f = sort(1:length(files),'descend')
for f=1:length(files)
    S=struct; % must clear this every time!
    disp(['file ' num2str(f) '/' num2str(length(files))]);
    S.dataset = fullfile(filepath,files(f).name);
    [pth nme ext] = fileparts(files(f).name);
    S.outfile = fullfile(outpath,[outprefix nme]);
    S.mode = dattype;
    S.timewin = timewin;
    S.inputformat = 'eeglab_set';
    
    % index of this EEG file within participant_data file
    subind=[];
    for i = 1:length(subjects)
        if ~isempty(strfind(allfiles{f},subjects{i}))
            subind=i;
        end
    end
    
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
        
        % flip channels right to left
        if ~isempty(cond2flip)
            load(path.chanlocs);
            flipidx = find(ismember(fnums,cond2flip));
            EEG.data(:,:,flipidx) = flipchan(EEG.data(:,:,flipidx),chanlocs);
            new_file = strrep(files(f).name,'.set','_flip.set');
            pop_saveset(EEG,new_file,filepath)
            S.dataset = fullfile(filepath,new_file);
        end
        
        % swap sides
        if strcmp(CRPSsides{subind},'R')
            new_fnums = nan(1,length(fnums));
            for nf = 1:length(fcond{1})
                new_fnums(fnums==fcond{1}(nf)) = fcond{2}(nf);
                new_fnums(fnums==fcond{2}(nf)) = fcond{1}(nf);
            end
        else
            new_fnums = fnums;
        end
    
        %spm_eeg_convert(S);
        S.conditionlabels = cellfun(@num2str, num2cell(new_fnums), 'UniformOutput', false);
        spm_eeg_convert_eeglab_epoched(S,1);  % faster version of spm_eeg_convert, only works with EEGLAB epoched data
        
        try
            if exist(fullfile(filepath,new_file),'file')
                delete(fullfile(filepath,new_file));
            end
        end
    end
    clear EEG;
end
