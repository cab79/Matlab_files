clear all
dbstop if error
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
%% SPECIFY DATA

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
% filepath = 'C:\Data\CORE\EEG\ana\prep\cleaned\part2'; 
filepath = 'C:\Data\CORE\EEG\ana\sim'; 
outpath = 'C:\Data\CORE\EEG\ana\spm\SPMdata'; 
outprefix = 'spm12_flip_CPavg_';
datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
fpref = '';
fmid = '';
fsuff = '_2_merged_cleaned_stats_BRR_all_chan_HGF_notrans_20190214T220715.set';
% fsuff = '_4_merged_cleaned.set';
% fsuff = '_2_cleaned_tm.set';

% set data type: 'epoched' or 'continuous'
dattype ='epoched';
overwrite = 1;

% set time window of interest in DATAPOINTS (advise whole epoch)
timewin = [-200 900];

cond2flip = []; % prior to averaging over conds
% cond2flip = [5:8 13:16 21:24]; % prior to averaging over conds
path.chanlocs = 'C:\Data\CORE\eeg\ana\prep\chanlocs.mat';
mark = {
    [1 9 17]; % left, Odd, DC1
    [2 10 18]; % left, Odd, DC3
    [3 11 19]; % left, stan, DC1
    [4 12 20]; % left, stan, DC3
    [5 13 21]; % right, Odd
    [6 14 22]; % right, Odd
    [7 15 23]; % right, stan
    [8 16 24]; % right, Stan
};
%% Sideavg
% mark = {
%     [1 5]; % CP10, Odd, DC1
%     [2 6]; % CP10, Odd, DC3
%     [3 7]; % CP10, Stan, DC1
%     [4 8]; % CP10, Stan, DC3
%     [9 13]; % CP30, Odd, DC1
%     [10 14]; % CP30, Odd, DC3
%     [11 15]; % CP30, Stan, DC1
%     [12 16]; % CP30, Stan, DC3
%     [17 21]; % CP50, Odd, DC1
%     [18 22]; % CP50, Odd, DC3
%     [19 23]; % CP50, Stan, DC1
%     [20 24]; % CP50, Stan, DC3
% };

% mark = {
%     [1 2]; % CP10, left, Odd
%     [3 4]; % CP10, left, stan
%     [5 6]; % CP10, right, Odd
%     [7 8]; % CP10, right, stan
%     [9 10]; % CP30, left, Odd
%     [11 12]; % CP30, left, stan
%     [13 14]; % CP30, right, Odd
%     [15 16]; % CP30, right, Stan
%     [17 18]; % CP50, left, Odd
%     [19 20]; % CP50, left, stan
%     [21 22]; % CP50, right, Odd
%     [23 24]; % CP50, right, Stan
% };


% conditions to swap
fcond = {[5:8 13:16 21:24],[1:4 9:12 17:20]}; 
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


        % for data recorded with EGI system and STIM/DIN markers
        [conds, tnums, fnums, bnums] = get_markers(EEG);
        
        % flip channels right to left
        if ~isempty(cond2flip)
            load(path.chanlocs);
            flipidx = find(ismember(conds,cond2flip));
            EEG.data(:,:,flipidx) = flipchan(EEG.data(:,:,flipidx),chanlocs);
            new_file = strrep(files(f).name,'.set','_flip.set');
            pop_saveset(EEG,new_file,filepath)
            S.dataset = fullfile(filepath,new_file);
        else
            S.dataset = fullfile(filepath,files(f).name);
        end
        
        % swap sides
        if strcmp(CRPSsides{subind},'R')
            new_conds = nan(1,length(conds));
            for nf = 1:length(fcond{1})
                new_conds(conds==fcond{1}(nf)) = fcond{2}(nf);
                new_conds(conds==fcond{2}(nf)) = fcond{1}(nf);
            end
            conds = new_conds;
        end
        
        cmark = nan(1,length(conds));
        for m = 1:length(mark)
            idx = ismember(conds,mark{m});
            cmark(idx) = m;
        end
        
        %spm_eeg_convert(S);
        S.conditionlabels = cellfun(@num2str, num2cell(cmark), 'UniformOutput', false);
        spm_eeg_convert_eeglab_epoched(S,1);  % faster version of spm_eeg_convert, only works with EEGLAB epoched data
%         spm_eeg_convert(S);  % SLOW
        
        try
            if exist(fullfile(filepath,new_file),'file')
                delete(fullfile(filepath,new_file));
            end
        end
    end
    clear EEG;
end
