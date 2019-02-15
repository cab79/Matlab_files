close all
clear all
dbstop if error
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
load('C:\Data\CORE\eeg\ana\prep\chanlocs.mat')
sample_res = 4;
xticks = -200:sample_res:900-sample_res;
topo_range = [min(xticks) max(xticks)];
eeglab_plot=0;
save_sim=1;

%% unpaired contrasts (cond numbers or subject indices to average over)
con.grp.i = {1:22, 23:44};
con.odd.i = {[1 2 5 6 9 10 13 14 17 18 21 22], [3 4 7 8 11 12 15 16 19 20 23 24]};
con.cp.i = {[1 2 5 6], [17 18 21 22]}; %oddballs only
con.dc.i = {1:4:19, 2:4:20}; %oddballs only

%% SET PATHS
%S.file.hgf = 'D_fit_r1_it7_pm3_rm4.mat';
%S.file.get_dt = 'CORE_fittedparameters_percmodel2_respmodel4_fractrain0_20190211T074650.mat';
% S.file.stats = 'stats_BRR_all_chan_cond_notrans_20190214T220713';
S.file.stats = 'stats_BRR_all_chan_HGF_notrans_20190214T220715';

% ERP
S.path.main = 'C:\Data\CORE\eeg\ana';
S.path.hgf = 'C:\Data\CORE\behaviour\hgf\fitted'; 
S.path.eeg = [S.path.main '\prep\cleaned\part2'];
S.path.stats = [S.path.main '\stats']; % folder to save outputs
S.path.design = ['C:\Data\CORE\design']; % 
S.path.datfile = 'C:\Data\CORE\Participants\Participant_data.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
S.path.chanlocs = 'C:\Data\CORE\eeg\ana\prep\chanlocs.mat';
S.path.GSNlocs = 'C:\Data\CORE\eeg\GSN-HydroCel-128-Flipmap.mat';
S.fname.parts = {'subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'set'}; 
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.load.suffixes = {'2_merged_cleaned'}; 

% load regression stats data and create simulation of EEG
%load(fullfile(S.path.hgf,S.file.hgf));
l = load(fullfile(S.path.stats,S.file.stats));
stats = l.stats;
for d = 1:length(stats.BRR.alldata.b)
    b{d,1} = stats.BRR.alldata.b{d,1};
    br{d,1} = reshape(b{d},size(b{d},1)*size(b{d},2),[]);
    X{d,1} = stats.BRR.alldata.pred{d,1};
    temp = br{d}*X{d}';
    Y{d,1} = reshape(temp,size(b{d},1),size(b{d},2),[]);
    
    % restore to original trial order
    [~,reverse_ind] = sort(stats.trialinfo{1}.idx{d});
    Y{d,1} = Y{d,1}(:,:,reverse_ind);
end

% load actual EEGlab data
S.load.suffix = S.load.suffixes(1);
S.path.file = S.path.eeg;
S = getfilelist(S);
for d = 1:length(S.select.subjects)
    
    % create output D
    D.subname = S.select.subjects{d};

    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,D.subname)))));

    %load
    filename = subfiles{1};
    S.load.suffix = S.load.suffixes{1};
    fprintf('\nImporting %s.\n\n', filename);
    EEG = pop_loadset(filename);

    % for data recorded with EGI system and STIM/DIN markers
    [conds, tnums, fnums, bnums] = get_markers(EEG);
    
    resampled=0;
    if eeglab_plot
        EEG = pop_resample(EEG,EEG.srate/sample_res);
        resampled=1;
        for ii=1:2
            con.odd.dat{d,ii} = mean(EEG.data(:,:,ismember(conds,con.odd.i{ii})),3);
            con.cp.dat{d,ii} = mean(EEG.data(:,:,ismember(conds,con.cp.i{ii})),3);
            con.dc.dat{d,ii} = mean(EEG.data(:,:,ismember(conds,con.dc.i{ii})),3);
        end
    else
        for ii=1:2
            con.odd.dat{d,ii} = mean(Y{d}(:,:,ismember(conds,con.odd.i{ii})),3);
            con.cp.dat{d,ii} = mean(Y{d}(:,:,ismember(conds,con.cp.i{ii})),3);
            con.dc.dat{d,ii} = mean(Y{d}(:,:,ismember(conds,con.dc.i{ii})),3);
        end
    end
    
    % save sim data in EEGlab format for further analysis
    if save_sim
        if ~resampled
            EEG = pop_resample(EEG,EEG.srate/sample_res);
        end
        sdir = fullfile(S.path.main,'sim');
        if ~exist(sdir,'dir')
            mkdir(sdir)
        end
        EEG.data = Y{d,1};
        sname = strrep(filename,'.set',['_' S.file.stats '.set']);
        pop_saveset(EEG,'filename',sname,'filepath',sdir); 
    end
end

%% GROUP contrast
for g = 1:2
    grpdat{g}=cat(3,con.odd.dat{con.grp.i{g},1}); % oddballs
    sizdat = size(grpdat{g});
    grpdat{g}=reshape(grpdat{g},size(grpdat{g},1)*size(grpdat{g},2),[])';
end
for s=1:size(grpdat{1},2)
    [h(s),p(s),~,stt] = ttest2(double(grpdat{1}(:,s)),double(grpdat{2}(:,s)));
    t(s)=stt.tstat;
end
% remove NaNs
t(isnan(t))=0;
h(isnan(h))=0;
p(isnan(p))=Inf;
% FDR correction
[~,fdr_mask] = fdr(p,0.05);
fdr_p=p.*double(fdr_mask);
fdr_t=t.*double(fdr_mask);
% reshape
h=reshape(h,sizdat(1),sizdat(2));
p=reshape(p,sizdat(1),sizdat(2));
t=reshape(t,sizdat(1),sizdat(2));
fdr_t=reshape(fdr_t,sizdat(1),sizdat(2));
fdr_p=reshape(fdr_p,sizdat(1),sizdat(2));
CORE_eeg_plot_stats('group effect',xticks,topo_range,t,fdr_t,chanlocs)

%% ODDBALL contrast
for g = 1:2
    grpdat{g}=cat(3,con.odd.dat{:,g});
    sizdat = size(grpdat{g});
    grpdat{g}=reshape(grpdat{g},size(grpdat{g},1)*size(grpdat{g},2),[])';
end
for s=1:size(grpdat{1},2)
    [h(s),p(s),~,stt] = ttest2(double(grpdat{1}(:,s)),double(grpdat{2}(:,s)));
    t(s)=stt.tstat;
end
% remove NaNs
t(isnan(t))=0;
h(isnan(h))=0;
p(isnan(p))=Inf;
% FDR correction
[~,fdr_mask] = fdr(p,0.05);
fdr_p=p.*double(fdr_mask);
fdr_t=t.*double(fdr_mask);
% reshape
h=reshape(h,sizdat(1),sizdat(2));
p=reshape(p,sizdat(1),sizdat(2));
t=reshape(t,sizdat(1),sizdat(2));
fdr_t=reshape(fdr_t,sizdat(1),sizdat(2));
fdr_p=reshape(fdr_p,sizdat(1),sizdat(2));
CORE_eeg_plot_stats('oddball effect',xticks,topo_range,t,fdr_t,chanlocs)

%% CP contrast
for g = 1:2
    grpdat{g}=cat(3,con.cp.dat{:,g});
    sizdat = size(grpdat{g});
    grpdat{g}=reshape(grpdat{g},size(grpdat{g},1)*size(grpdat{g},2),[])';
end
for s=1:size(grpdat{1},2)
    [h(s),p(s),~,stt] = ttest2(double(grpdat{1}(:,s)),double(grpdat{2}(:,s)));
    t(s)=stt.tstat;
end
% remove NaNs
t(isnan(t))=0;
h(isnan(h))=0;
p(isnan(p))=Inf;
% FDR correction
[~,fdr_mask] = fdr(p,0.05);
fdr_p=p.*double(fdr_mask);
fdr_t=t.*double(fdr_mask);
% reshape
h=reshape(h,sizdat(1),sizdat(2));
p=reshape(p,sizdat(1),sizdat(2));
t=reshape(t,sizdat(1),sizdat(2));
fdr_t=reshape(fdr_t,sizdat(1),sizdat(2));
fdr_p=reshape(fdr_p,sizdat(1),sizdat(2));
CORE_eeg_plot_stats('CP effect',xticks,topo_range,t,fdr_t,chanlocs)

%% DC contrast
for g = 1:2
    grpdat{g}=cat(3,con.dc.dat{:,g});
    sizdat = size(grpdat{g});
    grpdat{g}=reshape(grpdat{g},size(grpdat{g},1)*size(grpdat{g},2),[])';
end
for s=1:size(grpdat{1},2)
    [h(s),p(s),~,stt] = ttest2(double(grpdat{1}(:,s)),double(grpdat{2}(:,s)));
    t(s)=stt.tstat;
end
% remove NaNs
t(isnan(t))=0;
h(isnan(h))=0;
p(isnan(p))=Inf;
% FDR correction
[~,fdr_mask] = fdr(p,0.05);
fdr_p=p.*double(fdr_mask);
fdr_t=t.*double(fdr_mask);
% reshape
h=reshape(h,sizdat(1),sizdat(2));
p=reshape(p,sizdat(1),sizdat(2));
t=reshape(t,sizdat(1),sizdat(2));
fdr_t=reshape(fdr_t,sizdat(1),sizdat(2));
fdr_p=reshape(fdr_p,sizdat(1),sizdat(2));
CORE_eeg_plot_stats('DC effect',xticks,topo_range,t,fdr_t,chanlocs)
