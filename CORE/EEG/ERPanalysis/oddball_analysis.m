%% ODDBALL ANALYSIS
%% Settings

clear all
close all
%dbstop if error

% choose factor of interest
S.factcon = 'Odd';
S.runname = 'run2';

% time window to create temporal region of interest, e.g. decided from
% group statistics. Should be kept broad to accomodate all subjects
S.timewin = [33 73]; 
% central EEG channel to create spatial region of interest on scalp, e.g. decided from
% group statistics
S.centrechan= 'E93'; 
% number of neighbours to include in spatial region
S.nbr_levels= 2; 
% masking threshold to apply to individual subject ERP contrasts: percentage area of space * time       
S.mask_thresh=0.2; % e.g. if 0.2 = 20%
% baseline to correct the difference wave to
S.baseline = [-200 0];



%% SPECIFY DATA

addpath('C:\Data\Matlab\fieldtrip-20170113')

% Use previously saved data?
S.use_data=1;

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
S.data_path = 'C:\Data\CORE\Preprocessed_100Hz'; 
S.outpath = 'C:\Data\CORE\ERPs';

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
S.fpref = '';
S.fmid = '';
S.fsuff_bal = '4_cleaned_tm.set'; % use balanced data to analyse ERP contrasts
S.fsuff = '4_merged_cleaned.set'; % use unbalanced data to extract single-trials

% load .xlsx file containing 'Participant_ID', 'Group', and covariates
S.pdatfile = 'C:\Data\CORE\Participant_data.xlsx';
% names of headers in the above xls file:
    S.subhead = 'Subject';
    S.grphead = 'Group';
    S.inchead = 'Include';
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];

% conditions and factors
S.conds = 1:24;
S.useconds = 1:24;
S.flipside = 2;
S.factors = {'CP', 'Side', 'Odd', 'DC'}; % oly need to label those that 
S.factlev={{};{};{'Odd','Stan'};{}}; % label those that will be analyses
S.factor_matrix = [
              1 1 1 1
              1 1 1 2
              1 1 2 1
              1 1 2 2
              1 2 1 1
              1 2 1 2
              1 2 2 1
              1 2 2 2
              2 1 1 1
              2 1 1 2
              2 1 2 1
              2 1 2 2
              2 2 1 1
              2 2 1 2
              2 2 2 1
              2 2 2 2
              3 1 1 1
              3 1 1 2
              3 1 2 1
              3 1 2 2
              3 2 1 1
              3 2 1 2
              3 2 2 1
              3 2 2 2
              ];
          

%% get timelock data

% find data
files_bal = dir(fullfile(S.data_path,[S.fpref '*' S.fmid  '*' S.fsuff_bal]));
files = dir(fullfile(S.data_path,[S.fpref '*' S.fmid  '*' S.fsuff]));
cd(S.outpath)

if length(files) ~=length(files_bal)
    error('different number of full and balanced trials');
end

% CREATE FT CHANNEL LAYOUT STRUCTURE
% only need to ever run this once per experiment, and only if using ACSTP
if 0 
    % create FT channel layout
    cfglayout = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.sfp';
    cfgoutput = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.lay';
    cfg.layout = cfglayout;
    cfg.output = cfgoutput;
    layout = ft_prepare_layout(cfg);
    % define the neighborhood for channels
    cfg.senstype ='EEG';
    cfg.method = 'triangulation';
    cfg.elecfile=cfglayout;
    cfg.layout=cfgoutput;
    ft_nbrs = ft_prepare_neighbours(cfg);
    save(fullfile(S.outpath,'ft_nbrs.mat'), 'ft_nbrs');
end

rundir=fullfile(S.outpath,[S.runname '_timewin_' num2str(S.timewin(1)) '_' num2str(S.timewin(2)) '_basewin_' num2str(S.baseline(1)) '_' num2str(S.baseline(2)) '_centrechan_' S.centrechan '_nNeighbours_' num2str(S.nbr_levels) '_maskthresh_' num2str(S.mask_thresh)]);
if ~exist(rundir,'dir')
    mkdir(rundir)
end

%% create and save ERPs from balanced data
erp_sname = ['ERP_diffwave_flip' num2str(S.flipside) '_base_' num2str(S.baseline(1)) '_' num2str(S.baseline(2)) '.mat'];
% group data structure
GD=struct;
if ~exist(fullfile(S.outpath,erp_sname),'file') || ~S.use_data
    for f = 1:length(files_bal)
        
        % load BALANCED data
        data=fullfile(S.data_path,files_bal(f).name);
        EEG=pop_loadset(data);
        sname = strsplit(EEG.filename,'_');

        % whole epoch used if timewin not specified
        if isempty(S.timewin)
            S.timewin = [EEG.xmin EEG.xmax]*1000;
        end

        % obtain trial indices
        [conds, tnums, fnums, bnums] = get_markers(EEG);

        % create flipped version for RHS stimulation trials
        if S.flipside
            EEGflip = flipchan(EEG);
            sideind = find(strcmp(S.factors,'Side'));
            flipind = find(S.factor_matrix(:,sideind)==S.flipside);
            flipcond = S.conds(flipind);
            for i = flipcond
                trialind = find(conds==i);
                EEG.data(:,:,trialind)= EEGflip.data(:,:,trialind);
            end
            clear EEGflip
        end

        % reduce conditions to factors of interest
        factind = find(strcmp(S.factors,S.factcon));
        condana = S.factor_matrix(S.useconds,factind);
        targ = unique(condana)';

        % identify target indices
        ti = nan(size(conds'));
        for i = targ
            factcon{i} = find(condana==i);
            for c = 1:length(factcon{i})
                ti(find(conds==factcon{i}(c)),1)=i;
            end
        end
        filter_cond=~any(isnan(ti),2);
        ti=ti(filter_cond,:);


        %% create ERP difference wave and baseline correct

        % individual waves
        for i = 1:length(targ)
            GD(f).iw(:,:,targ(i)) = mean(EEG.data(:,:,ti==targ(i)),3);
        end

        % difference wave
        GD(f).dw = GD(f).iw(:,:,1)-GD(f).iw(:,:,2);

        % plot
        %figure
        %plot(mean(iw(elec,:,1),1),'r'); hold on;
        %plot(mean(iw(elec,:,2),1),'b'); hold on;
        %plot(mean(dw(elec,:),1),'k');
        %figure
        %topoplot(mean(dw(:,timewin(1):timewin(2)),2),EEG.chanlocs)

        % baseline correct dw
        basewin=dsearchn(EEG.times',S.baseline');
        base = mean(GD(f).dw(:,basewin(1):basewin(2)-1),2);
        GD(f).dw = GD(f).dw - repmat(base,1,size(GD(f).dw,2));

    end
    GD(1).chanlocs = EEG.chanlocs;
    GD(1).times = EEG.times;
    GD(1).S = S; % save settings for future reference
    save(fullfile(S.outpath,erp_sname),'GD');
else
    load(fullfile(S.outpath,erp_sname));
end

%% define channels for spatial mask
allchan = {GD(1).chanlocs.labels};
elec_centre = find(strcmp(S.centrechan,allchan));
load(fullfile(S.outpath,'ft_nbrs.mat'));
if S.nbr_levels==1
    elec = sort([elec_centre find(ismember({ft_nbrs(:).label},ft_nbrs(elec_centre).neighblabel))]);
elseif S.nbr_levels==2
    temp = sort([elec_centre find(ismember({ft_nbrs(:).label},ft_nbrs(elec_centre).neighblabel))]);
    temp2=[];
    for e = 1:length(temp)
        temp2 = [temp2 find(ismember({ft_nbrs(:).label},ft_nbrs(temp(e)).neighblabel))];
    end
    elec = unique(temp2);
end

%% create single-subject mask
MD = struct;
for f = 1:length(files_bal)
    %% create spatio-temporal mask to apply to each subject's ERP DW
    
    % spatial mask
    M(f).spatmask = zeros(size(GD(f).dw));
    M(f).spatmask(elec,:) = 1; 
    
    % temporal mask
    M(f).tempmask = zeros(size(GD(f).dw));
    timewin=dsearchn(GD(1).times',S.timewin');
    M(f).tempmask(:,timewin(1):timewin(2)) = 1; 

    % combine spatial and temporal
    M(f).erp_mask = M(f).tempmask .* M(f).spatmask;
    %figure
    %plot(mean(squeeze(mean(erp_mask,1)),3));
    %figure
    %topoplot(mean(squeeze(mean(erp_mask(:,timewin(1):timewin(2),:),2)),2),EEG.chanlocs)
    
    %% apply mask to ERP difference wave to identify individual subjects' spatiotemporal ROI
    % N data points in mask
    nd_mask = length(find(M(f).erp_mask==1));
    temp = reshape((GD(f).dw .* M(f).erp_mask),1,[]);
    M(f).mask_erp = zeros(size(temp));
    [~,idx] = sort(temp(:),'descend');
    erp_mask2 = idx(1:ceil(nd_mask*S.mask_thresh));
    M(f).mask_erp(erp_mask2)=1;
    M(f).mask_erp = reshape(M(f).mask_erp,size(GD(f).dw));
    %length(find(M(f).mask_erp))
    
end
dtstr = datestr(now,30);
mask_sname = ['ERP_mask_' dtstr];
save(fullfile(rundir,mask_sname),'M');

%% load single-trial data
st_sname = ['EEG_trials_flip' num2str(S.flipside) '.mat'];
% group single-trial structure
SS=struct;    

% for diagnostics
TS=struct;
Nplots=length(files);
Nx=ceil(sqrt(Nplots));
Ny=ceil(Nplots/Nx);
%fig=figure;

for f = 1:length(files)
    
    fname=strrep(files(f).name,S.fsuff,st_sname);
    %if ~exist(fullfile(S.outpath,fname),'file') || ~S.use_data
        
        % load UNBALANCED data
        data=fullfile(S.data_path,files(f).name);
        EEG=pop_loadset(data);
        %sname = strsplit(EEG.filename,'_');

        % obtain trial indices
        [conds, tnums, fnums, bnums] = get_markers(EEG);

        % create flipped version for RHS stimulation trials
        if S.flipside
            EEGflip = flipchan(EEG);
            sideind = find(strcmp(S.factors,'Side'));
            flipind = find(S.factor_matrix(:,sideind)==S.flipside);
            flipcond = S.conds(flipind);
            for i = flipcond
                trialind = find(conds==i);
                EEG.data(:,:,trialind)= EEGflip.data(:,:,trialind);
            end
            clear EEGflip
        end

        % reduce conditions to factors of interest
        factind = find(strcmp(S.factors,S.factcon));
        condana = S.factor_matrix(S.useconds,factind);
        targ = unique(condana)';

        % identify target indices
        SS.ti = nan(size(conds'));
        for i = targ
            factcon{i} = find(condana==i);
            for c = 1:length(factcon{i})
                SS.ti(find(conds==factcon{i}(c)),1)=i;
            end
        end
        filter_cond=~any(isnan(SS.ti),2);
        SS.ti=SS.ti(filter_cond,:);

        %% create single-trial difference waves and baseline correct

        % difference wave
        SS.dw = EEG.data-repmat(GD(f).iw(:,:,2),1,1,size(EEG.data,3));

        % baseline correct dw
        basewin=dsearchn(EEG.times',S.baseline');
        base = mean(SS.dw(:,basewin(1):basewin(2)-1,:),2);
        SS.dw = SS.dw - repmat(base,1,size(SS.dw,2),1);
        
        save(fullfile(S.outpath,fname),'SS','tnums');
    %else
    %    load(fullfile(S.outpath,fname));
    %    disp(['loading ' fname]);
    %end
    
    %% Apply single-subject mask to single-trial data to create deviant signal projection over time
    MM=struct;
    for t = 1:size(SS.dw,3)
        temp=SS.dw(:,:,t) .* M(f).mask_erp;
        MM.mm(t,1) = mean(temp(:));
    end
    
    mm_sname=strrep(files(f).name,S.fsuff,['mm_proj.mat']);
    save(fullfile(rundir,mm_sname),'MM','tnums');

    %% Diagnostics
    %for f = 1:length(files)
    % generate t stats
    [TS(f).h,TS(f).p,ci,stats] = ttest2(MM.mm(SS.ti==1),MM.mm(SS.ti==2))
    subplot(Nx,Ny,f); scatter(SS.ti,MM.mm'); hold on
    %end
    %plot(MM(f).mm)
    
end
%mm_sname = ['mm_proj_' dtstr];
%save(fullfile(rundir,mm_sname),'MM');
ts_sname = ['Tstat_' dtstr];
save(fullfile(rundir,ts_sname),'TS');


%figure
%topoplot(mean(mean(giw(:,timewin(1):timewin(2),1,:),2),3),EEG.chanlocs)
%figure
%topoplot(mean(mean(giw(:,timewin(1):timewin(2),2,:),2),3),EEG.chanlocs)
%figure
%topoplot(mean(mean(gdw(:,timewin(1):timewin(2),:),2),3),EEG.chanlocs)

%% plot example subject's data for figures
if 0
    clear all
    load('ERP_diffwave_flip2_base_-50_0.mat')
    load('CORE010_mm_proj.mat')
    plot(MM.mm(1:500))
    set(gca,'FontSize', 15);
    load('ERP_mask_20170821T091053.mat')
    load('chanlocs.mat')
    topoplot(any(M(10).mask_erp,2),chanlocs)
    area(GD(1).times,any(M(10).mask_erp,1))
    set(gca,'FontSize', 15);
end