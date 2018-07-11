%% MEEG time-lock searchlight
%
% Using a time-channel neighborhood, a searchlight map is computed
% indicating in time and space (channel) the pre-stimulation period can be
% distinguished from the peri/post stimulation perioS.
%
% The code presented here can be adapted for other MEEG analyses, but
% there are a few potential caveats:
% * assignment of targets (labels of conditions) is based here on
%   stimulation periods versus pre-stimulation periods. In typical
%   analyses the targets should be based on different trial conditions, for
%   example as set a FieldTrip .trialinfo fielS.
% * assignment of chunks (parts of the data that are assumed to be
%   independent) is based on a trial-by-trial basis. For cross-validation,
%   the number of chunks is reduced to two to speed up the analysis.
% * the current examples do not perform baseline corrections or signal
%   normalizations, which may reduce discriminatory power.
%
% Note: running this code requires FieldTrip.

%% Settings

clear all
close all
dbstop if error
plot_on=0;
%S.use_measure = 'correlation';
S.use_measure = 'crossvalidation';
S.parti='take-one-out';
%S.parti='splithalf'; 
%S.parti='oddeven';
%S.use_classifier = [];
S.use_classifier = 'LDA'; regularization = 0.5; matlab_lda = 0; logist=0; EEGLAB_logit = 0;% For ? = 1, the weight vector is proportional to ?2 ? ?1, while it is proportional to ?ˆ?1(?2??1) in the case without shrinkage (? = 0). http://www.sciencedirect.com/science/article/pii/S1053811910009067#bb0250
%S.use_classifier = 'Bayes';
%S.use_classifier = 'SVM';
%EEGLAB_logit = 0;
S.use_chunks = 'balance_targets';
    S.nchunks=4;
    
time_chan_SL=1;
time_SL =0;

% averaging: MAKES CALCULATIONS MUCH SLOWER
%average_train_count = 1; % average over this many trials
%average_train_resamplings = 1; % reuse trials this many times

% set time window of interest in ms (EMPTY: whole epoch)
%S.timewin = [0 50]; 
%S.timewin = [50 100]; 
%S.timewin = [100 150]; 
%S.timewin = [150 300];  
S.timewin = [0 150]; 

% choose factor of interest
S.factcon = 'Odd';

S.runname = 'LDA_part4_timechan';


%% SPECIFY DATA

addpath('C:\Data\Matlab\fieldtrip-20170113')
addpath('C:\Data\Matlab\CoSMoMVPA\mvpa')
addpath('C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions')
ft_defaults

% path of EEGLAB .set files after preprocessing, path of SPM outputs, and
% prefix of output files
S.data_path = 'C:\Data\CORE\Preprocessed_100Hz'; 
S.outpath = 'C:\Data\CORE\cosmo';

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
S.fpref = '';
S.fmid = '';
S.fsuff = '4_merged_cleaned.set'; balance_dataset_and_partitions =1;

% load .xlsx file containing 'Participant_ID', 'Group', and covariates
S.pdatfile = 'C:\Data\CORE\SPMstats\Participant_data.xlsx';
% names of headers in the above xls file:
    S.subhead = 'Subject';
    S.grphead = 'Group';
    S.inchead = 'Include';
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];

% conditions and factors
S.conds = 1:24;
S.useconds = 1:16;
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
          
% searchlight resolution          
time_radius=2; % e.g. if 2: 2*2+1=5 time bins

%% get timelock data in CoSMoMVPA format

% reset citation list
cosmo_check_external('-tic');

% find data
files = dir(fullfile(S.data_path,[S.fpref '*' S.fmid  '*' S.fsuff]));
cd(S.outpath)

% create FT channel layout
cfglayout = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.sfp';
cfgoutput = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.lay';
cfg.layout = cfglayout;
cfg.output = cfgoutput;
layout = ft_prepare_layout(cfg);

rundir=fullfile(S.outpath,[S.runname '_' num2str(S.timewin(1)) '_' num2str(S.timewin(2))]);
if ~exist(rundir,'dir')
    mkdir(rundir)
end

for f = 1:length(files)
    % load data
    data_fn=fullfile(S.data_path,files(f).name);
    EEG=pop_loadset(data_fn);
    sname = strsplit(EEG.filename,'_');
    
    if isempty(S.timewin)
        S.timewin = [EEG.xmin EEG.xmax]*1000;
    end
    
    % obtain trial indices
    [conds, tnums, fnums, bnums] = get_markers(EEG);
    
    % convert to cosmomvpa struct
    %ds_tl=cosmo_meeg_dataset(EEG,'data_field','data');
    ds_tl = eeglab2cosmo(EEG,S.timewin); % CAB custom function
    % create flipped version
    if S.flipside
        EEGflip = flipchan(EEG);
        ds_tl_flip = eeglab2cosmo(EEGflip,S.timewin); % CAB custom function
        sideind = find(strcmp(S.factors,'Side'));
        flipind = find(S.factor_matrix(:,sideind)==S.flipside);
        flipcond = S.conds(flipind);
        for i = flipcond
            trialind = find(ds_tl.sa.trialinfo(:,1)==i);
            ds_tl.samples(trialind,:)= ds_tl_flip.samples(trialind,:);
        end
        clear EEGflip ds_tl_flip
    end
    
    % assign full dataset
    ds_tl_full=ds_tl;

    % reduce conditions to factors of interest
    factind = find(strcmp(S.factors,S.factcon));
    condana = S.factor_matrix(S.useconds,factind);
    targ = unique(condana)';
    ti = nan(size(ds_tl.sa.trialinfo));
    for i = targ
        factcon{i} = find(condana==i);
        for c = 1:length(factcon{i})
            ti(find(ds_tl.sa.trialinfo(:,1)==factcon{i}(c)),1)=i;
        end
        condind = find(ti(:,1)==i);
        ti(condind,2) = 1:length(condind);
    end
    filter_cond=~any(isnan(ti),2);
    ti=ti(filter_cond,:);
    ds_tl.samples = ds_tl.samples(filter_cond,:);
    ds_tl.sa.trialinfo=ds_tl.sa.trialinfo(filter_cond,:);
    %ds_tl.sa.labels=ds_tl.sa.labels(filter_cond);
    %ds_tl.sa.trialinfo=ti;
    
    % initially, set the targets to be all conditions that need balancing across
    % chunks (later, targets will be simplified to the factor of interest only, ti)
    ds_tl.sa.targets=ds_tl.sa.trialinfo(:,1);

    % set the chunks (independent measurements)
    ds_tl.sa.chunks=[1:size(ds_tl.samples,1)]';
    
    if balance_dataset_and_partitions
        ds_tl_orig = ds_tl;
        ti_orig = ti;
        [ds_tl,idxs,classes]=cosmo_balance_dataset(ds_tl);
        ti = nan(size(ds_tl.sa.trialinfo));
        for i = targ
            factcon{i} = find(condana==i);
            for c = 1:length(factcon{i})
                ti(find(ds_tl.sa.trialinfo(:,1)==factcon{i}(c)),1)=i;
            end
            condind = find(ti(:,1)==i);
            ti(condind,2) = 1:length(condind);
        end
    end

    % get rid of features with at least one NaN value across samples
    fa_nan_mask=sum(isnan(ds_tl.samples),1)>0;
    fprintf('%d / %d features have NaN\n', ...
                sum(fa_nan_mask), numel(fa_nan_mask));
    ds_tl=cosmo_slice(ds_tl, ~fa_nan_mask, 2);


    %% set MVPA parameters
    fprintf('The input has feature dimensions %s\n', ...
                    cosmo_strjoin(ds_tl.a.fdim.labels,', '));

    if strcmp(S.use_chunks,'balance_targets') 
        ds_tl.sa.chunks=cosmo_chunkize(ds_tl,S.nchunks);
    end

    % set targets to be the factor of interest only
    ds_tl.sa.targets=ti(:,1);
    % in addition give a label to each trial
    ds_tl.sa.labels=cellfun(@(x)S.factlev{factind}(x),num2cell(ds_tl.sa.targets));

    % just to check everything is ok
    cosmo_check_dataset(ds_tl);

    % Set partition scheme. odd_even is fast; for publication-quality analysis
    % nfold_partitioner is recommended.
    % Alternatives are:
    % - cosmo_nfold_partitioner    (take-one-chunk-out crossvalidation)
    % - cosmo_nchoosek_partitioner (take-K-chunks-out  "             ").
    if strcmp(S.parti,'take-one-out')
        % do a take-one-fold out cross validation.
        % except when using a splithalf correlation measure it is important that
        % the partitions are *balanced*, i.e. each target (or class) is presented
        % equally often in each chunk
        partitions=cosmo_nchoosek_partitioner(ds_tl,1);
        if balance_dataset_and_partitions
            partitions=cosmo_balance_partitions(partitions, ds_tl);
        end
    elseif strcmp(S.parti,'splithalf')
        % split-half, if there are just two chunks
        % (when using a classifier, do not use 'half' but the number of chunks to
        % leave out for testing, e.g. 1).
        partitions= cosmo_nchoosek_partitioner(ds_tl,'half');
        partitions=cosmo_balance_partitions(partitions, ds_tl);
    elseif strcmp(S.parti,'oddeven')
        partitions = cosmo_oddeven_partitioner(ds_tl);
        partitions=cosmo_balance_partitions(partitions, ds_tl);
    end

    npartitions=numel(partitions);
    fprintf('There are %d partitions\n', numel(partitions.train_indices));
    fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            partitions.train_indices)));
    fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            partitions.test_indices)));

    %% define neighborhood parameters for each dimension
    
    % define the neighborhood for time
    time_nbrhood=cosmo_interval_neighborhood(ds_tl,'time',...
                                                'radius',time_radius);
                                            
    if time_SL
        nbrhood=time_nbrhood;
        nbrhood_nfeatures=cellfun(@numel,nbrhood.neighbors);
        center_ids=find(nbrhood_nfeatures>0);
    elseif time_chan_SL
 
        % define the neighborhood for channels
        cfg.senstype ='EEG';
        cfg.method = 'triangulation';
        cfg.elecfile=cfglayout;
        cfg.layout=cfgoutput;
        ft_nbrs = ft_prepare_neighbours(cfg);
        chan_nbrhood=cosmo_meeg_chan_neighborhood(ds_tl, ft_nbrs);
        
        % cross neighborhoods for chan-time searchlight
        nbrhood=cosmo_cross_neighborhood(ds_tl,{chan_nbrhood,...
                                            time_nbrhood});
        % print some info
        nbrhood_nfeatures=cellfun(@numel,nbrhood.neighbors);
        fprintf('Features have on average %.1f +/- %.1f neighbors\n', ...
                    mean(nbrhood_nfeatures), std(nbrhood_nfeatures));

        % only keep features with at least 10 neighbors
        center_ids=find(nbrhood_nfeatures>10);
    end

    measure_args=struct();
    if strcmp(S.use_measure,'crossvalidation')
        %% measure: LDA
        % Use the cosmo_cross_validation_measure and set its parameters
        % (classifier and partitions) in a measure_args struct.
        measure = @cosmo_crossvalidation_measure_CAB;
        h0_mean = 1/length(targ);
        try
            measure_args.average_train_count = average_train_count; 
            measure_args.average_train_resamplings =average_train_resamplings;
        end
        measure_args.priors = []; % leave blank to calculate from data
        if ~balance_dataset_and_partitions
            measure_args.check_partitions=0;
            %average_train_count = 1; % average over this many trials
            %average_train_resamplings = 1; % reuse trials this many times
        end
        %measure_args.priors = [sum(ti(:,1)==1),sum(ti(:,1)==2)];
    elseif strcmp(S.use_measure,'correlation')

        %% measure: correlation
        % for illustration purposes use the split-half measure because it is
        % relatively fast - but clasifiers can also be used
        measure=@cosmo_correlation_measure;
        measure_args.corr_type='Spearman';
        h0_mean = 0;
    end

    %if |isfield
    if strcmp(S.use_classifier,'LDA')
        % Define which classifier to use, using a function handle.
        % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes}
        measure_args.matlab_lda = matlab_lda;
        measure_args.logist = logist;
        measure_args.classifier = @cosmo_classify_lda_CAB;
        measure_args.regularization =regularization;
    elseif strcmp(S.use_classifier,'Bayes')
        measure_args.classifier=@cosmo_classify_naive_bayes;
    elseif strcmp(S.use_classifier,'SVM')
        measure_args.classifier=@cosmo_classify_svm;
    end

    measure_args.partitions=partitions;

    % print measure and arguments
    fprintf('Searchlight measure:\n');
    cosmo_disp(measure);
    fprintf('Searchlight measure arguments:\n');
    cosmo_disp(measure_args);

    %% average?
        % add the options to average samples to the measure arguments.
        % (if no averaging is desired, this step can be left out.)
       % average_train_args=average_train_args_cell{j};
       % measure_args=cosmo_structjoin(measure_args, average_train_args);

    %% run searchlight and save
    sl_tl_ds=cosmo_searchlight_CAB(ds_tl,nbrhood,measure,measure_args,...
                                          'center_ids',center_ids);


    %% visualize timeseries results
    % deduce layout from output
    %layout=cosmo_meeg_find_layout(sl_tl_ds);
    %fprintf('The output uses layout %s\n', layout.name);

    % map to FT struct for visualization 
    if isfield(sl_tl_ds,'weights') 
        weights=sl_tl_ds.weights; 
        offsets=sl_tl_ds.offsets; 
        prior=sl_tl_ds.prior; 
        sl_tl_ds = rmfield(sl_tl_ds,'weights'); 
        sl_tl_ds = rmfield(sl_tl_ds,'offsets'); 
        sl_tl_ds = rmfield(sl_tl_ds,'prior');
    end
    
    if time_SL
        time_values=sl_tl_ds.a.fdim.values{1}; % first dim (channels got nuked)
        plot(time_values,sl_tl_ds.samples);

        ylim([.4 .8])
        xlim([min(time_values),max(time_values)]);
        ylabel('classification accuracy (chance=.5)');
        xlabel('time');
    elseif time_chan_SL
    
        sl_tl_ft=cosmo_map2meeg(sl_tl_ds);

        % correct so that zero = chance accuracy
        %sl_tl_ft.avg = sl_tl_ft.avg-1/length(targ);

        if plot_on
            figure();
            cfg = [];
            cfg.interactive = 'yes';
            cfg.zlim=[-1 1];
            cfg.ylim ='maxmin';
            cfg.baseline='no';
            cfg.layout       = layout;

            % show figure with plots for each sensor
            ft_multiplotER(cfg, sl_tl_ft);

            %% visualize topology results
            % show figure with topology for 0 to 600ms after stimulus onset in bins of
            % 100 ms
            figure();
            cfg.xlim=0:0.005:0.1;
            cfg.zlim='maxmin';
            cfg.commentpos ='middlebottom';
            ft_topoplotER(cfg, sl_tl_ft);
        end

        %% for each feature and class, obtain all weights for that feature
        % obtained from analysis of neighbouring features
        if exist('weights','var') 
            weicell = cell(length(targ),size(ds_tl.samples,2));
            nh_centers=nbrhood.neighbors(center_ids);
            for w = 1:length(nh_centers) % for each feature
                for i = 1:length(nh_centers{w}) % for each neighbour of that feature
                    idx=nh_centers{w}(i); % find the feature index of that neighbour
                    if ~isempty(weights{w})
                        weicell{1,idx} = [weicell{1,idx},weights{w}(1,i)]; % add the weight for that neighbour to weicell
                        weicell{2,idx} = [weicell{2,idx},weights{w}(2,i)];
                    end
                end
                offmat(1,w) = offsets{w}(1); 
                offmat(2,w) = offsets{w}(2);
                primat(1,w) = prior{w}(1); 
                primat(2,w) = prior{w}(2);
            end
            for w = 1:length(weicell)
                weimat(:,w) = mean(vertcat(weicell{:,w}),2)/length(weicell{1,w}); % weights should get proportionally smaller with the number of features
            end
            offmat = mean(offmat,2);
            primat = mean(primat,2);
            weimat(find(all(isnan(weimat),2)),:) =[]; % remove NaN rows
            offmat(find(all(isnan(offmat),2)),:) =[]; % remove NaN rows
            primat(find(all(isnan(primat),2)),:) =[]; % remove NaN rows
            
            for t = 1:size(weimat,1)
                sl_wei{t}=sl_tl_ds;
                sl_wei{t}.samples=weimat(t,:);
            end
            if t==2
                sl_wei{t+1}=sl_wei{t};
                sl_wei{t+1}.samples = weimat(1,:)-weimat(2,:);
            end
            if plot_on
                for t = 1:size(sl_wei,2)
                    sl_wei_ft{t}=cosmo_map2meeg(sl_wei{t});
                    figure();
                    cfg = [];
                    cfg.interactive = 'yes';
                    cfg.zlim=[-1 1];
                    cfg.ylim ='maxmin';
                    cfg.baseline='no';
                    cfg.layout       = layout;

                    % show figure with plots for each sensor
                    ft_multiplotER(cfg, sl_wei_ft{t});

                    % visualize topology results
                    figure();
                    cfg.xlim=0:0.005:0.1;
                    cfg.zlim='maxmin';
                    cfg.commentpos ='middlebottom';
                    ft_topoplotER(cfg, sl_wei_ft{t});
                end
            end

        end
        %% 1 sample t-test
     %   opt=struct();
     %   opt.h0_mean=h0_mean;
     %   sl_tl_ds.sa.targets =1;
     %   sl_tl_ds.sa.chunks =1;
        % set the number of iterations.
        % At least 10000 is adviced for publication-quality analyses; because that
        % takes quite a while to compute, here we use 200
     %   opt.niter=200;
     %   cl_nh=cosmo_cluster_neighborhood(sl_tl_ds);
     %   tfce_sl_tl_ds=cosmo_montecarlo_cluster_stat(...
     %                                       sl_tl_ds,cl_nh,opt);
        save(fullfile(rundir,[sname{1} '_' sname{2} '_SLdata.mat']),'sl_tl_ds','weimat','offmat','primat','nbrhood');
    end
    
    %% create thresholded weight map
    if plot_on
        thresh = 0;
        wmt = weimat(1,:)-weimat(2,:);
        wmt(sl_tl_ds.samples(1,:)<thresh)=0;
        sl_wmt=sl_tl_ds;
        sl_wmt.samples=wmt;
        sl_wmt_ft=cosmo_map2meeg(sl_wmt);
        figure();
        cfg = [];
        cfg.interactive = 'yes';
        cfg.zlim=[-1 1];
        cfg.ylim ='maxmin';
        cfg.baseline='no';
        cfg.layout       = layout;

        % show figure with plots for each sensor
        ft_multiplotER(cfg, sl_wmt_ft);

        % visualize topology results
        figure();
        cfg.xlim=0:0.005:0.1;
        cfg.zlim='maxmin';%[-0.0005 0.0005];%
        cfg.commentpos ='middlebottom';
        ft_topoplotER(cfg, sl_wmt_ft);
    end
    
    % create mismatch estimate over trials of original data
    if ~exist('ds_tl_orig','var')
        ds_tl_orig=ds_tl;
        ti_orig=ti;
    end
    
    % calculate projection on all trials from factor of interest
    class_proj=bsxfun(@plus,offmat+primat,weimat*ds_tl_orig.samples');
    mm = class_proj(1,:)-class_proj(2,:);
    class1=targ(1)*single(mm>0);%CAB 
    class2=targ(2)*single(mm<=0);%CAB 
    predicted=class1+class2;%CAB 
    correct=bsxfun(@eq,ti_orig(:,1)',predicted);
    accuracy=[sum(correct)/length(ti_orig(:,1)) sum(ti_orig(:,1)==1)/length(ti_orig(:,1))]
    ratio = accuracy(1)/accuracy(2)
    
    % generate t stats
    [h,p,ci,stats] = ttest2(mm(ti_orig(:,1)==1),mm(ti_orig(:,1)==2))
    if plot_on
        figure
        plot(mm)
        figure
        scatter(ti_orig(:,1),mm');
    end
    
    % calculate projection on all trials including previously excluded
    % factors
    class_proj=bsxfun(@plus,offmat+primat,weimat*ds_tl_full.samples');
    mm = class_proj(1,:)-class_proj(2,:);
    
    % save
    
    save(fullfile(rundir,[sname{1} '_' sname{2} '_mm_projection.mat']),'mm','class_proj','ti_orig','tnums','predicted','targ');
    
    %
    if EEGLAB_logit
        fc = find(filter_cond);
        i1 = fc(ds_tl.sa.targets==1);
        i2 = fc(ds_tl.sa.targets==2);
        ALLEEG(1) = pop_select(EEG,'trial',i1);
        ALLEEG(2) = pop_select(EEG,'trial',i2);

        datasetlist = [1 2];
        chansubset=1:min(ALLEEG(datasetlist(1)).nbchan,ALLEEG(datasetlist(2)).nbchan);
        chansubset2=1:min(ALLEEG(datasetlist(1)).nbchan,ALLEEG(datasetlist(2)).nbchan);
        trainingwindowlength = 5;
        vinit = [];
        show=0;
        regularize =regularization>0;
        lambda = regularization;
        lambdasearch=1;
        eigvalratio = 0.01;
        dp = round(((S.timewin/1000)-EEG.xmin)*EEG.srate);
        trainingwindowoffset = dp(1):trainingwindowlength:dp(2)-(trainingwindowlength);
        LOO = 1;
        [ALLEEG, com] = pop_logisticregression(ALLEEG,datasetlist,chansubset,chansubset2,trainingwindowlength,trainingwindowoffset,regularize,lambda,lambdasearch,eigvalratio,vinit,show,LOO);
        figure;topoplot(ALLEEG(datasetlist(1)).icaweights(9,chansubset), EEG.chanlocs);
    end

end