%% MEEG time-lock searchlight
%
% Using a time-channel neighborhood, a searchlight map is computed
% indicating in time and space (channel) the pre-stimulation period can be
% distinguished from the peri/post stimulation period.
%
% The code presented here can be adapted for other MEEG analyses, but
% there are a few potential caveats:
% * assignment of targets (labels of conditions) is based here on
%   stimulation periods versus pre-stimulation periods. In typical
%   analyses the targets should be based on different trial conditions, for
%   example as set a FieldTrip .trialinfo field.
% * assignment of chunks (parts of the data that are assumed to be
%   independent) is based on a trial-by-trial basis. For cross-validation,
%   the number of chunks is reduced to two to speed up the analysis.
% * the current examples do not perform baseline corrections or signal
%   normalizations, which may reduce discriminatory power.
%
% Note: running this code requires FieldTrip.

%% Settings
clear all
%use_measure = 'correlation';
use_measure = 'crossvalidation';
    crossval='take-one-out';
    %crossval='splithalf'; % only if there are two chunks set for nchunks
    %crossval='oddeven';
use_classifier = 'LDA'; % options: 'LDA', 'Bayes'
use_chunks = 'balance_targets';
    nchunks=4;
    
%% file info
% set configuration
config=cosmo_config();
data_path=fullfile(config.tutorial_data_path,'meg_20hz');

%% get timelock data in CoSMoMVPA format


% show dataset information
readme_fn=fullfile(data_path,'README');
cosmo_type(readme_fn);

% reset citation list
cosmo_check_external('-tic');

% load data
data_fn=fullfile(data_path,'subj102_B01_20Hz_timelock.mat');
data_tl=load(data_fn);

% convert to cosmomvpa struct
ds_tl=cosmo_meeg_dataset(data_tl);

% set the target (trial condition)
ds_tl.sa.targets=ds_tl.sa.trialinfo(:,1); % 1=pre, 2=post
targ = unique(ds_tl.sa.targets);

% set the chunks (independent measurements)
% in this dataset, the first half of the samples (in order)
% are the post-trials;
% the second half the pre-trials
ds_tl.sa.chunks=[(1:145) (1:145)]';

% in addition give a label to each trial
index2label={'pre','post'}; % 1=pre, 2=peri/post
ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));

% just to check everything is ok
cosmo_check_dataset(ds_tl);

% get rid of features with at least one NaN value across samples
fa_nan_mask=sum(isnan(ds_tl.samples),1)>0;
fprintf('%d / %d features have NaN\n', ...
            sum(fa_nan_mask), numel(fa_nan_mask));
ds_tl=cosmo_slice(ds_tl, ~fa_nan_mask, 2);


%% set MVPA parameters
fprintf('The input has feature dimensions %s\n', ...
                cosmo_strjoin(ds_tl.a.fdim.labels,', '));
            
if strcmp(use_chunks,'balance_targets') 
    ds_tl.sa.chunks=cosmo_chunkize(ds_tl,nchunks);
end

% Set partition scheme. odd_even is fast; for publication-quality analysis
% nfold_partitioner is recommended.
% Alternatives are:
% - cosmo_nfold_partitioner    (take-one-chunk-out crossvalidation)
% - cosmo_nchoosek_partitioner (take-K-chunks-out  "             ").
if strcmp(crossval,'take-one-out')
    % do a take-one-fold out cross validation.
    % except when using a splithalf correlation measure it is important that
    % the partitions are *balanced*, i.e. each target (or class) is presented
    % equally often in each chunk
    partitions=cosmo_nchoosek_partitioner(ds_tl,1);
    partitions=cosmo_balance_partitions(partitions, ds_tl);
elseif strcmp(crossval,'splithalf')
    % split-half, if there are just two chunks
    % (when using a classifier, do not use 'half' but the number of chunks to
    % leave out for testing, e.g. 1).
    partitions= cosmo_nchoosek_partitioner(ds_tl,'half');
elseif strcmp(crossval,'oddeven')
    partitions = cosmo_oddeven_partitioner(ds_tl);
end

npartitions=numel(partitions);
fprintf('There are %d partitions\n', numel(partitions.train_indices));
fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                        partitions.train_indices)));
fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                        partitions.test_indices)));

%% define neighborhood parameters for each dimension

% channel neighborhood uses meg_combined_from_planar, which means that the
% input are planar channels but the output has combined-planar channels.
% to use the magnetometers, use 'meg_axial'
chan_type='meg_combined_from_planar';
chan_count=10;        % use 10 channel locations (relative to the combined
                      % planar channels)
                      % as we use meg_combined_from_planar there are
                      % 20 channels in each searchlight because
                      % gradiometers are paired
time_radius=2; % 2*2+1=5 time bines

% define the neighborhood for each dimensions
chan_nbrhood=cosmo_meeg_chan_neighborhood(ds_tl, 'count', chan_count, ...
                                                'chantype', chan_type);
time_nbrhood=cosmo_interval_neighborhood(ds_tl,'time',...
                                            'radius',time_radius);

% cross neighborhoods for chan-time searchlight
nbrhood=cosmo_cross_neighborhood(ds_tl,{chan_nbrhood,...
                                        time_nbrhood});

% print some info
nbrhood_nfeatures=cellfun(@numel,nbrhood.neighbors);
fprintf('Features have on average %.1f +/- %.1f neighbors\n', ...
            mean(nbrhood_nfeatures), std(nbrhood_nfeatures));

% only keep features with at least 10 neighbors
center_ids=find(nbrhood_nfeatures>10);

measure_args=struct();
if strcmp(use_measure,'crossvalidation')
    %% measure: LDA
    % Use the cosmo_cross_validation_measure and set its parameters
    % (classifier and partitions) in a measure_args struct.
    measure = @cosmo_crossvalidation_measure_CAB;
    
elseif strcmp(use_measure,'correlation')

    %% measure: correlation
    % for illustration purposes use the split-half measure because it is
    % relatively fast - but clasifiers can also be used
    measure=@cosmo_correlation_measure;
    measure_args.corr_type='Spearman';
end

if strcmp(use_classifier,'LDA')
    % Define which classifier to use, using a function handle.
    % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes}
    measure_args.classifier = @cosmo_classify_lda_CAB;
elseif strcmp(use_classifier,'Bayes')
    measure_args.classifier=@cosmo_classify_naive_bayes;
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
if isfield(sl_tl_ds,'weights') % CAB
    weights=sl_tl_ds.weights; % CAB
    sl_tl_ds = rmfield(sl_tl_ds,'weights'); % CAB
    weicell = cell(length(targ),size(ds_tl.samples,2));
    for w = 1:length(nbrhood.neighbors)
        for i = 1:length(nbrhood.neighbors{w})
            idx=nbrhood.neighbors{w}(i);
            weicell{1,idx} = [weicell{1,idx},weights{w}(1,i)];
            weicell{2,idx} = [weicell{2,idx},weights{w}(2,i)];
        end
    end
end

save(fullfile(data_path,['testrun_' datestr(now,30)]),'sl_tl_ds','weights','weicell','nbrhood');

%% visualize timeseries results
% deduce layout from output
layout=cosmo_meeg_find_layout(sl_tl_ds);
fprintf('The output uses layout %s\n', layout.name);

% map to FT struct for visualization
sl_tl_ft=cosmo_map2meeg(sl_tl_ds);

figure();
cfg = [];
cfg.interactive = 'yes';
cfg.zlim=[-1 1];
cfg.layout       = layout;

% show figure with plots for each sensor
ft_multiplotER(cfg, sl_tl_ft);

%% visualize topology results
% show figure with topology for 0 to 600ms after stimulus onset in bins of
% 100 ms
figure();
cfg.xlim=-0.1:0.1:0.5;
cfg.commentpos ='middlebottom';
ft_topoplotER(cfg, sl_tl_ft);

%% Show citation information
cosmo_check_external('-cite');

