function out = run_cosmo_machine(cos,S)
% multivariate classification between two groups
% using linear discriminant analysis

cos.sa.samples=round(cos.sa.samples,S.ndec);
rm=all(isnan(dcos.sa.samples),1) | all(diff(cos.sa.samples)==0); % remove nans and constants
cos.sa.samples(:,rm)=[];

if isempty(cos.sa.samples)
    out.empty=1;
    return
end

% balance trials
[cos,idxs,classes]=cosmo_balance_dataset(cos);

% get rid of features with at least one NaN value across samples
fa_nan_mask=sum(isnan(cos.samples),1)>0;
fprintf('%d / %d features have NaN\n', ...
            sum(fa_nan_mask), numel(fa_nan_mask));
cos=cosmo_slice(cos, ~fa_nan_mask, 2);

% balance targets over chunks
if strcmp(S.use_chunks,'balance_targets') 
    cos.sa.chunks=cosmo_chunkize(cos,S.nchunks);
end

% just to check everything is ok before analysis
cosmo_check_dataset(cos);
fprintf('The input has feature dimensions %s\n', ...
        cosmo_strjoin(cos.a.fdim.labels,', '));
    
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
    partitions=cosmo_nchoosek_partitioner(cos,1);
    if S.balance_dataset_and_partitions
        partitions=cosmo_balance_partitions(partitions, cos);
    end
elseif strcmp(S.parti,'splithalf')
    % split-half, if there are just two chunks
    % (when using a classifier, do not use 'half' but the number of chunks to
    % leave out for testing, e.g. 1).
    partitions= cosmo_nchoosek_partitioner(cos,'half');
    partitions=cosmo_balance_partitions(partitions, cos);
elseif strcmp(S.parti,'oddeven')
    partitions = cosmo_oddeven_partitioner(cos);
    partitions=cosmo_balance_partitions(partitions, cos);
elseif strcmp(S.parti,'nfold')
    partitions=cosmo_nchoosek_partitioner(cos,S.nchunks);
    if S.balance_dataset_and_partitions
        partitions=cosmo_balance_partitions(partitions, cos);
    end
end

npartitions=numel(partitions);
fprintf('There are %d partitions\n', numel(partitions.train_indices));
fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                        partitions.train_indices)));
fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                        partitions.test_indices)));
                                    
measure_args=struct();
if strcmp(S.use_measure,'crossvalidation')
    %% measure: LDA
    % Use the cosmo_cross_validation_measure and set its parameters
    % (classifier and partitions) in a measure_args struct.
    measure = @cosmo_crossvalidation_measure_CAB;
    %try
        measure_args.average_train_count = S.average_train_count; 
        measure_args.average_train_resamplings =S.average_train_resamplings;
    %end
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

% define neighborhood
time_nbrhood=cosmo_interval_neighborhood(cos,'time',...
                                        'radius',S.search_radius);
if strcmp(S.SL_type,'time')
    nbrhood=time_nbrhood;
    nbrhood_nfeatures=cellfun(@numel,nbrhood.neighbors);
    center_ids=find(nbrhood_nfeatures>0);
elseif strcmp(S.SL_type,'time_chan')

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


% run the searchlight using the measure, measure arguments, and
% neighborhood defined above.
% Note that while the input has both 'chan' and 'time' as feature
% dimensions, the output only has 'time' as the feature dimension

out=cosmo_searchlight_CAB(cos,nbrhood,measure,measure_args,...
                                          'center_ids',center_ids);
                                 
                          