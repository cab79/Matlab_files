function [S,D,stats] = CORE_eeg_trial_statistics(S,D)

%% set defaults
% if separate event groups or contrasts are not specified, pool all
% events
if ~isfield(S,'contrast_rows')
    S.contrast_rows = {};
end

% Mean subtraction setting - off by default
if ~isfield(S,'row_subtract')
    S.row_subtract = {};
end

%% find data
% get file names
% get filename using only first input type (filename suffix)
S.load.suffix = S.load.suffixes(1);

% GET FILE LIST
S.path.file = S.path.eeg;
S = getfilelist(S);

% change to the input directory
eval(sprintf('%s', ['cd(''' S.path.file ''')']));

% report if there are no such files
if isempty(S.filelist)
    error('No files found to import!\n');
end

if S.dsample
    S.total_samples = downsample(S.total_samples',S.dsample)';
    S.select_samples = downsample(S.select_samples',S.dsample)';
end

%% run though all files in a loop and get data
for d = 1:length(S.select.subjects)
    
    % create output D
    D.subname = S.select.subjects{d};
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,D.subname)))));

    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
    % add data to D structure which has a common format for any data type
    for f = 1:length(subfiles)
        filename = subfiles{f};
        S.load.suffix = S.load.suffixes{1};
        fprintf('\nImporting %s.\n\n', filename);

        switch S.fname.ext{:}
            case 'mat'
                temp = load(filename);
            case 'set'
                temp = pop_loadset(filename);
        end
        fename = fieldnames(temp(f));
        
        for fn = 1:length(fename)
            D.(fename{fn})(f).dat = temp.(fename{fn});
        end

        % get other filetypes with same name
        for fn = 1:length(S.load.suffixes)-1
            S.load.suffix = S.load.suffixes{fn+1};
            filename = strrep(filename,S.load.suffixes{1},S.load.suffix);

            temp = load(filename);
            fename = fieldnames(temp(f));
            for fn = 1:length(fename)
                if isstruct(temp.(fename{fn}))
                    D.(fename{fn})(f) = temp.(fename{fn});
                else
                    D.(fename{fn})(f).dat = temp.(fename{fn});
                end
            end
        end
    end

    switch S.fname.ext{:}
        case 'mat'
            topography = D.topography.dat;
            timecourse = D.timecourse.dat;
            eventType = D.eventType.dat;
            n_chans = D.n_chans.dat;
            n_samples = D.n_samples.dat;
            n_trials = D.n_trials.dat;
            tnums=1:length(eventType); % assumes all trials present
        case 'set'
            [eventType,tnums] = get_markers(temp);
            n_chans = D.nbchan.dat;
            n_samples = D.pnts.dat;
            n_trials = D.trials.dat;
            timecourse = reshape(D.data.dat,n_chans,[]);
    end
    

    switch S.data_type
        case {'chan','comp','comp_recon'} % process each chan or IC separately
            iter = 1:size(timecourse,1);
        case 'recon' % combine all ICs into a single reconstruction of the data
            iter = 1;
        case {'multicomp','all_chan'} % consider chans/ICs as multivariate data rather than separately
            iter = 1;
    end
    
    % load and format HGF data
    if isfield(S.path,'hgf')
        load(S.path.hgf);
        ind = find(strcmp({D_fit(:).subname},D.subname));
        D.HGF.dat = D_fit(ind).HGF;
        D.pred=[];
        D.pred_group=[];
        D.pred_label={};
        if strcmp(S.pred,'HGF')
            % for each traj group... reformat into predictor matrix
            for tg = 1:length(S.traj)
                for pm = 1:size(S.traj{tg})
                    pm_name = S.traj{tg}{pm,1};
                    for tr = 1:length(S.traj{tg}{pm,2})
                        tr_name = S.traj{tg}{pm,2}(tr);
                        tr_dat = D.HGF.dat.fit.traj.(pm_name{:}).(tr_name{:});
                        D.pred = [D.pred, tr_dat];
                        nvar=size(tr_dat,2);
                        D.pred_group = [D.pred_group, tg*ones(1,nvar)];
                        for nv = 1:nvar
                            D.pred_label = [D.pred_label, {[pm_name{:} '_' tr_name{:} num2str(nv)]}];
                        end
                    end
                end
            end
        elseif strcmp(S.pred,'RT')
            D.pred = D.HGF.dat.y(:,2);
            D.pred_group=1;
            D.pred_label={S.pred};
        end
        
        % remove trials not in EEG
        D.pred=D.pred(tnums,:);
        
        % remove all-NaN predictors
        ip = all(isnan(D.pred),1);
        D.pred(:,ip) = [];
        D.pred_group(ip)=[];
        D.pred_label(ip)=[];
        
    end
    
    % for each component/channel
    for c=iter
        
        disp(['file ' num2str(d) '/' num2str(length(S.select.subjects)) ' , comp/chan ' num2str(c) ' /' num2str(iter(end))])      
        
        switch S.data_type
            case {'chan','comp'}
                comps = c;
                data= timecourse(comps,:);  
            case 'comp_recon'
                comps = c;
                data= topography(:,comps) * timecourse(comps,:); 
            case 'recon'
                comps = 1:size(timecourse,1);
                data= topography(:,comps) * timecourse(comps,:);  
            case {'multicomp','all_chan'}
                data= timecourse;  
        end
        
        % smooth
        if S.smooth_samples
            disp('smoothing...')
            for i = 1:size(data,1)
                data(i,:) = smooth(data(i,:),S.smooth_samples,'moving');
            end
        end
        
        % downsample
        if S.dsample
            data = downsample(data',S.dsample)';
        end

        % make 3D
        data=reshape(data,size(data,1),[],n_trials);
            
        % select samples
        select = dsearchn(S.total_samples',[S.select_samples(1),S.select_samples(end)]');
        data = data(:,select(1):select(2),:);
        
        % create sets of trials 
        for set = 1:size(S.cond_idx,1)
            setidx{set} = find(ismember(eventType,[S.cond_idx{set}]));
            setData{set} = data(:,:,setidx{set});
        end
        
        % perform subtraction
        if ~isempty(S.row_subtract)
            for sb = 1:length(S.row_subtract{1})
                meandat = mean(setData{S.row_subtract{2}(sb)},3);
                subData{sb} = bsxfun(@minus,setData{S.row_subtract{1}(sb)},meandat);
                subidx{sb}=setidx{S.row_subtract{1}(sb)};
            end
            setData=subData;
            setidx = subidx;
        end
            
        % pool into contrasts
        conData={};
        idx={};
        if ~isempty(S.contrast_rows)
            for con = 1:length(S.contrast_rows)
                conData{con} = cat(3,setData{S.contrast_rows{con}});
                idx{con} = [setidx{S.contrast_rows{con}}];
            end
        else
            conData = {cat(3,setData{:})};
            idx = {[setidx{:}]};
        end
        
        % Data transformation
        if ~isempty(S.transform)
            for con = 1:length(conData)
                x=conData{con};
                if strcmp(S.transform,'arcsinh')
                    conData{con}=log(x+sqrt(x.^2+1));
                end
            end
        end
        
        % duplicate data over a number of runs
        if S.num_runs>1 && length(conData)==1
            for con = 1:S.num_runs
                conData{con} = conData{1};
                idx{con} = idx{1};
            end
        end
      
        for con = 1:length(conData)
            % split into training (encoding) and testing (decoding) fractions
            if S.balance_conds
                ucond = unique(eventType(idx{con}));
                S.trainidx{con} = [];
                for u = 1:length(ucond)
                    cond_idx = find(eventType(idx{con})==ucond(u));
                    rng(con); % seed random number generator for consistency
                    S.trainidx{con} = [S.trainidx{con} randsample(cond_idx,round(S.trainfrac*length(cond_idx)))];
                end
                S.trainidx{con} = sort(S.trainidx{con});
            else
                rng(1); % seed random number generator for consistency
                S.trainidx{con} = sort(randsample(length(idx{con}),round(S.trainfrac*length(idx{con}))));
            end
            if S.trainfrac<1
                S.testidx{con} = 1:length(idx{con});
                S.testidx{con}(S.trainidx{con}) = [];
            else
                S.testidx{con} = S.trainidx{con};
            end
        end
        
        % create con-specific predictors and update trainidx to remove NaN
        % predictor value indices
        % zscore normalisation for regression
        if c==1
            pred=D.pred;
        end
        D.pred={};
        for con = 1:length(conData)
            % idx is an index of eventType
            D.pred{con}=pred(idx{con},:);
            % remove NaNs
            D.pred_nonan{con}=find(~any(isnan(D.pred{con}),2));
            %D.pred{con}=D.pred{con}(D.pred_nonan{con},:);
            S.trainidx{con} = S.trainidx{con}(ismember(S.trainidx{con},D.pred_nonan{con}));
            S.testidx{con} = S.testidx{con}(ismember(S.testidx{con},D.pred_nonan{con}));
            
            % store info so we know which trials and how many were analysed
            stat.trialinfo.trainidx{con}=S.trainidx{con};
            stat.trialinfo.testidx{con}=S.testidx{con};
            stat.trialinfo.ntrials_traintest(d,c,:)= [length(S.trainidx{con}),length(S.testidx{con})];
            
        end

        for con = 1:length(conData)
            % reduce to Global Field Power
            if size(conData{con},1)>1
                gfpData{con} = squeeze(std(conData{con},{},1));
            else
                gfpData{con} = squeeze(conData{con});
            end
            stdgfpData{con} = std(gfpData{con},[],1);
        end
        
        disp('running stats...')
        
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'two_sample') 
            
            % non-parametric independent samples test
            if isfield(S,'ranksum_on') && S.ranksum_on
                for s = 1:length(S.select_samples)
                    [p,~,st]=ranksum(gfpData{1}(s,:),gfpData{2}(s,:));
                    stats.ranksum.all(d,c,s) = st.ranksum;
                    stats.ranksum.pval(d,c,s) = p;
                    stats.meandiff.gfp(d,c,s)=mean(gfpData{1}(s,:))-mean(gfpData{2}(s,:));
                end
                stats.ranksum.max(d,c) = max(stats.ranksum.all(d,c,:));
                stats.ranksum.min(d,c) = min(stats.ranksum.all(d,c,:));
                stats.ranksum.min_pval(d,c) = min(stats.ranksum.pval(d,c,:));
                stats.ranksum.stdgfp_pval(d,c)=ranksum(stdgfpData{1},stdgfpData{2});
                stats.meandiff.stdgfp(d,c)=mean(stdgfpData{1})-mean(stdgfpData{2});
            end
            
            % TFCE test
            if isfield(S,'tfce_on') && S.tfce_on
                shiftd=4-ndims(gfpData{1});
                img1=permute(shiftdim(gfpData{1},-shiftd),[3 1 2 4]);
                img2=permute(shiftdim(gfpData{2},-shiftd),[3 1 2 4]);
                varargout=matlab_tfce('independent',S.tfce_tail,img1,img2,'nperm',S.tfce_nperm);
                if ~iscell(varargout); varargout = {varargout}; end
                for v = 1:length(varargout)
                    stats.tfce(d,c,v) = min(varargout{v});
                end
            end
            
            % Cosmo TFCE
            if isfield(S,'cosmo_tfce_on') && S.cosmo_tfce_on
                % create cosmo data struct
                conds = nan(1,length(eventType));
                for cn = 1:length(idx)
                    conds(idx{cn}) = cn;
                end
                cos = eeglab2cosmo(data,S.select_samples,conds);

                % set the targets and chunks
                cos.sa.targets=cos.sa.trialinfo(:,1);
                cos.sa.chunks=[1:size(cos.samples,1)]'; % each trial is a chunk

                % set options
                opt=struct();
                %opt.h0_mean=h0_mean; % 1-sample test
                opt.niter=S.tfce_nperm;

                % create neighbourhood for time only
                cos_nh=cosmo_cluster_neighborhood(cos,'time',true,'chan',false);

                % create neighbourhood for time and chan (freq is another option)
                %cos_nh=cosmo_cluster_neighborhood(cos,'time',true,'chan',true);

                % display
                fprintf('Cluster neighborhood:\n');
                cosmo_disp(cluster_nbrhood);

                % run
                out=cosmo_montecarlo_cluster_stat(cos,cos_nh,opt);
            end
        end
        
        % non-parametric Spearman's correlation
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'SC')
            for con = 1:length(conData)
                
                if any(strcmp(S.data_form,'GFP'))
                    % GFP data
                    [stats.spear.GFP(con).rho{d,c}, stats.spear.GFP(con).p{d,c}] = corr(gfpData{con}(:,S.trainidx{con})',D.pred{con}(S.trainidx{con},:));

                    % FDR correction
                    [stats.spear.GFP(con).fdr_thres{d,c},stats.spear.GFP(con).fdr_masked{d,c}] = fdr(stats.spear.GFP(con).p{d,c},0.05);

                    % TFCE test
                    if isfield(S,'tfce_on') && S.tfce_on
                        shiftd=4-ndims(gfpData{con});
                        img=permute(shiftdim(gfpData{con},-shiftd),[3 1 2 4]);
                        for pr = 1:size(D.pred{con},2)
                            varargout=matlab_tfce('correlation',S.tfce_tail,img,[],D.pred{con}(S.trainidx{con},pr),S.tfce_nperm);
                            if ~iscell(varargout); varargout = {varargout}; end
                            % for each (pos and neg):
                            for v = 1:length(varargout)
                                stats.spear.GFP(con).tfce{d,c,v}(:,pr) = varargout{v};
                            end
                        end
                    end
                end
                
                if any(strcmp(S.data_form,'alldata'))
                    % all data
                    if ndims(conData{con})==3
                        alldata = reshape(conData{con},[],size(conData{con},3));
                        [rho, p] = corr(alldata(:,S.trainidx{con})',D.pred{con}(S.trainidx{con},:));
                        % FDR correction
                        [fdr_thresh,fdr_masked] = fdr(p,0.05);
                        % reshape
                        stats.spear.alldata(con).rho{d,c}=reshape(rho,size(conData{con},1),size(conData{con},2),[]);
                        stats.spear.alldata(con).p{d,c}=reshape(p,size(conData{con},1),size(conData{con},2),[]);
                        stats.spear.alldata(con).fdr_thresh{d,c}=fdr_thresh;
                        stats.spear.alldata(con).fdr_masked{d,c}=reshape(fdr_masked,size(conData{con},1),size(conData{con},2),[]);
                    end
                end
            end
        end
        
        % Multiple regression (linear, non-robust to outliers, non-robust to collinearity)
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'MR')
            
            for con = 1:length(conData)
                
                if size(D.pred{con}(S.trainidx{con},:),2)>1
                    % check predictor matrix is full rank
                    rk = rank(D.pred{con}(S.trainidx{con},:));
                    fprintf('Design matrix rank = %i\n',rk);
                    if rk<size(D.pred{con}(S.trainidx{con},:),2)
                        error('Design matrix is not full rank; predictors are collinear. Reduce the predictors or use a method that is robust to collinearity')
                    end
                end
                
                % zscore OVER TRIALS
                if S.zscore
                    D.pred{con}(S.trainidx{con},:) = zscore(D.pred{con}(S.trainidx{con},:));
                    conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
                    gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
                end
                
                D.pred{con}=[ones(size(D.pred{con},1),1) D.pred{con}];
                
                if any(strcmp(S.data_form,'GFP'))
                    % GFP data
                    for s = 1:length(S.select_samples)
                        [beta,~,~,~,stt] = regress(gfpData{con}(s,S.trainidx{con})',D.pred{con}(S.trainidx{con},:));
                        stats.MR.GFP(con).R2{d,c}(s) = stt(1);
                        stats.MR.GFP(con).F{d,c}(s) = stt(2);
                        stats.MR.GFP(con).p{d,c}(s) = stt(3);
                        stats.MR.GFP(con).beta{d,c}(s,:) = beta;
                    end

                    % FDR correction
                    [stats.MR.GFP(con).fdr_thres{d,c},stats.MR.GFP(con).fdr_masked{d,c}] = fdr(stats.MR.GFP(con).p{d,c},0.05);

                    % TFCE test
                    if isfield(S,'tfce_on') && S.tfce_on
                        shiftd=4-ndims(gfpData{con});
                        img=permute(shiftdim(gfpData{con},-shiftd),[3 1 2 4]);
                        varargout=matlab_tfce('regression',S.tfce_tail,img,[],D.pred{con}(S.trainidx{con},:), S.tfce_nperm,[],[],[],[],[],zeros(1,size(D.pred{con},2)));
                        if ~iscell(varargout); varargout = {varargout}; end
                        % for each (pos and neg):
                        for v = 1:length(varargout)
                            stats.MR.GFP(con).tfce{d,c,v} = varargout{v};
                        end
                    end
                end
                
                if any(strcmp(S.data_form,'alldata'))
                    % all data
                    if ndims(conData{con})==3
                        alldata = reshape(conData{con},[],size(conData{con},3));
                        for s = 1:size(alldata,1)
                            [beta,~,~,~,stt] = regress(alldata(s,S.trainidx{con})',D.pred{con}(S.trainidx{con},:));
                            stats.MR.alldata(con).R2{d,c}(s) = stt(1);
                            stats.MR.alldata(con).F{d,c}(s) = stt(2);
                            stats.MR.alldata(con).p{d,c}(s) = stt(3);
                            stats.MR.alldata(con).beta{d,c}(s,:) = beta;
                        end
                        % FDR correction
                        [fdr_thresh,fdr_masked] = fdr(stats.MR.alldata(con).p{d,c},0.05);
                        % reshape
                        stats.MR.alldata(con).R2{d,c}=reshape(stats.MR.alldata(con).R2{d,c},size(conData{con},1),size(conData{con},2),[]);
                        stats.MR.alldata(con).p{d,c}=reshape(stats.MR.alldata(con).p{d,c},size(conData{con},1),size(conData{con},2),[]);
                        stats.MR.alldata(con).F{d,c}=reshape(stats.MR.alldata(con).F{d,c},size(conData{con},1),size(conData{con},2),[]);
                        stats.MR.alldata(con).beta{d,c}=reshape(stats.MR.alldata(con).beta{d,c},size(conData{con},1),size(conData{con},2),[]);
                        stats.MR.alldata(con).fdr_thresh{d,c}=fdr_thresh;
                        stats.MR.alldata(con).fdr_masked{d,c}=reshape(fdr_masked,size(conData{con},1),size(conData{con},2),[]);
                    end
                end
            end
        end
        
        % Bayesian regression (linear, non-robust to outliers, non-robust to collinearity)
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'PEB')
            
            for con = 1:length(conData)
                
                if size(D.pred{con}(S.trainidx{con},:),2)>1
                    % check predictor matrix is full rank
                    rk = rank(D.pred{con}(S.trainidx{con},:));
                    fprintf('Design matrix rank = %i\n',rk);
                    if rk<size(D.pred{con}(S.trainidx{con},:),2)
                        error('Design matrix is not full rank; predictors are collinear. Reduce the predictors or use a method that is robust to collinearity')
                    end
                end
                
                % zscore OVER TRIALS
                if S.zscore
                    D.pred{con}(S.trainidx{con},:) = zscore(D.pred{con}(S.trainidx{con},:));
                    conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
                    gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
                end
                
                %D.pred{con}=[ones(size(D.pred{con},1),1) D.pred{con}];

                %enforces positively constraints on the covariance hyperparameters
                % by adopting a log-normal [flat] hyperprior.
                lognormhyper=0;
                
                % To implement non-hierarchical Bayes with priors on the parameters use
                % a two level model setting the second level design matrix to zeros, 
                % providing an unconstrained prior on the first-level parameters
                % which allows for single-level Bayesian inference
                P{1}.X = D.pred{con}(S.trainidx{con},:); %as your [n x p] design matrix
                P{1}.C{1}=eye(size(D.pred{con}(S.trainidx{con},:),1)); %as your covariance structure on observation errors
                P{2}.X=zeros(size(D.pred{con}(S.trainidx{con},:),2),1);
                P{2}.C{2}=eye(size(D.pred{con}(S.trainidx{con},:),2)); %as your prior covariance on parameters
                
                if any(strcmp(S.data_form,'GFP'))
                    % GFP data
                    for s = 1:length(S.select_samples)
                        [C,~,F] = spm_PEB(gfpData{con}(s,S.trainidx{con})',P,lognormhyper);
                        stats.PEB.GFP(con).F{d,c}(s) = F;
                        stats.PEB.GFP(con).cE{d,c}(s,:) = full(C{2}.E);
                    end
                end
                
                if any(strcmp(S.data_form,'alldata'))

                    % all data
                    if ndims(conData{con})==3
                        alldata = reshape(conData{con},[],size(conData{con},3));
                        for s = 1:size(alldata,1)
                            [C,~,F(s)] = spm_PEB(alldata(s,S.trainidx{con})',P,lognormhyper);
                            CE(s,:)=full(C{2}.E);
                        end

                        % reshape
                        stats.PEB.alldata(con).F{d,c}=reshape(F,size(conData{con},1),size(conData{con},2),[]);
                        stats.PEB.alldata(con).cE{d,c}=reshape(CE,size(conData{con},1),size(conData{con},2),[]);
                    end
                end
            end
        end
        
        % Ridge regression
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'RR')
            
            for con = 1:length(conData)
            
                % check predictor matrix is NOT full rank
                rk = rank(D.pred{con}(S.trainidx{con}));
                fprintf('Design matrix rank = %i\n',rk);
                if rk==size(D.pred{con}(S.trainidx{con}),2)
                    disp('Design matrix is full rank; can use standard linear regression')
                    %pause(5)
                end

                % zscore OVER TRIALS - occurs within subfunction
%                 if S.zscore
%                     D.pred{con}(S.trainidx{con},:) = zscore(D.pred{con}(S.trainidx{con},:));
%                     conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
%                     gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
%                 end

                % GFP data
                if any(strcmp(S.data_form,'GFP'))
                    %run ridge regression
                    [Beta_hat, Sigma_hat, X_mean, X_stand_de, Y_mean, Y_stand_de] = ...
                        cross_validation(S.rr.df_num, S.rr.folds, length(S.trainidx{con}), D.pred{con}(S.trainidx{con},:)', gfpData{con}(:,S.trainidx{con}), S.rr.z);
                    stats.RR.GFP(con).b{d,c} = Beta_hat;
                    stats.RR.GFP(con).s{d,c} = Sigma_hat;
                end
                
                % all data
                if any(strcmp(S.data_form,'alldata'))
                    if ndims(conData{con})==3
                        alldata = reshape(conData{con},[],size(conData{con},3));
                        %run ridge regression
                        [Beta_hat, Sigma_hat, X_mean, X_stand_de, Y_mean, Y_stand_de] = ...
                            cross_validation(S.rr.df_num, S.rr.folds, length(S.trainidx{con}), D.pred{con}(S.trainidx{con},:)', alldata(:,S.trainidx{con}), S.rr.z);
                        stats.RR.alldata(con).b{d,c}=reshape(Beta_hat,size(conData{con},1),size(conData{con},2),[]);
                        stats.RR.alldata(con).s{d,c} = Sigma_hat;
                    end
                end
            end
        end
        
        % Bayesian regularised regression
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'BRR')
            
            for con = 1:length(conData)
                
                % check predictor matrix is NOT full rank
                rk = rank(D.pred{con}(S.trainidx{con}));
                fprintf('Design matrix rank = %i\n',rk);
                if rk==size(D.pred{con}(S.trainidx{con}),2)
                    disp('Design matrix is full rank; can use standard linear regression')
                    pause(5)
                end

                % zscore OVER TRIALS - occurs within subfunction
%                 if S.zscore
%                     D.pred{con}(S.trainidx{con},:) = zscore(D.pred{con}(S.trainidx{con},:));
%                     conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
%                     gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
%                 end
                
                % GFP data
                if any(strcmp(S.data_form,'GFP'))
                    for s = 1:length(S.select_samples)
                        disp(['GFP sample ' num2str(s) '/' num2str(length(S.select_samples))])
                        out = bayesreg_crossval(D.pred{con}(S.trainidx{con},:),gfpData{con}(s,S.trainidx{con})',S,D.pred_group);
                        stats.BRR.GFP(con).b{d,c}(s,:) = out.muB';
                        stats.BRR.GFP(con).s{d,c}(s) = out.muSigma2;
                        stats.BRR.GFP(con).waic{d,c}(s) = out.waic;
                        stats.BRR.GFP(con).logl{d,c}(s) = out.logl;
                        stats.BRR.GFP(con).r2{d,c}(s) = out.r2;
                        stats.BRR.GFP(con).neglike{d,c}(s) = out.neglike;
                        stats.BRR.GFP(con).r2test{d,c}(s) = out.r2test;
                    end
                end
                
                % all data
                if any(strcmp(S.data_form,'alldata'))
                    if ndims(conData{con})==3
                        alldata = reshape(conData{con},[],size(conData{con},3));
                        for s = 1:size(alldata,1)
                            disp(['Data sample ' num2str(s) '/' num2str(size(alldata,1))])
                            out(s) = bayesreg_crossval(D.pred{con}(S.trainidx{con},:),alldata(s,S.trainidx{con})',S,D.pred_group);
                        end

                        % reshape
                        stats.BRR.alldata(con).b{d,c} = reshape([out(:).muB]',size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).s{d,c} = reshape([out(:).muSigma2],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).waic{d,c} = reshape([out(:).waic],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).logl{d,c} = reshape([out(:).logl],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).r2{d,c} = reshape([out(:).r2],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).neglike{d,c} = reshape([out(:).neglike],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).r2test{d,c} = reshape([out(:).r2test],size(conData{con},1),size(conData{con},2),[]);
                    end
                end
            end
        end
        
        
        %% DECODING METHODS

        % Inverse encoding: requires weights/error from an encoding model
        if isfield(S,'biem_on') && S.biem_on
            
             for con = 1:length(conData)
                 
                % prior: examples x features
                % prior = rand(100,size(D.pred,2));
                if strcmp(S.biem_prior,'subject_training')
                    prior = D.pred{con}(S.trainidx{con},:); 
                elseif strcmp(S.biem_prior,'uniform')
                    prior = 0;
                end
                
                % normalize prior
                prior = zscore(double(prior));
                
                % calculate the covariance of the prior
                % Prior is a multivariate Gaussian with zero mean
                R = cov(prior);
                if numel(R)>1
                    figure; imagesc(R);colormap('hot');
                end

                % run unimodal decoding
                % train_idx fed into this function to calculate mean/std of X/Y in
                % training set in order to normalise test set.
                % Input
                % X: brain data (trials x voxels)
                % Y: design matrix (trials x features)
                % trainidx: trials used for training
                % testidx: trials used for testing
                % B: encoding filters (regression weights)
                % Sigma: estimated residual variances used to define the model
                % R: prior covariance matrix (features x features)
                beta = stats.(S.analysis_type).alldata(con).b{d,c}';
                beta = reshape(beta',1,[]);
                sigma = stats.(S.analysis_type).alldata(con).s{d,c};
                if diff(size(sigma))~=0
                   % if sigma is not a square matrix with variance on the
                   % diagonal
                   sigma = diag(sigma);
                end
                % For decoding, use voxels contributing more to the prediction.
                % Eliminate betas if (normalised) variance more than 0.99
%                 for s = 1:size(sigma,1)
%                     if sigma(s,s)>0.99
%                         sigma(s,s)=1;
%                         beta(:,s)=0;
%                     end
%                 end
                alldata = reshape(conData{con},[],size(conData{con},3));
                [orig, recons] = decode_unimod(alldata', D.pred{con}, S.trainidx{con}, S.testidx{con}, beta, sigma, R);
                test_corr = corr(orig,recons,'type','Spearman');
                disp(['correlation: ' num2str(test_corr)])
                
                if S.num_runs>1
                    stats.biem.alldata(c).rho(d,con) = test_corr;
                else
                    stats.biem.alldata(con).rho(d,c) = test_corr;
                end
                
                % plot orig and reconstruction
                if 0
                    figure
                    for j=1:size(recons,2) 
                        scatter(orig(:,j),recons(:,j));
                        title(j);
                        pause(1); 
                    end
                end
             end
        end
        
        % MVPA: does not require an encoding model
        if isfield(S,'mvpa_on') && S.mvpa_on
            
            for con = 1:length(conData)
                disp(['running contrast' num2str(con)])
            
                % create cosmo data struct
                if iscell(S.trainidx{con}) % classification
                    conds = nan(1,length(S.trainidx{con}));
                    for cn = 1:length(S.trainidx{con})
                        conds(S.trainidx{con}{cn}) = cn;
                    end
                else % regression
                    conds = eventType(idx{con}(S.trainidx{con}));
                end
                cos = eeglab2cosmo(conData{con}(:,:,S.trainidx{con}),S.select_samples,conds);

                % set the targets and chunks
                if iscell(S.trainidx{con}) % classification
                    cos.sa.targets=cos.sa.trialinfo(:,1); 
                else % regression
                    cos.sa.targets=D.pred{con}(S.trainidx{con},:);
                end
                cos.sa.chunks=[1:size(cos.samples,1)]'; % each trial is a chunk

                if S.balance_dataset_and_partitions 
                    if c==1 % use same balancing for all components
                        S.balance_idx=[];
                    end
                end

                % run analysis
                [out,S] = run_cosmo_machine(cos,S);
                disp('MVPA complete')

                % prediction of test sample: dot product of weights with testdata
                testdata = reshape(conData{con}(:,:,S.testidx{con}),[],length(S.testidx{con}));
                out.testdata_pred = (out.weights * testdata) + out.offsets;
                out.testdata_corr = corr(out.testdata_pred',D.pred{con}(S.testidx{con}),'type','Spearman');
                disp(['correlation: ' num2str(out.testdata_corr)])

                if S.num_runs>1
                    stats.mvpa(c).alldata(d,con) = out;
                else
                    stats.mvpa(con).alldata(d,c) = out;
                end
                %stats.mvpa(con).alldata_cv_acc(d,c) = mean(stats.mvpa(d,c).samples);
            end
        end
    end
    % save stats
    if isfield(stats,'RR')
        stats.RR.alldata = rmfield(stats.RR.alldata,'s'); % sigma gets very large and will cause memory problems if accumulated
        stats.RR.GFP = rmfield(stats.RR.GFP,'s'); % sigma gets very large and will cause memory problems if accumulated
    end
    try
        save(fullfile(S.path.stats,['stats_' S.analysis_type '_' S.data_type '_' S.pred '_' S.transform '_' S.sname '.mat']),'stats','S');
    catch
        error('cannot save results')
    end
end

function out = bayesreg_crossval(X,y,S,groupvec);

[X, X_mean, X_stand_de] = zscore(X, [], 1);
[y, y_mean, Y_stand_de] = zscore(y);

N=size(X,1);
K=S.brr.folds;
if K
    x_val_in  = crossvalind('Kfold', N, K);
else
    % if no folds, run all trials for both training and testing
    x_val_in = ones(N,1);
    K=1;
end
for i = K : -1 : 1
    
    val_in    = x_val_in == i;
    es_in     = ~val_in;
    if ~any(es_in); es_in=val_in; end
    
    if S.brr.usegroups
        grps = unique(groupvec);
        for g = grps
            groups{g} = find(groupvec==g);
        end
        [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', S.brr.waic, 'groups',groups);
    else
        [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', S.brr.waic);
    end
    [pred, predstt(i)] = br_predict(X(val_in,:), beta, beta0, stt(i), 'ytest', y(val_in), 'CI', [2.5, 97.5], 'display', false);
    
    logl(i) = stt(i).modelstats.logl; 
    waic(i) = stt(i).modelstats.waic; 
    r2(i) = stt(i).modelstats.r2; 
end

out.muB = mean([stt(:).muB],2);
out.muSigma2 = mean([stt(:).muSigma2]);
out.waic = mean(waic);
out.logl = mean(logl);
out.r2 = mean(r2);
out.neglike = mean([predstt(:).neglike]);
out.r2test = mean([predstt(:).r2]);

