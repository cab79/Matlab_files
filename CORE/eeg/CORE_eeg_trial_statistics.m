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

if ~isfield(S,'pred_transform')
    S.pred_transform = 'notrans';
end

% if S.condor.on  
%     if ~exist(fullfile(S.path.eeg,'input_files'),'dir')
%         mkdir(fullfile(S.path.eeg,'input_files'))
%     end
% end

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
%     if ~any(ismember(S.select_samples,S.total_samples))
%         error('selected samples not suitable for the sampling rate')
%     end
end

%% run though all files in a loop and get data
for d = 1:length(S.select.subjects)
    
    % create output D
    D.subname = S.select.subjects{d};
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,D.subname)))));
    
    if isempty(subfiles)
        continue
    end

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
            if isfield(S.select,'freq')
                timecourse = squeeze(D.fdata.dat{1, 1}.powspctrm(:,:,S.select.freq,:));
                n_chans = length(D.fdata.dat{1, 1}.label);
                n_samples = length(D.fdata.dat{1, 1}.time{1});
                n_trials = length(D.fdata.dat{1, 1}.time);
                timecourse = reshape(permute(timecourse,[2 3 1]),n_chans,[]);
                eventType = D.fdata.dat{1, 1}.conds; tnums = D.fdata.dat{1, 1}.tnums; fnums = D.fdata.dat{1, 1}.fnums; bnums = D.fdata.dat{1, 1}.bnums;
            else
                try
                    topography = D.topography.dat;
                end
                timecourse = D.timecourse.dat;
                eventType = D.eventType.dat;
                n_chans = D.n_chans.dat;
                n_samples = D.n_samples.dat;
                n_trials = D.n_trials.dat;
                tnums=1:length(eventType); % assumes all trials present
            end
        case 'set'
            [eventType,tnums, fnums, bnums] = get_markers(temp);
            
            % make correction for CORE006, part2 data
            % block 2 tnums incorrectly restart from 0, when should start
            % from 452
            if strcmp(D.subname,'CORE006') && length(tnums)<1000 % i.e. it is part2 data
                tnums(bnums==2)=tnums(bnums==2)+451;
            end
            
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
    
    % create predictor variables
    S.pred=[];
    S.pred_group=[];
    S.pred_label={};
    S.pred_traintest={};
    if any(strcmp(S.pred_type,'null'))
        predtemp(1:n_trials,1)=1;
        S.pred = [S.pred predtemp];
        S.pred_group=[S.pred_group 1];
        S.pred_label=[S.pred_label 'null'];
        S.pred_traintest=[S.pred_traintest S.pred_type_traintest(strcmp(S.pred_type,'null'))];
        clear predtemp predgroup predlabel
    end
    if any(strcmp(S.pred_type,'cond'))
        % create sets of trials from S.cond_idx
        for con = 1:length(S.cond)
            condidx{con}=[];
            for row = 1:length(S.cond{con})
                condidx{con} = [condidx{con} find(ismember(eventType,S.cond_idx{S.cond{con}(row)}))];
            end
            % pool into conditions defined by S.cond
            predtemp(condidx{con},1)=con-1;
        end
        S.pred = [S.pred predtemp];
        S.pred_group=[S.pred_group 1];
        S.pred_label=[S.pred_label 'cond'];
        S.pred_traintest=[S.pred_traintest S.pred_type_traintest(strcmp(S.pred_type,'cond'))];
        clear predtemp predgroup predlabel
    end
    if any(strcmp(S.pred_type,'RT'))
        HGF=load(fullfile(S.path.hgf,S.file.hgf)); D_fit=HGF.D_fit; % prevents S loading.
        ind = find(strcmp({D_fit(:).subname},D.subname));
        D.HGF.dat = D_fit(ind).HGF;
        predtemp = D.HGF.dat.y(:,2);
        
        % remove trials not in EEG and put in EEG trial order
        if ~any(diff(tnums)>1)
            error('tnums is probably wrong')
        end
        predtemp=predtemp(tnums,:);
        
        S.pred = [S.pred predtemp];
        S.pred_group=[S.pred_group 2];
        S.pred_label=[S.pred_label 'RT'];
        S.pred_traintest=[S.pred_traintest S.pred_type_traintest(strcmp(S.pred_type,'RT'))];
        clear predtemp predgroup predlabel
    end
    if any(strcmp(S.pred_type,'HGF'))
        % load and format HGF data
        % for each traj group... reformat into predictor matrix
        HGF=load(fullfile(S.path.hgf,S.file.hgf)); D_fit=HGF.D_fit; % prevents S loading.
        if ~isfield(D_fit,'dt')
            dtfile = load(fullfile(S.path.hgf,S.file.get_dt));
            [D_fit(:).dt]=deal(dtfile.D_fit(:).dt);
            [D_fit(:).subname]=deal(dtfile.D_fit(:).subname);
        end
        ind = find(strcmp({D_fit(:).subname},D.subname));
        if isempty(ind)
            disp('USING ONE HGF FOR ALL SUBJECTS')
            ind=1;
        end
        D.HGF = D_fit(ind).HGF;
        [Stemp] = HGF_traj2mat(S,D);
        predtemp=Stemp.pred;
        predgroup=Stemp.pred_group;
        predlabel=Stemp.pred_label;
        predtt=repmat(S.pred_type_traintest(strcmp(S.pred_type,'HGF')),1,length(predlabel));
%         for tg = 1:length(S.traj)
%             for pm = 1:size(S.traj{tg})
%                 pm_name = S.traj{tg}{pm,1};
%                 for tr = 1:length(S.traj{tg}{pm,2})
%                     tr_name = S.traj{tg}{pm,2}(tr);
%                     tr_dat = D.HGF.dat.fit.traj.(pm_name{:}).(tr_name{:});
%                     if ~isempty(S.traj{tg}{pm,4}{tr})
%                         % select levels
%                         tr_dat = tr_dat(:,S.traj{tg}{pm,4}{tr});
%                     end
%                     if S.traj{tg}{pm,3}{1}(tr)>0
%                         if S.traj{tg}{pm,3}{1}(tr)==1
%                             tr_dat = abs(tr_dat); % absolute
%                         elseif S.traj{tg}{pm,3}{1}(tr)==2
%                             tr_dat(tr_dat<0) = 0; % positive rectification
%                         end
%                     end
%                     predtemp = [predtemp, tr_dat];
%                     nvar=size(tr_dat,2);
%                     predgroup = [predgroup, tg*ones(1,nvar)];
%                     for nv = 1:nvar
%                         try
%                             lev = S.traj{tg}{pm,4}{tr}(nv);
%                         catch
%                             lev=nv;
%                         end
%                         predlabel = [predlabel, {[pm_name{:} '_' tr_name{:} num2str(lev)]}];
%                     end
%                 end
%             end
%         end
        
        % remove trials not in EEG and put in order of EEG data trials
        predtemp=predtemp(tnums,:);
        
        S.pred = [S.pred predtemp];
        S.pred_group=[S.pred_group predgroup+length(S.pred_group)];
        S.pred_label=[S.pred_label predlabel];
        S.pred_traintest=[S.pred_traintest predtt];
        clear predtemp predgroup predlabel predtt
    end  
     
    % optionally plot the predictors (e.g. to compare two predictors)
    if 0 
        % normalise
        normpred = zscore(S.pred);
        figure; 
        plot(normpred); 
    end
    
    % remove all-NaN predictors
    ip = all(isnan(S.pred),1);
    S.pred(:,ip) = [];
    S.pred_group(ip)=[];
    S.pred_label(ip)=[];
    
    % Pred data transformation
    if ~isempty(S.pred_transform)
        for pr = 1:size(S.pred,2)
            if strcmp(S.pred_transform,'arcsinh') % need to modify this to apply to only selected predictors
                x=S.pred(:,pr);
                S.pred(:,pr)=log(x+sqrt(x.^2+1));
            
            elseif strcmp(S.pred_transform,'rank') % need to modify this to apply to only selected predictors
                x=S.pred(:,pr);
                [~,S.pred(:,pr)]=sort(x);
            end
        end
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
        
        % flip channels right to left
        if ~isempty(S.flipchan)
            load(S.path.chanlocs);
            
            flipidx=[];
            for row = 1:length(S.flipchan)
                flipidx = [flipidx find(ismember(eventType,S.cond_idx{S.flipchan(row)}))];
            end
            
            data(:,:,flipidx) = flipchan(data(:,:,flipidx),chanlocs,S);
        end
            
        % select samples
        select = dsearchn(S.total_samples',[S.select_samples(1),S.select_samples(end)]');
        data = data(:,select(1):select(2),:);
        
        % create sets of trials split into conditions, according to S.cond_idx
        % this enables subtraction of the mean of one condition from trials
        % of another condition (next step)
        for set = 1:size(S.cond_idx,1)
            setidx{set} = find(ismember(eventType,[S.cond_idx{set}])); % tnum indices organised in sets
            setData{set} = data(:,:,setidx{set}); % data ordered in sets
        end
        
        % perform subtraction of mean of one condition from trials of
        % another. Useful for isolating mismatch responses, i.e. mismatch
        % trials minus the mean of standard trials.
        if ~isempty(S.row_subtract)
            for sb = 1:length(S.row_subtract{1})
                meandat = mean(setData{S.row_subtract{2}(sb)},3);
                subData{sb} = bsxfun(@minus,setData{S.row_subtract{1}(sb)},meandat);
                subidx{sb}=setidx{S.row_subtract{1}(sb)};
            end
            setData=subData;
            setidx = subidx;
        end
            
        % pool sets into "contrasts" (con) according to S.contrast_rows
        % this reduces the number of sets down to only those of interest,
        % e.g. by removing the distinction between left and right-sided
        % stimuli by maintaining distinction between mismatch and standard
        % trials (if standards not subtracted already).
        conData={};
        idx={};
        if ~isempty(S.contrast_rows)
            for con = 1:length(S.contrast_rows)
                conData{con} = cat(3,setData{S.contrast_rows{con}}); % data still ordered according to original sub-sets
                idx{con} = [setidx{S.contrast_rows{con}}]; % indices of tnum for each contrast, so we know what trial order the data is in
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
      
        switch S.traintest
            case 'random'
                for con = 1:length(conData)
                    % Split into training (encoding) and testing (decoding) fractions
                    % Produces S.trainidx and S.testidx. These are indices
                    % of conData-ordered trials that will be used for training and
                    % testing.
                    if S.balance_conds
                        set_events = eventType(idx{con}); % eventypes organised into sets within this contrast (current conData trial order)
                        ucond = unique(set_events);
                        S.trainidx{con} = [];
                        for u = 1:length(ucond)
                            cond_idx = find(set_events==ucond(u)); % for each condition, it's index within conData
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
            case 'cond'
                % split into training (encoding) and testing (decoding)
                % trials by trial type (condition/contrast)
                if length(idx)~=2
                    error('wrong number of contrasts: needs 2')
                end
                conData{1} = cat(3,conData{1},conData{2});
                conData(2)=[]; 
                %S.pred = cat(1,S.pred(idx{1},:),S.pred(idx{2},:));
                S.trainidx{1}= 1:length(idx{1});
                S.testidx{1} = length(idx{1})+1 : length(idx{1})+length(idx{2});
                idx{1} = cat(2,idx{1},idx{2});
                idx(2)=[];
                
        end
        
        % each pred need train/test elements
        % 
        
        % create con-specific predictors and update trainidx to remove NaN
        % predictor value indices
        if c==1
            pred=S.pred;
        end
        S.pred_train={};
        S.pred_test={};
        for con = 1:length(conData)
            % idx is an index of tnum derived from eventType
            S.pred_train{con}=pred(idx{con},strcmp(S.pred_traintest,'train'));
            S.pred_test{con}=pred(idx{con},strcmp(S.pred_traintest,'test'));
            if isempty(S.pred_test{con})
                S.pred_test{con}=S.pred_train{con};
            end
            % remove NaNs
            pred_nonan_train{con}=find(~any(isnan(S.pred_train{con}),2));
            pred_nonan_test{con}=find(~any(isnan(S.pred_test{con}),2));
            %S.pred{con}=S.pred{con}(D.pred_nonan{con},:);
            S.trainidx{con} = S.trainidx{con}(ismember(S.trainidx{con},pred_nonan_train{con}));
            S.testidx{con} = S.testidx{con}(ismember(S.testidx{con},pred_nonan_test{con}));% if no need to remove nans from testidx, this can be commented out
            
            % store info so we know which trials and how many were analysed
            stats.subname{d,1}=D.subname;
            stats.trialinfo{con}.tnums{d,c}=tnums;
            stats.trialinfo{con}.idx{d,c}=idx{con};
            stats.trialinfo{con}.trainidx{d,c}=S.trainidx{con};
            stats.trialinfo{con}.testidx{d,c}=S.testidx{con};
            stats.trialinfo{con}.ntrials_traintest(d,c,:)= [length(S.trainidx{con}),length(S.testidx{con})];
            stats.trialinfo{con}.pred_train{d,c}= S.pred_train{con};
            stats.trialinfo{con}.pred_test{d,c}= S.pred_test{con};
            
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
                    [stats.spear.GFP(con).rho{d,c}, stats.spear.GFP(con).p{d,c}] = corr(gfpData{con}(:,S.trainidx{con})',S.pred_train{con}(S.trainidx{con},:));

                    % FDR correction
                    [stats.spear.GFP(con).fdr_thres{d,c},stats.spear.GFP(con).fdr_masked{d,c}] = fdr(stats.spear.GFP(con).p{d,c},0.05);

                    % TFCE test
                    if isfield(S,'tfce_on') && S.tfce_on
                        shiftd=4-ndims(gfpData{con});
                        img=permute(shiftdim(gfpData{con},-shiftd),[3 1 2 4]);
                        for pr = 1:size(S.pred_train{con},2)
                            varargout=matlab_tfce('correlation',S.tfce_tail,img,[],S.pred_train{con}(S.trainidx{con},pr),S.tfce_nperm);
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
                        [rho, p] = corr(alldata(:,S.trainidx{con})',S.pred_train{con}(S.trainidx{con},:));
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
                
                if size(S.pred_train{con}(S.trainidx{con},:),2)>1
                    % check predictor matrix is full rank
                    rk = rank(S.pred_train{con}(S.trainidx{con},:));
                    fprintf('Design matrix rank = %i\n',rk);
                    if rk<size(S.pred_train{con}(S.trainidx{con},:),2)
                        error('Design matrix is not full rank; predictors are collinear. Reduce the predictors or use a method that is robust to collinearity')
                    end
                end
                
                % zscore OVER TRIALS
                if S.zscore
%                     S.pred_train{con}(S.trainidx{con},:) = zscore(S.pred_train{con}(S.trainidx{con},:));
%                     conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
%                     gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
                    predX=zscore(S.pred_train{con}(S.trainidx{con},:));
                    conX=zscore(conData{con}(:,:,S.trainidx{con}),[],3);
                    gfpX=zscore(gfpData{con}(:,S.trainidx{con}),[],2);
                else
                    predX=S.pred_train{con}(S.trainidx{con},:);
                    conX=conData{con}(:,:,S.trainidx{con});
                    gfpX=gfpData{con}(:,S.trainidx{con});
                end
                
                predX=[ones(size(predX,1),1) predX];
                %S.pred_label = ['constant' S.pred_label];
                %S.pred_group = [0 S.pred_group];
                
                if any(strcmp(S.data_form,'GFP'))
                    % GFP data
                    for s = 1:length(S.select_samples)
                        [beta,~,~,~,stt] = regress(gfpX(s,:)',predX);
                        stats.MR.GFP(con).R2{d,c}(s) = stt(1);
                        stats.MR.GFP(con).F{d,c}(s) = stt(2);
                        stats.MR.GFP(con).p{d,c}(s) = stt(3);
                        stats.MR.GFP(con).s{d,c}(s) = stt(4);
                        stats.MR.GFP(con).b{d,c}(s,:) = beta(2:end);
                    end

                    % FDR correction
                    [stats.MR.GFP(con).fdr_thres{d,c},stats.MR.GFP(con).fdr_masked{d,c}] = fdr(stats.MR.GFP(con).p{d,c},0.05);

                    % TFCE test
                    if isfield(S,'tfce_on') && S.tfce_on
                        shiftd=4-ndims(gfpX);
                        img=permute(shiftdim(gfpX,-shiftd),[3 1 2 4]);
                        varargout=matlab_tfce('regression',S.tfce_tail,img,[],predX, S.tfce_nperm,[],[],[],[],[],zeros(1,size(predX,2)));
                        if ~iscell(varargout); varargout = {varargout}; end
                        % for each (pos and neg):
                        for v = 1:length(varargout)
                            stats.MR.GFP(con).tfce{d,c,v} = varargout{v};
                        end
                    end
                end
                
                if any(strcmp(S.data_form,'alldata'))
                    % all data
                    if ndims(conX)==3
                        alldata = reshape(conX,[],size(conX,3));
                        for s = 1:size(alldata,1)
                            [beta,~,resid,~,stt] = regress(alldata(s,:)',predX);
                            stats.MR.alldata(con).R2{d,c}(s) = stt(1);
                            stats.MR.alldata(con).F{d,c}(s) = stt(2);
                            stats.MR.alldata(con).p{d,c}(s) = stt(3);
                            stats.MR.alldata(con).s{d,c}(s) = stt(4);
                            stats.MR.alldata(con).b{d,c}(s,:) = beta(2:end); % remove constant
                            stats.MR.alldata(con).hnorm{d,c}(s)=kstest(resid);
                            stats.MR.alldata(con).skew{d,c}(s)=skewness(resid);
                            stats.MR.alldata(con).kurt{d,c}(s)=kurtosis(resid);
                            if S.save_residuals
                                stats.MR.alldata(con).resid{d,c}(s,:) = resid;
                            end

                        end
                        % FDR correction
                        [fdr_thresh,fdr_masked] = fdr(stats.MR.alldata(con).p{d,c},0.05);
                        % reshape
                        stats.MR.alldata(con).R2{d,c}=reshape(stats.MR.alldata(con).R2{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).p{d,c}=reshape(stats.MR.alldata(con).p{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).F{d,c}=reshape(stats.MR.alldata(con).F{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).b{d,c}=reshape(stats.MR.alldata(con).b{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).s{d,c}=reshape(stats.MR.alldata(con).s{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).skew{d,c}=reshape(stats.MR.alldata(con).skew{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).kurt{d,c}=reshape(stats.MR.alldata(con).kurt{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).hnorm{d,c}=reshape(stats.MR.alldata(con).hnorm{d,c},size(conX,1),size(conX,2),[]);
                        stats.MR.alldata(con).fdr_thresh{d,c}=fdr_thresh;
                        stats.MR.alldata(con).fdr_masked{d,c}=reshape(fdr_masked,size(conX,1),size(conX,2),[]);
                        if S.save_residuals
                            stats.MR.alldata(con).resid{d,c} = reshape(stats.MR.alldata(con).resid{d,c},size(conX,1),size(conX,2),[]);
                        end
                        
                    end
                end
            end
        end
        
        % Bayesian regression (linear, non-robust to outliers, non-robust to collinearity)
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'PEB')
            
            for con = 1:length(conData)
                
                if size(S.pred_train{con}(S.trainidx{con},:),2)>1
                    % check predictor matrix is full rank
                    rk = rank(S.pred_train{con}(S.trainidx{con},:));
                    fprintf('Design matrix rank = %i\n',rk);
                    if rk<size(S.pred_train{con}(S.trainidx{con},:),2)
                        error('Design matrix is not full rank; predictors are collinear. Reduce the predictors or use a method that is robust to collinearity')
                    end
                end
                
                % zscore OVER TRIALS
                if S.zscore
                    S.pred_train{con}(S.trainidx{con},:) = zscore(S.pred_train{con}(S.trainidx{con},:));
                    conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
                    gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
                end
                
                %S.pred_train{con}=[ones(size(S.pred_train{con},1),1) S.pred_train{con}];

                %enforces positively constraints on the covariance hyperparameters
                % by adopting a log-normal [flat] hyperprior.
                lognormhyper=0;
                
                % To implement non-hierarchical Bayes with priors on the parameters use
                % a two level model setting the second level design matrix to zeros, 
                % providing an unconstrained prior on the first-level parameters
                % which allows for single-level Bayesian inference
                P{1}.X = S.pred_train{con}(S.trainidx{con},:); %as your [n x p] design matrix
                P{1}.C{1}=eye(size(S.pred_train{con}(S.trainidx{con},:),1)); %as your covariance structure on observation errors
                P{2}.X=zeros(size(S.pred_train{con}(S.trainidx{con},:),2),1);
                P{2}.C{2}=eye(size(S.pred_train{con}(S.trainidx{con},:),2)); %as your prior covariance on parameters
                
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
                        
                        
                        checkp = gcp('nocreate')
                        if S.parallel
                            if isempty(checkp)
                                myPool = parpool
                            end
                            parforArg = Inf;
                        else
                            parforArg = 0;
                        end
                        
                        parfor (s = 1:size(alldata,1),parforArg)
                        %for s = 1:size(alldata,1)
                        
                            disp(['Data sample ' num2str(s) '/' num2str(size(alldata,1))])
                            [C,~,F(s)] = spm_PEB(alldata(s,S.trainidx{con})',P,lognormhyper);
                            CE(s,:)=full(C{2}.E);
                        end
                        
%                         for s = 1:size(alldata,1)
%                             [C,~,F(s)] = spm_PEB(alldata(s,S.trainidx{con})',P,lognormhyper);
%                             CE(s,:)=full(C{2}.E);
%                         end

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
                rk = rank(S.pred_train{con}(S.trainidx{con}));
                fprintf('Design matrix rank = %i\n',rk);
                if rk==size(S.pred_train{con}(S.trainidx{con}),2)
                    disp('Design matrix is full rank; can use standard linear regression')
                    %pause(5)
                end

                % zscore OVER TRIALS - occurs within subfunction
%                 if S.zscore
%                     S.pred_train{con}(S.trainidx{con},:) = zscore(S.pred_train{con}(S.trainidx{con},:));
%                     conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
%                     gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
%                 end

                % GFP data
                if any(strcmp(S.data_form,'GFP'))
                    %run ridge regression
                    [Beta_hat, Sigma_hat, X_mean, X_stand_de, Y_mean, Y_stand_de] = ...
                        cross_validation(S.rr.df_num, S.rr.folds, length(S.trainidx{con}), S.pred_train{con}(S.trainidx{con},:)', gfpData{con}(:,S.trainidx{con}), S.rr.z);
                    stats.RR.GFP(con).b{d,c} = Beta_hat;
                    stats.RR.GFP(con).s{d,c} = Sigma_hat;
                end
                
                % all data
                if any(strcmp(S.data_form,'alldata'))
                    if ndims(conData{con})==3
                        alldata = reshape(conData{con},[],size(conData{con},3));
                        %run ridge regression
                        [Beta_hat, Sigma_hat, X_mean, X_stand_de, Y_mean, Y_stand_de] = ...
                            cross_validation(S.rr.df_num, S.rr.folds, length(S.trainidx{con}), S.pred_train{con}(S.trainidx{con},:)', alldata(:,S.trainidx{con}), S.rr.z);
                        stats.RR.alldata(con).b{d,c}=reshape(Beta_hat,size(conData{con},1),size(conData{con},2),[]);
                        stats.RR.alldata(con).s{d,c} = diag(Sigma_hat);
                    end
                end
            end
        end
        
        % Bayesian regularised regression
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'BRR')
            
            for con = 1:length(conData)
                
                % check predictor matrix is NOT full rank
                rk = rank(S.pred_train{con}(S.trainidx{con}));
                fprintf('Design matrix rank = %i\n',rk);
                if rk==size(S.pred_train{con}(S.trainidx{con}),2)
                    disp('Design matrix is full rank; can use standard linear regression')
                    %pause(5)
                end

                % zscore OVER TRIALS - occurs within subfunction
%                 if S.zscore
%                     S.pred_train{con}(S.trainidx{con},:) = zscore(S.pred_train{con}(S.trainidx{con},:));
%                     conData{con}(:,:,S.trainidx{con}) = zscore(conData{con}(:,:,S.trainidx{con}),[],3);
%                     gfpData{con}(:,S.trainidx{con}) = zscore(gfpData{con}(:,S.trainidx{con}),[],2);
%                 end
                
                % GFP data
                if any(strcmp(S.data_form,'GFP'))
                    for s = 1:length(S.select_samples)
                        disp(['GFP sample ' num2str(s) '/' num2str(length(S.select_samples))])
                        out = bayesreg_crossval(S.pred_train{con}(S.trainidx{con},:),gfpData{con}(s,S.trainidx{con})',S,S.pred_group);
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
                        
                        checkp = gcp('nocreate')
                        if S.parallel
                            if isempty(checkp)
                                myPool = parpool
                            end
                            parforArg = Inf;
                        else
                            parforArg = 0;
                        end
                        
                        if ~S.condor.on
                            parfor (s = 1:size(alldata,1),parforArg)
                            %for s = 1:size(alldata,1)

                                disp(['Data sample ' num2str(s) '/' num2str(size(alldata,1))])
                                out(s) = bayesreg_crossval(S.pred_train{con}(S.trainidx{con},:),alldata(s,S.trainidx{con})',S,S.pred_group);

                                % residuals
                                yhat = S.pred_train{con}(S.trainidx{con},:)*out(s).muB; % Predicted responses at each data point.
                                resid = alldata(s,S.trainidx{con})'-yhat; % Residuals.
                                try
                                    hnorm(s)=kstest(resid);
                                    skew(s)=skewness(resid);
                                    kurt(s)=kurtosis(resid);
                                    if S.save_residuals
                                        resids(s,:) = resid;
                                    end
                                catch
                                    hnorm(s)=NaN;
                                    skew(s)=NaN;
                                    kurt(s)=NaN;
                                end
                                
                                % std
                                sd(s) = std(alldata(s,S.trainidx{con}));
                            end
                        else
                            % chunk the data for Condor - ideal run time is
                            % about 10 mins. runs at about 250 samples per
                            % minutes, per processor, so 2500 sample is
                            % ideal.
                            n_chunks = ceil(size(alldata,1)/S.condor.chunksize);
                            chunksize = ceil(size(alldata,1)/n_chunks);
                            
                            for nc = 1:n_chunks
                                chunk_index = chunksize*(nc-1)+1 : min(chunksize*nc,size(alldata,1));
                                
                                S.data_info.S=size(alldata,1);
                                S.data_info.D=length(S.select.subjects);
                                S.data_info.C=length(iter);
                                S.data_info.CON=length(conData);
                                S.data_info.n_chunks = n_chunks;
                                S.data_info.chunk_index=chunk_index;
                                S.data_info.d=d;
                                S.data_info.c=c;
                                S.data_info.con=con;
                                S.data_info.dim=size(conData{con});
                                X = S.pred_train{con}(S.trainidx{con},:);
                                Y = alldata(chunk_index,S.trainidx{con})';

                                % index
                                condor_index = (d-1)*length(iter)*length(conData)*n_chunks +(c-1)*length(conData)*n_chunks +(con-1)*n_chunks +nc;
                                index_length = length(S.select.subjects)*length(iter)*length(conData)*n_chunks;

                                % save
                                disp(['creating input file ' num2str(condor_index) '/' num2str(index_length)])
                                save(fullfile(S.path.stats,['input' num2str(condor_index-1) '.mat']),'X','Y','S','stats');
                            end
                            continue
                        end

                        % reshape
                        stats.BRR.alldata(con).b{d,c} = reshape([out(:).muB]',size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).s{d,c} = reshape([out(:).muSigma2],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).waic{d,c} = reshape([out(:).waic],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).logl{d,c} = reshape([out(:).logl],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).r2{d,c} = reshape([out(:).r2],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).neglike{d,c} = reshape([out(:).neglike],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).r2test{d,c} = reshape([out(:).r2test],size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).skew{d,c}=reshape(skew,size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).kurt{d,c}=reshape(kurt,size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).hnorm{d,c}=reshape(hnorm,size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).sd{d,c}=reshape(sd,size(conData{con},1),size(conData{con},2),[]);
                        stats.BRR.alldata(con).pred{d,c}=S.pred_train{con}(S.trainidx{con},:);
                        if S.save_residuals
                            timecourse = reshape(resids,size(conData{con},1),size(conData{con},2),[]);
                            [n_chans,n_samples,n_trials] = size(timecourse);
                            resid_sname = strrep(filename,['.' S.fname.ext{:}],['_resid_comp' num2str(c) '_con' num2str(con) '_' S.sname '.mat']);
                            save(fullfile([S.path.stats '\residuals'], resid_sname),'timecourse','eventType','n_chans','n_samples','n_trials')
                        end
                        
                        clear out hnorm skew kurt resids;
                    end
                end
            end
        end
        
        
        %% DECODING METHODS

        % Bayesian inverse encoding model: requires weights/error from an encoding model
        if isfield(S,'biem_on') && S.biem_on
            
             for con = 1:length(conData)
                 
                 % update according to whether data should be averaged
                 if ~isempty(S.biem_trialavg_pre)
                     for ii = 1:length(S.biem_trialavg_pre)
                         trialind = find(ismember(eventType,S.biem_trialavg_pre{ii}));
                        datatemp(:,:,ii) = mean(conData{con}(:,:,trialind),3);
                        pred_train_temp(ii,:) = mean(S.pred_train{con}(trialind,:),1);
                        pred_test_temp(ii,:) = mean(S.pred_test{con}(trialind,:),1);
                     end
                         
                    conData{con} = datatemp;
                    S.pred_train{con} = pred_train_temp;
                    S.trainidx{con} = 1:length(S.biem_trialavg_pre);
                    S.pred_test{con} = pred_test_temp;
                    S.testidx{con} = 1:length(S.biem_trialavg_pre);
                 end
                 
                % prior: examples x features
                % prior = rand(100,size(S.pred_train,2));
                if strcmp(S.biem_prior,'subject_training')
                    prior = S.pred_train{con}(S.trainidx{con},:); 
                elseif strcmp(S.biem_prior,'uniform')
                    prior = nan; % unfinished
                end
                
                % normalize prior
                prior_nonan = ~isnan(prior);
                prior(prior_nonan) = zscore(double(prior(prior_nonan)));
                
                % calculate the covariance of the prior
                % Prior is a multivariate Gaussian with zero mean
                R = cov(prior(prior_nonan));
                if 0%numel(R)>1
                    figure; imagesc(R);colormap('hot');
                end
                
                for dtype = 1:length(S.data_form)
                    datatype = S.data_form{dtype};

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
                    
                    % load group weights if specified
                    grpweights=1;
                    if ~isempty(S.biem_groupweights)
                        grpstats=load(fullfile(S.path.stats,S.biem_groupweights));
                        try
                            grpbeta = grpstats.stats.(S.biem_group_analysis_type).(datatype).b_mean;
                            grpsigma = grpstats.stats.(S.biem_group_analysis_type).(datatype).s_mean;
                        catch
                            disp('no group weights!!')
                            grpweights=0;
                        end
                    else
                        grpweights=0;
                    end
                    
                    % get subject weights
                    switch datatype 
                        case 'alldata'
                            data = reshape(conData{con},[],size(conData{con},3));
                            beta = stats.(S.analysis_type).(datatype)(con).b{d,c};
                            sigma = stats.(S.analysis_type).(datatype)(con).s{d,c};
                            beta = reshape(beta,[],size(S.pred_train{con},2))';
                            if grpweights; grpbeta = reshape(grpbeta,1,[]);end
                        case 'GFP'
                            data = gfpData{con};
                            beta = stats.(S.analysis_type).(datatype)(con).b{d,c}';
                            sigma = stats.(S.analysis_type).(datatype)(con).s{d,c};
                    end
                    
                    % if sigma is not a square matrix with variance on the diagonal
                    if diff(size(sigma))~=0
                       sigma = reshape(sigma,1,[]);
                       sigma = diag(sigma);
                    end
                    if grpweights
                        if diff(size(grpsigma))~=0
                           grpsigma = reshape(grpsigma,1,[]);
                           grpsigma = diag(grpsigma);
                        end
                    end
                    
                    % decode
                    [orig, recons] = decode_unimod(data', S.pred_train{con}, S.trainidx{con}, S.testidx{con}, beta, sigma, R);
                    if grpweights
                        [grporig, grprecons] = decode_unimod(data', S.pred_train{con}, S.trainidx{con}, S.testidx{con}, grpbeta, grpsigma, R);
                    end
                    
                    % post-decoding averaging
                     if ~isempty(S.biem_trialavg_post)
                         for ii = 1:length(S.biem_trialavg_post)
                             trialind = find(ismember(eventType,S.biem_trialavg_post{ii}));
                            origtemp(ii,:) = mean(orig(trialind,:),1);
                            recontemp(ii,:) = mean(recons(trialind,:),1);
                            pred_train_temp(ii,:) = mean(S.pred_train{con}(trialind,:),1);
                            pred_test_temp(ii,:) = mean(S.pred_test{con}(trialind,:),1);
                            if grpweights
                                grporigtemp(ii,:) = mean(grporig(trialind,:),1);
                                grprecontemp(ii,:) = mean(grprecons(trialind,:),1);
                            end
                         end

                        orig = origtemp;
                        recons = recontemp;
                        S.trainidx{con} = 1:length(S.biem_trialavg_post);
                        S.testidx{con} = 1:length(S.biem_trialavg_post);
                        
                        if grpweights
                            grporig = grporigtemp;
                            grprecons = grprecontemp;
                        end
                     end
                    
                    % test correlations
                    clear testvar
                    if isequal(S.pred_train{con},S.pred_test{con})
                        testvar=orig;
                        nonan = ~isnan(testvar(:,S.biem_pred));
                    elseif size(S.pred_test{con}(S.testidx{con}),2)==1
                        testvar=S.pred_test{con}(S.testidx{con});
                        nonan = ~isnan(testvar(:,1));
                    else
                        error('check input variables')
                    end
                    
                    test_corr=[];test_grpcorr=[];
                    if size(recons,1)>1
                        for pn = 1:min(length(S.biem_pred),size(recons,2))
                            test_corr(pn) = corr(testvar(nonan),recons(nonan,S.biem_pred(pn)),'type','Spearman');
                            disp(['correlation: ' num2str(test_corr(pn))])
                        end
                        if grpweights
                            for pn = 1:length(S.biem_pred)
                                test_grpcorr(pn) = corr(testvar(nonan),grprecons(nonan,S.biem_pred(pn)),'type','Spearman');
                                disp(['correlation from grp weights: ' num2str(test_grpcorr(pn))])
                            end
                        end
                    end
                    
                    % store results
                    if S.num_runs>1
                        stats.biem.(datatype)(c).rho{d,con} = test_corr;
                        stats.biem.(datatype)(c).recons{d,con} = recons;
                        stats.biem.(datatype)(c).orig{d,con} = testvar;
                    else
                        stats.biem.(datatype)(con).rho{d,c} = test_corr;
                        stats.biem.(datatype)(con).recons{d,c} = recons;
                        stats.biem.(datatype)(con).orig{d,c} = testvar;
                    end
                    if grpweights
                        if S.num_runs>1
                            stats.biem.(datatype)(c).grprho{d,con} = test_grpcorr;
                        else
                            stats.biem.(datatype)(con).grprho{d,c} = test_grpcorr;
                            stats.biem.(datatype)(con).grprecons{d,c} = grprecons;
                        end
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
        end
        
        % MVPA: does not require an encoding model
        if isfield(S,'mvpa_on') && S.mvpa_on
            
            %if strcmp(S.mvpa_type,'class')
            %    runcon = 1;
            %else
                runcon = 1:length(conData);
            %end
            
            for con = 1:length(conData)
                disp(['running contrast' num2str(con)])
            
                % create cosmo data struct
                if strcmp(S.mvpa_type,'class') % classification
                    if ~isempty(S.contrast_rows)
                        conds = nan(1,length([S.trainidx{:}]));
                        conds=[];
                        trainidx=[];
                        testidx=[];
                        for cn = 1:length(conData)
                            conds = [conds, cn*ones(1,length(conData{cn}))];
                            if cn-1
                                trainidx = [trainidx, S.trainidx{cn} + size(conData{cn-1},3)];
                                testidx = [testidx, S.testidx{cn} + size(conData{cn-1},3)];
                            else
                                trainidx = [trainidx, S.trainidx{cn}];
                                testidx = [testidx, S.testidx{cn}];
                            end
                        end
                        conData = cat(3,conData{:});
                        cos = eeglab2cosmo(conData(:,:,trainidx),S.select_samples,conds(trainidx));
                        testdata = reshape(conData(:,:,testidx),[],length(testidx));
                    else
                        conds=S.pred_train{con}+1;
                        cos = eeglab2cosmo(conData{con}(:,:,S.trainidx{con}),S.select_samples,conds(S.trainidx{con}));
                        testdata = reshape(conData{con}(:,:,S.testidx{con}),[],length(S.testidx{con}));
                        testidx=S.testidx{con};
                    end
                else % regression
                    conds = eventType(idx{con}(S.trainidx{con}));
                    cos = eeglab2cosmo(conData{con}(:,:,S.trainidx{con}),S.select_samples,conds);
                    testdata = reshape(conData{con}(:,:,S.testidx{con}),[],length(S.testidx{con}));
                end

                % set the targets and chunks
                if strcmp(S.mvpa_type,'class') % classification
                    cos.sa.targets=cos.sa.trialinfo(:,1); 
                else % regression
                    cos.sa.targets=S.pred_train{con}(S.trainidx{con},:);
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
                if strcmp(S.mvpa_type,'class') % classification
                    
                    % variable to predict
                    predicted=S.pred_test{con}(testidx)+1;
                    out.predicted=predicted;
                    
                    % weight predictions
                    out.testdata_pred = (out.weights(1,:) * testdata);
                    if length(unique(predicted))==2 % if predicted consists of two values
                        out.testdata_predbin=nan(1,length(out.testdata_pred));
                        out.testdata_predbin(out.testdata_pred>0)=1;
                        out.testdata_predbin(out.testdata_pred<0)=2; % this is opposite to intuition for GPC
                        out.testdata_corr = -corr(out.testdata_pred',predicted,'type','Spearman');
                        disp(['correlation: ' num2str(out.testdata_corr)])
                        out.testdata_correct = sum(out.testdata_predbin==predicted')/length(out.testdata_pred);
                        disp(['correct: ' num2str(out.testdata_correct)])
                    else
                        out.testdata_corr = corr(out.testdata_pred',predicted,'type','Spearman');
                        disp(['correlation: ' num2str(out.testdata_corr)])
                    end
                    % transweight predictions
                    out.testdata_predtrans = (out.transweights(1,:) * testdata);
                    if length(unique(predicted))==2 % if predicted consists of two values
                        out.testdata_predbintrans=nan(1,length(out.testdata_predtrans));
                        out.testdata_predbintrans(out.testdata_predtrans>0)=1;
                        out.testdata_predbintrans(out.testdata_predtrans<0)=2; % this is opposite to intuition for GPC
                        out.testdata_corrtrans = -corr(out.testdata_predtrans',predicted,'type','Spearman');
                        disp(['correlation trans: ' num2str(out.testdata_corrtrans)])
                        out.testdata_correcttrans = sum(out.testdata_predbintrans==predicted')/length(out.testdata_predtrans);
                        disp(['correct trans: ' num2str(out.testdata_correcttrans)])
                    else
                        out.testdata_corrtrans = corr(out.testdata_predtrans',predicted,'type','Spearman');
                        disp(['correlation trans: ' num2str(out.testdata_corrtrans)])
                    end
                else % regression
                    % weight predictions
                    out.testdata_pred = (out.weights * testdata) + out.offsets;
                    out.testdata_corr = corr(out.testdata_pred',S.pred_test{con}(S.testidx{con}),'type','Spearman');
                    disp(['correlation weigh: ' num2str(out.testdata_corr)])
                    % transweight predictions
                    out.testdata_predtrans = (out.transweights * testdata) + out.offsets;
                    out.testdata_corrtrans = corr(out.testdata_predtrans',S.pred_test{con}(S.testidx{con}),'type','Spearman');
                    disp(['correlation trans: ' num2str(out.testdata_corrtrans)])
                end
                
                % reshape weights to 2D
                sizdat = size(conData{con});
                out.weights = reshape(out.weights(1,:),sizdat(1),sizdat(2));
                out.transweights = reshape(out.transweights(1,:),sizdat(1),sizdat(2));

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
%     if isfield(stats,'RR')
%         if isfield(stats.RR,'alldata')
%             stats.RR.alldata = rmfield(stats.RR.alldata,'s'); % sigma gets very large and will cause memory problems if accumulated
%         end
%         if isfield(stats.RR,'GFP')
%             stats.RR.GFP = rmfield(stats.RR.GFP,'s'); % sigma gets very large and will cause memory problems if accumulated
%         end
%     end
    if ~S.condor.on
        try
            save(fullfile(S.path.stats,['stats_' S.analysis_type '_' S.data_type '_' S.pred_type{:} '_' S.transform '_' S.sname '.mat']),'stats','S');
        catch
            error('cannot save results')
        end
    end
end
if S.condor.on
    create_job_submission_file(S.path.stats,index_length)
    quit
end

function create_job_submission_file(pth,nf)

disp('creating job submission file')
A = {
    'executable=CORE_condor_EEG_job.exe'
    'indexed_input_files=input.mat'
    'indexed_output_files=output.mat'
    'indexed_stdout=CORE_condor_fit_job.out'
    'indexed_stderr=CORE_condor_fit_job.err'
    'indexed_log=CORE_condor_fit_job.log'
    'max_run_time=60'
    ['total_jobs=' num2str(nf)]
};

fid = fopen(fullfile(pth, ['CORE_condor_EEG_job_run.sub']),'w');
for i = 1:size(A,1)
    fprintf(fid,'%s\n',A{i});
end
fclose(fid);
