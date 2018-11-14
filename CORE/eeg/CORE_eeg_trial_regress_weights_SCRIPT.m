% Representational similiarity between encoding models.
% Takes regression weight matrices from encoding models and regresses
% these weights. Performs model comparison to answer questions such as:
% Which neural representations (e.g. from a cognitive model) are most
% similar to representations predicting response time?

% preliminaries
clear all
close all
dbstop if error
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

% path and data
S.file = mfilename;
S.spath = 'C:\Data\CORE\eeg\ana\stats';
S.iv_files = {
    %'stats_BRR_all_chan_HGF_arcsinh_20180831T160611.mat'; % real PE
%     'stats_BRR_all_chan_HGF_arcsinh_20180901T124507.mat' % real Dau epsi2'
%     'stats_BRR_all_chan_HGF_arcsinh_20180901T124507.mat' % real Dau epsi2'
%     'stats_BRR_all_chan_HGF_arcsinh_20180901T124507.mat' % real Dau epsi2'
    %'stats_BRR_all_chan_HGF_arcsinh_20180831T205032.mat' % real Dau
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    'stats_BRR_all_chan_condHGF_arcsinh_20180831T131157.mat'
    %'stats_BRR_all_chan_cond_arcsinh_20180830T145121.mat';
    };
S.dv_file = 'stats_BRR_all_chan_RT_arcsinh_20180813T211745.mat';
S.use_iv = {[1],[2:7],[2],[3],[4],[5],[6],[7]}; % which IVs to include in regression?
S.nperm=0; % number of permutations
S.plot_scatter=0;

% fields within data structures
S.encoding_type = 'BRR';
S.statfield = 'b';

% Bayesian regularised regression (BRR) settings
S.brr.folds = 0;            % number of folds in the traindata. Set to 0 to not conduct predictive cross-validation.
S.brr.model = 'gaussian';   % error distribution - string, one of {'gaussian','laplace','t','binomial'}
S.brr.prior = 'g';        %- string, one of {'g','ridge','lasso','horseshoe','horseshoe+'}
S.brr.nsamples = 1000;   %- number of posterior MCMC samples (Default: 1000)  
S.brr.burnin = 1000;     %- number of burnin MCMC samples (Default: 1000)
S.brr.thin = 5;       %- level of thinning (Default: 5)
S.brr.catvars = [];    %- vector of variables (as column numbers of X) to treat
%                       as categorical variables, with appropriate expansion. 
%                       See examples\br_example5 (Default: none)
S.brr.nogrouping = false; %- stop automatic grouping of categorical predictors
%                       that is enabled with the 'catvars' options. (Default: false)
S.brr.usegroups = 0;     % ****Specified by S.traj cells**** - create groups of predictors. Grouping of variables
%                       works only with HS, HS+ and lasso prior
%                       distributions. The same variable can appear in
%                       multiple groups. See examples\br_example[9,10,11,12]  (Default: { [] } )  
S.brr.waic = true;       %- whether to calculate the WAIC -- disabling can lead
%                       to large speed-ups, especially for Gaussian models with large n
%                       (default: true)
S.zscore=1;

for ii = 1:length(S.iv_files)

    % load data
    temp=load(fullfile(S.spath,S.iv_files{ii}));
    iv_stat = temp.stats;
    temp=load(fullfile(S.spath,S.dv_file));
    dv_stat = temp.stats;
    clear temp

    % check data and produce subject indices
    if ~isequal(dv_stat.subname,iv_stat.subname)
        disp('data has different subjects - attempting correction')
        iv_ind = ismember(iv_stat.subname,dv_stat.subname);
        dv_ind = ismember(dv_stat.subname,iv_stat.subname);

        if ~isequal(dv_stat.subname(dv_ind),iv_stat.subname(iv_ind))
            error('subject orders do not match')
        end
    else
        iv_ind = ones(length(iv_stat.subname),1);
        dv_ind = ones(length(dv_stat.subname),1);
    end

    % initiate figures
    if S.plot_scatter
        for i = 1:length(S.use_iv{ii})
            fig(i)=figure('Name',['IV ' num2str(i)]);
        end
    end

    % cycle through subjects and regress
    iv_all = iv_stat.(S.encoding_type).alldata.(S.statfield)(iv_ind);
    dv_all = dv_stat.(S.encoding_type).alldata.(S.statfield)(dv_ind);
    nsub=length(iv_all);
    spdim = ceil(sqrt(nsub)); % subplot dims
    for ns = 1:nsub
        disp(['running subject ' num2str(ns) ' ...'])
        % get subject data
        dv = dv_all{ns};
        iv = iv_all{ns}(:,:,S.use_iv{ii});

        % reshape to 1D
        datsize = size(dv);
        dv=reshape(dv,prod(datsize),1);
        iv=squeeze(reshape(iv,prod(datsize),1,length(S.use_iv{ii})));

        % multiple regression
    %     dv=zscore(dv);
    %     iv=zscore(iv,[],1);
    %     iv = [ones(length(iv),1),iv]; % add columns of 1s to IV
    %     [beta,~,~,~,stt] = regress(dv,iv);
    %     stats(ns).R2 = stt(1);
    %     stats(ns).F = stt(2);
    %     stats(ns).p = stt(3);
    %     stats(ns).s = stt(4);
    %     stats(ns).b = beta(2:end); % remove constant

        BRR = bayesreg_crossval(iv,dv,S);
        stats(ns,ii).R2 = BRR.r2;
        stats(ns,ii).waic = BRR.waic;
        stats(ns,ii).b = BRR.muB;

        % plot
        if S.plot_scatter
            for i = 1:length(S.use_iv{ii})
                figure(fig(i));
                subplot(spdim,spdim,ns);
                scatter(iv(:,i),dv);
                title(num2str(ns));
            end
        end

        % perform randomization to obtain statistical significance for partial coefficients
        if S.nperm
            beta_perm = nan(S.nperm,length(S.use_iv{ii}));
            for n=1:S.nperm
                disp(['running subject ' num2str(ns) ', permutation ' num2str(n)])

                % vector of indices (a random ordering of the integers between
                % 1 and n where n is the number of data points)
                iperm = randperm(length(dv));

                % compute regression between randomly shuffled IVs and the DVs
                [beta,~,~,~,stt] = regress(dv,iv(iperm,:));
                stats(ns,ii).b_perm(n,:) = beta(2:end);
            end

            % calculate permuted p value
            for bn = 1:length(S.use_iv{ii})
                stats(ns,ii).p_perm(bn) = sum(abs(stats(ns,ii).b_perm(:,bn)) > abs(stats(ns,ii).b(bn))) / length(stats(ns,ii).b_perm(:,bn));
            end
        end

    end

    % betas
    b_all = cat(2,stats(:,ii).b);
    r2_all = cat(2,stats(:,ii).R2);

    % one sample t-test on betas
    clear p h
    for bi = 1:size(b_all,1)
        [h(bi),p(bi),~,~] = ttest(double(b_all(bi,:)));
    end

    % plots
    x=repmat([1:size(r2_all,1)]',1,size(r2_all,2));
    figure;scatter(x(:),r2_all(:))
    title(['R2 for all subjects'])
    figure;boxplot(r2_all(:))
    title(['R2 for all subjects'])

    x=repmat([1:size(b_all,1)]',1,size(b_all,2));
    figure;scatter(x(:),b_all(:))
    title(['beta weights for each predictor'])
end

sname = ['stat_regress_weights_' datestr(now,30)];
save(fullfile(S.spath,sname),'stats');

if ii>1
    LME=reshape(cat(1,stats.waic),[],ii)';
    [bmc.posterior,bmc.out] = VBA_groupBMC(LME)
end
