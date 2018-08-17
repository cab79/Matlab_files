dbstop if error
clear all
close all

% y column to analyse
% 0 = u col
% 1 = Choice
% 2 = RTs
% 3 = decoded EEG variables
ycol = 2;
S.use_group_recons=0; % for decoded EEG only

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'part2'; %
S.path.seq = ['C:\Data\CORE\design']; % unprocessed data in original format
S.path.raw = ['C:\Data\CORE\behaviour\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\CORE\behaviour\processed']; % folder to save processed data
S.path.hgf = ['C:\Data\CORE\behaviour\hgf']; % folder to save processed data
S.path.design = ['C:\Data\CORE\design']; % 
S.fname.parts = {'prefix','subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.groups = {};
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
%save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

% data import
S.load.prefixes = {'RT','dt'};
S.load.suffix = {'*'};
[S,D]=SCIn_data_import(S);

% response data processing
S.fitsim=1; % is the analysis of recorded data (1) or simulated data (2)?
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 1;
S.RT.on = 1;
S.RT.min = 0.2; % min RT to consider
S.save.tables = 0;
[S,D_prep]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

% Import decoded EEG
if ycol==3
    S.path.stats = 'C:\Data\CORE\eeg\ana\stats';
    S.file = 'stats_MR_all_chan_cond_arcsinh_20180804T144006.mat';
    [S,D_prep]=CORE_decEEG_import(S,D_prep,3);  
end

% which subjects?
run_subs=1:length(D_prep)
%run_subs(6) = []; % remove outlier
D_prep=D_prep(run_subs);

% Bayesian regularised regression (BRR) settings
S.brr.folds = 0;            % number of folds in the traindata. Set to 0 to not conduct predictive cross-validation.
S.brr.model = 'gaussian';   % error distribution - string, one of {'gaussian','laplace','t','binomial'}
S.brr.prior = 'horseshoe+';        %- string, one of {'g','ridge','lasso','horseshoe','horseshoe+'}
S.brr.nsamples = 1000;   %- number of posterior MCMC samples (Default: 1000)  
S.brr.burnin = 1000;     %- number of burnin MCMC samples (Default: 1000)
S.brr.thin = 5;       %- level of thinning (Default: 5)
S.brr.catvars = [];    %- vector of variables (as column numbers of X) to treat
%                       as categorical variables, with appropriate expansion. 
%                       See examples\br_example5 (Default: none)
S.brr.nogrouping = false; %- stop automatic grouping of categorical predictors
%                       that is enabled with the 'catvars' options. (Default: false)
S.brr.usegroups = 1;     % ****Specified by S.traj cells**** - create groups of predictors. Grouping of variables
%                       works only with HS, HS+ and lasso prior
%                       distributions. The same variable can appear in
%                       multiple groups. See examples\br_example[9,10,11,12]  (Default: { [] } )  
S.brr.waic = true;       %- whether to calculate the WAIC -- disabling can lead
%                       to large speed-ups, especially for Gaussian models with large n
%                       (default: true)

% Load Bayes optimal results
fname = 'CORE_fittedparameters_percmodel10_bayesopt_20180805T093044.mat';
ls=load(fullfile(S.path.hgf,'fitted',fname));
D_fit=ls.D_fit; 

%% HGF trajectory models
% 1st col: model type
% 2nd col: traj names 
% 3rd col: 0=no transformation, 1=absolute, 2=rectified, 3=exp
% 4th col: select levels

S.trajmodels(1).traj{1} = {
     {'PL'},{'sa','sahat','pv','inferv','surp'},{[0 0 0 0 0]},{[],[],[],[],[]};
     }; 
S.trajmodels(1).traj{2} = {
     {'PL'},{'dau','da','epsi'},{[0 0 0]},{[],[],[]};
     }; 
 
S.trajmodels(2).traj{1} = {
     {'PL'},{'sa','sahat','pv','inferv','surp'},{[0 0 0 0 0]},{[],[],[],[],[]};
     }; 
S.trajmodels(2).traj{2} = {
     {'PL'},{'dau','da','epsi'},{[1 1 1]},{[],[],[]};
     }; 
 
S.trajmodels(3).traj{1} = {
     {'PL'},{'sa','sahat','pv','inferv','surp'},{[0 0 0 0 0]},{[],[],[],[],[]};
     }; 
S.trajmodels(3).traj{2} = {
     {'PL'},{'dau','da','epsi'},{[2 2 2]},{[],[],[]};
     }; 
 

S.trajmodels(4).traj{1} = {
     {'PL'},{'surp','dau','da','epsi'},{[0 0 0 0]},{[],[],[1:2],[]};
     }; 
 
S.trajmodels(5).traj{1} = {
     {'PL'},{'surp','dau','da','epsi'},{[0 1 1 1]},{[],[],[1:2],[]};
     }; 
 
S.trajmodels(6).traj{1} = {
     {'PL'},{'surp','dau','da','epsi'},{[0 2 2 2]},{[],[],[1:2],[]};
     }; 
 
% beliefs
% S.trajmodels(1).traj{1} = {
%      {'PL'},{'mu','muhat'},{[0 0]},{[],[]};
%      }; 
 
% belief variances
% S.trajmodels(2).traj{1} = {
%      {'PL'},{'sa',},{[0]},{[]};
%      }; 
%  
% % belief variances
% S.trajmodels(3).traj{1} = {
%      {'PL'},{'sahat',},{[0]},{[]};
%      }; 
 
% % prediction errors
% S.trajmodels(4).traj{1} = {
%      {'PL'},{'dau','da'},{[0 0]},{[],[1:2]};
%      }; 
%  
% % prediction errors - abs
% S.trajmodels(5).traj{1} = {
%      {'PL'},{'dau','da'},{[1 1]},{[],[1:2]};
%      }; 
%  
% % prediction errors - rectified
% S.trajmodels(6).traj{1} = {
%      {'PL'},{'dau','da'},{[2 2]},{[],[1:2]};
%      }; 
%  
% % prediction errors (precision-weighted)
% S.trajmodels(7).traj{1} = {
%      {'PL'},{'epsi'},{[0]},{[1:2]};
%      }; 
%  
% % prediction errors - abs (precision-weighted)
% S.trajmodels(8).traj{1} = {
%      {'PL'},{'epsi'},{[1]},{[1:2]};
%      }; 
%  
% % prediction errors - rectified (precision-weighted)
% S.trajmodels(9).traj{1} = {
%      {'PL'},{'epsi'},{[2]},{[1:2]};
%      }; 
%  
% % combined
% S.trajmodels(10).traj{1} = {
%      {'PL'},{'mu','sa','muhat','sahat'},{[0 0 0 0]},{[],[],[],[]};
%      }; 
% S.trajmodels(10).traj{2} = {
%      {'PL'},{'da','dau','epsi'},{[0 0 0 0]},{[],[],[],[]};
%      }; 
% 
% % combined (abs)
% S.trajmodels(11).traj{1} = {
%      {'PL'},{'mu','sa','muhat','sahat'},{[0 0 0 0]},{[],[],[],[]};
%      }; 
% S.trajmodels(11).traj{2} = {
%      {'PL'},{'da','dau','epsi'},{[1 1 1 1]},{[],[],[],[]};
%      }; 
%  
% % combined (rect)
% S.trajmodels(12).traj{1} = {
%      {'PL'},{'mu','sa','muhat','sahat'},{[0 0 0 0]},{[],[],[],[]};
%      }; 
% S.trajmodels(12).traj{2} = {
%      {'PL'},{'da','dau','epsi'},{[2 2 2 2]},{[],[],[],[]};
%      }; 
%  
% % benchmark uncertainty model 1
% S.trajmodels(13).traj{1} = {
%      {'PL'},{'mu','sa','surp','sahat'},{[3 0 0 0]},{[3],[2],[],[1]}; % exp(mu3), sa2, surprise, sa1hat, 
%      }; 
% 
% % benchmark uncertainty model 2
% S.trajmodels(14).traj{1} = {
%      {'PL'},{'pv','inferv','surp','sahat'},{[0 0 0 0]},{[],[],[],[1]}; % pv, inferv, surprise, sa1hat, 
%      }; 
%  
% % beliefs and their variances
% S.trajmodels(15).traj{1} = {
%      {'PL'},{'mu','sa','muhat','sahat'},{[0 0 0 0]},{[],[],[],[]};
%      }; 
%  
% % beliefs and their variances
% S.trajmodels(16).traj{1} = {
%      {'PL'},{'sa','dau','da'},{[0 0 0]},{[],[],[1:2]};
%      }; 
%  
% % comb
% S.trajmodels(17).traj{1} = {
%      {'PL'},{'sa','dau','da'},{[0 1 1]},{[],[],[1:2]};
%      }; 
%  
% % combined
% S.trajmodels(18).traj{1} = {
%      {'PL'},{'sa','dau','da'},{[0 2 2]},{[],[],[1:2]};
%      }; 

% load and format HGF data
% for each traj group... reformat into predictor matrix
S.pred_transform = 'notrans'; % arcsinh, rank or notrans
S.zscore = 1;
S.trajmodelselect = [1 2 3 4 5 6];

for tm = S.trajmodelselect
    S.traj = S.trajmodels(tm).traj;
    [S] = HGF_traj2mat(S,D_fit(1));

    % remove nan columns from S.pred
    rmnan = all(isnan(S.pred),1);
    S.pred(:,rmnan)=[];
    S.pred_label(rmnan)=[];
    S.pred_group(rmnan)=[];

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

    predcat=[];
    Vcat = [];

    for d = 1:length(D_prep)
        % add Bayes opt fit to D_prep containing the y struct
        D_prep(d).HGF.fit=D_fit.HGF.fit;

        % V: variable of interest to predict
        if ycol
            V=D_prep(d).HGF.y(:,ycol);
        else
            V=D_prep(d).HGF.u(:,1);
        end

        % V matrix over subjects (for averaging later)
        Vmat(:,d)=V;

        % remove trials without V
        rmtrial = isnan(V) | any(isnan(S.pred),2);
        V(rmtrial,:)=[];
        pred = S.pred;
        pred(rmtrial,:)=[];

        % Regression
        stats(tm).BRR(d) = bayesreg_crossval(pred,V,S,S.pred_group);
        disp(['model ' num2str(tm) ',  subject ' num2str(d) ', r2_brr = ' num2str(stats(tm).BRR(d).r2)])

        % concatentate
%         predcat = cat(1,predcat,pred);
%         Vcat = cat(1,Vcat,zscore(V));
    end
    
    % separate LME into groups
    grplist = S.designmat(2:end,strcmp(S.designmat(1,:),'groups'));
    grpuni = unique(grplist,'stable');
    LME=[stats(tm).BRR.waic];
    for g = 1:length(grpuni)
        LMEgrp{1,g}(tm,:) = LME(strcmp(grplist,grpuni{g}));
    end

    % put outputs into a matrix for plotting
    all_mu=cat(2,stats(tm).BRR(:).muB)';

    % average over subjects, remove nans, and regress
    Vavg = nanmean(Vmat,2);
    rmtrial = isnan(Vavg) | any(isnan(S.pred),2);
    Vavg(rmtrial,:)=[];
    predavg = S.pred;
    predavg(rmtrial,:)=[];
    stats(tm).BRRgrpavg = bayesreg_crossval(predavg,Vavg,S,S.pred_group);
    r2_brr_groupavg=stats(tm).BRRgrpavg.r2

    % identify params with greatest contribution (e.g. 75%)
    [sorted,si]=sort(abs(stats(tm).BRRgrpavg.muB),'descend');
    sortedcum = cumsum(sorted);
    numpara = find(sortedcum>.75*max(sortedcum),1,'first');
    si_max = sort(si(1:numpara));
    max_muB=stats(tm).BRRgrpavg.muB;
    max_muB(~ismember(stats(tm).BRRgrpavg.muB,stats(tm).BRRgrpavg.muB(si_max)))=0;

    % plot
    figure('name',['model ' num2str(tm) ' predictor weights'])
    subplot(1,3,1)
    hold on
    bar(mean(all_mu,1))
    errorbar(mean(all_mu,1),std(all_mu,[],1),'.')
    set(gca,'xtick',1:length(S.pred_label),'XTickLabel',S.pred_label,'XTickLabelRotation',90);
    hold off
    subplot(1,3,2)
    x=repmat(1:length(S.pred_label),size(all_mu,1),1);
    y=all_mu;
    scatter(x(:),y(:))
    set(gca,'xtick',1:length(S.pred_label),'XTickLabel',S.pred_label,'XTickLabelRotation',90);
    subplot(1,3,3)
    hold on
    bar(stats(tm).BRRgrpavg.muB,'b')
    bar(max_muB,'r')
    set(gca,'xtick',1:length(S.pred_label),'XTickLabel',S.pred_label,'XTickLabelRotation',90);
    hold off

    % plot timecourse of response
    xticks=1:size(Vmat,1);
    Vavg = nanmean(Vmat,2);
    Vstd = nanstd(Vmat,[],2);
    nsub = size(Vmat,2); 
    SEM = Vstd/sqrt(nsub);
    tscore = -tinv(0.025,nsub-1);
    CI = tscore*SEM;
    upper = Vavg+CI;
    lower = Vavg-CI;
    nonan = ~isnan(lower);% remove nan
    figure('name',['model ' num2str(tm) ' timecourse']); hold on
    fill([xticks(nonan), fliplr(xticks(nonan))], [upper(nonan)', fliplr(lower(nonan)')], ...
    'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    plot(xticks(nonan),Vavg(nonan),'b'); 
    title('timecourse of actual response')

    % scatter plots of predictors vs. predicted
    figure('name',['model ' num2str(tm) ' scatter pred'])
    spdim = ceil(sqrt(length(si_max)));
    for p = 1:length(si_max)
        subplot(spdim,spdim,p)
        scatter(predavg(:,si_max(p)),Vavg(~rmtrial))
        title(S.pred_label{si_max(p)})
        ylabel('predicted variable')
        lsline
    end

    % standard linear regression
    rmtrial = isnan(Vavg) | any(isnan(S.pred),2);
    Y=Vavg(~rmtrial);
    X=S.pred(~rmtrial,si_max);
    [b_all,~,~,~,st] = regress(Y, [ones(length(X),1) X]);
    r2_MR_groupavg=st(1)
    
    % plot against residuals from standard linear regression
    b=[];
    figure('name',['model ' num2str(tm) ' residuals'])
    spdim = ceil(sqrt(length(si_max)));
    for p = 1:length(si_max)
        not_V = si_max; not_V(p)=[];
        X=predavg(:,not_V);
        Y=Vavg(~rmtrial);
        b(:,p) = regress(Y, [ones(length(X),1) X]); 
        Yhat = [ones(length(X),1) X]*b(:,p);
        resid2 = Y - Yhat;

        subplot(spdim,spdim,p)
        scatter(predavg(:,si_max(p)),resid2)
        title(['marginal ' S.pred_label{si_max(p)}])
        ylabel('predicted variable')
        lsline
    end
    
    % % save
    % save(fullfile(S.path.stats(tm),['stats(tm)_models' num2str(testEEGmodels) '_' S.pred_transform '_zscore' num2str(S.zscore) '_' S.sname '.mat']),'stats(tm)','bmc','S');
end

% group model comparison
if length(S.trajmodelselect)>1
    [bmc.gposterior,bmc.gout] = VBA_groupBMC_btwGroups(LMEgrp)
end