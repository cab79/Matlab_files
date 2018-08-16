%% prediction of behaviour using regression
dbstop if error
clear all
close all

% y column: response to predict
% 0 = u col
% 1 = Choice
% 2 = RTs
% 3 = decoded EEG variables
ycol = 2;
%S.use_group_recons=0; % for decoded EEG only
S.mismatch_trials_only=1;
S.use_grpavg_weights=1; 
S.brr_var = 0.75; % 
S.y_predictor = 3; % additional ycol to use as predictor for comparison
S.bayesopt = 1; % use Bayesopt traj as predictors

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

% FOLDER AND FILENAME DEFINITIONS
S.path.main = 'C:\Data\CORE\eeg\ana';
S.path.stats = [S.path.main '\stats']; % folder to save outputs
S.path.hgf = 'C:\Data\CORE\behaviour\hgf'; 

% filenames
fnames = {
    %'CORE_fittedparameters_percmodel10_respmodel16_fractrain0_20180802T180446.mat'
    %'CORE_fittedparameters_percmodel10_respmodel22_fractrain0_20180806T171509.mat'
    %'CORE_fittedparameters_percmodel10_respmodel26_fractrain0_20180806T171509.mat'
    'CORE_fittedparameters_percmodel10_respmodel28_fractrain0_20180806T171509.mat'
    }
fname_bayesopt = 'CORE_fittedparameters_percmodel10_bayesopt_20180805T093044.mat';

% settings for HGF trajectory predictors
S.traj{1} = {
     {'PL'},{'mu','sa','muhat','sahat'},{[0 0 0 0]},{[],[],[],[]};
     }; 
S.traj{2} = {
     {'PL'},{'da','dau','psi','epsi'},{[0 0 0 0]},{[],[],[],[]};
     }; 
S.pred_transform = 'arcsinh'; % arcsinh, rank or notrans
S.zscore = 1;

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
S.brr.usegroups = 0;     % ****Specified by S.traj cells**** - create groups of predictors. Grouping of variables
%                       works only with HS, HS+ and lasso prior
%                       distributions. The same variable can appear in
%                       multiple groups. See examples\br_example[9,10,11,12]  (Default: { [] } )  
S.brr.waic = true;       %- whether to calculate the WAIC -- disabling can lead
%                       to large speed-ups, especially for Gaussian models with large n
%                       (default: true)

if S.bayesopt
    bs=load(fullfile(S.path.hgf,'fitted',fname_bayesopt));
    D_bayesopt=bs.D_fit;
end

for f=1:length(fnames)
    ls=load(fullfile(S.path.hgf,'fitted',fnames{f}));
    D_fit=ls.D_fit;
    
    Vmat=[];
    Vpmat=[];
    pred_all={};
    for d = 1:length(D_fit)
        
        %% predictor variables from S.traj
        
        % load and format HGF data
        % for each traj group... reformat into predictor matrix
        if S.bayesopt
            D_traj=D_bayesopt;
        else
            D_traj=D_fit(d);
        end
        [S] = HGF_traj2mat(S,D_traj);
        
        % remove nan columns from S.pred
        rmnan = any(isnan(S.pred),1);
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
        
        %% predicted variable V
        
        % V: variable of interest to predict
        MM=D_fit(d).HGF.u(:,1);
        if ycol
            V=D_fit(d).HGF.y(:,ycol);
            if S.mismatch_trials_only
                V=V(find(MM));
                S.pred=S.pred(find(MM),:);
            end
        else
            V=MM;
        end
        if S.y_predictor
            Vp=D_fit(d).HGF.y(:,S.y_predictor);
            if S.mismatch_trials_only
                Vp=Vp(find(MM));
            end
            Vpmat(:,d)=Vp;
        end
        
        % V matrix over subjects (for averaging later)
        Vmat(:,d)=V;
        % add predictors to cell array over subjects
        pred_all{d}=S.pred;

        % remove trials without V for this regression
        rmtrial = isnan(V);
        V(rmtrial,:)=[];
        S.pred(rmtrial,:)=[];
        
        %% Regression
        disp(['BRR: model ' num2str(f) ',  subject ' num2str(d)])
        stats(f).BRR(d) = bayesreg_crossval(S.pred,V,S,S.pred_group);
        

    end
    
    
    %% separate LME into groups
    grplist = ls.S.designmat(2:end,strcmp(ls.S.designmat(1,:),'groups'));
    grpuni = unique(grplist,'stable');
    LME=[stats(f).BRR.waic];
    for g = 1:length(grpuni)
        LMEgrp{1,g}(f,:) = LME(strcmp(grplist,grpuni{g}));
    end
    
    % put outputs into a matrix for plotting
    all_mu=cat(2,stats(f).BRR(:).muB)';

    % average over subjects, remove nans, and regress
    Vavg = nanmean(Vmat,2);
    rmtrial = isnan(Vavg);
    Vavg(rmtrial,:)=[];
    predavg = mean(cat(3,pred_all{:}),3);
    predavg(rmtrial,:)=[];
    stats(f).BRRgrpavg = bayesreg_crossval(predavg,Vavg,S,S.pred_group);
    r2_brr=stats(f).BRRgrpavg.r2

    % identify params with greatest contribution (e.g. 75%)
    [sorted,si]=sort(abs(stats(f).BRRgrpavg.muB),'descend');
    sortedcum = cumsum(sorted);
    numpara = find(sortedcum>S.brr_var*max(sortedcum),1,'first');
    si_max = sort(si(1:numpara));
    max_muB=stats(f).BRRgrpavg.muB;
    max_muB(~ismember(stats(f).BRRgrpavg.muB,stats(f).BRRgrpavg.muB(si_max)))=0;

    % plot
    figure
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
    bar(stats(f).BRRgrpavg.muB,'b')
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
    figure; hold on
    fill([xticks(nonan), fliplr(xticks(nonan))], [upper(nonan)', fliplr(lower(nonan)')], ...
    'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    plot(xticks(nonan),Vavg(nonan),'b'); 
    title('timecourse of actual response')

    % scatter plots of predictors vs. predicted
    figure
    spdim = ceil(sqrt(length(si_max)));
    for p = 1:length(si_max)
        subplot(spdim,spdim,p)
        scatter(predavg(:,si_max(p)),Vavg(~rmtrial))
        title(S.pred_label{si_max(p)})
        ylabel('predicted variable')
        lsline
    end
    
    % plot S.y_predictor vs. response
    if S.y_predictor
        Vpavg=nanmean(Vpmat,2);
        rmtrial = isnan(Vpavg) | isnan(Vavg);
        figure
        cc_avg(f) = corr(Vpavg(~rmtrial),Vavg(~rmtrial),'type','Spearman');
        scatter(Vpavg(~rmtrial),Vavg(~rmtrial))
        title(['y column ' num2str(S.y_predictor) ' predictor of mean response over subjects; rho = ' num2str(cc_avg)])
        ylabel('predicted variable')
        lsline
    end

    % standard linear regression on selected predictors
    X=predavg(:,si_max);
    rmtrial = isnan(Vavg);
    Y=Vavg(~rmtrial);
    [b_all,~,~,~,st] = regress(Y, [ones(length(X),1) X]);
    r2_MR=st(1)
    
    % plot against residuals from standard linear regression
    b=[];
    figure
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
    
    % loop though subjects and simulate V over trials from S.pred
    Vsim=nan(size(Vmat,1),length(pred_all));
    for d = 1:length(pred_all)
        
        % standard linear regression on selected predictors
        X=[ones(length(pred_all{d}),1) pred_all{d}(:,si_max)];
        
        if ~S.use_grpavg_weights
            
            % remove nan trials 
            Vactd=Vmat(:,d);
            rmtrial = isnan(Vactd);
            Vactd(rmtrial,:)=[];
            X(rmtrial,:)=[];
            
            b_all = regress(Vactd,X);
        else
            rmtrial=zeros(1,length(Vsim));
        end
        
        % simulate
        Vsim(~rmtrial,d) = b_all'*X';
        
    end
    
    % plot timecourse of predicted response
    xticks=1:size(Vsim,1);
    Vsimavg = nanmean(Vsim,2);
    Vsimstd = nanstd(Vsim,[],2);
    nsub = size(Vsim,2); 
    SEM = Vsimstd/sqrt(nsub);
    tscore = -tinv(0.025,nsub-1);
    CI = tscore*SEM;
    upper = Vsimavg+CI;
    lower = Vsimavg-CI;
    nonan = ~isnan(lower);% remove nan
    figure; hold on
    fill([xticks(nonan), fliplr(xticks(nonan))], [upper(nonan)', fliplr(lower(nonan)')], ...
    'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    plot(xticks(nonan),Vsimavg(nonan),'b'); 
    title('timecourse of predicted response')
    
    % scatter actual and predicted
    figure
    scatter(Vsimavg(~rmtrial),Vavg(~rmtrial))
    xlabel('predicted')
    ylabel('actual')
    
    % correlations for each subject
    for d = 1:length(pred_all)
        
        % remove nan trials 
        Vactd=Vmat(:,d);
        Vsimd=Vsim(:,d);
        rmtrial = isnan(Vactd);
        Vactd(rmtrial,:)=[];
        Vsimd(rmtrial,:)=[];

        %correlate
        cc(d,f) = corr(Vsimd,Vactd,'type','Spearman');
        
        % plot S.y_predictor vs. response
        if S.y_predictor
            Vp=Vpmat(:,d);
            rmtrial = isnan(Vp) | isnan(Vmat(:,d));
            cc_vp(d,f) = corr(Vp(~rmtrial),Vmat(~rmtrial,d),'type','Spearman');
        end
    end
    
    % separate predicted responses into conds
    if isfield(D_fit,'dt')
        condnum = [1 2 5 6 9 10 13 14 17 18 21 22]; % change (mismatch) conds only
        for d = 1:length(D_fit)
            Vsimd=Vsim(:,d);
            for cn = 1:length(condnum)
                cond_alltrials = D_fit(d).dt.design(2,:)==condnum(cn);
                cond_MMtrials = cond_alltrials(find(MM));
                cond_res(f,d,cn) = nanmean(Vsimd(cond_MMtrials));
                if S.y_predictor
                    cond_res_vp(f,d,cn) = nanmean(Vpmat(cond_MMtrials,d));
                end
            end
        end
    end
    
    % plot condition effects
    figure('Name',['condition effects: model ' num2str(f)]);
    hold on
    m=squeeze(nanmean(cond_res(f,:,:),2));
    sd=squeeze(nanstd(cond_res(f,:,:),[],2));
    bar(condnum,m);
    errorbar(condnum,m,sd,'.');
    
    if S.y_predictor
        % plot condition effects
        figure('Name',['condition effects for ycol predictor: model ' num2str(f)]);
        hold on
        m=squeeze(nanmean(cond_res_vp(f,:,:),2));
        sd=squeeze(nanstd(cond_res_vp(f,:,:),[],2));
        bar(condnum,m);
        errorbar(condnum,m,sd,'.');
    end
    
    % predict magnitude of DC effects across subjects for RTs only
    if ycol == 2;
        [num,txt,raw]=xlsread('C:\Data\CORE\behaviour\processed\condition effects\logreactiontime_aff_20180702T065408.xlsx')
        actRT = num(:,4:end);
        subind=ismember(txt(2:end,1),{D_fit.subname});
        actRTdiff=nanmean(actRT(subind,[1 3 5 7 9 11])-actRT(subind,[2 4 6 8 10 12]),2);
        figure('Name',['DC effect: model ' num2str(f)]);
        hold on
        m=squeeze(cond_res(f,:,:));
        simRTdiff=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2);
        scatter(simRTdiff,actRTdiff);
        rho=corr(simRTdiff,actRTdiff,'type','Spearman');
        title(['DC effect over subjects; rho = ' num2str(rho)]);
        if S.y_predictor
            figure('Name',['DC effect for ycol predictor: model ' num2str(f)]);
            hold on
            m=squeeze(cond_res_vp(f,:,:));
            simRTdiff=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2);
            scatter(simRTdiff,actRTdiff);
            rho=corr(simRTdiff,actRTdiff,'type','Spearman');
            title(['DC effect over subjects for ycol predictor; rho = ' num2str(rho)]);
        end
    end
    
end

% plot correlation coefficients
figure
subplot(1,2,1)
hold on
bar(1:f,nanmean(cc,1))
errorbar(1:f,nanmean(cc,1),nanstd(cc,[],1),'.')
hold off
subplot(1,2,2)
x=repmat(1:f,size(cc,1),1);
y=cc;
scatter(x(:),y(:))

if S.y_predictor
    figure
    subplot(1,2,1)
    hold on
    bar(1:f,nanmean(cc_vp,1))
    errorbar(1:f,nanmean(cc_vp,1),nanstd(cc_vp,[],1),'.')
    hold off
    subplot(1,2,2)
    x=repmat(1:f,size(cc_vp,1),1);
    y=cc_vp;
    scatter(x(:),y(:))
end

% group model comparison
if f>1
    [bmc.gposterior,bmc.gout] = VBA_groupBMC_btwGroups(LMEgrp)
end

