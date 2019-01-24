function varargout=HGF_test_model_predictions(D_act,D_sim,S) 
% outputs summary stats of model predictions of behaviour over trials

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
S.zscore = 1;

switch S.resp_modelspec.responses{:}
    case 'Ch' % choices
        for d = 1:length(D_act)
            actual_choices = D_act(d).HGF.y(:,1);
            for rep = 1:S.numsimrep
                fitted_choices = D_sim(d).HGF(rep).sim.y(:,1);
                fitted_correct(rep) = sum(fitted_choices == actual_choices)/length(actual_choices);

                % to create benchmark for predictions, randomise simulated
                % responses over trials separately for each stimulus type (0 and 1)
                u0i = find(D_sim(d).HGF(rep).sim.u(:,1)==0);
                u1i = find(D_sim(d).HGF(rep).sim.u(:,1)==1);
                rand_choices = nan(length(fitted_choices),1);
                rand_choices(u0i) = D_sim(d).HGF(rep).sim.y(u0i(randperm(length(u0i))),1);
                rand_choices(u1i) = D_sim(d).HGF(rep).sim.y(u1i(randperm(length(u1i))),1);
                rand_correct(rep) = sum(rand_choices == actual_choices)/length(actual_choices);
            end
            means_fitted(d,1)=mean(fitted_correct);
            stds_fitted(d,1)=std(fitted_correct);
            % random-corrected
            means_randcorr_fitted(d,1)=mean(fitted_correct-rand_correct);
            stds_randcorr_fitted(d,1)=std(fitted_correct-rand_correct);
        end
        varargout = {means_fitted,stds_fitted,means_randcorr_fitted,stds_randcorr_fitted};

    case 'RT' % response times
        for d = 1:length(D_act)
            actual_rt = D_act(d).HGF.y(:,2);
            for rep = 1:S.numsimrep
                fitted_rt = D_sim(d).HGF(rep).sim.y(:,2);
                
                %only consider trials in which both fitted and actual
                %responses both occur
                ii = find(~isnan(actual_rt) .* ~isnan(fitted_rt));
                
                %correlate
                cc(rep) = corr(fitted_rt(ii),actual_rt(ii),'type','Spearman');
                
                % bayes reg
                brr(d,rep) = bayesreg_crossval(fitted_rt(ii),actual_rt(ii),S,0);
                
            end
            cc_mean(d,1)=mean(cc);
        end
        varargout = {cc_mean,brr};
        
        
end