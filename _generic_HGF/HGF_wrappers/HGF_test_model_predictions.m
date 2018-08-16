function varargout=HGF_test_model_predictions(D_act,D_sim,S) 
% outputs summary stats of model predictions of behaviour over trials

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
            end
            cc_mean(d,1)=mean(cc);
        end
        varargout = {cc_mean};
        
        
end