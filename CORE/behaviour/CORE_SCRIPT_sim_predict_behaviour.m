%% Analysis: Perceptual model simulation and prediction of behaviour
dbstop if error
clear all
close all
% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
addpath('C:\Data\Matlab\Violinplot-Matlab-master');
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

% FOLDER AND FILENAME DEFINITIONS
S.path.hgf = 'C:\Data\CORE\behaviour\hgf'; 

% filename prefix and extension
%fname_pref = 'CORE_fittedparameters_percmodel10_respmodel';
%fname_ext = '_fractrain0.5_20180818T070202.mat';

%sname='1_20190114T061324_it'; 
S.sname='0_20190112T081048'; 
%sname='0.5_20180818T070202'; % 50% training
S.fname_pref= 'CORE_fittedparameters';
S.fname_ext= ['_fractrain' S.sname];

%     fname_pref = 'CORE_fittedparameters_percmodel';
%     fname_ext = '_respmodel10_20180719T131912.mat';

% which models?
S.it=[0];
%S.perc_model_fnames=[1205 1210 1215 1220 1225 1225 1230 1235 1240 1245 1250];
S.perc_models=[1:3];
S.resp_models = [1:7];
% S.perc_models=[1 3 9 10 11 12];
% S.resp_models = [2];
% S.perc_models=[3];
% S.resp_models = [2 21];


tm=0;testmodels=[];
for ii = 1:length(S.it)
    for pm=1:length(S.perc_models)
    %for pm=1:length(S.perc_models)

             S.perc_model = S.perc_models(pm);
        %     load(fullfile(S.path.hgf,'fitted',[fname_pref num2str(S.perc_model) fname_ext])); 

        for rm=1:length(S.resp_models)
            %ls=load(fullfile(S.path.hgf,'fitted',[fname_pref num2str(testmodels(rm)) fname_ext]));
            if S.it(ii)~=0
                fname=[S.fname_pref '_percmodel' num2str(S.perc_models(pm)) '_respmodel' num2str(S.resp_models(rm)) S.fname_ext num2str(S.it(ii)) '.mat'];
            else
                fname=[S.fname_pref '_percmodel' num2str(S.perc_models(pm)) '_respmodel' num2str(S.resp_models(rm)) S.fname_ext '.mat'];
            end
            ls = load(fullfile(S.path.hgf,'fitted',fname));
            D_fit=ls.D_fit;
            tm=tm+1; testmodels=[testmodels tm];

            %S.perc_model = S.perc_models(pm); 
            S.resp_model = S.resp_models(rm); 
            S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
            S=CORE_perceptual_models(S);
            S=CORE_response_models(S);

            % split data into training and testing sets (if we want to test for prediction of behaviour)
            if ~isfield(S,'frac_train')
                S.frac_train = 0; % set to 0 to include all data in training set AND test set
            end

            D_test=D_fit;
            if S.frac_train>0
                for d=1:length(D_fit)
                    if ~isfield(S,'train_idx')
                        % just taje a random sample
                        cond = D_prep(d).dt.design(2,:); % get conditions
                        ucond = unique(cond);
                        % random indices from each cond
                        S.train_idx = [];
                        for u = 1:length(ucond)
                            cond_idx = find(cond==ucond(u));
                            S.train_idx = [S.train_idx randsample(cond_idx,S.frac_train*length(cond_idx))];
                        end
                    end
                    %D_test(d).HGF.u = D_prep(d).HGF.u(sort(S.test_idx),:);
                    %D_test(d).HGF.y = D_prep(d).HGF.y(sort(S.test_idx),:);
                    D_test(d).HGF.y(S.train_idx) = nan;
                end
            end

            % add parameters to test data
        %     if S.frac_train>0
        %         for d=1:length(D_fit)
        %             D_test(d).HGF.fit.p_prc = D_fit(d).HGF.fit.p_prc;
        %         end
        %     else
        %         D_test=D_fit;
        %     end

            % Simulations
            %S.resp_model = 12;S=CORE_response_models(S);
            S.numsimrep = 1; % number of simulations to run per parameter combination
            S.sim=[];
            [D_sim,S] = HGF_sim(D_test,S); 
            switch S.resp_modelspec.responses{:}
                case 'Ch'
                    [means_fitted(:,tm),stds_fitted(:,tm),means_randcorr_fitted(:,tm),stds_randcorr_fitted(:,tm)]=HGF_test_model_predictions(D_test,D_sim,S);
                    %save(fullfile(S.path.hgf,'fitted',['CORE_fittedpredictions_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'means_fitted','stds_fitted');
                case 'RT'
                    [cc(:,tm),brr(tm).sub]=HGF_test_model_predictions(D_test,D_sim,S);
                    %save(fullfile(S.path.hgf,'fitted',['CORE_fittedpredictions_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'cc');

                    % separate into conds
                    if isfield(D_sim,'dt')
                        condnum = [1 2 5 6 9 10 13 14 17 18 21 22]; % change (mismatch) conds only
                        for d = 1:length(D_sim)
                            for rep = 1:S.numsimrep
                                fitted_rt = D_sim(d).HGF(rep).sim.y(:,2);
                                for cn = 1:length(condnum)
                                    cond_rt(tm,d,cn,rep) = nanmean(fitted_rt(D_sim(1).dt.design(2,:)==condnum(cn)));
                                end
                            end
                        end
                        cond_rt=nanmean(cond_rt,4);
                    end

                    % extract model parameters
    %                 for d = 1:length(D_sim)
    %                     al = D_sim(d).HGF(1).fit.p_prc.like_al0; 
    %                     al = mean(al([2,4])) - mean(al([1,3])); % only work for 4 alpha model!
    %                     params(tm).val(d,1:length(al)) = al;
    %                     params(tm).name(d,1:length(al)) = deal({'al'});
    %                     om = D_sim(d).HGF(1).fit.p_prc.PL_om(2:3);
    %                     params(tm).val(d,length(al)+1:length(al)+length(om)) = om;
    %                     params(tm).name(d,length(al)+1:length(al)+length(om)) = deal({'om'});
    %                 end

        %         case 'EEG'
        %             [cc(:,tm)]=HGF_test_model_predictions(D_test,D_sim,S);
            end

        end
    end
end


%S.resp_model = 12; S=CORE_response_models(S);
switch S.resp_modelspec.responses{:}
    case 'Ch'
        % plot fraction predicted correctly
        figure
        hold on
        bar(testmodels,mean(means_fitted,1))
        errorbar(testmodels,mean(means_fitted,1),std(means_fitted,[],1),'.')
        %plot(xlim,mean(means_bayesopt)*[1 1], 'k--')

        % plot random corrected predictions
        figure
        hold on
        bar(testmodels,mean(means_randcorr_fitted,1))
        errorbar(testmodels,mean(means_randcorr_fitted,1),std(means_randcorr_fitted,[],1),'.')
    case 'RT'
        
        % plot correlation coefficients: bar chart
        figure
        hold on
        bar(testmodels,nanmean(cc,1))
        errorbar(testmodels,nanmean(cc,1),nanstd(cc,[],1),'.')
        title('Spearmans')
        
        % plot correlation coefficients: violin plot
        figure
        %x=repmat(testmodels,size(cc,1),1);
        %y=cc;
        %scatter(x(:),y(:))
        violinplot(cc,1:tm,testmodels)
        title('Spearmans')
        
        % plot WAIC: violin plot
        figure
        for tm = 1:length(brr)
            waic(:,tm) = [brr(tm).sub(:).waic];
        end
        violinplot(waic,1:tm,testmodels)
        title('WAIC')
        
        % model comparison
        %pep_flag=1;
        %[bmc.gposterior,bmc.gout] = VBA_groupBMC_cab(-waic',[],pep_flag);
        if tm>1
            S.fname_ext= ['_fractrain' S.sname '.mat'];
            S.designmat = ls.S.designmat;
            S.pep_flag = 1; % display PEP instead of EP
            [~,~,bmc.pm_out,bmc.rm_out]=HGF_group_model_comparison(S,-waic');
        end
        
        % plot condition effects, averaged over side
        marg_cond_rt = mean(cat(4,cond_rt(:,:,[1 2 5 6 9 10]),cond_rt(:,:,[3 4 7 8 11 12])),4);
        S.x_pos = [1:6];
        for tm=1:length(testmodels)
            
            % bars
            figure('Name',['condition effects: model ' num2str(testmodels(tm))]);
            hold on
            
            m=squeeze(mean(marg_cond_rt(tm,:,:),2));
            sd=squeeze(std(marg_cond_rt(tm,:,:),[],2));
            bar(S.x_pos,m);
            errorbar(S.x_pos,m,sd,'.');
            
            % violins
            S.colour_code = repmat([1 2],1,3);
            ucol = unique(S.colour_code);
            [cb] = cbrewer('qual', 'Set1', length(ucol), 'pchip');
            plotdat = squeeze(marg_cond_rt(tm,:,:));
            xlab = [];
            figure('Name',['condition effects: model ' num2str(testmodels(tm))]);
            violinplot(plotdat,xlab,S.x_pos,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6)
            title(['model ' num2str(tm)])
            
        end
        
        % PREDICT ACTUAL BEHAVIOUR
        [num,txt,raw]=xlsread('C:\Data\CORE\behaviour\processed\condition effects\logreactiontime_aff_20180702T065408.xlsx')
        actRT = num(:,4:end);
        subind=ismember(txt(2:end,1),{D_test.subname});
        [num,txt,raw]=xlsread('C:\Data\CORE\behaviour\processed\condition effects\accuracy_aff_20180702T065408.xlsx')
        actACC = num(:,4:end);
        
        % violins of actual RT
        plotdat = mean(cat(3,actRT(:,[1 2 5 6 9 10]),actRT(:,[3 4 7 8 11 12])),3);
        xlab = [];
        figure('Name',['condition effects: RT']);
        violinplot(plotdat,xlab,S.x_pos,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.8 0.8 0.8],'MedianColor',[0 0 0])
        
        % violins of actual Acc
        plotdat = mean(cat(3,actACC(:,[1 2 5 6 9 10]),actACC(:,[3 4 7 8 11 12])),3);
        xlab = [];
        figure('Name',['condition effects: Acc']);
        violinplot(plotdat,xlab,S.x_pos,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.8 0.8 0.8],'MedianColor',[0 0 0])
        
        % predict magnitude of DC effects across subjects
        actRTdiff=nanmean(actRT(subind,[1 3 5 7 9 11])-actRT(subind,[2 4 6 8 10 12]),2);
        for tm=1:length(testmodels)
            figure('Name',['DC abs effect: model ' num2str(testmodels(tm))]);
            hold on
            m=squeeze(cond_rt(tm,:,:));
            simRTdiff=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2);
            scatter(simRTdiff,actRTdiff);
            rho=corr(simRTdiff,actRTdiff,'type','Spearman');
            title(['rho = ' num2str(rho)]);
        end
        
        % predict z-scores of DC effects across subjects (check for
        % consistency with actual magnitudes, to validate use of z-scores
        % for CP condition effects below).
%         actRTz=nanmean(actRT(subind,[1 3 5 7 9 11])-actRT(subind,[2 4 6 8 10 12]),2)./nanstd(actRT(subind,1:12),[],2);
%         for tm=1:length(testmodels)
%             figure('Name',['DC z effect: model ' num2str(testmodels(tm))]);
%             hold on
%             m=squeeze(cond_rt(tm,:,:));
%             simRTz=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2)./nanstd(m(:,1:12),[],2);
%             scatter(simRTz,actRTz);
%             rho=corr(simRTz,actRTz,'type','Spearman');
%             title(['rho = ' num2str(rho)]);
%         end
        
        % predict CP effects (10 vs. 50) across subjects
        actRTdiff=(nanmean(actRT(subind,[1:4]),2)-nanmean(actRT(subind,[9:12]),2));
        for tm=1:length(testmodels)
            figure('Name',['CP 10v50 effect: model ' num2str(testmodels(tm))]);
            hold on
            m=squeeze(cond_rt(tm,:,:));
            simRTdiff=nanmean(m(:,[1:4]),2)-nanmean(m(:,[9:12]),2);
            scatter(simRTdiff,actRTdiff);
            rho=corr(simRTdiff,actRTdiff,'type','Spearman');
            title(['rho = ' num2str(rho)]);
        end
        
        % predict CP effects (30 vs. 50) across subjects
%         actRTdiff=(nanmean(actRT(subind,[5:8]),2)-nanmean(actRT(subind,[9:12]),2));
%         for tm=1:length(testmodels)
%             figure('Name',['CP 30v50 effect: model ' num2str(testmodels(tm))]);
%             hold on
%             m=squeeze(cond_rt(tm,:,:));
%             simRTdiff=nanmean(m(:,[5:8]),2)-nanmean(m(:,[9:12]),2);
%             scatter(simRTdiff,actRTdiff);
%             rho=corr(simRTdiff,actRTdiff,'type','Spearman');
%             title(['rho = ' num2str(rho)]);
%         end
        
        % predict actual DC effects from model parameters
%         actRTdiff=nanmean(actRT(subind,[1 3 5 7 9 11])-actRT(subind,[2 4 6 8 10 12]),2);
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['DC effect: param ' params(1).name{1,pi}]);
%             for tm=1:length(testmodels)
%                 prm = params(tm).val(subind,pi);
%                 [rho(1,tm),pval(1,tm)]=corr(prm,actRTdiff,'type','Spearman');
%             end
%             bar(testmodels,rho);
%             title(['rho = ' num2str(rho) ', pval = ' num2str(pval)]);
%         end
        
        % predict predicted DC effects from model parameters
%         simRTdiff=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2);
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['DC effect, sim: param ' params(1).name{1,pi}]);
%             for tm=1:length(testmodels)
%                 prm = params(tm).val(subind,pi);
%                 [rho(1,tm),pval(1,tm)]=corr(prm,simRTdiff,'type','Spearman');
%             end
%             bar(testmodels,rho);
%             title(['rho = ' num2str(rho) ', pval = ' num2str(pval)]);
%         end
        
        % predict actual CP effects from model parameters
%         actRTdiff=(nanmean(actRT(subind,[1:4]),2)-nanmean(actRT(subind,[9:12]),2));
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['CP 10v50 effect: param ' params(1).name{1,pi}]);
%             for tm=1:length(testmodels)
%                 prm = params(tm).val(subind,pi);
%                 [rho(1,tm),pval(1,tm)]=corr(prm,actRTdiff,'type','Spearman');
%             end
%             bar(testmodels,rho);
%             title(['rho = ' num2str(rho) ', pval = ' num2str(pval)]);
%         end
        
        % predict actual DC effects on ACCURACY from model parameters
%         actACCdiff=nanmean(actACC(subind,[1 3 5 7 9 11])-actACC(subind,[2 4 6 8 10 12]),2);
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['DC effect on Acc: param ' params(1).name{1,pi}]);
%             for tm=1:length(testmodels)
%                 prm = params(tm).val(subind,pi);
%                 [rho(1,tm),pval(1,tm)]=corr(prm,actACCdiff,'type','Spearman');
%             end
%             bar(testmodels,rho);
%             title(['rho = ' num2str(rho) ', pval = ' num2str(pval)]);
%         end
        
        % predict actual CP effects on ACCURACY from model parameters
%         actACCdiff=(nanmean(actACC(subind,[1:4]),2)-nanmean(actACC(subind,[9:12]),2));
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['CP 10v50 effect on Acc: param ' params(1).name{1,pi}]);
%             for tm=1:length(testmodels)
%                 prm = params(tm).val(subind,pi);
%                 [rho(1,tm),pval(1,tm)]=corr(prm,actACCdiff,'type','Spearman');
%             end
%             bar(testmodels,rho);
%             title(['rho = ' num2str(rho) ', pval = ' num2str(pval)]);
%         end
        
end

