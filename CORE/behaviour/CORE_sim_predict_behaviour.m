%% Analysis: Perceptual model simulation and prediction of behaviour
function CORE_sim_predict_behaviour(varargin)
% use: CORE_sim_predict_behaviour(S,D_fit)
% or: run as script

if isempty(varargin)
    dbstop if error
    clear all
    close all
    clear S
    % add toolbox paths
    run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
    
    % FOLDER AND FILENAME DEFINITIONS
    S.path.hgf = ['C:\Data\CORE\behaviour\hgf']; 
    
    % flilename prefix and extension
    
    % 100%
    % CORE_fittedparameters_percmodel6_respmodel10_20180720T195113.mat
    % CORE_fittedparameters_percmodel3_respmodel11_20180722T095006.mat
    % CORE_fittedpredictions_percmodel2_respmodel12_20180722T095040.mat
    
    % 50%
    % CORE_fittedparameters_percmodel8_respmodel10_20180719T131912.mat
    % CORE_fittedparameters_percmodel8_respmodel11_20180726T074654.mat
    % CORE_fittedparameters_percmodel8_respmodel12_20180726T220524.mat
    
%     fname_pref = 'CORE_fittedparameters_percmodel10_respmodel';
%     fname_ext = '_fractrain0_20180731T224814.mat';
    fname_pref = 'CORE_fittedparameters_percmodel';
    fname_ext = '_respmodel10_20180719T131912.mat';
    % which models?
    S.resp_model=[10];
    testmodels = [1:8];
%     S.resp_model=[10];
%     testmodels = [12:17];
else
    S=varargin{1};
    D_fit = varargin{2};
end

for pm=1:length(testmodels)
    S.perc_model = testmodels(pm);
    load(fullfile(S.path.hgf,'fitted',[fname_pref num2str(S.perc_model) fname_ext])); 
% for pm=1:length(testmodels)
%     S.resp_model = testmodels(pm); 
%     load(fullfile(S.path.hgf,'fitted',[fname_pref num2str(S.resp_model) fname_ext]));
    
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
    if S.frac_train>0
        for d=1:length(D_fit)
            D_test(d).HGF.fit.p_prc = D_fit(d).HGF.fit.p_prc;
        end
    else
        D_test=D_fit;
    end
    
    % Simulations
    %S.resp_model = 12;S=CORE_response_models(S);
    S.numsimrep = 100; % number of simulations to run per parameter combination
    S.sim=[];
    [D_sim,S] = HGF_sim(D_test,S); 
    switch S.resp_modelspec.responses{:}
        case 'Ch'
            [means_fitted(:,pm),stds_fitted(:,pm),means_randcorr_fitted(:,pm),stds_randcorr_fitted(:,pm)]=HGF_test_model_predictions(D_test,D_sim,S);
            %save(fullfile(S.path.hgf,'fitted',['CORE_fittedpredictions_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'means_fitted','stds_fitted');
        case 'RT'
            [cc(:,pm)]=HGF_test_model_predictions(D_test,D_sim,S);
            %save(fullfile(S.path.hgf,'fitted',['CORE_fittedpredictions_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'cc');
            
            % separate into conds
            if isfield(D_sim,'dt')
                condnum = [1 2 5 6 9 10 13 14 17 18 21 22]; % change (mismatch) conds only
                for d = 1:length(D_sim)
                    for rep = 1:S.numsimrep
                        fitted_rt = D_sim(d).HGF(rep).sim.y(:,2);
                        for cn = 1:length(condnum)
                            cond_rt(pm,d,cn,rep) = mean(fitted_rt(D_sim(1).dt.design(2,:)==condnum(cn)));
                        end
                    end
                end
                cond_rt=mean(cond_rt,4);
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
        % plot correlation coefficients
        figure
        hold on
        bar(testmodels,mean(cc,1))
        errorbar(testmodels,mean(cc,1),std(cc,[],1),'.')
        % plot condition effects
        figure
        x=repmat(testmodels,size(cc,1),1);
        y=cc;
        scatter(x(:),y(:))
        for pm=1:length(testmodels)
            figure('Name',['condition effects: model ' num2str(testmodels(pm))]);
            hold on
            m=squeeze(mean(cond_rt(pm,:,:),2));
            sd=squeeze(std(cond_rt(pm,:,:),[],2));
            bar(condnum,m);
            errorbar(condnum,m,sd,'.');
        end
        % predict magnitude of DC effects across subjects
        [num,txt,raw]=xlsread('C:\Data\CORE\behaviour\processed\condition effects\logreactiontime_aff_20180702T065408.xlsx')
        actRT = num(:,4:end);
        subind=ismember(txt(2:end,1),{D_test.subname});
        actRTdiff=nanmean(actRT(subind,[1 3 5 7 9 11])-actRT(subind,[2 4 6 8 10 12]),2);
        for pm=1:length(testmodels)
            figure('Name',['DC effect: model ' num2str(testmodels(pm))]);
            hold on
            m=squeeze(cond_rt(pm,:,:));
            simRTdiff=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2);
            scatter(simRTdiff,actRTdiff);
            rho=corr(simRTdiff,actRTdiff,'type','Spearman');
            title(['rho = ' num2str(rho)]);
        end
end

