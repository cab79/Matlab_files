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
get_dt = 'CORE_fittedparameters_percmodel2_respmodel4_fractrain0_20190211T074650.mat';

% filename prefix and extension
% S.sname='0_20190211T074650'; 
% S.fname_pref= 'CORE_fittedparameters';
% S.fname_ext= ['_fractrain' S.sname];
S.sname=''; 
% S.fname_pref= 'D_fit'; test_frac=0;
S.fname_pref= 'D_fit_05'; test_frac=1;
S.fname_ext= '';
fontsize=14; % use 14 for singel plots, 10 for multiplots

% which models?
S.it=[0];
%S.perc_model_fnames=[1205 1210 1215 1220 1225 1225 1230 1235 1240 1245 1250];
S.perc_models=[3];
S.resp_models = [4];
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
            try
                if S.it(ii)~=0
                    fname=[S.fname_pref '_percmodel' num2str(S.perc_models(pm)) '_respmodel' num2str(S.resp_models(rm)) S.fname_ext num2str(S.it(ii)) '.mat'];
                else
                    fname=[S.fname_pref '_percmodel' num2str(S.perc_models(pm)) '_respmodel' num2str(S.resp_models(rm)) S.fname_ext '.mat'];
                end
                ls = load(fullfile(S.path.hgf,'fitted',fname));
            catch
                if S.it(ii)~=0
                    fname=[S.fname_pref '_pm' num2str(S.perc_models(pm)) '_rm' num2str(S.resp_models(rm)) S.fname_ext num2str(S.it(ii)) '.mat'];
                else
                    fname=[S.fname_pref '_pm' num2str(S.perc_models(pm)) '_rm' num2str(S.resp_models(rm)) S.fname_ext '.mat'];
                end
                ls = load(fullfile(S.path.hgf,'fitted',fname));
            end
            if ~isfield(ls,'S') && isfield(ls,'GS')
                ls.S=ls.GS;
            end
            D_fit=ls.D_fit;
            if ~isfield(D_fit,'dt')
                dtfile = load(fullfile(S.path.hgf,'fitted',get_dt));
                [D_fit(:).dt]=deal(dtfile.D_fit(:).dt);
                [D_fit(:).subname]=deal(dtfile.D_fit(:).subname);
            end
            tm=tm+1; testmodels=[testmodels tm];

            %S.perc_model = S.perc_models(pm); 
            S.resp_model = S.resp_models(rm); 
            S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
            S=CORE_perceptual_models(S);
            S=CORE_response_models(S);

            % split data into training and testing sets (if we want to test for prediction of behaviour)
            if ~isfield(ls.S,'train_idx')
                S.train_idx = []; % set to 0 to include all data in training set AND test set
            else
                S.train_idx = ls.S.train_idx; 
            end

            D_test=D_fit;
            if ~isempty(S.train_idx) && test_frac==1
                for d=1:length(D_fit)
%                     if ~isfield(S,'train_idx')
%                         % just taje a random sample
%                         cond = D_fit(d).dt.design(2,:); % get conditions
%                         ucond = unique(cond);
%                         % random indices from each cond
%                         S.train_idx = [];
%                         for u = 1:length(ucond)
%                             cond_idx = find(cond==ucond(u));
%                             S.train_idx = [S.train_idx randsample(cond_idx,S.frac_train*length(cond_idx))];
%                         end
%                     end
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
            if ~isempty(S.train_idx) && test_frac==1
                for d=1:length(D_sim)
                    D_sim(d).HGF.sim.y(S.train_idx,:) = nan; 
                end
            end
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
        S.colour_code = repmat([1:6],1,3);
        ucol = unique(S.colour_code);
        [cb] = cbrewer('qual', 'Paired', max(length(ucol),max(ucol)), 'pchip');
        violinplot(cc,1:tm,testmodels,'DotSize',36,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.8 0.8 0.8],'MedianColor',[1 1 1])
        ylabel('Spearman''s Rho'); xlabel('model')
        title('correlation between simulated and actual log(RT)')
        set(gca,'fontsize',fontsize)
        
        % plot WAIC: violin plot
        figure
        for tm = 1:length(brr)
            waic(:,tm) = [brr(tm).sub(:).waic];
        end
        violinplot(-waic,1:tm,testmodels)
        title('-WAIC')
        
        % model comparison
        %pep_flag=1;
        %[bmc.gposterior,bmc.gout] = VBA_groupBMC_cab(-waic',[],pep_flag);
        if tm>1
            S.fname_ext= ['_fractrain' S.sname '.mat'];
            S.designmat = ls.S.designmat;
            S.family_on=1;
            S.pep_flag = 1; % display PEP instead of EP
            [~,~,bmc.pm_out,bmc.rm_out]=HGF_group_model_comparison(S,-waic');
        end
        
        % plot condition effects, averaged over side
        marg_cond_rt = mean(cat(4,cond_rt(:,:,[1 2 5 6 9 10]),cond_rt(:,:,[3 4 7 8 11 12])),4);
        S.x_pos = [1:6];
%         h1=figure('Name','condition effects')
        h2=figure('Name','condition effects')
        plot_order = testmodels;%reshape(reshape(1:length(testmodels),length(S.resp_models),length(S.perc_models))',[],1);
        
        % paired colours
%         S.colour_code = [9 10 7 8 3 4];
%         ucol = unique(S.colour_code);
%         [cb] = cbrewer('qual', 'Paired', max(length(ucol),max(ucol)), 'pchip');
%         [cb2] = cbrewer('qual', 'Paired', 6, 'pchip');
        
        % scale colours
        S.colour_code = [1 2 3 4 5 6];
        ucol = unique(S.colour_code);
        cb([1 3 5],:) = [0.25 0.25 1; 0 0 0.75; 0 0 0.5]; % violets
        cb([2 4 6],:) = [1 0.25 0.25; 0.75 0 0; 0.5 0 0]; % reds
        
        for tm=1:length(testmodels)
            
            % bars
            %figure('Name',['condition effects: model ' num2str(testmodels(tm))]);
            %hold on
%             figure(h1)
%             subplot(length(S.perc_models),length(S.resp_models),plot_order(tm))
%             m=squeeze(mean(marg_cond_rt(tm,:,:),2));
%             sd=squeeze(std(marg_cond_rt(tm,:,:),[],2));
%             bar(S.x_pos,m);
%             hold on
%             errorbar(S.x_pos,m,sd,'.');
%             title(['model ' num2str(tm)])
            
            % violins
            figure(h2)
            subplot(length(S.perc_models),length(S.resp_models),plot_order(tm))
            plotdat = squeeze(marg_cond_rt(tm,:,:));
            if length(testmodels)>1
                xlab = [];
                violins = violinplot(plotdat,xlab,S.x_pos,'DotSize',10,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
                title(['model ' num2str(tm)])
            else
                xlab = {'CP10, CD1','CP10, CD3','CP30, CD1','CP30, CD3','CP50, CD1','CP50, CD3'};
                violins = violinplot(plotdat,xlab,S.x_pos,'DotSize',72,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
              
                title(['model ' num2str(tm)])
                set(gca,'fontsize',fontsize)
                xtickangle(45)
                ylabel('log(RT)')
                %title('response times')
            end
            lineplot = {[1 2],[3 4],[5 6]};
            for lp = 1:length(lineplot)
                XData=[];YData=[];
                for x=1:length(lineplot{lp})
                    plotind=lineplot{lp}(x);
                    ii=1:size(plotdat,1);
                    ni=find(isnan(plotdat(:,plotind)));
                    ii(ni)=[];
                    XData(x,ii)=violins(plotind).ScatterPlot.XData;
                    YData(x,ii)=violins(plotind).ScatterPlot.YData;
                    XData(x,ni)=mean(violins(plotind).ScatterPlot.XData);
                    YData(x,ni)=mean(violins(plotind).ScatterPlot.YData);
                end
                hold on; hb=plot(XData,YData,'Color',[0.8 0.8 0.8]); 
                uistack(hb,'bottom')
            end
            
        end
        
        % PREDICT ACTUAL BEHAVIOUR
        [numRT,txtRT,rawRT]=xlsread('C:\Data\CORE\behaviour\processed\condition effects\logreactiontime_aff_20180702T065408.xlsx')
        actRT = numRT(:,4:end);
        subind=ismember(txtRT(2:end,1),{D_test.subname});
        [numACC,txtACC,rawACC]=xlsread('C:\Data\CORE\behaviour\processed\condition effects\accuracy_aff_20180702T065408.xlsx')
        actACC = numACC(:,4:end);
        [numDIT,txtDIT,rawDIT]=xlsread('C:\Data\CORE\behaviour\processed\condition effects\DIT_RT.xlsx')
        actDITRT = numDIT(:,2:3);
        
        % violins of actual RT
        plotdat = mean(cat(3,actRT(:,[1 2 5 6 9 10]),actRT(:,[3 4 7 8 11 12])),3);
        xlab = {'CP10, CD1','CP10, CD3','CP30, CD1','CP30, CD3','CP50, CD1','CP50, CD3'};
        h3=figure('Name',['condition effects: RT']);
        %plot(plotdat','Color',[0.8 0.8 0.8]); hold on
        violins = violinplot(plotdat,xlab,S.x_pos,'DotSize',72,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
        lineplot = {[1 2],[3 4],[5 6]};
        for lp = 1:length(lineplot)
            XData=[];YData=[];
            for x=1:length(lineplot{lp})
                plotind=lineplot{lp}(x);
                ii=1:size(plotdat,1);
                ni=find(isnan(plotdat(:,plotind)));
                ii(ni)=[];
                XData(x,ii)=violins(plotind).ScatterPlot.XData;
                YData(x,ii)=violins(plotind).ScatterPlot.YData;
                XData(x,ni)=mean(violins(plotind).ScatterPlot.XData);
                YData(x,ni)=mean(violins(plotind).ScatterPlot.YData);
            end
            hold on; hb=plot(XData,YData,'Color',[0.8 0.8 0.8]); 
            uistack(hb,'bottom')
        end
        set(gca,'fontsize',fontsize)
        xtickangle(45)
        ylabel('log(RT)')
        %title('response times')
        
        % violins of actual Acc
        plotdat = mean(cat(3,actACC(:,[1 2 5 6 9 10]),actACC(:,[3 4 7 8 11 12])),3);
        h4=figure('Name',['condition effects: Acc']);
        violins = violinplot(plotdat,xlab,S.x_pos,'DotSize',72,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
        lineplot = {[1 2],[3 4],[5 6]};
        for lp = 1:length(lineplot)
            XData=[];YData=[];
            for x=1:length(lineplot{lp})
                plotind=lineplot{lp}(x);
                ii=1:size(plotdat,1);
                ni=find(isnan(plotdat(:,plotind)));
                ii(ni)=[];
                XData(x,ii)=violins(plotind).ScatterPlot.XData;
                YData(x,ii)=violins(plotind).ScatterPlot.YData;
                XData(x,ni)=mean(violins(plotind).ScatterPlot.XData);
                YData(x,ni)=mean(violins(plotind).ScatterPlot.YData);
            end
            hold on; hb=plot(XData,YData,'Color',[0.8 0.8 0.8]); 
            uistack(hb,'bottom') 
        end
        set(gca,'fontsize',fontsize)
        xtickangle(45)
        ylabel('% accuracy')
        %title('accuracy')
        
        % grp effect on actual RT
        grps = rawRT(2:end,strcmp(rawRT(1,:),'Group'));
        ugrps = unique(grps);
        plotdat=struct;
        plotdat.(ugrps{1}) = mean(actRT(ismember(grps,ugrps(1)),1:12),2);
        plotdat.(ugrps{2}) = mean(actRT(ismember(grps,ugrps(2)),1:12),2);
        xlab = grps;
        h3b=figure('Name',['group effects: RT']);
        violins = violinplot(plotdat,xlab,[1 2],'DotSize',72,'ViolinColor',cb([6 3],:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
        set(gca,'fontsize',fontsize)
%         xtickangle(45)
        ylabel('log(RT)')
        %title('response times')
        clear plotdat
        
        % grp * hand effect on actual RT
        plotdat=struct;
        plotdat.([ugrps{1} 'aff']) = mean(actRT(ismember(grps,ugrps(1)),[1 2 5 6 9 10]),2);
        plotdat.([ugrps{1} 'unaff']) = mean(actRT(ismember(grps,ugrps(1)),[3 4 7 8 11 12]),2);
        plotdat.([ugrps{2} 'aff']) = mean(actRT(ismember(grps,ugrps(2)),[1 2 5 6 9 10]),2);
        plotdat.([ugrps{2} 'unaff']) = mean(actRT(ismember(grps,ugrps(2)),[3 4 7 8 11 12]),2);
        xlab = {'CPRS-A','CRPS-U','HC-A','HC-U'};
        h3c=figure('Name',['group effects: RT']);
        violins = violinplot(plotdat,xlab,[1 2],'DotSize',72,'ViolinColor',cb([6 2 5 1],:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
        fn = fieldnames(plotdat);
        lineplot = {[1 2],[3 4]};
        for lp = 1:length(lineplot)
            XData=[];YData=[];
            for x=1:length(lineplot{lp})
                plotind=lineplot{lp}(x);
                ii=1:size(plotdat.(fn{plotind}),1);
                ni=find(isnan(plotdat.(fn{plotind})));
                ii(ni)=[];
                XData(x,ii)=violins(plotind).ScatterPlot.XData;
                YData(x,ii)=violins(plotind).ScatterPlot.YData;
                XData(x,ni)=mean(violins(plotind).ScatterPlot.XData);
                YData(x,ni)=mean(violins(plotind).ScatterPlot.YData);
            end
            hold on; hb=plot(XData,YData,'Color',[0.8 0.8 0.8]); 
            uistack(hb,'bottom')
        end
        set(gca,'fontsize',fontsize)
        ylabel('log(RT)')
%         xtickangle(45)
        %title('response times')
        clear plotdat
        
        
        % grp * hand effect on DIT RT
        plotdat=struct;
        plotdat.([ugrps{1} 'aff']) = actDITRT(ismember(grps,ugrps(1)),[1]);
        plotdat.([ugrps{1} 'unaff']) = actDITRT(ismember(grps,ugrps(1)),[2]);
        plotdat.([ugrps{2} 'aff']) = actDITRT(ismember(grps,ugrps(2)),[1]);
        plotdat.([ugrps{2} 'unaff']) = actDITRT(ismember(grps,ugrps(2)),[2]);
        xlab = {'CPRS-A','CRPS-U','HC-A','HC-U'};
        h3c=figure('Name',['group effects: DIT RT']);
        violins = violinplot(plotdat,xlab,[1 2],'DotSize',72,'ViolinColor',cb([6 2 5 1],:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
        fn = fieldnames(plotdat);
        lineplot = {[1 2],[3 4]};
        for lp = 1:length(lineplot)
            XData=[];YData=[];
            for x=1:length(lineplot{lp})
                plotind=lineplot{lp}(x);
                ii=1:size(plotdat.(fn{plotind}),1);
                ni=find(isnan(plotdat.(fn{plotind})));
                ii(ni)=[];
                XData(x,ii)=violins(plotind).ScatterPlot.XData;
                YData(x,ii)=violins(plotind).ScatterPlot.YData;
                XData(x,ni)=mean(violins(plotind).ScatterPlot.XData);
                YData(x,ni)=mean(violins(plotind).ScatterPlot.YData);
            end
            hold on; hb=plot(XData,YData,'Color',[0.8 0.8 0.8]); 
            uistack(hb,'bottom')
        end
        set(gca,'fontsize',fontsize)
%         xtickangle(45)
        ylabel('completion time (s)')
        %title('completion time')
        clear plotdat
        % grp * hand effect on DIT Acc
        plotdat=struct;
        plotdat.([ugrps{1} 'aff']) = actDITRT(ismember(grps,ugrps(1)),[1]);
        plotdat.([ugrps{1} 'unaff']) = actDITRT(ismember(grps,ugrps(1)),[2]);
        plotdat.([ugrps{2} 'aff']) = actDITRT(ismember(grps,ugrps(2)),[1]);
        plotdat.([ugrps{2} 'unaff']) = actDITRT(ismember(grps,ugrps(2)),[2]);
        xlab = {'CPRS-A','CRPS-U','HC-A','HC-U'};
        h3c=figure('Name',['group effects: DIT RT']);
        violins = violinplot(plotdat,xlab,[1 2],'DotSize',72,'ViolinColor',cb([6 2 5 1],:),'ViolinAlpha',0.6,'EdgeColor',[0.5 0.5 0.5],'MedianColor',[1 1 1])
        fn = fieldnames(plotdat);
        lineplot = {[1 2],[3 4]};
        for lp = 1:length(lineplot)
            XData=[];YData=[];
            for x=1:length(lineplot{lp})
                plotind=lineplot{lp}(x);
                ii=1:size(plotdat.(fn{plotind}),1);
                ni=find(isnan(plotdat.(fn{plotind})));
                ii(ni)=[];
                XData(x,ii)=violins(plotind).ScatterPlot.XData;
                YData(x,ii)=violins(plotind).ScatterPlot.YData;
                XData(x,ni)=mean(violins(plotind).ScatterPlot.XData);
                YData(x,ni)=mean(violins(plotind).ScatterPlot.YData);
            end
            hold on; hb=plot(XData,YData,'Color',[0.8 0.8 0.8]); 
            uistack(hb,'bottom')
        end
        set(gca,'fontsize',fontsize)
        ylabel('completion time (s)')
%         xtickangle(45)
        %title('completion time')
        clear plotdat
        
        h5=figure('Name','CD effect')
        % predict magnitude of CD effects across subjects
        actRTdiff=nanmean(actRT(subind,[1 3 5 7 9 11])-actRT(subind,[2 4 6 8 10 12]),2);
        for tm=1:length(testmodels)
            figure(h5)
            subplot(length(S.perc_models),length(S.resp_models),plot_order(tm))
            m=squeeze(cond_rt(tm,:,:));
            simRTdiff=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2);
            rho=corr(simRTdiff,actRTdiff,'type','Spearman');
            scatter(simRTdiff,actRTdiff,50,(1-rho)*ones(1,3),'filled');
            rho=num2str(rho);rho=rho(1:4);
            if length(testmodels)==1
                set(gca,'fontsize',fontsize*2)
            else
                set(gca,'fontsize',fontsize)
                title(['model ' num2str(tm) ', rho = ' rho])
            end
        end
        
        h6=figure('Name','CP effect')
        % predict CP effects (10 vs. 50) across subjects
        actRTdiff=(nanmean(actRT(subind,[1:4]),2)-nanmean(actRT(subind,[9:12]),2));
        for tm=1:length(testmodels)
            figure(h6)
            subplot(length(S.perc_models),length(S.resp_models),plot_order(tm))
            m=squeeze(cond_rt(tm,:,:));
            simRTdiff=nanmean(m(:,[1:4]),2)-nanmean(m(:,[9:12]),2);
            rho=corr(simRTdiff,actRTdiff,'type','Spearman');
            scatter(simRTdiff,actRTdiff,50,(1-abs(rho))*ones(1,3),'filled');
            rho=num2str(rho);rho=rho(1:4);
            if length(testmodels)==1
                set(gca,'fontsize',fontsize*2)
            else
                set(gca,'fontsize',fontsize)
                title(['model ' num2str(tm) ', rho = ' rho])
            end
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
        
        % predict actual CD effects from model parameters
%         actRTdiff=nanmean(actRT(subind,[1 3 5 7 9 11])-actRT(subind,[2 4 6 8 10 12]),2);
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['CD effect: param ' params(1).name{1,pi}]);
%             for tm=1:length(testmodels)
%                 prm = params(tm).val(subind,pi);
%                 [rho(1,tm),pval(1,tm)]=corr(prm,actRTdiff,'type','Spearman');
%             end
%             bar(testmodels,rho);
%             title(['rho = ' num2str(rho) ', pval = ' num2str(pval)]);
%         end
        
        % predict predicted CD effects from model parameters
%         simRTdiff=nanmean(m(:,[1 3 5 7 9 11])-m(:,[2 4 6 8 10 12]),2);
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['CD effect, sim: param ' params(1).name{1,pi}]);
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
        
        % predict actual CD effects on ACCURACY from model parameters
%         actACCdiff=nanmean(actACC(subind,[1 3 5 7 9 11])-actACC(subind,[2 4 6 8 10 12]),2);
%         for pi = 1:size(params(1).val,2)
%             figure('Name',['CD effect on Acc: param ' params(1).name{1,pi}]);
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

