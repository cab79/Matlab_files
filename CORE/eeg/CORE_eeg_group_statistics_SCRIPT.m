clear all
close all
dbstop if error
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
spath = 'C:\Data\CORE\eeg\ana\stats';
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

sfiles = {
    
   
%     'stats_BRR_all_chan_HGF_arcsinh_20190115T215155.mat' % t dist, Dau, real, g prior, CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it14.mat
%     'stats_BRR_all_chan_HGF_arcsinh_20190117T075855.mat' % Da1
%     'stats_BRR_all_chan_HGF_arcsinh_20190116T145135.mat' % Da2
%     'stats_BRR_all_chan_HGF_arcsinh_20190115T215243.mat' % t dist, Epsi1, real, g prior, CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it14.mat
%     'stats_BRR_all_chan_HGF_arcsinh_20190116T145202.mat' % epsi2
%     'stats_BRR_all_chan_HGF_arcsinh_20190117T075334.mat' % epsi3
%      'stats_BRR_all_chan_HGF_arcsinh_20190116T231951.mat' % t dist, Dau/Da1-2, ridge
%     'stats_BRR_all_chan_HGF_arcsinh_20190116T232149.mat' % t dist, epsi1-3, ridge
%     'stats_BRR_all_chan_condHGF_arcsinh_20190117T180646.mat' % t dist, Dau/Da1-2, ridge
%     'stats_BRR_all_chan_condHGF_arcsinh_20190117T154534.mat' % t dist, epsi1-3, ridge
%     'stats_BRR_all_chan_HGF_arcsinh_20190117T234711.mat' % t dist, Dau/Da1-2, g
%     'stats_BRR_all_chan_HGF_arcsinh_20190118T042700.mat' % t dist, epsi1-3, g

% 'stats_BRR_all_chan_condHGF_arcsinh_20190118T145255.mat' % t dist, Dau/Da1-2, hs+
% 'stats_BRR_all_chan_HGF_arcsinh_20190118T092602.mat' % t dist, Dau/Da1-2, hs+
%'stats_BRR_all_chan_HGF_arcsinh_20190118T185327.mat' % t dist, Mu1 Dau/Da1-2, hs+

% 'stats_BRR_all_chan_HGF_arcsinh_20190119T093237.mat' % t dist, Mu, ridge
% 'stats_BRR_all_chan_HGF_arcsinh_20190119T093021.mat' % t dist, Ud, ridge

% S.transform
% S.pred_transform  
% S.traj{1, 1}{1, 2}
% S.brr.model

% '' % no, no, lp, dau/da
% 'stats_BRR_all_chan_HGF_notrans_20190126T101431.mat' % no, no, lp, epsi
% '' % no, no, t, dau/da
% 'stats_BRR_all_chan_HGF_notrans_20190126T080704.mat' % no, no, t, epsi
% '' % no, arc, lp, dau/da
% 'stats_BRR_all_chan_HGF_notrans_20190126T140224.mat' % no, arc, lp, epsi
% '' % no, arc, t, dau/da
% 'stats_BRR_all_chan_HGF_notrans_20190126T120000.mat' % no, arc, t, epsi *** BEST IN TERMS OF MODEL EVIDENCE
% '' % arc, no, lp, dau/da
% 'stats_BRR_all_chan_HGF_arcsinh_20190126T101431.mat' % arc, no, lp, epsi
% '' % arc, no, t, dau/da
% 'stats_BRR_all_chan_HGF_arcsinh_20190126T080844.mat' % arc, no, t, epsi
% '' % arc, arc, lp, dau/da
% 'stats_BRR_all_chan_HGF_arcsinh_20190126T140224.mat' % arc, arc, lp, epsi
% '' % arc, arc, t, dau/da
% 'stats_BRR_all_chan_HGF_arcsinh_20190126T120002.mat' % arc, arc, t, epsi
    
% 'stats_BRR_all_chan_HGF_notrans_20190126T160735.mat' % BRR, no, arc, epsi1, 20 downsample
% 'stats_BRR_all_chan_HGF_notrans_20190126T180844.mat' % BRR, no, arc, epsi2, 20 downsample
% 'stats_BRR_all_chan_HGF_notrans_20190126T200615.mat' % BRR, no, arc, epsi3, 20 downsample

% 'stats_PEB_all_chan_null_arcsinh_20190125T213550.mat' %PEB, null

% 'stats_PEB_all_chan_null_notrans_20190126T201046.mat' %PEB null
% 'stats_PEB_all_chan_HGF_notrans_20190126T184044.mat' %PEB, no, arc, epsi3, 20 downsample
% 'stats_PEB_all_chan_HGF_notrans_20190126T171529.mat' %PEB, no, arc, epsi2, 20 downsample
% 'stats_PEB_all_chan_HGF_notrans_20190126T154728.mat' %PEB, no, arc, epsi1, 20 downsample
% 'stats_PEB_all_chan_HGF_notrans_20190129T184033.mat' % ud3
% 'stats_PEB_all_chan_HGF_notrans_20190129T165435.mat' % ud2
% 'stats_PEB_all_chan_HGF_notrans_20190129T165351.mat' % ud1

% 0 - 780ms
% 'stats_BRR_all_chan_HGF_notrans_20190126T225906.mat' % ep1
% 'stats_BRR_all_chan_HGF_notrans_20190127T013657.mat' % ep2
% 'stats_BRR_all_chan_HGF_notrans_20190127T041903.mat' % ep3


% 'stats_PEB_all_chan_condHGF_notrans_20190127T214945.mat' % condHGF with epsi x3
% 'stats_PEB_all_chan_HGF_notrans_20190128T033817.mat' % minus cond
% 'stats_PEB_all_chan_condHGF_notrans_20190127T215347.mat' % minus ep1
% 'stats_PEB_all_chan_condHGF_notrans_20190128T030529.mat' % minus ep2
% 'stats_PEB_all_chan_condHGF_notrans_20190128T075910.mat' % minus ep3

% PEB  highres
% 'stats_PEB_all_chan_HGF_notrans_20190127T095132.mat' %PEB, no, arc, epsi3
% 'stats_PEB_all_chan_HGF_notrans_20190128T092123.mat' %PEB, no, arc, epsi2
% 'stats_PEB_all_chan_HGF_notrans_20190127T095014.mat' %PEB, no, arc, epsi1
% 'stats_PEB_all_chan_cond_notrans_20190128T092155.mat' %PEB, no, arc, cond
% 'stats_PEB_all_chan_HGF_notrans_20190128T231850.mat' %Dau

% BRR with residuals
%'stats_BRR_all_chan_condHGF_notrans_20190130T220808.mat' % cond, mu1-3.
% BRR on above residual data
% 'stats_BRR_all_chan_HGF_notrans_20190131T102539.mat' %ep3
% 'stats_BRR_all_chan_HGF_notrans_20190131T084414.mat' %ep2
% 'stats_BRR_all_chan_HGF_notrans_20190131T084057.mat' %ep1

% 'stats_BRR_all_chan_HGF_notrans_20190131T142456.mat'% BRR with full 15 traj HGF (psi), ridge t
% 'stats_BRR_all_chan_HGF_notrans_20190131T175758.mat'% BRR with full 15 traj HGF (ud), ridge t
% 'stats_BRR_all_chan_condHGF_notrans_20190131T215250.mat'
% 'stats_BRR_all_chan_cond_notrans_20190131T170005.mat'%cond
% 'stats_BRR_all_chan_HGF_notrans_20190201T194803.mat' % Dau
% 'stats_BRR_all_chan_HGF_notrans_20190202T085955.mat' % dau abs
% 'stats_BRR_all_chan_HGF_notrans_20190202T105718.mat' % dau rect
% 'stats_BRR_all_chan_HGF_notrans_20190202T090539.mat' % epsi2 abs
% 'stats_BRR_all_chan_HGF_notrans_20190202T110013.mat' %epsi2 rect 

% downsample 4, 800ms
% 'stats_BRR_all_chan_HGF_notrans_20190202T220500.mat'% BRR with full 15 traj HGF (ud), ridge t
% 'stats_BRR_all_chan_condHGF_notrans_20190203T201944.mat' % 15 traj HGF (ud), ridge t
% 'stats_BRR_all_chan_condHGF_notrans_20190203T133856.mat' % cond, 12 traj HGF (NO ud)
% 'stats_BRR_all_chan_HGF_notrans_20190203T133828.mat' % 12 traj HGF (NO ud)

% 'stats_BRR_all_chan_condHGF_notrans_20190206T220218.mat' % Dau/da
% 'stats_BRR_all_chan_HGF_notrans_20190206T220406.mat' % Dau/da
% 'stats_BRR_all_chan_condHGF_notrans_20190204T222014.mat' % cond, epsi1-3
% 'stats_BRR_all_chan_HGF_notrans_20190204T222045.mat' % epsi1-3
% 'stats_BRR_all_chan_condHGF_notrans_20190206T220526.mat' % All PE
% 'stats_BRR_all_chan_HGF_notrans_20190206T220644.mat' % All PE

% 'stats_BRR_all_chan_HGF_notrans_20190205T181917.mat' % epsi1-3 Lasso
% 'stats_BRR_all_chan_condHGF_notrans_20190205T182003.mat' % t lasso, cond, epsi1-3
% 'stats_BRR_all_chan_condHGF_notrans_20190206T093958.mat' % t hs+, cond, epsi1-3
% 'stats_BRR_all_chan_condHGF_notrans_20190206T165012.mat' % gaus ridge, cond, epsi1-3 it14
% 'stats_BRR_all_chan_condHGF_notrans_20190206T175108.mat' % gaus ridge, cond, epsi1-3 it7
% 'stats_BRR_all_chan_condHGF_notrans_20190206T135013.mat' % gaus lasso, cond, epsi1-3
% 'stats_BRR_all_chan_condHGF_notrans_20190206T094015.mat' % gaus hs+, cond, epsi1-3

% 'stats_BRR_all_chan_cond_notrans_20190207T093201.mat' % no catvars
% 'stats_BRR_all_chan_cond_notrans_20190207T093200.mat' % catvars

% 'stats_BRR_all_chan_condHGF_notrans_20190208T113720.mat' % freq1
% 'stats_BRR_all_chan_condHGF_notrans_20190207T181656.mat' % freq2
% 'stats_BRR_all_chan_condHGF_notrans_20190207T181535.mat' % freq3
% 'stats_BRR_all_chan_condHGF_notrans_20190207T181734.mat' % freq4

% updated models; odd only
% 'stats_BRR_all_chan_HGF_notrans_20190213T111930.mat' % all PE, ridge, t
% 'stats_BRR_all_chan_HGF_notrans_20190213T135853.mat' % epsis, ridge, t

% contain pred
'stats_BRR_all_chan_condHGF_notrans_20190216T072318' % epsi, ridge
% 'stats_BRR_all_chan_condHGF_notrans_20190216T072519' % epsi, hs+
    };

subtract = [];
encoding_types = {'spear','MR','PEB','RR','BRR','mvpa'}; % only these will be analysed
% statfield = {'cE','beta','b','rho','weights','transweights'}; % for t-test stats
statfield = {'b'}; % for t-test stats
recon = []; % multiple beta by predictors specified here
% statfield = {'waic','F',};
% statfield = {'logl'}; % for PEB
%statfield = {'beta','b','s','rho','weights','transweights','kurt','skew'}; % for t-test stats
%statfield = {'s','r2'};
% xticks = 0:20:600;
% xticks = 0:20:780;
xticks = 0:4:796;
% xticks = 0:8:792;
% xticks = 0:4:600;
topo_range = [min(xticks) max(xticks)];
plotsmooth = 10;
no_plot_ele = [];
% param={1}; % multiple param's betas are summed within each cell
% param={[1],[2],[3],[4],[5],[6]}; % multiple param's betas are summed within each cell
% param_legend = {'Dau','Da1','Da2','Ep1','Ep2','Ep3'};
% param={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13]}; % multiple param's betas are summed within each cell
% param_legend = {'Cond','Dau','Da1','Da2','Ep1','Ep2','Ep3','Mu1','Mu2','Mu3','Sa1','Sa2','Sa3'};
% param={[1],[2],[3],[4]}; % multiple param's betas are summed within each cell
% param_legend = {};%'cond','Ep1','Ep2','Ep3'};
param={[4],[3],[2],[1]}; % multiple param's betas are summed within each cell
param_legend = {'epsi3','epsi2','epsi1','oddball'};
colormaps = {
    cbrewer('seq', 'Greens', 100, 'pchip')
    cbrewer('seq', 'Blues', 100, 'pchip')
    cbrewer('seq', 'Purples', 100, 'pchip')
    cbrewer('seq', 'YlOrBr', 100, 'pchip')
    };
load('C:\Data\CORE\eeg\ana\prep\chanlocs.mat')
cmap=colormap('parula'); close all;
chan_summary = 'std'; % mean or std
test_models = [];%2:length(sfiles); % for comparing waic
null_model = [];%1; % if more than one model (file), a null is created from the mean of them
% test_models = 1:length(sfiles); 
% null_model = 1:length(sfiles); % if more than one model (file), a null is created from the mean of them
model_comparison = 'sum'; % none, sum, time, or chan_time
options.DisplayWin = 0;
F_smooth = 0; % variance of Gaussian filter
save_stats = 0;
grp_effect = 0;

% get data
allstat={};
i=0; % index for fig handles
if length(sfiles)>1 
    for tm=1:length(test_models)
        h0m(tm)=figure('Name','overlap');
    end
end
% for each model/file
for f = 1:length(sfiles)
    load(fullfile(spath,sfiles{f}));
    allstat{f} = stats;
    
    % get group info
    grplist = S.designmat(2:end,strcmp(S.designmat(1,:),'groups'));
    [grpuni,iA,iB] = unique(grplist,'stable');
    
    % groups stats for each analysis
    an = fieldnames(allstat{f});
    an = an(ismember(an,encoding_types));
    for a = 1:length(an)
        da = fieldnames(allstat{f}.(an{a}));
        sf_names=fieldnames(allstat{f}.(an{a}).(da{1}));
        statfield = statfield(ismember(statfield,sf_names));
        
        for sf = 1:length(statfield)
        
            % if there are multiple runs, average the outputs
            for d = 1:length(da)
                nsub = size(allstat{f}.(an{a}).(da{d}),1);
                nrun = size(allstat{f}.(an{a}).(da{d}),2);
                if nrun>1
                    runcell={};meandat={};
                    for ns = 1:nsub
                        for nr = 1:nrun
                            runcell{ns,nr} = allstat{f}.(an{a}).(da{d})(ns,nr).(statfield{sf});
                        end
                        meandat{ns,1} = mean(cat(3,runcell{ns,:}),3);
                    end
                    allstat{f}.(an{a}).(da{d})(1).(statfield{sf}) = meandat;
                end
            end
        end
        allstat{f}.(an{a}).(da{d})(2:end) = [];
        
        for sf = 1:length(statfield)
            
            % number of columns
            nc = size(allstat{f}.(an{a}).(da{1}).(statfield{sf}),2);
            for c = 1:nc 
                    
                for d = 1:length(da)
                    disp(['file' num2str(f) ', analysis: ' an{a} ', data type: ' da{d}, ', comp: ' num2str(c)])
                    
                    % for each parameter
                    if length(param)>1
                        h0=figure('Name','overlap');
                        h0f=figure('Name','overlap');
                    end
                    for np = 1:length(param)
                        h1=figure;  
                        h2=figure('Name',['file_' num2str(f) '_' an{a} ', comp: ' num2str(c)]);
                        if grp_effect
                            h3=figure('Name',['file_' num2str(f) '_' an{a} ', comp: ' num2str(c)]);pl3 = 0; % plot index
                            h4=figure('Name',['file_' num2str(f) '_' an{a} ', comp: ' num2str(c)]);pl4 = 0; % plot index
                        end

                        % number of rows (e.g. subjects)
                        nr = size(allstat{f}.(an{a}).(da{d}).(statfield{sf}),1);

                        clear alldat h p ci
                        for r=1:nr
                            dat=allstat{f}.(an{a}).(da{d}).(statfield{sf}){r,c};
                            if strcmp(statfield{sf},'waic')
                                dat = -dat;
                            end
                            if ndims(dat)==3
                                dat=sum(dat(:,:,param{np}),3);
                            elseif ndims(dat)==2 && strcmp(da{d},'GFP')
                                dat=sum(dat(:,param{np}),3);
    %                         else
    %                             dat = dat;
                            end
                            sizdat=size(dat);
                            if ~any(sizdat==1) % if a matrix
                                dat = reshape(dat,prod(sizdat),[]);
                            end
                            if ~isempty(recon)
                                pred=[stats.trialinfo{1}.pred_train{r}(:,recon)];
                                dat = mean(pred*repmat(dat,1,length(recon))',1)';
                            end
                            alldat(r,:)=dat;
                        end
                        grpmeandat{f} = mean(alldat,1);

                        %% one-sample t-test (vs. zero) on beta weights
                        for s=1:size(alldat,2)
                            try
                                [h(s),p(s),~,stt] = ttest(double(alldat(:,s)));
                                t{np}(s)=stt.tstat;
                            catch
                                h(s)=NaN;
                                p(s)=NaN;
                                t{np}(s)=NaN;
                            end
                        end

                        % remove NaNs
                        t{np}(isnan(t{np}))=0;
                        h(isnan(h))=0;
                        p(isnan(p))=Inf;

                        % FDR correction
                        [~,fdr_mask] = fdr(p,0.05);
                        fdr_p=p.*double(fdr_mask);
                        fdr_t{np}=t{np}.*double(fdr_mask);
                        % reshape
                        if ~any(sizdat==1) 
                            h=reshape(h,sizdat(1),sizdat(2));
                            p=reshape(p,sizdat(1),sizdat(2));
                            t{np}=reshape(t{np},sizdat(1),sizdat(2));
                            fdr_t{np}=reshape(fdr_t{np},sizdat(1),sizdat(2));
                            fdr_p=reshape(fdr_p,sizdat(1),sizdat(2));
                            grpmeandat{f}=reshape(grpmeandat{f},sizdat(1),sizdat(2));
                        end
                        % output
                        allstat{f}.(an{a}).(da{d}).([statfield{sf} '_mean'])=reshape(mean(alldat,1)',sizdat(1),sizdat(2));
                        allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest']).h=h;
                        allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest']).p=p;
                        allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest']).t=t{np};
                        allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest']).fdr_t=fdr_t{np};
                        allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest']).fdr_p=fdr_p;

                        trange = dsearchn(xticks',[topo_range(1);topo_range(2)]);

                        %% plot group mean of data

                        figure(h1)
                        clim=[min(grpmeandat{f}(:)) max(grpmeandat{f}(:))];
                        if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
                        subplot(length(sfiles),1,f)
                        imagesc(xticks,[],grpmeandat{f},clim); 
%                         [~,mi]=max(mean(abs(grpmeandat{f}(:,trange(1):trange(2))),1));
%                         mi=mi+trange(1)-1;
%                         line(xticks(mi*[1 1]),[1 92],'color','k','linewidth',2)
%                         title([statfield{sf} ', group mean values']);
                        colorbar

%                         pl=pl+1;
%                         if ~any(sizdat==1) 
%                             subplot(3,2,3)
%                             topoplot(grpmeandat{f}(:,mi),chanlocs,'maplimits',clim,'electrodes','off','style','map','intrad',0.55,'colormap',cmap);
%                         end
                        
%                         subplot(3,2,[5:6])
%                         plot(xticks,mean(grpmeandat{f},1),'b'); 
%                         title([statfield{sf} ', ' chan_summary ' over channels of group mean values']);
%                         ylabel(chan_summary)
%                         xlabel('post-stimulus time (ms)')

                        %% plot t stats
                        figure(h2)
                        clim=[min(t{np}(:)) max(t{np}(:))];
                        if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
                        npeaks = 5;
                        subplot(length(da)*3,10,[1:5])
                        imagesc(xticks,[],t{np},clim); 
                       % meanabs = mean(abs(t{np}(:,trange(1):trange(2))),1);
                       if strcmp(chan_summary,'std')
                            tval_summ(np,:) = std(t{np}(:,trange(1):trange(2)),[],1);
                       elseif strcmp(chan_summary,'mean')
                            tval_summ(np,:) = mean(t{np}(:,trange(1):trange(2)),1);
                       elseif strcmp(chan_summary,'absmean')
                            tval_summ(np,:) = mean(abs(t{np}(:,trange(1):trange(2))),1);
                       end
                        tval_summ(np,:) = smooth(tval_summ(np,:),plotsmooth,'moving');
                        [ma,mi]=findpeaks(tval_summ(np,:),'MinPeakWidth',2);
                        mi=mi+trange(1)-1;
                        % highest peaks, plot as lines
                        [maxma,maxmi]=sort(ma,'descend');
                        maxmi = maxmi(1:min(npeaks,length(maxmi)));
                        peaks = sort(mi(maxmi));
                        hold on
                        for mii=1:length(peaks)
                            plot(xticks(mi(maxmi(mii))*[1 1]),[1 92],'k:','linewidth',1)
                        end
                        title([statfield{sf} ', t values for one-sample t-test']);
                        ylabel('channel')
                        %colorbar

                        subplot(length(da)*3,10,[11:15])
                        plot(xticks,tval_summ(np,:),'b'); 
                        yl = ylim;
                        hold on
                        for mii=1:length(peaks)
                            plot(xticks(mi(maxmi(mii))*[1 1]),[yl(1) ma(maxmi(mii))],'k:','linewidth',0.5)
                        end
                        title([statfield{sf} ', ' chan_summary ' over channels of t values']);
                        ylabel(chan_summary)
                        xlabel('post-stimulus time (ms)')

                        % plot 4 highest peaks
                        if ~any(sizdat==1) 
                            for mii = 1:npeaks
                                subplot(length(da)*3,10,20+mii)
                                try
                                    topoplot(t{np}(:,peaks(mii)),chanlocs,'maplimits',clim,'electrodes','off','style','map','intrad',0.55,'colormap',cmap);
                                    title([num2str(xticks(peaks(mii))) ' ms'])
                                end
                            end
                        end

                        clim=[min(fdr_t{np}(:)) max(fdr_t{np}(:))];
                        if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
                        if ~any(clim)
                            clim=[-5 5];
                        end

                        subplot(length(da)*3,10,[6:10])
                        imagesc(xticks,[],fdr_t{np},clim); 
                       if strcmp(chan_summary,'std')
                            fdr_tval_summ(np,:) = std(fdr_t{np}(:,trange(1):trange(2)),[],1);
                       elseif strcmp(chan_summary,'mean')
                            fdr_tval_summ(np,:) = mean(fdr_t{np}(:,trange(1):trange(2)),1);
                       elseif strcmp(chan_summary,'absmean')
                            fdr_tval_summ(np,:) = mean(abs(fdr_t{np}(:,trange(1):trange(2))),1);
                       end
                        fdr_tval_summ(np,:) = smooth(fdr_tval_summ(np,:),plotsmooth,'moving');
                        [ma,mi]=findpeaks(fdr_tval_summ(np,:),'MinPeakWidth',2);
                        mi=mi+trange(1)-1;
                        % 4 highest peaks, plot as lines
                        [maxma,maxmi]=sort(ma,'descend');
                        maxmi = maxmi(1:min(npeaks,length(maxmi)));
                        peaks = sort(mi(maxmi));
                        hold on
                        for mii=1:length(peaks)
                            plot(xticks(mi(maxmi(mii))*[1 1]),[1 92],'k:','linewidth',1)
                        end
                        title([statfield{sf} ', fdr thresholded t values']);
                        ylabel('channel')
                        %colorbar

                        subplot(length(da)*3,10,[16:20])
                        plot(xticks,fdr_tval_summ(np,:),'b'); 
                        yl=ylim;
                        hold on
                        for mii=1:length(peaks)
                            plot(xticks(mi(maxmi(mii))*[1 1]),[yl(1) ma(maxmi(mii))],'k:','linewidth',0.5)
                        end
                        title([statfield{sf} ', ' chan_summary ' over channels of fdr-thresholded t values']);
                        ylabel(chan_summary)
                        xlabel('post-stimulus time (ms)')

                        if ~any(sizdat==1) 
                            for mii = 1:npeaks
                                subplot(length(da)*3,10,25+mii)
                                try
                                    topoplot(fdr_t{np}(:,peaks(mii)),chanlocs,'maplimits',clim,'electrodes','off','style','map','intrad',0.55,'colormap',cmap);
                                    title([num2str(xticks(peaks(mii))) ' ms'])
                                end
                            end
                        end

                        if grp_effect
                            %% two-sample t-test
                            clear t h p stt
                            for g = 1:length(grpuni)
                                grpdat{g}=alldat(strcmp(grplist,grpuni{g}),:);
                            end
                            for s=1:size(alldat,2)
                                [h(s),p(s),~,stt] = ttest2(double(grpdat{1}(:,s)),double(grpdat{2}(:,s)));
                                t{np}(s)=stt.tstat;
                            end

                            % remove NaNs
                            t{np}(isnan(t{np}))=0;
                            h(isnan(h))=0;
                            p(isnan(p))=Inf;

                            % FDR correction
                            [~,fdr_mask] = fdr(p,0.05);
                            fdr_p=p.*double(fdr_mask);
                            fdr_t{np}=t{np}.*double(fdr_mask);
                            % reshape
                            if ~any(sizdat==1) 
                                h=reshape(h,sizdat(1),sizdat(2));
                                p=reshape(p,sizdat(1),sizdat(2));
                                t{np}=reshape(t{np},sizdat(1),sizdat(2));
                                fdr_t{np}=reshape(fdr_t{np},sizdat(1),sizdat(2));
                                fdr_p=reshape(fdr_p,sizdat(1),sizdat(2));
                            end
                            % output
                            allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest2']).h=h;
                            allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest2']).p=p;
                            allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest2']).t=t{np};
                            allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest2']).fdr_t=fdr_t{np};
                            allstat{f}.(an{a}).(da{d}).([statfield{sf} '_grpttest2']).fdr_p=fdr_p;

                            trange = dsearchn(xticks',[topo_range(1);topo_range(2)]);

                            clim=[min(t{np}(:)) max(t{np}(:))];

                            figure(h3)
                            hold on
                            pl3=pl3+1;
                            subplot(length(da)*2,2,pl3)
                            imagesc(xticks,[],t{np},clim); 
                            [~,mi]=max(mean(abs(t{np}(:,trange(1):trange(2))),1));
                            mi=mi+trange(1)-1;
                            line(xticks(mi*[1 1]),[1 92],'color','k','linewidth',2)
                            title([statfield{sf} ', t values for two-sample t-test']);
                            colorbar

                            pl3=pl3+1;
                            if ~any(sizdat==1) 
                                subplot(length(da)*2,2,pl3)
                                topoplot(t{np}(:,mi),chanlocs,'maplimits',clim);
                            end

                            clim=[min(fdr_t{np}(:)) max(fdr_t{np}(:))];
                            if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
                            if ~any(clim)
                                clim=[-5 5];
                            end

                            pl3=pl3+1;
                            subplot(length(da)*2,2,pl3)
                            imagesc(xticks,[],fdr_t,clim); 
                            [~,mi]=max(mean(abs(fdr_t{np}(:,trange(1):trange(2))),1));
                            mi=mi+trange(1)-1;
                            line(xticks(mi*[1 1]),[1 92],'color','k','linewidth',2)
                            title([statfield{sf} ', fdr thresholded t values']);
                            colorbar

                            pl3=pl3+1;
                            if ~any(sizdat==1) 
                                subplot(length(da)*2,2,pl3)
                                topoplot(fdr_t{np}(:,mi),chanlocs,'maplimits',clim);
                            end
% 

                            figure(h4)
                            pl4=pl4+1;
                            subplot(length(da),1,pl4)
                            scatdat=[];
                            for g = 1:length(grpdat)
                                dat=reshape(grpdat{g},[],sizdat(1),sizdat(2));
                                scatdat=[scatdat;mean(abs(dat(:,:,mi)),2)];
                            end
                            scatter(iB,scatdat);
                            set(gca,'xtick',1:2,'XTickLabel',grpuni);
                            title('group difference in mean absolute EEG at time of greatest group difference')
    % 
                            
                        end
                        % model comparison
                        if strcmp(an{a},'BRR')
                            % separate LME into groups for model comparison
%                             LME=[stats.BRR.alldata.logl];
%                             for g = 1:length(grpuni)
%                                 LMEall=LME(strcmp(grplist,grpuni{g}));
%                                 LMEmat = cat(3,LMEall{:});
%                                 LMEmean =double(squeeze(mean(mean(LMEmat,1),2)));
%                                 LMEgrp{1,g}(f,:) = LMEmean;
%                             end

                            % separate WAIC into groups for model comparison
                            LME=[stats.BRR.alldata.waic];
                            for g = 1:length(grpuni)
                                LMEall=LME(strcmp(grplist,grpuni{g}));
                                LMEmat = cat(3,LMEall{:});
                                LMEmean =double(squeeze(mean(mean(LMEmat,1),2)));
                                LMEgrp{1,g}(f,:) = -LMEmean;
                                LMEgrpmat{1,g}(f,:,:,:) = -LMEmat;
                            end
                        elseif strcmp(an{a},'PEB')
                            % separate WAIC into groups for model comparison
                            LME=[stats.PEB.alldata.F];
                            for g = 1:length(grpuni)
                                LMEall=LME(strcmp(grplist,grpuni{g}));
                                LMEmat = cat(3,LMEall{:});

                                if F_smooth
                                    for s = 1:size(LMEmat,3)
                                        LMEmat(:,:,s) = imgaussfilt(LMEmat(:,:,s),F_smooth);
                                    end
                                end

                                LMEmean =double(squeeze(mean(mean(LMEmat,1),2)));
                                LMEgrp{1,g}(f,:) = LMEmean;
                                LMEgrpmat{1,g}(f,:,:,:) = LMEmat;
                            end
                        end
                    end

                    % overlapping summary stats over channels for each
                    % param
                    if exist('h0','var')
                        nparam=size(tval_summ,1);
                        
                        % no FDR
                        fig=figure(h0)
                        set(fig, 'Units', 'normalized', 'Position', [0,0,0.2,0.5]);
                        for np = 1:nparam
                            ax(np)=subplot(nparam+1,1,np)
                            clim=[min(t{np}(:)) max(t{np}(:))];
                            if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
                            if ~any(clim)
                                clim=[-5 5];
                            end
                            imagesc(xticks,[],t{np},clim); 
                            ylabel(param_legend{np},'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                            set(ax(np),'YColor',colormaps{np}(90,:),'XTick',[]);
                            cl(:,np) = caxis;
                            colorbar
                        end
                        set(ax, 'CLim', [mean(cl(1,:)), mean(cl(2,(cl(2,:)>1)))]);
                        pause(1)
                        pos=get(ax,'Position'); 
                        pos=cat(1,pos{:}); 
                        pos(:,1) = pos(:,1)*1.2;
                        ax(np+1)=subplot(nparam+1,1,np+1);
                        for np = 1:nparam
                            set(ax(np),'Position',pos(np,:))
                            hold on
                            plot(xticks,tval_summ(np,:),'color',colormaps{np}(70,:),'Linewidth',2); 
                        end
                        ylabel({'standard';'deviation';'over';'channels'},'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                        xlabel('post-stimulus time (ms)')  
                        try
                            [~, hobj, ~, ~] = legend(param_legend);
                        catch
                            [~, hobj, ~, ~] = legend(cellfun(@num2str,param,'UniformOutput',0));
                        end
                        lh = findobj(hobj,'type','line');
                        set(lh,'LineWidth',2);
                        pos1 = get(ax(1),'Position');
                        pos2 = get(ax(np+1),'Position');
                        pos2([1 3]) = pos1([1 3]);
                        set(ax(np+1),'Position',pos2)

                        % FDR
                        fig=figure(h0f)
                        set(fig, 'Units', 'normalized', 'Position', [0,0,0.2,0.5]);
                        clear ax
                        for np = 1:nparam
                            ax(np)=subplot(nparam+1,1,np)
                            clim=[min(fdr_t{np}(:)) max(fdr_t{np}(:))];
                            if ~any(clim)
                                clim=[-5 5];
                            end
                            if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
                            imagesc(xticks,[],fdr_t{np},clim); 
                            ylabel(param_legend{np},'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                            set(ax(np),'YColor',colormaps{np}(90,:),'XTick',[]);
                            cl(:,np) = caxis;
                            colorbar
                        end
                        set(ax, 'CLim', [mean(cl(1,:)), mean(cl(2,(cl(2,:)>1)))]);
                        pause(1)
                        pos=get(ax,'Position'); 
                        pos=cat(1,pos{:}); 
                        pos(:,1) = pos(:,1)*1.2;
                        ax(np+1)=subplot(nparam+1,1,np+1)
                        for np = 1:nparam
                            set(ax(np),'Position',pos(np,:))
                            hold on
                            plot(xticks,fdr_tval_summ(np,:),'color',colormaps{np}(70,:),'Linewidth',2); 
                        end
                        ylabel({'standard';'deviation';'over';'channels'},'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                        xlabel('post-stimulus time (ms)')
                        try
                            [~, hobj, ~, ~] = legend(param_legend);
                        catch
                            [~, hobj, ~, ~] = legend(cellfun(@num2str,param,'UniformOutput',0));
                        end
                        lh = findobj(hobj,'type','line');
                        set(lh,'LineWidth',2);
                        pos1 = get(ax(1),'Position');
                        pos2 = get(ax(np+1),'Position');
                        pos2([1 3]) = pos1([1 3]);
                        set(ax(np+1),'Position',pos2)
                    end
                    
                    
                end
            end
        end
        
        
        % plot residuals
        if isfield(allstat{f}.(an{a}).(da{1}),'resid')
            resid=allstat{f}.(an{a}).(da{1}).resid;
            % subject means
            for i = 1:length(resid)
                resid{i,1} = mean(resid{i,1},3);
            end
            % group mean
            resid_mean = mean(cat(3,resid{:}),3);
            % group std
            resid_gfp = std(resid_mean,[],1);
            
            clim=[min(resid_mean(:)) max(resid_mean(:))];
            if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
            figure
            subplot(2,2,[1 2])
            imagesc(xticks,[],resid_mean,clim); 
            title(['residuals, group mean']);
            subplot(2,2,[3 4])
            plot(xticks,resid_gfp,'b'); 
            title(['residuals, group mean global field power']);
            xlabel('post-stimulus time (ms)')
        end
        
        % compare groups regarding predictions
        try
            recons = allstat{f}.biem.alldata.recons;
        catch
            try
                for i = 1:length(allstat{f}.mvpa.alldata)
                    recons{i} = mean(cat(3,{allstat{f}.mvpa.alldata(i,:).testdata_pred}),3);
                end
            end
        end
        try
            dat=cellfun(@mean,recons,'uniformoutput',0);
            subidx=ismember(stats.subname,S.designmat(2:end,2));
            dat=dat(subidx);
            p_rs = ranksum(dat(iB==1),dat(iB==2))
            figure; scatter(iB,dat);
            set(gca,'xtick',1:2,'XTickLabel',grpuni);
            title(['group difference in mean reconstructed predictor variable: p = ' num2str(p_rs)])
        end
        
    end
    
    if save_stats
        stats=allstat{f};
        sname = strrep(sfiles{f},'stats','stats_grp');
        save(fullfile(spath,sname),'stats');
    end
    
end

% overlapping means
if exist('h0m','var')
    null = mean(cat(3,grpmeandat{null_model}),3);
    test = grpmeandat(test_models);
    for tm = 1:length(test)
        if ~isempty(null)
            grpmeandat_tm = test{tm} - null;
        else
            grpmeandat_tm = test{tm};
        end
        clim=[min(grpmeandat_tm(:)) max(grpmeandat_tm(:))];
        if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
        figure(h0m(tm))
        hold on
        subplot(2,2,[1 2])
        imagesc(xticks,[],grpmeandat_tm,clim); 
        [~,mi]=max(mean(abs(grpmeandat_tm(:,trange(1):trange(2))),1));
        mi=mi+trange(1)-1;
        line(xticks(mi*[1 1]),[1 92],'color','k','linewidth',2)
        title([statfield{sf} ', group mean values']);
        colorbar

        subplot(2,2,[3 4])
        plot(xticks,mean(grpmeandat_tm,1),'b'); 
        title([statfield{sf} ', ' chan_summary ' over channels of group mean values']);
        ylabel(chan_summary)
        xlabel('post-stimulus time (ms)')
    end
end

% subtract R2 between datasets
if ~isempty(subtract) && length(allstat)==2
    for r=1:nr 
        %datsub{1}=allstat{1}.MR.alldata.R2{r,1};
        %datsub{2}=allstat{2}.MR.alldata.R2{r,1};
        datsub{1}=allstat{1}.BRR.alldata.r2{r,1};
        datsub{2}=allstat{2}.BRR.alldata.r2{r,1};
        subdat{r,1} = subtract(1)*datsub{1} + subtract(2)*datsub{2};
    end
    stats = allstat{1};
    %stats.MR.alldata.R2 = subdat;
    stats.BRR.alldata.r2 = subdat;
    sname = strrep(sfiles{f},'stats','stats_subtractR2');
    save(fullfile(spath,sname),'stats','S');
end
    
switch model_comparison
    case 'sum'
        %group model comparison
        if exist('LMEgrp','var')
            runvar=LMEgrp;
            if exist('LMEgrp','var') && size(runvar{1},1)>1
                pep_flag=1;
                [~,bmc.gposterior,bmc.gout] = VBA_groupBMC_btwGroups_CAB(runvar,options,pep_flag)
            end
        end
    case 'chan_time'
        % group model comparison over chans/times
        if exist('LMEgrp','var')
            if exist('LMEgrpmat','var') && size(LMEgrpmat{1},1)>1
                sz=size(LMEgrpmat{1});
                ndat = sz(2)*sz(3);
                pep_flag=1;
                options.DisplayWin = 0;
                for g = 1:length(LMEgrpmat)
                    temp = reshape(LMEgrpmat{g},sz(1),ndat,sz(4));
                    for n = 1:ndat
                        runvar{n}{g}=squeeze(temp(:,n,:));
                    end
                end
                for n = 1:ndat
                    L = [runvar{n}{1} runvar{n}{2}]; % for now, don't separate groups
                    [bmc(n).gposterior,bmc(n).gout] = VBA_groupBMC_cab(L,options,pep_flag);
                end
                for n = 1:ndat
                    ep(n,:) = bmc(n).gout.ep;
                    pep(n,:) = bmc(n).gout.pep;
                end
                for f = 1:size(ep,2)
                    all_ep{f} = reshape(ep(:,f),sz(2),sz(3));
                    all_pep{f} = reshape(pep(:,f),sz(2),sz(3));
                end
                thresh = 0.95;
                figure
                for f = 1:length(all_ep)
                    subplot(length(all_ep),1,f)
                    imagesc(xticks,[],all_ep{f},[thresh 1]); 
                    set(gca,'fontsize',12)
                end
                figure
                for f = 1:length(all_pep)
                    subplot(length(all_pep),1,f)
                    imagesc(xticks,[],all_pep{f},[thresh 1]); 
                    set(gca,'fontsize',12)
                end
            end
        end
     case 'time'
        % group model comparison over time
        if exist('LMEgrp','var')
            if exist('LMEgrpmat','var') && size(LMEgrpmat{1},1)>1
                sz=size(LMEgrpmat{1});
                ndat = sz(3);
                pep_flag=1;
                options.DisplayWin = 0;
                for g = 1:length(LMEgrpmat)
                    temp = squeeze(mean(LMEgrpmat{g},2));
                    for n = 1:ndat
                        runvar{n}{g}=squeeze(temp(:,n,:));
                    end
                end
                for n = 1:ndat
                    L = [runvar{n}{1} runvar{n}{2}]; % for now, don't separate groups
                    [bmc(n).gposterior,bmc(n).gout] = VBA_groupBMC_cab(L,options,pep_flag);
                end
                for n = 1:ndat
                    ep(n,:) = bmc(n).gout.ep;
                    pep(n,:) = bmc(n).gout.pep;
                end
                for f = 1:size(ep,2)
                    all_ep{f} = ep(:,f);
                    all_pep{f} = pep(:,f);
                end
                thresh = 0.95;
                figure
                for f = 1:length(all_ep)
                    hold on
                    plot(xticks,all_ep{f}); 
                end
                figure
                for f = 1:length(all_pep)
                    hold on
                    plot(xticks,all_pep{f}); 
                end
            end
        end
end
% best predictors
if 0
    for f = 1:length(allstat)
        
        % absolute betas (probably not very useful!)
        meanb=squeeze(mean(mean(abs(cat(4,allstat{f}.BRR.alldata.b{:})),2),1));
        x=repmat([1:size(meanb,1)]',1,size(meanb,2));
        figure;scatter(x(:),meanb(:))
        title(['mean b for file ' num2str(f)])
    end
end

% plot R2
clear meanr2
if 0
    for f = 1:length(allstat)
        meanr2{f}=squeeze(mean(mean(cat(3,allstat{f}.BRR.alldata.r2{:}),2),1));
        x=repmat([1:size(meanr2{f},1)]',1,size(meanr2{f},2));
        figure;scatter(x(:),meanr2{f}(:))
        title(['r2 for file ' num2str(f)])
    end
    if f==2
        submeanr2= subtract(1)*meanr2{1} + subtract(2)*meanr2{2};
        x=repmat([1:size(submeanr2,1)]',1,size(submeanr2,2));
        figure;scatter(x(:),submeanr2(:))
        title(['subtracted r2'])
        pc_r2 = submeanr2*100./meanr2{1};
        figure;scatter(x(:),pc_r2(:))
        title(['percent increase in r2'])
    end
end

