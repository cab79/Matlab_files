clear all
dbstop if error

addpath('C:\Data\Matlab\CoSMoMVPA\mvpa')
addpath('C:\Data\Matlab\Matlab_files\CORE\EEG\cosmo\modified functions')
addpath(genpath('C:\Data\Matlab\PRoNTo_dev-2.0.1'));
addpath(genpath('C:\Data\Matlab\TFCE'));
S.icaDir = 'C:\Data\CORE\eeg\ana\groupICA';
%S.compname = {'CORE_part2_ica_c','-1.mat'};
S.compname = {'CORE','_2_merged_cleaned_grp-ica_c.mat'};
%S.infoname = {'fileinfo_','.mat'};
S.infoname = {'CORE','_2_merged_cleaned_grp-fileinfo_.mat'};
S.cond_idx = {
%     [1 2 9 10 17 18] %left hand, mismatch
%     [3 4 11 12 19 20] %left hand, standard
%     [5 6 13 14 21 22] %right hand, mismatch
%     [7 8 15 16 23 24]} %right hand, standard

    [1 2 9 10] %left hand, mismatch
    [3 4 11 12] %left hand, standard
    [5 6 13 14] %right hand, mismatch
    [7 8 15 16]}; %right hand, standard
S.contrast_rows = {[1 3],[2 4]}; % row of above cond_idx to contrast
S.total_samples = -200:799;
S.select_samples = 0:600;
S.smooth_samples = 1;
S.dsample = 4;
S.analysis_type='multicomp'; % comp is considerably faster for LDA than comp_recon (latter reconstructs all channels)
S.ndec=8;

% non-parametric independent samples rank test
S.ranksum_on=0;

% TFCE
S.tfce_on = 0;
S.tfce_nperm=100; % for initial testing, between 100 and 1000. 5000 is most robust, but slow

% MVPA settings
S.mvpa_on=1;
S.nRandSamp = 1; % currently not functional
S.SL_type = 'time';
S.search_radius = Inf;
S.use_measure = 'crossvalidation';
S.balance_dataset_and_partitions =1;
S.parti='nchunks'; % 'take-one-out', 'splithalf', 'oddeven', 'nchunks'
%S.use_classifier = 'LDA'; S.regularization = 0.5; S.matlab_lda = 0; S.logist = 0; S.output_weights = 1;
S.use_classifier = 'GP'; 
S.use_chunks = 'balance_targets';
S.nchunks=4;
S.average_train_count = 1;
S.average_train_resamplings = 1;

sname=datestr(now,30);

% get file names
compfiles = dir(fullfile(S.icaDir,[S.compname{1} '*' S.compname{2}]));
if S.dsample
    S.total_samples = downsample(S.total_samples',S.dsample)';
    S.select_samples = downsample(S.select_samples',S.dsample)';
end
for f = 1:length(compfiles)
    currentFile = fullfile(S.icaDir,compfiles(f).name);
    currentInfo = fullfile(S.icaDir,strrep(strrep(compfiles(f).name,S.compname{1},S.infoname{1}),S.compname{2},S.infoname{2}));
    
    fprintf('\n');
    disp(['Loading ', currentFile]);
    disp(['Loading ', currentInfo]);

    load(currentFile);
    load(currentInfo);
    
%     fulldata= topography * timecourse; 
%     fulldata=reshape(fulldata,size(fulldata,1),[],n_trials);
%     ERP = mean(fulldata,3);

    switch S.analysis_type
        case {'comp','comp_recon'} % process each ICA component separately
            iter = 1:size(timecourse,1);
        case 'recon' % combine all ICs into a single reconstruction of the data
            iter = 1;
        case 'multicomp' % consider ICs as multivariate data rather than separately
            iter = 1;
    end
    
    for c=iter
        
        disp(['file ' num2str(f) '/' num2str(length(compfiles)) ' , comp ' num2str(c) ' /' num2str(size(timecourse,1))])      
        
        switch S.analysis_type
            case 'comp'
                comps = c;
                data= timecourse(comps,:);  
            case 'comp_recon'
                comps = c;
                data= topography(:,comps) * timecourse(comps,:); 
            case 'recon'
                comps = 1:size(timecourse,1);
                data= topography(:,comps) * timecourse(comps,:);  
            case 'multicomp'
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
            
        % select samples
        select = dsearchn(S.total_samples',[S.select_samples(1),S.select_samples(end)]');
        data = data(:,select(1):select(2),:);

        % create two sets of trials to contrast
        for con = 1:length(S.contrast_rows)
            idx{con} = find(ismember(eventType,[S.cond_idx{S.contrast_rows{con},:}]));
            conData{con} = data(:,:,idx{con});

            % reduce to Global Field Power
            gfpData{con} = squeeze(std(conData{con},{},1));
            stdgfpData{con} = std(gfpData{con},[],1);
        end

        % non-parametric independent paired test
        if S.ranksum_on
            for s = 1:length(S.select_samples)
                [p,~,st]=ranksum(gfpData{1}(s,:),gfpData{2}(s,:));
                stats.ranksum.all(f,c,s) = st.ranksum;
                stats.ranksum.pval(f,c,s) = p;
                stats.meandiff.gfp(f,c,s)=mean(gfpData{1}(s,:))-mean(gfpData{2}(s,:));
            end
            stats.ranksum.max(f,c) = max(stats.ranksum.all(f,c,:));
            stats.ranksum.min(f,c) = min(stats.ranksum.all(f,c,:));
            stats.ranksum.min_pval(f,c) = min(stats.ranksum.pval(f,c,:));
            stats.ranksum.stdgfp_pval(f,c)=ranksum(stdgfpData{1},stdgfpData{2});
            stats.meandiff.stdgfp(f,c)=mean(stdgfpData{1})-mean(stdgfpData{2});
        end
        
        % TFCE test
        if S.tfce_on
            shiftd=4-ndims(gfpData{1});
            img1=permute(shiftdim(gfpData{1},-shiftd),[3 1 2 4]);
            img2=permute(shiftdim(gfpData{2},-shiftd),[3 1 2 4]);
            pvals=matlab_tfce('independent',1,img1,img2,'nperm',tfce_nperm);
            stats.tfce(f,c) = min(pvals);
        end

        if S.mvpa_on
            % create cosmo data struct
            conds = nan(1,length(eventType));
            for cn = 1:length(idx)
                conds(idx{cn}) = cn;
            end
            cos = eeglab2cosmo(data,S.select_samples,conds);
            
            % set the targets 
            cos.sa.targets=cos.sa.trialinfo(:,1);

            % set the chunks (independent measurements)
            cos.sa.chunks=[1:size(cos.samples,1)]'; % each trial is a chunk - CORRECT?
            
            if S.balance_dataset_and_partitions 
                if c==1 % use same balancing for all components
                    S.balance_idx=[];
                end
            end
            
            % run analysis
            %try
                [out,S] = run_cosmo_machine(cos,S);
                stats.mvpa(f,c) = out;
                stats.mvpa_cv_acc(f,c) = mean(stats.mvpa(f,c).samples);
                disp('MVPA complete')
            %catch
            %    disp('MVPA failed!')
            %end
        end
    end
    save(fullfile(S.icaDir,['stats_' sname '.mat']),'stats'); % temporary - will be overwritten. Allows re-starting a failed analysis.
end

if S.ranksum_on
    figure
    hold on
    bar(1:c,mean(stats.ranksum.max,1))
    errorbar(1:c,mean(stats.ranksum.max,1),std(stats.ranksum.max,[],1),'.')
    figure
    hold on
    bar(1:c,mean(stats.ranksum.min,1))
    errorbar(1:c,mean(stats.ranksum.min,1),std(stats.ranksum.min,[],1),'.')
end
if S.tfce_on
    figure
    hold on
    bar(1:c,mean(stats.tfce,1))
    errorbar(1:c,mean(stats.tfce,1),std(stats.tfce,[],1),'.')
end
if S.mvpa_on
    figure
    imagesc(stats.mvpa_cv_acc,[0.4 0.6])
    colormap(parula(10))
    colorbar
    title([S.icaDir ' ' S.analysis_type])
    figure
    hold on
    bar(1:c,mean(stats.mvpa_cv_acc,1))
    errorbar(1:c,mean(stats.mvpa_cv_acc,1),std(stats.mvpa_cv_acc,[],1),'.')
    plot(xlim,[0.6 0.6], 'k--')
    ylim([0.5 0.6])
    figure
    hold on
    bar(1:f,max(stats.mvpa_cv_acc,[],2))
    plot(xlim,[0.6 0.6], 'k--')
    ylim([0.5 0.7])
end

% group-level stats
if S.ranksum_on
    for c = 1:size(timecourse,1)
        stats.signrank.stdgfp(c,1)=signrank(stats.meandiff.stdgfp(:,c));
    end
end