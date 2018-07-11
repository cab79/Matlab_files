clear all
dbstop if error
addpath(genpath('C:\Data\Matlab\TFCE'));
icaDir = 'C:\Data\CORE\eeg\ana\groupICA';
compname = {'CORE_part2_ica_c','-1.mat'};
infoname = {'fileinfo_','.mat'};
cond_idx = {
%     [1 2 9 10 17 18] %left hand, mismatch
%     [3 4 11 12 19 20] %left hand, standard
%     [5 6 13 14 21 22] %right hand, mismatch
%     [7 8 15 16 23 24]} %right hand, standard

    [1 2 9 10] %left hand, mismatch
    [3 4 11 12] %left hand, standard
    [5 6 13 14] %right hand, mismatch
    [7 8 15 16]} %right hand, standard
contrast_rows = {[1 3],[2 4]}; % row of above cond_idx to contrast
total_samples = -200:799;
select_samples = 0:600;
smooth_samples = 0;
dsample = 1;
analysis_type='multicomp'; % comp is considerably faster for LDA than comp_recon (latter reconstructs all channels)

% non-parametric independent samples rank test
ranksum_on=0;

% TFCE
tfce_on = 0;
tfce_nperm=100; % for initial testing, between 100 and 1000. 5000 is most robust, but slow

% LDA settings
lda_on=1;
nRandSamp = 1;
nfold=10;
ndec=12;

sname=datestr(now,30);

% get file names
compfiles = dir(fullfile(icaDir,[compname{1} '*' compname{2}]));
if dsample
    total_samples = downsample(total_samples',dsample)';
    select_samples = downsample(select_samples',dsample)';
end
for f = 1:length(compfiles)
    currentFile = fullfile(icaDir,compfiles(f).name);
    currentInfo = fullfile(icaDir,strrep(strrep(compfiles(f).name,compname{1},infoname{1}),compname{2},infoname{2}));
    
    fprintf('\n');
    disp(['Loading ', currentFile]);
    disp(['Loading ', currentInfo]);

    % Use EEGLAB function to load set file
    load(currentFile);
    load(currentInfo);
    
%     fulldata= topography * timecourse; 
%     fulldata=reshape(fulldata,size(fulldata,1),[],n_trials);
%     ERP = mean(fulldata,3);

    switch analysis_type
        case {'comp','comp_recon'} % process each ICA component separately
            iter = 1:size(timecourse,1);
        case 'all_recon' % combine all ICs into a single reconstruction of the data
            iter = 1;
        case 'multicomp' % consider ICs as multivariate data rather than separately
            iter = 1;
    end
    
    for c=iter
        
        disp(['file ' num2str(f) '/' num2str(length(compfiles)) ' , comp ' num2str(c) ' /' num2str(size(timecourse,1))])      
        
        switch analysis_type
            case 'comp'
                comps = c;
                data= timecourse(comps,:);  
            case 'comp_recon'
                comps = c;
                data= topography(:,comps) * timecourse(comps,:); 
            case 'all_recon'
                comps = 1:size(timecourse,1);
                data= topography(:,comps) * timecourse(comps,:);  
            case 'multicomp'
                data= timecourse;  
        end
        
        % smooth
        if smooth_samples
            disp('smoothing...')
            for i = 1:size(data,1)
                data(i,:) = smooth(data(i,:),smooth_samples,'moving');
            end
        end
        
        % downsample
        if dsample
            data = downsample(data',dsample)';
        end

        % make 3D
        data=reshape(data,size(data,1),[],n_trials);
            
        % select samples
        select = dsearchn(total_samples',[select_samples(1),select_samples(end)]');
        data = data(:,select(1):select(2),:);

        % create two sets of trials to contrast
        for con = 1:length(contrast_rows)
            idx = find(ismember(eventType,[cond_idx{contrast_rows{con},:}]));
            conData{con} = data(:,:,idx);

            % reduce to Global Field Power
            gfpData{con} = squeeze(std(conData{con},{},1));
            stdgfpData{con} = std(gfpData{con},[],1);
        end

        % non-parametric independent paired test
        if ranksum_on
            for s = 1:length(select_samples)
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
        if tfce_on
            shiftd=4-ndims(gfpData{1});
            img1=permute(shiftdim(gfpData{1},-shiftd),[3 1 2 4]);
            img2=permute(shiftdim(gfpData{2},-shiftd),[3 1 2 4]);
            pvals=matlab_tfce('independent',1,img1,img2,'nperm',tfce_nperm);
            stats.tfce(f,c) = min(pvals);
        end

        if lda_on
            % Linear Discriminant Analysis on balanced trials
            nT = cell2mat(cellfun(@size,conData,'UniformOutput',0)');
            [minT,minI] = min(nT(:,3));
            [maxT,maxI] = max(nT(:,3));

            randData = conData;
            for n = 1:nRandSamp
                disp(['rep ' num2str(n) ' /' num2str(nRandSamp)])
                randData{maxI} = conData{maxI}(:,:,randsample(maxT,minT));
                datmat = cat(2,reshape(randData{1},[],minT),reshape(randData{2},[],minT))';
                groups = [ones(minT,1);2*ones(minT,1)];
                try
                    out = lda_class(datmat,groups,nfold,ndec);
                catch
                    out.empty=1;
                    disp('LDA failed!')
                end
                if isfield(out,'empty')
                    continue
                else
                    stats.lda(f,c,n) = out;
                    disp('LDA complete')
                end
                stats.lda_cv_error(f,c) = mean(stats.lda(f,c).ldaCVErr);
                clear out datmat
            end
            clear data randData conData
        end
    end
    save(fullfile(icaDir,['stats_' sname '.mat']),'stats'); % temporary - will be overwritten. Allows re-starting a failed analysis.
end

if ranksum_on
    figure
    hold on
    bar(1:c,mean(stats.ranksum.max,1))
    errorbar(1:c,mean(stats.ranksum.max,1),std(stats.ranksum.max,[],1),'.')
    figure
    hold on
    bar(1:c,mean(stats.ranksum.min,1))
    errorbar(1:c,mean(stats.ranksum.min,1),std(stats.ranksum.min,[],1),'.')
end
if tfce_on
    figure
    hold on
    bar(1:c,mean(stats.tfce,1))
    errorbar(1:c,mean(stats.tfce,1),std(stats.tfce,[],1),'.')
end
if lda_on
    figure
    imagesc(stats.lda_cv_error,[0.4 0.6])
    colormap(parula(10))
    colorbar
end

% group-level stats
if ranksum_on
    for c = 1:size(timecourse,1)
        stats.signrank.stdgfp(c,1)=signrank(stats.meandiff.stdgfp(:,c));
    end
end