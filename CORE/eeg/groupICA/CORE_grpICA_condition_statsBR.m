clear all
dbstop if error
addpath(genpath('C:\Data\Matlab\TFCE'));
icaDir = 'C:\Data\CORE\eeg\ana\groupICA';
compname = {'CORE_part2_ica_br','.mat'};
infoname = {'fileinfo_','.mat'};
cond_idx = {
    [1 2 9 10 17 18] %left hand, mismatch
    [3 4 11 12 19 20] %left hand, standard
    [5 6 13 14 21 22] %right hand, mismatch
    [7 8 15 16 23 24]} %right hand, standard
contrast_rows = {[1 3],[2 4]}; % row of above cond_idx to contrast
total_samples = -200:799;
select_samples = 0:600;
smooth_samples = 0;
dsample = 4;

% TFCE
tfce_on = 0;
tfce_nperm=100; % no difference between 100 and 500. 5000 is most robust

% LDA settings
lda_on=0;
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
    
    timecourse=compSet.timecourse;
    topography=compSet.topography;
    
    for c=1:size(timecourse,1)
        
        disp(['file ' num2str(f) '/' num2str(length(compfiles)) ' , comp ' num2str(c) ' /' num2str(size(timecourse,1))])
                
        % reconstruct
        data= topography(:,c) * timecourse(c,:);  
        
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
        end

        % non-parametric independent paired test
        for s = 1:length(select_samples)
            [p,~,st]=ranksum(gfpData{1}(s,:),gfpData{2}(s,:));
            stats.ranksum.all(f,c,s) = st.ranksum;
        end
        stats.ranksum.max(f,c) = max(stats.ranksum.all(f,c,:));
        stats.ranksum.min(f,c) = min(stats.ranksum.all(f,c,:));
        
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

figure
hold on
bar(1:c,mean(stats.ranksum.max,1))
errorbar(1:c,mean(stats.ranksum.max,1),std(stats.ranksum.max,[],1),'.')
figure
hold on
bar(1:c,mean(stats.ranksum.min,1))
errorbar(1:c,mean(stats.ranksum.min,1),std(stats.ranksum.min,[],1),'.')
if tfce_on
    figure
    hold on
    bar(1:c,mean(stats.tfce,1))
    errorbar(1:c,mean(stats.tfce,1),std(stats.tfce,[],1),'.')
end