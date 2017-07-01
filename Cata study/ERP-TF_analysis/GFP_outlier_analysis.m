%% uses GFP to identify subject outliers, allowing further data cleaning, removal, or flagging of data prior to parametric stats
% currently averages data over all conditions - can be modified in the
% future to consider each cell of the analysis (e.g. each event type)
% separately.
close all
clear all

%% -------SETTINGS------- %%
% path to cleaned EEGLAB data
Dpath = 'C:\Data\Catastrophising study\Preprocessed';
% name of file containing ERPs from all subjects, produced by "ERP_TF_analysis.m"
erpfilename = 'ERP_1_2_3_4_5_6_7_8_data.mat';
% list all event types in the file
eventtypes = {'c1','c2','c3','c4','c5','c6','c7','c8'}; 
% select which event types to include this this analysis
use_etype = [1 2 3 4 5 6 7 8]; 
% only include conditions if they have more than min_trials number of trials
min_trials =5; % minimum no. of trials per condition
% data window (peristimulus ms)
plotwin = [-5500 -2500];
% baseline correct
baseline = [-5500 -5000];


%% -------RUN------- %%
cd(Dpath)
load(erpfilename)

base = find(ismember(times,baseline));
ptimes = find(ismember(times,plotwin));

% put subject data into a structure array (FUTURE VERSION: RE-WRITE ERP_TF_ANALYSIS TO DO THIS) 
SUB = struct;
sg=0;
for g = 1:length(subjects)
    for s = 1:length(subjects{g,1})
        sg=sg+1;
        SUB(sg).ID = subjects{g,1}{s,1};
        SUB(sg).group = 'g';%['g' num2str(g)]; % OR as separate groups
        for c = 1:length(use_etype)
            SUB(sg).data(c).eventname = eventtypes{use_etype(c)};
            SUB(sg).data(c).ntrials = DATs{g,1}{s,1}{c,4};
            if SUB(sg).data(c).ntrials < min_trials
                SUB(sg).data(c).erp = NaN(el,length(times));
                SUB(sg).data(c).gfp = NaN(size(times));
            else
                % baseline correct
                tempdata = DATs{g,1}{s,1}{c,1};
                tempdata = tempdata-repmat(mean(tempdata(:,base(1):base(2)),2),1,size(tempdata,2));
                SUB(sg).data(c).erp = tempdata;
                SUB(sg).data(c).gfp = std(tempdata,[],1);
            end
        end
    end
end
   
% put stat data into a separate structure
STAT=struct;
grp = unique({SUB(:).group});
Ngrp = length(grp);
for g = 1:Ngrp
    STAT(g).name = grp{g};
    subind = find(strcmp({SUB(:).group},grp{g}));
    for c = 1:length(SUB(subind(1)).data) % for each event type
        STAT(g).data(c).eventname = SUB(subind(1)).data(c).eventname;
        tempdata = arrayfun(@(x) x.data(c).erp, SUB(subind), 'unif', false);
        tempdata = cat(3,tempdata{:});
        sizet = size(tempdata);
        STAT(g).data(c).gavg = nanmean(tempdata,3);
        STAT(g).data(c).gstd = nanstd(tempdata,1,3);
        STAT(g).data(c).zscore = nanzscore(tempdata,1,3);
        tempdata = arrayfun(@(x) x.data(c).gfp, SUB(subind), 'unif', false);
        tempdata = cat(1,tempdata{~cellfun(@isempty,tempdata)})';
        STAT(g).data(c).gfp_gavg = nanmean(tempdata,2);
        STAT(g).data(c).gfp_gstd = nanstd(tempdata,1,2);
        STAT(g).data(c).gfp_zscore = nanzscore(tempdata,1,2);
    end
    tempdata = arrayfun(@(x) mean(cat(3,x.data(:).erp),3), SUB(subind), 'unif', false);
    tempdata = cat(3,tempdata{:});
    STAT(g).gavg_all = nanmean(tempdata,3);
    STAT(g).std_all = nanstd(tempdata,1,3);
    STAT(g).zs_all = nanzscore(tempdata,1,3);
    tempdata = arrayfun(@(x) nanmean(cat(1,x.data(:).gfp),1), SUB(subind), 'unif', false);
    tempdata = cat(1,tempdata{~cellfun(@isempty,tempdata)})';
    STAT(g).gfp_gavg_all = nanmean(tempdata,2);
    STAT(g).gfp_std_all = nanstd(tempdata,1,2);
    STAT(g).gfp_zs_all = nanzscore(tempdata,1,2);
end

%plot grand averages for all electrodes
close all
for g = 1:Ngrp
    figure
    plot(times(ptimes(1):ptimes(2)),STAT(g).gavg_all(:,ptimes(1):ptimes(2)));
end

%plot zscores of GFPs averaged over conditions
close all
for g = 1:Ngrp
    figure
    plot(times(ptimes(1):ptimes(2)),STAT(g).gfp_zs_all(ptimes(1):ptimes(2),:));
end

% find most extreme subjects from 4 different metrics
sub_order={};
[~,~,sub_order{1}]=unique(nanmean(squeeze(max(abs(STAT(g).zs_all(:,ptimes(1):ptimes(2),:)),[],1)),1));
[~,~,sub_order{2}]=unique(max(squeeze(max(abs(STAT(g).zs_all(:,ptimes(1):ptimes(2),:)),[],1)),[],1));
[~,~,sub_order{3}]=unique(nanmean(abs(STAT(g).gfp_zs_all(ptimes(1):ptimes(2),:)),1));
[~,~,sub_order{4}]=unique(max(abs(STAT(g).gfp_zs_all(ptimes(1):ptimes(2),:)),[],1));
[~,sub_order{5}]=sort(mean(cell2mat(sub_order),2),'descend');
allsub={SUB(:).ID};
allsub(sub_order{5})

outliers={};
%plot outliers based on max zscores over electrodes of ERPs averaged over conditions
outZ=5;
outliers{1} = find(nanmean(squeeze(max(abs(STAT(g).zs_all(:,ptimes(1):ptimes(2),:)),[],1)),1)>outZ); % average over time > outZ
outliers(2) = find(any(squeeze(max(abs(STAT(g).zs_all(:,ptimes(1):ptimes(2),:)),[],1))>outZ)); % any data point over time > outZ
close all
for g = 1:Ngrp
    figure
    plot(times(ptimes(1):ptimes(2)),squeeze(max(abs(STAT(g).zs_all(:,ptimes(1):ptimes(2),outliers{1})),[],1)));
    leg = {SUB(outliers).ID};
    legend(leg)
end

%plot outliers from  zscores of GFPs averaged over conditions
outZ=1;
outliers{3} = find(nanmean(STAT(g).gfp_zs_all(ptimes(1):ptimes(2),:),1)>outZ); % average over time > outZ
outliers{4} = find(any(STAT(g).gfp_zs_all(ptimes(1):ptimes(2),:)>outZ)); % any data point over time > outZ
close all
for g = 1:Ngrp
    figure
    plot(times(ptimes(1):ptimes(2)),STAT(g).gfp_zs_all(ptimes(1):ptimes(2),outliers{3}));
    leg = {SUB(outliers).ID};
    legend(leg)
end

