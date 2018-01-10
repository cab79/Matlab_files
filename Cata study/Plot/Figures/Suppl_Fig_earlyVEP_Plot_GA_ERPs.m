clear all
close all

load('C:\Data\Catastrophising study\Preprocessed\ERP_1_3_5_7_a_a_a_a_b_b_b_b_data_SPN.mat')
% load .xlsx file containing 'Participant_ID', 'Group', and covariates
pdatfile = 'C:\Data\Catastrophising study\Behavioural\Participant_data_nocodes.xlsx';
load('C:\Data\Catastrophising study\Orig\chanlocs.mat');

min_trials_gavg=3;
anatype = {'evoked', 'induced','itc'}; anaDAT=1; % 1:evoked, 2:induced, 3: ITC
lat_ms = [-2266]; %ms
sigtimes =[-2346 -2266]; % grey background at cluster extent. Can use find_cluster_extent.m
%sigtimes = [104 131]; %268ms peak
%sigtimes = [80 95]; %132ms peak
%sigtimes = [79 96]; %124ms peak
no_plot_ele=[];
pos_neg_peak = 1; %1 or -1
num_ele = 3;
lat = (lat_ms+3000)/2;
no_cond=9;
av_grps = {2,3};
plotidx=[1 2];
col = {'r','b'};
%col = {'b','b--','r','r--'};
plottimewin = [-3000 -1500];

% read participant data
[~,~,pdata] = xlsread(pdatfile);
grp_col = find(strcmp(pdata(1,:),'Group'));
sub_col = find(strcmp(pdata(1,:),'Subject'));
inc_col = find(strcmp(pdata(1,:),'Include'));
Ngrp = max([pdata{2:end,grp_col}]);
SubInd = cell(Ngrp,1);
Subs = [];
for g = 1:Ngrp
    inc_idx = find(cellfun(@(x) x==1, pdata(2:end,inc_col), 'UniformOutput', 1));
    grp_idx = find(cellfun(@(x) x==g, pdata(2:end,grp_col), 'UniformOutput', 1));
    SubInd{g,1} = intersect(inc_idx,grp_idx);
    Nsub(g,1) = length(SubInd{g,1});
end

% plot grand average
gavg = cell(4,1);
for s = 1:length(subjects)
    tempgavg = zeros(el,length(times));
    for s2 = 1:length(subjects{s,1}) 
        for i = use_etype
            if DATs{s,1}{s2,1}{i,6}>=min_trials_gavg
                tempgavg = cat(3,tempgavg,DATs{s,1}{s2,1}{i,anaDAT});
            end
        end
        tempgavg = mean(tempgavg,3);
        gavg{s,1} = cat(3,gavg{s,1},tempgavg);
    end
end

% average over groups
gavgplot = cell(length(av_grps),2);
for s = 1:length(av_grps)
    for g = 1:length(av_grps{s})
        gavgplot{s,1} = cat(4,gavgplot{s,1},gavg{av_grps{s}(g),1});
    end
    gavgplot{s,1} = mean(gavgplot{s,1},4);
    gavgplot{s,2} = squeeze(std(permute(gavgplot{s,1},[3 1 2])));
    gavgplot{s,1} = mean(gavgplot{s,1},3);
end

%% plot ERPs

% identify peak electrodes
mdata=[];
for ii = plotidx; 
    mdata = cat(3,mdata,gavgplot{ii,1});
end
mdata = mean(mdata,4);
ERP = mean(mdata(:,lat),2);
if pos_neg_peak==1
    [~,tps] = sort(ERP,'descend');
    tp = tps(1:num_ele);
elseif pos_neg_peak==-1
    [~,tps] = sort(ERP,'ascend');
    tp = tps(1:num_ele);
end
    
%% plot topo
figure
for ii = 1:length(gavgplot);
    peakdata = gavgplot{ii,1};
    plotchans=1:length(chanlocs);
    plotchans(no_plot_ele)=[];
    [~,markchans] = intersect(plotchans,tp);
    subplot(4,1,ii); topoplot(mean(peakdata(:,lat),2), chanlocs,'maplimits',[-4 4],'electrodes','on','plotchans',plotchans,'emarker2',{markchans,'o','w',7,1}); colorbar
end

% plot data
figure
ptimes = find(((EEG.times >= plottimewin(1)).*(EEG.times <= plottimewin(2)))==1);
ts=EEG.times(ptimes);
hold all
st=sigtimes;
fill([st, fliplr(st)], [ones(1,length(st))*10, fliplr(ones(1,length(st))*-10)], ...
[0.5 0.5 0.5], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
plot([lat_ms lat_ms],[-10 10],'k');
lowest=Inf;
highest=-Inf;
for ii = plotidx; 
    peakdata = gavgplot{ii,1};
    vardata = gavgplot{ii,2};
    ERP = mean(peakdata(tp,:),1);
    VAR = mean(vardata(tp,:),1);
    nsub = Nsub(ii); 
    SEM = VAR/sqrt(nsub);               % Standard Error
    tscore = -tinv(0.025,nsub-1);      % T-Score
    CI = tscore*SEM;                      % Confidence Intervals
    upper = ERP(ptimes)+CI(ptimes);
    lower = ERP(ptimes)-CI(ptimes);
    lowest = min(min(lower),lowest);
    highest = max(max(upper),highest);
    fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
    col{ii}, 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    plot(ts,ERP(ptimes),col{ii}); 
end
%legend('Healthy controls', 'healthy unaff', 'CRPS patients','patient unaff')
%legend('Healthy controls', 'CRPS patients','patient unaff');
ylabel('Amplitude, uV') % label for y axis
xlabel('Time (ms)') % label for x axis
set(gca,'FontSize',15);
set(findall(gcf,'type','text'),'FontSize',15);
ylim([lowest highest]);

%figure
%upper = ERP(ptimes)+1;
%lower=ERP(ptimes)-1;
%ts=EEG.times(ptimes);
%fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
%'r', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
%hold on
%plot(EEG.times(ptimes),ERP(ptimes),col{ii});
