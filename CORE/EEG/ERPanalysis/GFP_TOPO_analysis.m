%---select data to analyse---%
%select = 'TF'; TFmethod = '-FT'; % or, '-EL'
select = 'ERP'; TFmethod = '';
savenametype = '_CNUM_4';
load([select savenametype '_data.mat'])
use_etype = 1:24;
eventtypes = cellfun(@num2str, num2cell(use_etype), 'UniformOutput', false);
no_cond = length(eventtypes);
%-------------------------%

%---settings---%
anatype = {'evoked', 'induced','itc'}; anaDAT=1; % 1:evoked, 2:induced, 3: ITC
singtrial=0;
GFPnormtype = {'nonorm', 'GFPnorm', 'GFPpriornorm1','GFPpriornorm2','GFPbase','GFPp1norm'}; GFPnorm=1;%
%avcond = {[1:2]}; avcondname = 'avcond'; % average all conditions together
avcond = {[1],[2]}; avcondname = 'eachcond'; % all conditions separately
min_trials_gavg =10; % minimum no. of trials per condition
output_ntrials=0;
data_max=0;
%--------------%

gavg = cell(length(subjects)*2,2); % for left and right separately
for g = 1:length(subjects)
    gavg{g,1} = zeros(el,length(times),nf);
    gavg{g,2} = 0;
    for s = 1:length(subjects{g,1}) 
        for i = 1:no_cond
            if DATs{g,1}{s,1}{i,4}>=min_trials_gavg
                gavg{g,1} = gavg{g,1} + DATs{g,1}{s,1}{i,anaDAT};
                gavg{g,2} = gavg{g,2}+1;
            end
        end
    end
    gavg{g,1} = gavg{g,1}/gavg{g,2};
end

%for i = 1:size(gavg,1)
%    gavg_all(:,:,:) = gavg{i,1};
%end
%gavg_all = nanmean(gavg_all,4);

% normalise by common average and GFP mean over time, for each frequency
for g = 1:length(subjects)
    for f = 1:size(gavg{g,1},3)
        for e = 1:size(gavg{g,1},1)
            gavg{g,1}(e,:,f) = gavg{g,1}(e,:,f) - nanmean(gavg{g,1}(:,:,f),1);
        end
        gavg{g,1}(:,:,f) = gavg{g,1}(:,:,f)/nanmean(std(gavg{g,1}(:,:,f),1));
    end
end

if strcmp(select,'ERP')
    
    no_freqs=1;
    freq_idx={1};
    freqs_nme = {select};
    freqs_limits=cell(1);
    times=EEG.times;
    posnegpeak = [];
    
    % smooth gavg for peak identification
    gavg_sm = gavg;
    for g = 1:length(subjects)
        for f = 1:size(gavg{g,1},3)
            for e = 1:size(gavg{g,1},1)
                gavg_sm{g,1}(e,:,f) = smooth(squeeze(gavg{g,1}(e,:,f)),30,'lowess');
            end
        end
    end
    
    %gavg_sm = gavg;
    
    % plot grand average
    for ii = 1:size(gavg_sm,1);
        plotdata = gavg_sm{ii,1};
        [maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
        [~, maxmaxidx] = max(maxval);
        plottime = times(maxidx(maxmaxidx));
        if plottime == times(end)
            plottime = times(end-1);
        end
        figure('Name','grand average');
        timtopo(plotdata,chanlocs,...
            'limits',[times(1) times(end)],...
            'plottimes',plottime);
        set(gcf,'Color','white');
    end
    
elseif strcmp(select,'TF')
    
    % define frequency limits for data reduction
    no_freqs=10;
    [freqs_nme,freqs_limits,freq_idx] = find_freq_limits(gavg,no_freqs);

    % plot grand-grand average
    close all
    for f = 1:length(freq_idx)
        chanlocs = EEG.chanlocs;
        plotdata = mean(gavg_all(:,:,freq_idx{f}),3);
        [maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
        [x, maxmaxidx] = max(maxval);
        plottime = times(maxidx(maxmaxidx));
        if plottime == times(end)
            plottime = times(end-1);
        end
        figure('Name',num2str(freqs(freq_idx{f})));
        timtopo(plotdata(:,:),chanlocs,'limits',[times(1) times(end)],'plottimes',plottime);
        set(gcf,'Color','white');
    end
    
    % smooth gavg for peak identification
    gavg_sm = gavg;
    for g = 1:length(subjects)
        for f = 1:size(gavg{g,1},3)
            for e = 1:size(gavg{g,1},1)
                gavg_sm{g,1}(e,:,f) = smooth(squeeze(gavg{g,1}(e,:,f)),30,'lowess');
            end
        end
    end
    figure;
    imagesc(times,freqs,squeeze(mean(gavg_sm{g,1}([6 62],:,:),1))'); 
    colormap jet;
    set(gca,'YDir','normal')
    
    
   
end

% data for each subject
%subdat_all=[];
%for s = 1:length(subjects)
%    for s2 = 1:length(subjects{s,1}) 
%        subdat = [];
%        for i = 1:size(DATs{s,1}{s2,1})
%            subdat(:,:,:,i) = DATs{s,1}{s2,1}{i,1};
%        end
%        subdat_all(:,:,:,s,s2) = permute(mean(subdat,4),[3 2 1]);
%    end
%end

% identify time bins
limits=cell(no_freqs,1);
limits_all=cell(no_freqs,1);
bins_all=cell(no_freqs,1);
basetimesGFP = intersect(find(times<basebin(2)*1000), find(times>=basebin(1)*1000));
stimtimesGFP = find(times>basebin(2)*1000);
timesRMSE = times(2:end);
basetimesRMSE = intersect(find(timesRMSE<basebin(2)*1000), find(timesRMSE>=basebin(1)*1000));
stimtimesRMSE = find(timesRMSE>basebin(2)*1000);

close all
for f = 1:no_freqs
    
    % GFP and CIs across subjects
    GFPg=[];
    GFP_sub=[];
    for g = 1:length(subjects)
        GFPg(:,g) = std(mean(gavg_sm{g,1}(:,:,freq_idx{f}),3),1);
        %for s2 = 1:length(subjects{s,1}) 
        %    GFP_sub(:,s,s2) = squeeze(std(squeeze(mean(subdat_all(:,:,freq_idx{f},s,s2),3)),1));
        %end
    end
    %GFP_CI=[];
    %for t = 1:length(times)
    %    GFPg(t,:) = mean(squeeze(GFP_sub(t,:,:)),2);
    %    GFP_CI(t,:) = 2.575*std(squeeze(GFP_sub(t,:,:))');
    %end
    GFPg=mean(GFPg,2);
    GFP_CI = repmat(mean(GFPg(basetimesGFP))+std(GFPg(basetimesGFP)),length(GFPg),1);
    GFP_base = repmat(mean(GFPg(basetimesGFP)),length(GFPg),1);
    figure
    yyaxis left
    plot(times,GFPg);
    hold on
    plot(times,GFP_CI,'b');
    hold on
    plot(times,GFP_base,'b--');
    
    % GFP peaks
    [peaks_GFP peaklocs_GFP] = findpeaks_Rbase(GFPg,'MinPeakHeight',unique(GFP_CI),'MinPeakProminence',unique(GFP_CI)-unique(GFP_base));
    peakloctimes_GFP = times(peaklocs_GFP);
    hold on;
    plot(peakloctimes_GFP,peaks_GFP,'v','MarkerFaceColor','b');
    % GFP troughs
    [tro_GFP trolocs_GFP] = findpeaks_Lbase(-GFPg,'MinPeakProminence',unique(GFP_CI)-unique(GFP_base));
    troloctimes_GFP = times(trolocs_GFP);
    tro_GFP = -tro_GFP;
    hold on;
    plot(troloctimes_GFP,tro_GFP,'s','MarkerFaceColor','b');
    
    % RMSE
    RMSE = [];
    RMSE_sub = [];
    for g = 1:length(subjects)
        for t = 2:size(gavg_sm{g,1},2)
            y = mean(gavg_sm{g,1}(:,t,freq_idx{f}),3);
            yhat = mean(gavg_sm{g,1}(:,t-1,freq_idx{f}),3);
            RMSE(t-1,g) = sqrt(mean((y - yhat).^2));
            %for s2 = 1:length(subjects{s,1}) 
            %    y = squeeze(mean(subdat_all(:,t,freq_idx{f},s,s2),3));
            %    yhat = squeeze(mean(subdat_all(:,t-1,freq_idx{f},s,s2),3));
            %    RMSE_sub(t-1,s,s2) = sqrt(mean((y - yhat).^2));
            %end
        end
    end
    %for t = 1:length(times)-1
    %    RMSE(t,:) = mean(squeeze(RMSE_sub(t,:,:)),2);
    %    RMSE_CI(t,:) = 2.575*std(squeeze(RMSE_sub(t,:,:))');
    %end
    RMSEav = mean(RMSE,2);
    RMSE_CI = repmat(mean(RMSEav(basetimesRMSE))+std(RMSEav(basetimesRMSE)),length(RMSEav),1);
    RMSE_base = repmat(mean(RMSEav(basetimesRMSE)),length(RMSEav),1);
    hold on
    yyaxis right
    plot(times(2:end),RMSEav,'r');
    hold on
    plot(times(2:end),RMSE_CI,'r');
    hold on
    plot(times(2:end),RMSE_base,'r--');
    
    % RMSE peaks
    [peaks_RMSE peaklocs_RMSE] = findpeaks_Rbase(RMSEav,'MinPeakHeight',unique(RMSE_CI),'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    peaklocs_RMSE = peaklocs_RMSE+1;
    hold on;
    plot(times(peaklocs_RMSE),peaks_RMSE,'v','MarkerFaceColor','r');
    % RMSE troughs
    [tro_RMSE trolocs_RMSE] = findpeaks_Lbase(-RMSEav,'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    trolocs_RMSE = trolocs_RMSE+1;
    tro_RMSE = -tro_RMSE;
    hold on;
    plot(times(trolocs_RMSE),tro_RMSE,'s','MarkerFaceColor','r');
    
    % find peak limits
    peaklocs_RMSE=peaklocs_RMSE(peaklocs_RMSE>stimtimesRMSE(1));
    limit1=[];
    limit2=[];
    for i = 1:length(peaklocs_RMSE)
        limit1=[limit1 peaklocs_RMSE(i)+1];
        if i<length(peaklocs_RMSE)
            next_tros = trolocs_RMSE(intersect(find(trolocs_RMSE>peaklocs_RMSE(i)),find(trolocs_RMSE<peaklocs_RMSE(i+1))));
        else
            next_tros = trolocs_RMSE(find(trolocs_RMSE>peaklocs_RMSE(i)));
        end
        if ~isempty(next_tros)
            limit2=[limit2 next_tros(end)];
        else
            limit1 = limit1(1:end-1);
            continue
        end
        limits_all{f}{i}=[limit1(end):limit2(end)];
    end
    
    limits_all{f}=limits_all{f}(~cellfun('isempty',limits_all{f}))  
    
    if isempty(limits_all{f})
        no_peaks_temp=0;
    else
        no_peaks_temp=length(limits_all{f});
    end
    
    %TOPOSIM
    if no_peaks_temp>0
        TOPO = [];
        TOPOSIM = [];
        for g = 1:length(subjects)
            for n = 1:no_peaks_temp
                TOPO(:,g,n) = mean(mean(gavg{g,1}(:,limits_all{f}{n},freq_idx{f}),3),2);
                %subplot(length(subjects),no_peaks,(s-1)*no_peaks+n); topoplot(TOPO(:,s,n), EEG.chanlocs);
            end
            TOPOSIM(g,:) = pdist(squeeze(TOPO(:,g,:))','cosine');
        end
        TOPOSIMav=mean(TOPOSIM,1);
        TOPOSIMav=squareform(TOPOSIMav);

        for thresh=0.2;
            test=TOPOSIMav.*(TOPOSIMav<thresh);
            test = tril(test);
            for t = 1:size(test,2)
                %tz = test(:,t)==0;
                %tz1 = tz(2:end)-tz(1:end-1);
                %tz2 = find(tz1==1);
                %if ~isempty(tz2)
                %    tz3 = tz2(1)+1;
                %    test(tz3:end,t)=0; 
                %end
                if t<size(test,2)-1
                    test(t+2:end,t)=0; % only consider adjacent connections
                end
            end
            test = test'+test;
            bins = conncomp(graph(test));
            test2 = bins(2:end)-bins(1:end-1);
            %if all(test2>=0)
            %    break
            %end
        end

        no_peaks=max(bins);
        for n = 1:no_peaks
            bin_idx = find(bins==n);
            no_temp_peaks = length(bin_idx);
            limits_temp_np = limits_all{f}(bin_idx);
            limits{f}{n}=[];
            for tp = 1:no_temp_peaks
                limits{f}{n}= [limits{f}{n} limits_temp_np{tp}];
            end
        end
        figure
        for g = 1:length(subjects)
            for n = 1:no_peaks
                subplot(length(subjects),no_peaks,(g-1)*no_peaks+n); topoplot(mean(mean(gavg{g,1}(:,limits{f}{n},freq_idx{f}),3),2), EEG.chanlocs);
            end
        end
        
        % plot grand average
        for ii = 1:size(gavg_sm,1);
            plotdata = gavg_sm{ii,1};
            figure('Name','grand average');
            plot(times,plotdata');
            cm=colormap(jet(no_peaks));
            hold on
            for p = 1:no_peaks
                binidx = find(bins==p);
                for n = 1:length(binidx)
                    line([times(limits_all{f}{binidx(n)}(1)) times(limits_all{f}{binidx(n)}(end))],[0 0],'color',cm(p,:),'LineWidth',8);
                end
            end
            hold off
        end
        
    else
        no_peaks=0;
    end
    bins_all{f}=bins;
end
timefreq_limits = struct;
timefreq_limits.freqs_limits = freqs_limits;
timefreq_limits.limits_all = limits_all;
timefreq_limits.bins = bins_all;
timefreq_limits.limits = limits;
save(['timefreq_limits_' select '_' anatype{anaDAT} savenametype],'timefreq_limits');
