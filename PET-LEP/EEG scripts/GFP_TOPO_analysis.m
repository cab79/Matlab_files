%---select data to load---%
%select = 'TF'; TFmethod = '-FT'; % or, '-EL'
close all
clear all
select = 'ERP'; TFmethod = '';
eventtypes = {'S  1','S  2'};
hand_nme={'right'};
use_etype = 1:2;
savenametype='';
for i = use_etype
    savenametype = [savenametype '_' eventtypes{1,i}(end)];
end
load([select TFmethod savenametype '_data.mat'])

%---settings---%
anatype = {'evoked', 'induced','itc'}; anaDAT=1; % 1:evoked, 2:induced, 3: ITC
handtype = {'R'}; hand_ana=1; % 1:LR, 2:AU
singtrial=0;
GFPnormtype = {'nonorm', 'GFPnorm', 'GFPpriornorm1','GFPpriornorm2','GFPbase','GFPp1norm'}; GFPnorm=1;%
%avcond = {[1:2]}; avcondname = 'avcond'; % average all conditions together
avcond = {[1],[2]}; avcondname = 'eachcond'; % keep conditions separate
min_trials_gavg =10; % minimum no. of trials per condition
output_ntrials=0;
data_max=0;
%--------------%

gavg = cell(length(subjects),2);
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
    no_freqs=8;
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

        for thresh=0;
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

% identify GFPs for each subject and condition
%LATgfp = NaN('double');
%GFP = NaN('double');
ugrplist = unique(sublist_grp','stable');
grplistall = num2cell(reshape(repmat(1:length(ugrplist),size(sublist,2),1),size(sublist,2)*length(ugrplist),1));
usublist = unique(sublist','stable');
usub = reshape(usublist,length(usublist)/length(ugrplist),length(ugrplist));
GFP = struct('subject',usublist,'group',grplistall);
%LAT = GFP;
POW = GFP;
GFPvar = GFP;
TOPdist = GFP;
TOPTdist = GFP;
TOPdistZ = GFP;
TOPdistC = GFP;
TOPind = GFP;
TOPindorig = GFP;
TOPmean = GFP;
TOPstd = GFP;
TOPmeanorig = GFP;
TOPstdorig = GFP;
TOPit = GFP;
W = [];

elec_radius = 0.35;
elec_select = [EEG.chanlocs.radius]<elec_radius;
topstd_thresh = 2;
meansim = 0.001;
it_max = 8;

if singtrial
    % Single trial similarity to the condition mean for selecting trials above
    % a lower thresh of std from that mean
    for g = 1:length(ugrplist)
        for s = 1:size(sublist,2) 
            s_name = usub{s,g};
            for f = 1:no_freqs
                f_name = freqs_nme{f};
                no_peaks = size(limits{f},2);
                for p = 1:no_peaks
                    p_name = ['peak' num2str(p)];
                    lim = limits{f}(:,p); % in datapoints
                    for h = 1:length(hands)
                        h_name = hands{h};
                        load([s_name 'singdat_' h_name '.mat']);
                        for ac = 1:length(avcond)
                            top=[];
                            c_name = ['c' num2str(ac)];
                            SDAT=[];
                            for i = 1:length(avcond{ac}) % for selected conditions
                               SDAT(:,:,:,:,i) = SINGDAT{avcond{ac}(i),1};
                            end
                            if size(SDAT,5)>1
                                SDAT = permute(SDAT,[1 2 3 5 4]);
                                SDAT = reshape(SDAT,size(SDAT,1),size(SDAT,2),size(SDAT,3),size(SDAT,4)*size(SDAT,5));
                                SDAT = permute(SDAT,[1 2 4 3]);
                            end
                            ntrial = size(SDAT,3);
                            iThresh = find(ones(ntrial,1));
                            iThresh_prev = find(ones(ntrial,1));
                            meantop = zeros(length(find(elec_select)),1);
                            meantod = 0;
                            stdtod = 0;
                            for i = 1:it_max
                                trialind = zeros(ntrial,1);
                                trialind(iThresh_prev(iThresh)) = ones(length(iThresh),1);
                                iTlen = length(iThresh);
                                meantop_prev = meantop;
                                alltop = squeeze(mean(mean(SDAT(elec_select,find(times==lim(1)):find(times==lim(2)),find(trialind),freq_idx{f}),4),2));
                                meantop = mean(alltop,2);
                                meandist = pdist([meantop meantop_prev]','cosine');
                                if meandist<meansim
                                    break
                                end
                                top(:,2:iTlen+1) = alltop;
                                top(:,1) = mean(top(:,2:end),2);
                                topdist = pdist(top','cosine');
                                tdfmean=topdist(1:iTlen); % topo dist from mean topo
                                meantod_prev = meantod;
                                stdtod_prev = stdtod;
                                meantod = mean(tdfmean);
                                stdtod = std(tdfmean);
                                if i==1
                                    orig_stdtod = stdtod;
                                    orig_meantod = meantod;
                                    orig_ind = find(trialind);
                                end
                                if orig_stdtod<stdtod || orig_meantod<meantod
                                    meantod =meantod_prev;
                                    stdtod =stdtod_prev;
                                    break
                                end
                                iThresh_prev=iThresh;
                                iThresh = find(tdfmean<meantod+(stdtod*topstd_thresh)); % index of above thresh trials
                            end
                            trialind = zeros(ntrial,1);
                            trialind(iThresh_prev(iThresh)) = ones(length(iThresh),1);
                            fieldname = [f_name '_' p_name '_' h_name '_' c_name];
                            TOPind(strcmp(s_name,{TOPind.subject})).(fieldname) = find(trialind); % trial indices of kept trials (topographically above std thresh)
                            TOPindorig(strcmp(s_name,{TOPindorig.subject})).(fieldname) = orig_ind; % trial indices of kept trials (topographically above std thresh)
                            TOPmean(strcmp(s_name,{TOPmean.subject})).(fieldname) = meantod;% meantod
                            TOPmeanorig(strcmp(s_name,{TOPmeanorig.subject})).(fieldname) = orig_meantod;% meantod
                            TOPstd(strcmp(s_name,{TOPstd.subject})).(fieldname) = stdtod;% meantod
                            TOPstdorig(strcmp(s_name,{TOPstdorig.subject})).(fieldname) = orig_stdtod;% meantod
                            TOPit(strcmp(s_name,{TOPit.subject})).(fieldname) = i;% meantod
                        end
                    end
                end
            end
        end
    end
    T_temp = struct2table(TOPstdorig);
    writetable(T_temp,['TOPstdorig_results_' select '_' anatype{anaDAT} savenametype '.xlsx']);
    T_temp = struct2table(TOPmeanorig);
    writetable(T_temp,['TOPmeanorig_results_' select '_' anatype{anaDAT} savenametype '.xlsx']);
end
   


for g = 1:length(ugrplist)
    for s = 1:size(sublist,2) 
        s_name = usub{s,g};
        for f = 1:no_freqs
            f_name = freqs_nme{f};
            no_peaks = length(limits{f});
            for p = 1:no_peaks
                p_name = ['peak' num2str(p)];
                hands = hand_nme;
                for h = 1:length(hands)
                    h_name = hands{h};
                    top=[];
                    topt=[];
                    GFPrvar=[];
                    for ac = 1:length(avcond)
                        c_name = ['c' num2str(ac)];
                        DAT = [];
                        ntrials=[];
                        for i = 1:length(avcond{ac}) % for selected conditions
                            if ~isempty(DATs{(g-1)*2+1,1}{s,1}{avcond{ac}(i),1})
                                DAT(:,:,:,i) = DATs{(g-1)*2+1,1}{s,1}{avcond{ac}(i),1};
                                ntrials(i) = DATs{(g-1)*2+1,1}{s,1}{avcond{ac}(i),4};
                            else
                                DAT(:,:,:,i) = NaN;
                            end
                        end
                        ntrials =sum(ntrials);
                        if data_max
                            [~,maxcond] = max(mean(squeeze(std(squeeze(mean(DAT(:,:,freq_idx{f},:),3)),1))));
                            DAT = DAT(:,:,:,maxcond);
                        else
                            DAT = nanmean(DAT,4);
                        end
                        if isnan(DAT); continue; end;
                        
                        GFPr = std(mean(DAT(:,limits{f}{p},freq_idx{f}),3),1);
                        [~, GFPmaxLOC] = max(GFPr);
                        GFPr = sort(GFPr(:), 'descend');
                        if GFPnorm==2
                            GFPr=GFPr/mean(squeeze(std(mean(DAT(:,:,freq_idx{f}),3),1)));
                        elseif GFPnorm==3
                            GFPr=GFPr/mean(squeeze(std(mean(DAT(:,1:limits{f}{p}(1)-1,freq_idx{f}),3),1)));
                        elseif GFPnorm==4
                            if p==1; 
                                GFPr=GFPr/mean(squeeze(std(mean(DAT(:,1:limits{f}{p}(1)-1,freq_idx{f}),3),1)));
                            else
                                GFPr=GFPr/mean(squeeze(std(mean(DAT(:,limits{f}{p-1},freq_idx{f}),3),1)));
                            end
                        elseif GFPnorm==5
                            GFPr=GFPr-mean(squeeze(std(mean(DAT(:,basetimesGFP,freq_idx{f}),3),1)));
                        elseif GFPnorm==6
                            GFPr=GFPr/mean(squeeze(std(mean(DAT(:,1:limits{f}{1},freq_idx{f}),3),1)));
                        end

                        POWr = mean(DAT(:,limits{f}{p},freq_idx{f}),3);
                        POWr = sort(POWr(:), 'descend');
                        %if posnegpeak(p)==1
                        %    POWr = sort(POWr(:), 'descend');
                        %elseif posnegpeak(p)==2
                        %    POWr = sort(POWr(:), 'ascend');
                        %end
                        if GFPnorm==2
                            POWr=POWr/mean(squeeze(std(mean(DAT(:,:,freq_idx{f}),3),1)));
                        elseif GFPnorm==3
                            POWr=POWr/mean(squeeze(std(mean(DAT(:,1:limits{f}{p}(1)-1,freq_idx{f}),3),1)));
                        elseif GFPnorm==4
                            if p==1; 
                                POWr=POWr/mean(squeeze(std(mean(DAT(:,1:limits{f}{p}(1)-1,freq_idx{f}),3),1)));
                            else
                                POWr=POWr/mean(squeeze(std(mean(DAT(:,limits{f}{p-1},freq_idx{f}),3),1)));
                            end
                        elseif GFPnorm==5
                            POWr=POWr-mean(squeeze(std(mean(DAT(:,basetimesGFP,freq_idx{f}),3),1)));
                        elseif GFPnorm==6
                            POWr=POWr/mean(squeeze(std(mean(DAT(:,1:limits{f}{1},freq_idx{f}),3),1)));
                        end

                        fieldname = [f_name '_' p_name '_' h_name '_' c_name];
                        percent_mean = 1;
                        GFP(strcmp(s_name,{GFP.subject})).(fieldname) = mean(GFPr(1:ceil(length(GFPr)*percent_mean)));%mean(GFPr);
                        %LAT(strcmp(s_name,{LAT.subject})).(fieldname) = times(find(times==lim(1))+GFPmaxLOC-1);
                        POW(strcmp(s_name,{POW.subject})).(fieldname) = mean(POWr(1:floor(length(POWr)*percent_mean)));
                        
                        %iThresh = TOPind(strcmp(s_name,{TOPind.subject})).(fieldname);
                        topdat=mean(mean(DAT(:,limits{f}{p},freq_idx{f}),3),2);
                       
                        top(:,ac)=topdat(elec_select);
                        
                        %topt(:,:) = mean(DAT(elec_select,find(times==lim(1)):find(times==lim(2)),freq_idx{f}),3);
                        %TOPTdist(strcmp(s_name,{TOPdist.subject})).(fieldname) = pdist(topt','cosine');
                        
                        GFPrvar(ac) = GFP(strcmp(s_name,{GFP.subject})).(fieldname);
                        
                        if output_ntrials
                            fieldname_trials = [h_name '_' c_name '_ntrials'];
                            GFP(strcmp(s_name,{GFP.subject})).(fieldname_trials) = ntrials;%mean(GFPr);
                        end
                        
                        if g==1 && s==1; W = [W; f,p,h,ac];end
                    end
                    fieldname = [f_name '_' p_name '_' h_name];
                    TOPdist(strcmp(s_name,{TOPdist.subject})).(fieldname) = pdist(top','cosine');
                    GFPvar(strcmp(s_name,{GFPvar.subject})).(fieldname) = std(GFPrvar);
                    
                end
            end
        end
    end
end

T_GFP = struct2table(GFP);
T_GFPvar = struct2table(GFPvar);
%T_LAT = struct2table(LAT);
T_POW = struct2table(POW);

writetable(T_GFP,['gfp_results_' select '_' anatype{anaDAT} savenametype '_' avcondname '_' GFPnormtype{GFPnorm} '.xlsx']);
%writetable(T_GFPvar,['gfp_var_results_' select '_' anatype{anaDAT} savenametype '_' handtype{hand_ana} '_' avcondname '_' GFPnormtype{GFPnorm} '.xlsx']);
%writetable(T_LAT,['lat_results_' select '_' anatype{anaDAT} savenametype '_' handtype{hand_ana} '_' avcondname '_' GFPnormtype{GFPnorm} '.xlsx']);
%writetable(T_POW,['pow_results_' select '_' anatype{anaDAT} savenametype '_' handtype{hand_ana} '_' avcondname '_' GFPnormtype{GFPnorm} '.xlsx']);

% modify table for RM ANOVA 
%T_GFP.group = categorical(T_GFP.group);
%T_LAT.group = categorical(T_LAT.group);
%T_POW.group = categorical(T_POW.group);
%ys = {};
%for ysi = 1:length(T_GFP.Properties.VariableNames(3:end))
%    ys{ysi} = ['y' num2str(ysi)];
%end
%T_GFP.Properties.VariableNames(3:end) = ys;
%T_LAT.Properties.VariableNames(3:end) = ys;
%T_POW.Properties.VariableNames(3:end) = ys;


% specify and fit model
%w1u = unique(W(:,1))';
%for f = w1u
%    w1Ind = find(W(:,1)==f);
%    w2u = unique(W(w1Ind,2))';
%    for p = w2u
%        w2Ind = find(W(:,1)==f & W(:,2)==p);
%        within = array2table(W(w2Ind,3:4));
%        between = T_GFP(:,[1:2 2+w2Ind']);
%        modelspec = [ys{w2Ind(1)} '-' ys{w2Ind(end)} ' ~ group'];
%        GFP_rm = fitrm(between,modelspec,'WithinDesign',within);
%        [ranovatbl] = ranova(GFP_rm,'WithinModel','Var1*Var2');
%        if any(ranovatbl.pValue([2 5 8 11])<0.05)
%            indx = find(ranovatbl.pValue([2 5 8 11])<0.05);
%            sprintf('%d %d %d',f,p,indx);
%        end
%    end
%end

% analyse topo cosine distance
top_dist_ana = {[1 5 8 9], [3 4 7]}; % 1-digit differences vs. 3/4 digit differences

for f = 1:no_freqs
    f_name = freqs_nme{f};
    no_peaks = size(limits{f},2);
    for p = 1:no_peaks
        p_name = ['peak' num2str(p)];
        lim = limits{f}(:,p); % in datapoints
        for h = 1:length(hand_nme)
            h_name = hand_nme{h};
            fieldname = [f_name '_' p_name '_' h_name];
            for n = 1:26
                topnum = TOPdist(n).(fieldname);
                meanall = mean(topnum);
                stdall = std(topnum);
                TOPdistZ(n).(fieldname) = mean(topnum(top_dist_ana{2}))-mean(topnum(top_dist_ana{1}));%(mean(topnum(top_dist_ana{2}))-meanall)/std(topnum) - (mean(topnum(top_dist_ana{1}))-meanall)/std(topnum); 
                for c = 1:length(top_dist_ana)
                    c_name = ['c' num2str(c)];
                    fieldnameC = [f_name '_' p_name '_' h_name '_' c_name];
                    TOPdistC(n).(fieldnameC) = mean(topnum(top_dist_ana{c}));%(mean(topnum(top_dist_ana{c}))-meanall)/std(topnum);
                end
            end
            for n = 1:13
                topmatH(n,:) = TOPdist(n).(fieldname);
            end
            for n = 14:26
                topmatP(n-13,:) = TOPdist(n).(fieldname);
            end
            H = squareform(mean(topmatH,1))
            P = squareform(mean(topmatP,1))
        end
    end
end
T = struct2table(TOPdistZ);
writetable(T,['topdistZ_results_' select '_' anatype{anaDAT} savenametype '.xlsx']);
T = struct2table(TOPdistC);
writetable(T,['topdistC_results_' select '_' anatype{anaDAT} savenametype '.xlsx']);

%Fn = fieldnames(GFP);
%Fn = Fn(3:2:end-1);
%GFPc = struct2cell(GFP)'; 
%GFP_L=cell2mat(GFPc(:,3:2:end-1));
%GFP_R=cell2mat(GFPc(:,4:2:end));
%GFP_LL=combvec(GFP_L,GFP_L);
%GFP_RR=combvec(GFP_R,GFP_R);
%GFP_LL = GFP_LL(1:26,:)./GFP_LL(27:52,:);
%GFP_RR = GFP_RR(1:26,:)./GFP_RR(27:52,:);
%xlswrite('GFP_LL.xls',GFP_LL);
%xlswrite('GFP_RR.xls',GFP_RR);
%Fn_ind = combvec(1:length(Fn),1:length(Fn));
%ii = [8 9 11 105 106]; %exp1
%%ii = [84 104 105 111 112 172 176 190]; %exp2
%for i = 1:length(ii)
%    one = Fn(Fn_ind(1,ii(i)))
%    two = Fn(Fn_ind(2,ii(i)))
%end
