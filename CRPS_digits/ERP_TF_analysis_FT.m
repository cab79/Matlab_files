clear all
close all

select = 'TF';
%select = 'ERP';

grplist = [1 2 29 30]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\';%Exp2
%grplist = [35 36 37 38]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\';%Exp1 left v right
%grplist = [51 52 53 54]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\';%Exp1 left v right


run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
cd(filepath)
ele_left = [];
ele_right = [];
%ele_left = [29 23 24 28 30 33 34];
%ele_right = [78 66 70 77 79 83 84];

aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1];
%aff_side = [1 1 1 1 1 1 1 1 1 1 1 1 1];

cond_nme = {'little','ring','middle','index','thumb'};
hand_nme = {'L','R'};
el = 92;
no_cond = length(cond_nme); % no of conditions per data file (arm)
subjects = subjlists(grplist);
DATs = subjects;

if strcmp(select,'ERP')
    dp = 250;
    nf = 1;
elseif strcmp(select,'TF')
    dp = 251;
    freqs = [4:2:40];
    nf = length(freqs);
    times = -0.2:0.004:0.8;
end
    
gavgL1 = zeros(el,dp,nf);
gavgR1 = zeros(el,dp,nf);
gavgL2 = zeros(el,dp,nf);
gavgR2 = zeros(el,dp,nf);
countavgL1 = 0;
countavgR1 = 0;
countavgL2 = 0;
countavgR2 = 0;

% create DATs and grand average
sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        sublist{s,s2} = subj(1:3);
        EEG = pop_loadset([subj '.set'],filepath);
        chanlocs = EEG.chanlocs;
        dsize = size(EEG.data);
        EEG.data = reshape(EEG.data,dsize(1),dsize(2),dsize(3));
        fevents = zeros(1,size(EEG.data,3));
        cevents = zeros(1,size(EEG.data,3));
        %j=0;
        for i = 1:size(EEG.data,3)
            if sum(strcmp(EEG.epoch(i).eventtype(1,:),'STIM'))>0
                ind = find(strcmp(EEG.epoch(i).eventtype(1,:),'STIM'));
                findfnum=[];
                %findcnum=[];
                if iscell(EEG.epoch(i).eventinit_index)
                    findfnum = find(strcmp('FNUM',EEG.epoch(i).eventcodes{1,ind(1)}(:,1)));
                    %findcnum = find(strcmp('CNUM',EEG.epoch(i).eventcodes{1,ind(1)}(:,1)));
                    fevents(1,i) = EEG.epoch(i).eventcodes{1,ind(1)}{findfnum(end),2};
                    %cevents(1,i) = EEG.epoch(i).eventcodes{1,ind(1)}{findcnum(end),2};
                else 
                    findfnum = find(strcmp('FNUM',EEG.epoch(i).eventcodes(:,ind(1))));
                    %findcnum = find(strcmp('CNUM',EEG.epoch(i).eventcodes(:,ind(1)))); 
                    fevents(1,i) = EEG.epoch(i).eventcodes{findfnum(end),2};
                    %cevents(1,i) = EEG.epoch(i).eventcodes{findcnum(end),2};
                end
            end
        end
        
        cond = min(fevents):max(fevents);
        DATs{s,1}{s2,1} = cell(length(cond),1);
        if strcmp(select,'TF')
            EEG=eeglab2fieldtrip(EEG,'preprocessing','none');
            cfg = [];
            cfg.output     = 'pow';
            cfg.foi     = freqs;
            %cfg.method     = 'mtmconvol';
            %cfg.taper       = 'hanning';
            %cfg.t_ftimwin    = 1./cfg.foi;  % cycles per time window
            cfg.method     = 'wavelet';                
            cfg.width      = 1;
            cfg.toi          = times;
            %cfg.pad        = 'nextpow2';
            cfg.keeptrials = 'yes';
            %cfg.precision = 'single';
            EEG = ft_freqanalysis(cfg,EEG); % returns power (auto-spectra) and shared activity (cross-spectra) between all channels
            eegdatasing = permute(EEG.powspctrm,[2 4 1 3]);
        end
        for i = 1:no_cond
            eegdata = [];
            if strcmp(select,'ERP')
                eegdata = double(squeeze(mean(EEG.data(:,:,find(fevents==cond(i))),3)));
            elseif strcmp(select,'TF')
                eegdata = squeeze(nanmean(eegdatasing(:,:,find(fevents==cond(i)),:),3));
            end
            if strcmp(sublist_side{s},'L') && s==1
                ord = i;
                DATs{s,1}{s2,1}{ord,1} = eegdata;
                %DATs{s,1}{s2,1}{ord,2} = eegdatasing;
                gavgL1 = gavgL1 + DATs{s,1}{s2,1}{ord,1};
                countavgL1 = countavgL1+1;
            elseif strcmp(sublist_side{s},'L') && s==3
                ord = i;
                DATs{s,1}{s2,1}{ord,1} = eegdata;
                %DATs{s,1}{s2,1}{ord,2} = eegdatasing;
                gavgL2 = gavgL2 + DATs{s,1}{s2,1}{ord,1};
                countavgL2 = countavgL2+1;
            elseif  strcmp(sublist_side{s},'R') && s==2
                ord = (no_cond+1)-i;
                DATs{s,1}{s2,1}{ord,1} = eegdata;
                %DATs{s,1}{s2,1}{ord,2} = eegdatasing;
                gavgR1 = gavgR1 + DATs{s,1}{s2,1}{ord,1};
                countavgR1 = countavgR1+1;
            elseif strcmp(sublist_side{s},'R') && s==4
                ord = (no_cond+1)-i;
                DATs{s,1}{s2,1}{ord,1} = eegdata;
                %DATs{s,1}{s2,1}{ord,2} = eegdatasing;
                gavgR2 = gavgR2 + DATs{s,1}{s2,1}{ord,1};
                countavgR2 = countavgR2+1;
            end
        end
    end
end
gavgL1 = gavgL1/countavgL1;
gavgR1 = gavgR1/countavgR1;
gavgL2 = gavgL2/countavgL2;
gavgR2 = gavgR2/countavgR2;

gavg = cell(4,1);
gavg{1,1} = gavgL1;
gavg{2,1} = gavgR1;
gavg{3,1} = gavgL2;
gavg{4,1} = gavgR2;

gavg_all = (gavg{1,1}+gavg{2,1}+gavg{3,1}+gavg{4,1})/4;

save([select '_dataFT.mat']);


% normalise by common average and GFP mean over time, for each frequency
for s = 1:length(subjects)
    for f = 1:size(gavg{s,1},3)
        for e = 1:size(gavg{s,1},1)
            gavg{s,1}(e,:,f) = gavg{s,1}(e,:,f) - nanmean(gavg{s,1}(:,:,f),1);
        end
        gavg{s,1}(:,:,f) = gavg{s,1}(:,:,f)/nanmean(nanstd(gavg{s,1}(:,:,f),1));
    end
end

if strcmp(select,'ERP')
    no_freqs=1;
    freq_idx={1};
    freqs_nme = {select};
    freqs_limits=cell(1);
    times=EEG.times;
    posnegpeak = [2 1 1 2 2 1 1 1 1 1 1 1 1 1];
    gavg_sm = gavg;
elseif strcmp(select,'TF')
    
    % cross-correlations
    cormat=[];
    pcamat=[];
    for s = 1:length(subjects)
        for e = 1:size(gavg{s,1},1)
            cormat(e,:,:,s) = corrcoef(squeeze(gavg{s,1}(e,:,:)),'rows','complete');
            %pcamat(e,:,:,s) = pca(squeeze(gavg{s,1}(e,:,:)));
        end
    end
    cormatm = nanmean(cormat,4);
    corav = squeeze(nanmean(cormatm,1));

    % Freq bin identification
    no_freqs=6;
    for threshcor=0.8:0.001:1;
        corth=corav.*(corav>threshcor);
        bins = conncomp(graph(corth));
        if length(unique(bins))==no_freqs
            break
        end
    end

    freq_idx={};
    for f = 1:no_freqs
        freq_idx{f}=find(bins==f);
    end
    %for cc=1:size(corth,2)
    %    if any(corth(:,cc))==0
    %        freq_idx{length(freq_idx)+1}=cc;
    %    elseif any(corth(:,cc))
    %        tempi=find(corth(:,cc)==1);
    %        if length(freq_idx)>1
    %            [te fr] = intersect(tempi,1:max(freq_idx{cc-1}));
    %            tempi(fr)=[];
    %        end
    %        freq_idx{length(freq_idx)+1}=tempi;
    %        if any(tempi==size(corth,1))
    %            break
    %        end
    %    end
    %end

    freqs_nme = cell(no_freqs,1);
    freqs_limits = cell(no_freqs,1);
    for f = 1:no_freqs
        fr = freqs(freq_idx{f});
        C1 = strsplit(num2str(fr(1)),'.');
        if length(fr)>1
            C2 = strsplit(num2str(fr(end)),'.');
            freqs_nme{f} = ['f' C1{1} '_f' C2{1}];
        else
            freqs_nme{f} = ['f' C1{1}];
        end
        freqs_limits{f} = [fr(1) fr(end)];
    end

    % plot grand-grand average
    %close all
    %for f = 1:length(freq_idx)
    %    chanlocs = EEG.chanlocs;
    %    plotdata = mean(gavg_all(:,:,freq_idx{f}),3);
    %    [maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
    %    [x, maxmaxidx] = max(maxval);
    %    plottime = times(maxidx(maxmaxidx));
    %    if plottime == times(end)
    %        plottime = times(end-1);
    %    end
    %    figure('Name',num2str(freqs(freq_idx{f})));
    %    timtopo(plotdata(:,:),chanlocs,...
    %        'limits',[times(1) times(end)],...
    %        'plottimes',plottime);
    %    set(gcf,'Color','white');
    %end
    
    % smooth gavg for peak identification
    gavg_sm = gavg;
    %for s = 1:length(subjects)
    %    for f = 1:size(gavg{s,1},3)
    %        for e = 1:size(gavg{s,1},1)
    %            gavg_sm{s,1}(e,:,f) = smooth(squeeze(gavg{s,1}(e,:,f)),30,'lowess');
    %        end
    %    end
    %end
    figure;
    imagesc(squeeze(gavg_sm{s,1}(7,:,:))'); 
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
stimtimesGFP = find(times>0);
timesRMSE = times(2:end);
stimtimesRMSE = find(timesRMSE>0);

close all
for f = 3:no_freqs
datex = squeeze(gavg_sm{s,1}(7,:,f))';
basetimesGFP = find(times<0 & ~isnan(datex(1:end))');
basetimesRMSE = find(timesRMSE<0 & ~isnan(datex(2:end))');
   
    % GFP and CIs across subjects
    GFPg=[];
    GFP_sub=[];
    for s = 1:length(subjects)
        GFPg(:,s) = nanstd(nanmean(gavg_sm{s,1}(:,:,freq_idx{f}),3),1);
        %for s2 = 1:length(subjects{s,1}) 
        %    GFP_sub(:,s,s2) = squeeze(std(squeeze(mean(subdat_all(:,:,freq_idx{f},s,s2),3)),1));
        %end
    end
    %GFP_CI=[];
    %for t = 1:length(times)
    %    GFPg(t,:) = mean(squeeze(GFP_sub(t,:,:)),2);
    %    GFP_CI(t,:) = 2.575*std(squeeze(GFP_sub(t,:,:))');
    %end
    GFPg=nanmean(GFPg,2);
    GFP_CI = repmat(nanmean(GFPg(basetimesGFP))+nanstd(GFPg(basetimesGFP)),length(GFPg),1);
    GFP_base = repmat(nanmean(GFPg(basetimesGFP)),length(GFPg),1);
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
    for s = 1:length(subjects)
        for t = 2:size(gavg_sm{s,1},2)
            y = nanmean(gavg_sm{s,1}(:,t,freq_idx{f}),3);
            yhat = nanmean(gavg_sm{s,1}(:,t-1,freq_idx{f}),3);
            RMSE(t-1,s) = sqrt(nanmean((y - yhat).^2));
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
    RMSEav = nanmean(RMSE,2);
    RMSE_CI = repmat(nanmean(RMSEav(basetimesRMSE))+nanstd(RMSEav(basetimesRMSE)),length(RMSEav),1);
    RMSE_base = repmat(nanmean(RMSEav(basetimesRMSE)),length(RMSEav),1);
    hold on
    yyaxis right
    plot(times(2:end),RMSEav,'r');
    hold on
    plot(times(2:end),RMSE_CI,'r');
    hold on
    plot(times(2:end),RMSE_base,'r--');
    
    % RMSE peaks
    [peaks_RMSE peaklocs_RMSE] = findpeaks_Rbase(RMSEav,'MinPeakHeight',unique(RMSE_CI),'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    peakloctimes_RMSE = times(peaklocs_RMSE+1);
    hold on;
    plot(peakloctimes_RMSE,peaks_RMSE,'v','MarkerFaceColor','r');
    % RMSE troughs
    [tro_RMSE trolocs_RMSE] = findpeaks_Lbase(-RMSEav,'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    troloctimes_RMSE = times(trolocs_RMSE+1);
    tro_RMSE = -tro_RMSE;
    hold on;
    plot(troloctimes_RMSE,tro_RMSE,'s','MarkerFaceColor','r');
    
    % find peak limits
    peakloctimes_RMSE=peakloctimes_RMSE(peakloctimes_RMSE>timesRMSE(stimtimesRMSE(1)));
    limit1=[];
    limit2=[];
    for i = 1:length(peakloctimes_RMSE)
        limit1=[limit1 peakloctimes_RMSE(i)];
        next_tros = troloctimes_RMSE(find(troloctimes_RMSE>peakloctimes_RMSE(i)));
        limit2=[limit2 next_tros(1)];
    end
    limits{f}=[limit1;limit2];
    no_peaks=size(limits{f},2);
    
    figure
    for s = 1:length(subjects)
        for n = 1:no_peaks
            subplot(length(subjects),no_peaks,(s-1)*no_peaks+n); topoplot(mean(mean(gavg{s,1}(:,find(times==limits{f}(1,n)):find(times==limits{f}(2,n)),freq_idx{f}),3),2), chanlocs);
        end
    end
end
timefreq_limits = horzcat(freqs_limits,limits);
save(['timefreq_limitsFT_' select],'timefreq_limits');

% identify GFPs for each subject and condition
%LATgfp = NaN('double');
%GFP = NaN('double');
ugrplist = unique(sublist_grp','stable');
grplistall = reshape(repmat(ugrplist,1,size(sublist,2))',size(sublist,2)*length(ugrplist),1);
usublist = unique(sublist','stable');
usub = reshape(usublist,length(usublist)/length(ugrplist),length(ugrplist));
GFP = struct('subject',usublist,'group',grplistall);
LAT = GFP;
POW = GFP;
TOPdist = GFP;
TOPTdist = GFP;
W = [];

%avcond = {[1 5],[2 3 4]};
avcond = {1,2,3,4,5};
   
for g = 1:length(ugrplist)
    for s = 1:size(sublist,2) 
        s_name = usub{s,g};
        for f = 1:no_freqs
            f_name = freqs_nme{f};
            no_peaks = size(limits{f},2);
            for p = 1:no_peaks
                p_name = ['peak' num2str(p)];
                lim = limits{f}(:,p); % in datapoints
                for h = 1:length(hand_nme)
                    h_name = hand_nme{h};
                    top=[];
                    topt=[];
                    for ac = 1:length(avcond)
                        c_name = ['c' num2str(ac)];
                        for i = 1:length(avcond{ac}) % for selected conditions
                            DAT(:,:,:,i) = DATs{(g-1)*2+h,1}{s,1}{avcond{ac}(i),1};
                        end
                        DAT = mean(DAT,4);

                        GFPr = std(mean(DAT(:,find(times==lim(1)):find(times==lim(2)),freq_idx{f}),3),1);
                        [~, GFPmaxLOC] = max(GFPr);
                        GFPr = sort(GFPr(:), 'descend');

                        POWr = mean(DAT(:,find(times==lim(1)):find(times==lim(2)),freq_idx{f}),3);
                        POWr = sort(POWr(:), 'descend');

                        fieldname = [f_name '_' p_name '_' h_name '_' c_name];
                        percent_mean = 1;
                        GFP(strcmp(s_name,{GFP.subject})).(fieldname) = mean(GFPr(1:ceil(length(GFPr)*percent_mean)));%mean(GFPr);
                        LAT(strcmp(s_name,{LAT.subject})).(fieldname) = times(find(times==lim(1))+GFPmaxLOC-1);
                        POW(strcmp(s_name,{POW.subject})).(fieldname) = mean(POWr(1:floor(length(POWr)*percent_mean)));
                        
                        elec_radius = 0.35;
                        elec_select = [EEG.chanlocs.radius]<elec_radius;
                        top(:,ac) = mean(mean(DAT(elec_select,find(times==lim(1)):find(times==lim(2)),freq_idx{f}),3),2);
                        topt(:,:) = mean(DAT(elec_select,find(times==lim(1)):find(times==lim(2)),freq_idx{f}),3);
                        TOPTdist(strcmp(s_name,{TOPdist.subject})).(fieldname) = pdist(topt','cosine');
                        
                        if g==1 && s==1; W = [W; f,p,h,ac];end
                    end
                    fieldname = [f_name '_' p_name '_' h_name];
                    TOPdist(strcmp(s_name,{TOPdist.subject})).(fieldname) = pdist(top','cosine');
                end
            end
        end
    end
end

T_GFP = struct2table(GFP);
T_LAT = struct2table(LAT);
T_POW = struct2table(POW);

writetable(T_GFP,['gfp_resultsFT_' select '.xlsx']);
writetable(T_LAT,['lat_resultsFT_' select '.xlsx']);
writetable(T_POW,['pow_resultsFT_' select '.xlsx']);

% modify table for RM ANOVA 
T_GFP.group = categorical(T_GFP.group);
T_LAT.group = categorical(T_LAT.group);
T_POW.group = categorical(T_POW.group);
ys = {};
for ysi = 1:length(T_GFP.Properties.VariableNames(3:end))
    ys{ysi} = ['y' num2str(ysi)];
end
T_GFP.Properties.VariableNames(3:end) = ys;
T_LAT.Properties.VariableNames(3:end) = ys;
T_POW.Properties.VariableNames(3:end) = ys;


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
for f = 1:no_freqs
    f_name = freqs_nme{f};
    no_peaks = size(limits{f},2);
    for p = 1:no_peaks
        p_name = ['peak' num2str(p)];
        lim = limits{f}(:,p); % in datapoints
        for h = 1:length(hand_nme)
            h_name = hand_nme{h};
            fieldname = [f_name '_' p_name '_' h_name]
            for n = 1:13
                topmatH(n,:) = TOPdist(n).(fieldname);
            end
            for n = 14:26
                topmatP(n,:) = TOPdist(n).(fieldname);
            end
            H = squareform(mean(topmatH,1))
            P = squareform(mean(topmatP,1))
        end
    end
end

