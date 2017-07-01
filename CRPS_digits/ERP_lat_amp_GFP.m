clear all
close all
%grplist = [1 3 10 12]; sublist_side = {'L','R','L','R'}; %flipped
%grplist = [33 34 31 32]; sublist_side = {'L','R','L','R'}; %Left, right, OR mostly left, mostly right %unflipped
%grplist = [1 2 29 30]; sublist_side = {'L','R','L','R'}; %Exp2
%grplist = [35 36 37 38]; sublist_side = {'L','R','L','R'}; %Exp1 left v right
grplist = [39 40 41 42]; sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
filepath = 'W:\Data\CRPS_Digit_Perception_exp1\correcttrials\';
cd(filepath)
ele_left = [];
ele_right = [];
%ele_left = [29 23 24 28 30 33 34];
%ele_right = [78 66 70 77 79 83 84];

aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1];
%aff_side = [1 1 1 1 1 1 1 1 1 1 1 1 1];
peakpol = [2 1 2 1 1 1];
range_fig=0;

peaks_nme = {'1','2','3','4'};
%peaks_wind = {20,48;
%              48,88;
%              88,132;
%              168,228;
%              280,340;
%              552,748;};
peaks_wind = {20,48;
              48,88;
              88,128;
              128,272;
              272,352;}; % for exp2: 252,500; for exp1: 252,352;
 peaks_lat = {40,20;
              88,20;
              132,40;
              268,40;}; % for exp2: 252,500; for exp1: 252,352;
cond_nme = {'little','ring','middle','index','thumb'};
elegrp_nme = {'Lpos','Rpos','Lneg','Rneg'};
hand_nme = unique(sublist_side, 'stable');

%sublist = [1 2];
%sublist_side = {'L','R'}; %Left, right
dp = 250;
el = 92;
latvar = 20;
no_ele = 5; % number of electrodes to use to identify latencies
no_cond = length(cond_nme); % no of conditions per data file (arm)
no_peaks = length(peaks_nme);
no_elegrp = length(elegrp_nme);
end_search = 360; % time (ms) to use as end limit for final peak 


loadsubj
subjects = subjlists(grplist);
ERPs = subjects;
gavgL1 = zeros(el,dp);
gavgR1 = zeros(el,dp);
countavgL1 = 0;
countavgR1 = 0;
gavgL2 = zeros(el,dp);
gavgR2 = zeros(el,dp);
countavgL2 = 0;
countavgR2 = 0;

% create ERPs and grand average
sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        sublist{s,s2} = subj(1:3);
        EEG = pop_loadset([subj '.set'],filepath);
        dsize = size(EEG.data);
        %EEG = pop_eegfilt(EEG,0,25,0,0);
        EEG.data = reshape(EEG.data,dsize(1),dsize(2),dsize(3));
        %EEG.data = iirfilt(EEG.data,250,0,25);
        events = zeros(1,size(EEG.data,3));
        %j=0;
        for i = 1:size(EEG.data,3)
            %if strcmp(EEG.event(i).type,'STIM')
                %j = j+1;
             if strcmp(filepath,'W:\Data\CRPS_Digit_Perception\');
                  events(1,i) = EEG.epoch(i).eventcodes{2,2};
             else
                 if length(EEG.epoch(i).eventcodes)==1 || length(EEG.epoch(i).eventcodes)==3
                    events(1,i) = EEG.epoch(i).eventcodes{1,1}{2,2};
                 elseif length(EEG.epoch(i).eventcodes)==2 || length(EEG.epoch(i).eventcodes)==4
                    events(1,i) = EEG.epoch(i).eventcodes{1,2}{2,2};
                 end
            end
        end
        cond = min(events):max(events);
        ERPs{s,1}{s2,1} = cell(length(cond),1);
        for i = 1:no_cond
            if strcmp(sublist_side{s},'L') && s==1
                ord = i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgL1 = gavgL1 + ERPs{s,1}{s2,1}{ord,1};
                countavgL1 = countavgL1+1;
            elseif strcmp(sublist_side{s},'L') && s==3
                ord = i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgL2 = gavgL2 + ERPs{s,1}{s2,1}{ord,1};
                countavgL2 = countavgL2+1;
            elseif  strcmp(sublist_side{s},'R') && s==2
                ord = (no_cond+1)-i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgR1 = gavgR1 + ERPs{s,1}{s2,1}{ord,1};
                countavgR1 = countavgR1+1;
            elseif strcmp(sublist_side{s},'R') && s==4
                ord = (no_cond+1)-i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgR2 = gavgR2 + ERPs{s,1}{s2,1}{ord,1};
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

gavg_all = (gavgL1+gavgL2+gavgR1+gavgR2)/4;

% plot grand average
chanlocs = EEG.chanlocs;
for ii = 1:length(gavg);
    plotdata = gavg{ii,1};
    [maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
    [~, maxmaxidx] = max(maxval);
    plottime = EEG.times(maxidx(maxmaxidx));
    if plottime == EEG.times(end)
        plottime = EEG.times(end-1);
    end
    figure('Name','grand average');
    timtopo(plotdata,chanlocs,...
        'limits',[EEG.times(1) EEG.times(end)],...
        'plottimes',plottime);
    set(gcf,'Color','white');
end

% identify limits for each peak
for ii = 1:length(gavg);
   %gavgtmp = gavg_all;
   gavgtmp = gavg{ii,1};

    mm_gavg = std(gavgtmp,1);

    [Maxima,MaxIdx] = findpeaks(mm_gavg);
    DataInv = 1.01*max(mm_gavg) - mm_gavg;
    [Minima,MinIdx] = findpeaks(DataInv);
    Minima = mm_gavg(MinIdx);
    poststim=EEG.times(MinIdx)>0;
    MinIdx = MinIdx(poststim);
    Minima = Minima(poststim);
    %figure
    %plot(EEG.times',mm_gavg);
    %hold on;
    %scatter(EEG.times(MinIdx),Minima);
    
    figure
    loc = zeros(1,no_peaks-1);
    peak = loc;
    for i=1:no_peaks-1
        peakdata = -mm_gavg(:,find(EEG.times==peaks_wind{i,1}):find(EEG.times==peaks_wind{i,2}));
        [peaks,locs] = max(peakdata);
        loc(1,i) = EEG.times(find(EEG.times==peaks_wind{i,1})-1+locs);
        peak(1,i) = peaks;
    end
    plot(EEG.times',mm_gavg);
    hold on;
    scatter(loc,-peak);

    if loc(1)~=0; loc = [0 loc]; end
    if loc(end)~=EEG.times(end);loc = [loc EEG.times(end)]; end
    for i = 1:length(loc)-1
        limit{i,1} = [loc(i) loc(i+1)];
    end
    limits{ii,1} = limit;
end

% identify GFP peaks within each limit

for ii = 1:length(gavg);
    limit = limits{ii,1};
    clear peaks locs
    peakdata = gavg{ii,1};
    
    for i = 1:no_peaks
        lim = limit{i,1};
        range = find(EEG.times==lim(1)):find(EEG.times==lim(2));
        GFP = std(peakdata(:,range),1);
        peakGFP(ii,i) = range(find(GFP==max(GFP)));
    end
end

% plot electrodes for each peak
%for ii = 1:length(gavg);
%    peakdata = gavg{ii,1};
%    figure
%    subplot(1,6,1); topoplot(peakdata(:,peakGFP(ii,1)), chanlocs,'maplimits',[-1 1]);
%    subplot(1,6,2); topoplot(peakdata(:,peakGFP(ii,2)), chanlocs,'maplimits',[-5 5]);
%    subplot(1,6,3); topoplot(peakdata(:,peakGFP(ii,3)), chanlocs,'maplimits',[-5 5]);
%    subplot(1,6,4); topoplot(peakdata(:,peakGFP(ii,4)), chanlocs,'maplimits',[-3 3]);
%    subplot(1,6,5); topoplot(peakdata(:,peakGFP(ii,5)), chanlocs,'maplimits',[-7 7]);
%    subplot(1,6,6); topoplot(peakdata(:,peakGFP(ii,6)), chanlocs,'maplimits',[-4 4]);
%end

% identify GFPs for each subject and condition
LAT = NaN('double');
LATgfp = NaN('double');
AMP = NaN('double');
GFP = NaN('double');
ranges = NaN('double');

for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        affs = aff_side(s2);
        for i = 1:length(ERPs{s,1}{s2,1}) % for 5 conditions
            ERP = ERPs{s,1}{s2,1}{i,1};
            plot(EEG.times',std(ERP(:,:),1));
            pause
            %figure
            %peakdata = ERP;
            %plot(EEG.times',peakdata);hold on;
            for p = 1:no_peaks % peak
                %ppol = peakpol(p);
                lim = limits{s,1}{p,1};
                srange = latvar / (1000/dp); %ceil(length(find(EEG.times==lim(1)):find(EEG.times==lim(2)))/8);
                range = find(EEG.times==peaks_lat{p,1})-(peaks_lat{p,2}/(1000/dp)):find(EEG.times==peaks_lat{p,1})+(peaks_lat{p,2}/(1000/dp));%peakGFP(s,p)-srange:peakGFP(s,p)+srange; % CHANGED TO USE FIXED LATENCIES!!!
                GFPr = std(ERP(:,range),1);
                GFPmax = max(GFPr);
                GFPdesc = sort(GFPr,'descend');
                GFP(s,s2,i,p) = mean(GFPdesc(1:floor(length(GFPdesc)*0.5))); % select top 50% to average
                LATgfp(s,s2,i,p) = EEG.times(range(find(GFPr==GFPmax)));
            end
  
        end
    end
end


% create results tables for export to excel / SPSS
usublist = unique(sublist','stable');
nrow = length(usublist);
ncol = 2*size(GFP,3)*size(GFP,4); % hand, condition, peak
nrowhead = 1;
ncolhead = 1;
gfp_results = cell(nrow+ncolhead,ncol+nrowhead);
gfp_results(ncolhead+1:end,1) = usublist;
for p = 1:no_peaks
    for c = 1:no_cond
        for h = 1:2
            ind = ((p-1)*no_cond*2 + (h-1)*no_cond + c);
            gfp_results{1,nrowhead+ind} = [peaks_nme{p} '_' hand_nme{h} '_' cond_nme{c}];
        end
    end
end
gfp_lat_results = gfp_results;
GFP2 = reshape(GFP,size(GFP,1)/2,size(GFP,1)/2,size(GFP,2),no_cond,no_peaks); % hand, group, subject, cond, peaks
GFP3 = permute(GFP2,[3 2 4 1 5]); % subject, group, cond, hand,  peaks 
GFP4 = reshape(GFP3,size(GFP3,1)*size(GFP3,2),2*no_cond*no_peaks);
gfp_results(ncolhead+1:end,nrowhead+1:end) = num2cell(GFP4);

LATgfp3 = reshape(LATgfp,size(LATgfp,1)/2,size(LATgfp,1)/2,size(LATgfp,2),no_cond,no_peaks); % hand, group, subject, cond, peaks
LATgfp4 = permute(LATgfp3,[3 2 4 1 5]); % subject, group, peaks,cond, hand,
LATgfp5 = reshape(LATgfp4,size(LATgfp4,1)*size(LATgfp4,2),2*no_cond*no_peaks);
gfp_lat_results(ncolhead+1:end,nrowhead+1:end) = num2cell(LATgfp5);

save gfp_results gfp_results
save gfp_lat_results gfp_lat_results
