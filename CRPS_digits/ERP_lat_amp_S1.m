clear all
close all
%grplist = [1 3 10 12]; sublist_side = {'L','R','L','R'}; %flipped
%grplist = [33 34 31 32]; sublist_side = {'L','R','L','R'}; %Left, right, OR mostly left, mostly right %unflipped
%grplist = [1 29 2 30]; sublist_side = {'L','L','R','R'}; %sources
grplist = [35 36 37 38]; sublist_side = {'L','R','L','R'}; %Exp1 left v right
%grplist = [39 40 41 42]; sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
filepath = 'W:\Data\CRPS_Digit_Perception_exp1\';
ele_left = [29 23 24 28 30 33 34];
ele_right = [78 66 70 77 79 83 84];

peaks_nme = {'1','2','3','4','5','6','7'};
%peaks_wind = {20,48;
%              48,88;
%              88,132;
%              168,228;
%              280,340;
%              552,748;};
peaks_wind = {20,48;
              48,88;
              88,120;
              152,200;
              252,352; % for exp2: 200,352; for exp1: 252,352;
              552,748;}; 
cond_nme = {'little','ring','middle','index','thumb'};
elegrp_nme = {'Lpos','Rpos','Lneg','Rneg'};
hand_nme = unique(sublist_side, 'stable');

%sublist = [1 2];
%sublist_side = {'L','R'}; %Left, right
dp = 250;
el = 92;
no_ele = 7; % number of electrodes to use to identify latencies
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
        events = zeros(1,size(EEG.data,3));
        %j=0;
        for i = 1:size(EEG.data,3)
            %if strcmp(EEG.event(i).type,'STIM')
                %j = j+1;
             if strcmp(filepath,'W:\Data\CRPS_Digit_Perception\');
                  events(1,i) = EEG.epoch(i).eventcodes{2,2};
             elseif strcmp(filepath,'W:\Data\CRPS_Digit_Perception_exp1\');
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

% plot grand average
chanlocs = EEG.chanlocs;
for ii = 1:length(gavg);
    plotdata = gavg{ii,1};
    [maxval, maxidx] = min(abs(plotdata(ele_right,130:230)),[],2);
    [~, maxmaxidx] = min(maxval);
    plottime = EEG.times(maxidx(maxmaxidx));
    if plottime == EEG.times(end)
        plottime = EEG.times(end-1);
    end
    figure('Name','grand average');
    timtopo(plotdata(ele_right,:),chanlocs(ele_right),...
        'limits',[EEG.times(1) EEG.times(end)],...
        'plottimes',plottime);
    set(gcf,'Color','white');
end

% identify limits for each peak
for ii = 1:length(gavg);
    gavgtmp = gavg{ii,1};
    %start_time = 0;
    %end_time = 200;

    for i = 1:dp
        mm_gavg(1,i) = max(gavgtmp(:,i)) - min(gavgtmp(:,i));
    end

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

% identify top x electrodes for each peak from grand average, for each
% polarity


for ii = 1:length(gavg);
    limit = limits{ii,1};
    clear peaks locs
    peakdata = gavg{ii,1};

    for i = 1:no_peaks
        lim = limit{i,1};
        range = find(EEG.times==lim(1)):find(EEG.times==lim(2));
        
        if ~isempty(ele_left) && ~isempty(ele_right)
            e_sel = ele_left;
            ele{ii,1}{i,1}=e_sel;
            [M I] = max(median(peakdata(e_sel,range),1));
            ele{ii,1}{i,2}=range(I);
            ele{ii,1}{i,3}=chanlocs(e_sel); 
            ele{ii,3}{i,1}=e_sel;
            [M I] = min(median(peakdata(e_sel,range),1));
            ele{ii,3}{i,2}=range(I);
            ele{ii,3}{i,3}=chanlocs(e_sel);
            
            e_sel = ele_right;
            ele{ii,2}{i,1}=e_sel;
            [M I] = max(median(peakdata(e_sel,range),1));
            ele{ii,2}{i,2}=range(I);
            ele{ii,2}{i,3}=chanlocs(e_sel);
            ele{ii,4}{i,1}=e_sel;
            [M I] = min(median(peakdata(e_sel,range),1));
            ele{ii,4}{i,2}=range(I);
            ele{ii,4}{i,3}=chanlocs(e_sel);
        else
            e_peaks = zeros(el,1); 
            for e = 1:el
                [pospeak] = max(peakdata(e,range)); % greatest amplitude for each electrode, i.e. peak
                [negpeak] = min(peakdata(e,range)); % lowest amplitude for each electrode, i.e. trough
                if abs(pospeak) > abs(negpeak); e_peaks(e,1) = pospeak; else e_peaks(e,1) = negpeak;end
                e_peaks(e,2) = find(peakdata(e,:)==e_peaks(e,1));
            end

            [e_peaks_max ei] = sort(e_peaks(:,1),'descend'); % 
            e_locs_max = e_peaks(:,2);
            e_locs_max = e_locs_max(ei);
            e_max = e_peaks_max(1:no_ele);
            e_locs_max = e_locs_max(1:no_ele);
            [e_peaks_min ei] = sort(e_peaks(:,1),'ascend');
            e_locs_min = e_peaks(:,2);
            e_locs_min = e_locs_min(ei);
            e_min = e_peaks_min(1:no_ele);
            e_locs_min = e_locs_min(1:no_ele);

            e_sel=[];

            for j = 1:length(e_max)
                e_sel(1,j) = find(e_peaks(:,1)==e_max(j));
                ele{ii,1}{i,2}=median(e_locs_max);
                ele{ii,1}{i,4}=median(e_max);
            end
            ele{ii,1}{i,1}=e_sel;
            ele{ii,1}{i,3}=chanlocs(e_sel);

            for j = 1:length(e_min)
                e_sel(1,j) = find(e_peaks(:,1)==e_min(j));
                ele{ii,2}{i,2}=median(e_locs_min);
                ele{ii,2}{i,4}=median(e_min);
            end
            ele{ii,2}{i,1}=e_sel;
            ele{ii,2}{i,3}=chanlocs(e_sel);



        end
    end
end


% plot electrodes for each peak
for ii = 1:length(gavg);
    peakdata = gavg{ii,1};
    figure
    subplot(4,7,1); topoplot(zeros(1,no_ele), ele{ii,3}{1,3});
    subplot(4,7,2); topoplot(zeros(1,no_ele), ele{ii,3}{2,3});
    subplot(4,7,3); topoplot(zeros(1,no_ele), ele{ii,3}{3,3});
    subplot(4,7,4); topoplot(zeros(1,no_ele), ele{ii,3}{4,3});
    subplot(4,7,5); topoplot(zeros(1,no_ele), ele{ii,3}{5,3});
    subplot(4,7,6); topoplot(zeros(1,no_ele), ele{ii,3}{6,3});
    subplot(4,7,7); topoplot(zeros(1,no_ele), ele{ii,3}{7,3});
    subplot(4,7,8); topoplot(peakdata(:,ele{ii,3}{1,2}), chanlocs,'maplimits',[-1 1]);
    subplot(4,7,9); topoplot(peakdata(:,ele{ii,3}{2,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,10); topoplot(peakdata(:,ele{ii,3}{3,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,11); topoplot(peakdata(:,ele{ii,3}{4,2}), chanlocs,'maplimits',[-3 3]);
    subplot(4,7,12); topoplot(peakdata(:,ele{ii,3}{5,2}), chanlocs,'maplimits',[-7 7]);
    subplot(4,7,13); topoplot(peakdata(:,ele{ii,3}{6,2}), chanlocs,'maplimits',[-4 4]);
    subplot(4,7,14); topoplot(peakdata(:,ele{ii,3}{7,2}), chanlocs,'maplimits',[-7 7]);
    subplot(4,7,15); topoplot(zeros(1,no_ele), ele{ii,4}{1,3});
    subplot(4,7,16); topoplot(zeros(1,no_ele), ele{ii,4}{2,3});
    subplot(4,7,17); topoplot(zeros(1,no_ele), ele{ii,4}{3,3});
    subplot(4,7,18); topoplot(zeros(1,no_ele), ele{ii,4}{4,3});
    subplot(4,7,19); topoplot(zeros(1,no_ele), ele{ii,4}{5,3});
    subplot(4,7,20); topoplot(zeros(1,no_ele), ele{ii,4}{6,3});
    subplot(4,7,21); topoplot(zeros(1,no_ele), ele{ii,4}{7,3});
    subplot(4,7,22); topoplot(peakdata(:,ele{ii,4}{1,2}), chanlocs,'maplimits',[-1 1]);
    subplot(4,7,23); topoplot(peakdata(:,ele{ii,4}{2,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,24); topoplot(peakdata(:,ele{ii,4}{3,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,25); topoplot(peakdata(:,ele{ii,4}{4,2}), chanlocs,'maplimits',[-3 3]);
    subplot(4,7,26); topoplot(peakdata(:,ele{ii,4}{5,2}), chanlocs,'maplimits',[-7 7]);
    subplot(4,7,27); topoplot(peakdata(:,ele{ii,4}{6,2}), chanlocs,'maplimits',[-4 4]);
    subplot(4,7,28); topoplot(peakdata(:,ele{ii,4}{7,2}), chanlocs,'maplimits',[-7 7]);
end

peakselect = 6;
elegrpselect = 1;
figure
col = {'b','b--','r','r--'};
same_elegrp = 2;
times = find(((EEG.times >= -100).*(EEG.times <= 500))==1);
for ii = [1 3]; 
    peakdata = gavg{ii,1};
    if same_elegrp>0
        ERP = mean(peakdata(ele{same_elegrp,elegrpselect}{peakselect,1},:),1);
    else
        ERP = mean(peakdata(ele{ii,elegrpselect}{peakselect,1},:),1);
    end
    plot(EEG.times(times),ERP(times),col{ii}); hold on;
end
%legend('Healthy controls', 'healthy unaff', 'CRPS patients','patient unaff')
legend('Healthy controls', 'CRPS patients','patient unaff');
ylabel('Amplitude, uV') % label for y axis
xlabel('Time (ms)') % label for x axis
set(gca,'FontSize',15);
set(findall(gcf,'type','text'),'FontSize',15);

% identify peaks for each subject and condition
LAT = NaN('double');
AMP = NaN('double');
GFP = NaN('double');

for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        for i = 1:length(ERPs{s,1}{s2,1}) % for 5 conditions
            ERP = ERPs{s,1}{s2,1}{i,1};
            %figure
            %peakdata = ERP;
            %plot(EEG.times',peakdata);hold on;
            for p = 1:length(limits{1,1}) % peak
                lim = limits{s,1}{p,1};
                range = find(EEG.times==lim(1)):find(EEG.times==lim(2)); 
                for iv=1:2 % left or right elegrp
                    for iii = 1:2 % pos or neg elegrp
                        ii = (2*iii-2)+iv;
                        ERPer = mean(ERP(ele{s,ii}{p,1},range),1);
                        GFPer = std(ERP(:,range),1);
                        if iii==1
                            [pospeak locs] = max(ERPer);
                            actpeak = pospeak; %list, subject, condition, peak, EleGroup
                        else
                            [negpeak locs] = min(ERPer);
                            actpeak = negpeak;
                        end
                        [gfppeak locs] = max(GFPer);
                        LAT(s,s2,i,p,iv,iii) = EEG.times(range(1)-1+find(ERPer==actpeak));  % iii = pos/neg; iv = left (1,3)/right(2,4) ele_group
                        %AMP(s,s2,i,p,iv,iii) = squeeze(mean(ERPer(1,max(1,(find(ERPer==actpeak)-floor(((lim(2)-lim(1))/4)/8))):min(length(ERPer),(find(ERPer==actpeak)+floor(((lim(2)-lim(1))/4)/8)))),2)); % range divided by 8
                        AMP(s,s2,i,p,iv,iii) = squeeze(mean(ERPer(1,max(1,(find(GFPer==gfppeak)-floor(12/4))):min(length(ERPer),(find(GFPer==gfppeak)+floor(12/4)))),2)); % range divided by 8
                        GFP(s,s2,i,p) = squeeze(mean(GFPer(1,max(1,(find(GFPer==gfppeak)-floor(12/4))):min(length(GFPer),(find(GFPer==gfppeak)+floor(12/4)))),2)); % range divided by 8
                        
                        %plot(EEG.times, mean(ERP(ele{ii,iii}{p,1},:),1)); 
                        %hold on;
                        %scatter(LAT(s,s2,i,p,iv,iii),AMP(s,s2,i,p,iv,iii));
                        %pause;
                        %close all;
                        
                        %topoplot(ERP(:,find(EEG.times==LAT(s,s2,i,p,iv,iii))), chanlocs);
                        %title(['hand ' num2str(iv) '; peak ' num2str(p)]);
                        %figure
                        %topoplot([0 0 0 0 0], ele{iv,iii}{p,3});
                        %datap = find(EEG.times==LAT(s,s2,i,p,iv,iii));
                        %scatter(LAT(s,s2,i,p,iv,iii)*ones(size(ERP,1),1),ERP(:,datap));
                    end
                    %hold off;             
            
                end
            end
            
        end
        %topoplot(ERP(:,find(EEG.times==LAT(s,s2,i,p,iv,iii))), chanlocs);
        %title(['hand ' num2str(iv) '; peak ' num2str(p)]);
        %pause
        %close all;
    end
end

figure; 
ss= [2 4];
for s2 = 1:length(subjects{s,1}) 
    ERP_right = zeros(92,250);
    ERP_left = zeros(92,250);
    countright = 0;
    countleft = 0;
    for s =ss
        
        
        for i = 1:length(ERPs{s,1}{s2,1})
            ERP = ERPs{s,1}{s2,1}{i,1};
            for p = 4 % peak
                
                for iv=1:2 % hand
                    if strcmp(sublist_side{s},'L'); 
                        ii=1; 
                        lim = limits{ii,1}{p,1};
                        range_left = find(EEG.times==lim(1)):find(EEG.times==lim(2)); 
                        ERP_left = ERP_left+ERP;
                        countleft=countleft+1;
                    elseif strcmp(sublist_side{s},'R') 
                        ii=2; 
                        lim = limits{ii,1}{p,1};
                        range_right = find(EEG.times==lim(1)):find(EEG.times==lim(2)); 
                        ERP_right = ERP_right+ERP;
                        countright=countright+1;
                    end % which side is it?
                    
                end
            end
            
        end
       
    end
    ERP_left = ERP_left/countleft;
    ERP_right = ERP_right/countright;
    ind = s2*2-1;
    subplot(13,2,ind); topoplot(mean(ERP_left(:,range_left),2), chanlocs);
    subplot(13,2,ind+1); topoplot(mean(ERP_right(:,range_right),2), chanlocs);
end

% create results tables for export to excel / SPSS
usublist = unique(sublist','stable');
nrow = length(usublist);
AMP2 = reshape(AMP,size(AMP,1),size(AMP,2),no_cond,no_peaks,no_elegrp);
LAT2 = reshape(LAT,size(LAT,1),size(LAT,2),no_cond,no_peaks,no_elegrp);
ncol = 2*size(AMP2,3)*size(AMP2,4)*size(AMP2,5); % hand, condition, peak, EleGroup
nrowhead = 1;
ncolhead = 1;
amp_results = cell(nrow+ncolhead,ncol+nrowhead);
amp_results(ncolhead+1:end,1) = usublist;
for p = 1:no_peaks
    for e = 1:no_elegrp
        for c = 1:no_cond
            for h = 1:2
                ind = ((p-1)*no_cond*2*no_elegrp + (e-1)*2*no_cond + (h-1)*no_cond + c);
                amp_results{1,nrowhead+ind} = [peaks_nme{p} '_' elegrp_nme{e} '_' hand_nme{h} '_' cond_nme{c}];
            end
        end
     end
end

lat_results = amp_results;

ncol = 2*size(AMP2,3)*size(AMP2,4); % hand, condition, peak
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

AMP3 = reshape(AMP2,size(AMP2,1)/2,size(AMP2,1)/2,size(AMP2,2),no_cond,no_peaks,no_elegrp); % hand, group, subject, cond, peaks, elegrp
AMP4 = permute(AMP3,[3 2 4 1 6 5]); % subject, group, cond, hand, elegrp, peaks 
AMP5 = reshape(AMP4,size(AMP4,1)*size(AMP4,2),2*no_cond*no_peaks*no_elegrp);
amp_results(ncolhead+1:end,nrowhead+1:end) = num2cell(AMP5);

LAT3 = reshape(LAT2,size(LAT2,1)/2,size(LAT2,1)/2,size(LAT2,2),no_cond,no_peaks,no_elegrp);
LAT4 = permute(LAT3,[3 2 6 5 4 1]);
LAT5 = reshape(LAT4,size(LAT4,1)*size(LAT4,2),2*no_cond*no_peaks*no_elegrp);
lat_results(ncolhead+1:end,nrowhead+1:end) = num2cell(LAT5);

GFP2 = reshape(GFP,size(GFP,1)/2,size(GFP,1)/2,size(AMP2,2),no_cond,no_peaks); % hand, group, subject, cond, peaks
GFP3 = permute(GFP2,[3 2 4 1 5]); % subject, group, cond, hand,  peaks 
GFP4 = reshape(GFP3,size(GFP3,1)*size(GFP3,2),2*no_cond*no_peaks);
gfp_results(ncolhead+1:end,nrowhead+1:end) = num2cell(GFP4);


% create latency table for SPM source analysis

latencies = squeeze(mean(LAT4,3)); %mean elegrp
latencies = squeeze(mean(latencies,4)); %mean cond
latencies = permute(latencies,[1 4 2 3]); % subject, group, peaks, hand
latencies = reshape(latencies,size(latencies,1)*size(latencies,2)*2,no_peaks); % group * hand * subject, condition * peak
save latencies_s1.dat latencies -ascii

% create matrices with meaned values

%AMP_allfingers = squeeze(mean(AMP4,5));
%AMP_allfingers_peak1 = squeeze(AMP_allfingers(:,:,:,1,:));
%AMP_allfingers_peak2 = squeeze(AMP_allfingers(:,:,:,2,:));
%AMP_allfingers_peak3 = squeeze(AMP_allfingers(:,:,:,3,:));
%AMP_allfingers_peak4 = squeeze(AMP_allfingers(:,:,:,4,:));%
%AMP_allfingers_peak5 = squeeze(AMP_allfingers(:,:,:,5,:));
%AMP_allfingers_peak1 = reshape(AMP_allfingers_peak1,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak2 = reshape(AMP_allfingers_peak2,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak3 = reshape(AMP_allfingers_peak3,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak4 = reshape(AMP_allfingers_peak4,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak5 = reshape(AMP_allfingers_peak5,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peaks = cat(2,AMP_allfingers_peak1,AMP_allfingers_peak2,AMP_allfingers_peak3,AMP_allfingers_peak4,AMP_allfingers_peak5);
%AMP_peak_group = squeeze(mean(AMP4,3));
%AMP_peak_group = squeeze(mean(AMP_peak_group,4));
%AMP_peak_group = squeeze(mean(AMP_peak_group,4));

AMP_peaks_elegrp = squeeze(mean(AMP4,3));
AMP_peaks_elegrp = squeeze(mean(AMP_peaks_elegrp,3));
AMP_peaks_elegrp = reshape(AMP_peaks_elegrp,size(AMP4,1)*size(AMP4,2),no_peaks*no_elegrp);

AMP_peaks_elegrp_hand = squeeze(mean(AMP4,3));
AMP_peaks_elegrp_hand = reshape(AMP_peaks_elegrp_hand,size(AMP4,1)*size(AMP4,2),size(AMP4,4)*no_peaks*no_elegrp);
ncol = size(AMP_peaks_elegrp_hand,2); % subject*group, hand, elegrp, peaks 
nrowhead = 1;
ncolhead = 1;
amp_results_nocond = cell(nrow+ncolhead,ncol+nrowhead);
amp_results_nocond(ncolhead+1:end,1) = usublist;
for p = 1:no_peaks
    for e = 1:no_elegrp
        %for c = 1:no_cond
            for h = 1:2
                ind = ((p-1)*2*no_elegrp + (e-1)*2 + h);
                amp_results_nocond{1,nrowhead+ind} = [peaks_nme{p} '_' elegrp_nme{e} '_' hand_nme{h}];
            end
        %end
     end
end
amp_results_nocond(ncolhead+1:end,nrowhead+1:end) = num2cell(AMP_peaks_elegrp_hand);


GFP_peaks_hand = squeeze(mean(GFP3,3));
GFP_peaks_hand = reshape(GFP_peaks_hand,size(GFP3,1)*size(GFP3,2),2*no_peaks);
ncol = size(GFP_peaks_hand,2); % subject*group, hand * peaks 
nrowhead = 1;
ncolhead = 1;
gfp_results_nocond = cell(nrow+ncolhead,ncol+nrowhead);
gfp_results_nocond(ncolhead+1:end,1) = usublist;
for p = 1:no_peaks
    %for e = 1:no_elegrp
        %for c = 1:no_cond
            for h = 1:2
                ind = ((p-1)*2 + h);
                gfp_results_nocond{1,nrowhead+ind} = [peaks_nme{p} '_' hand_nme{h}];
            end
        %end
     %end
end
gfp_results_nocond(ncolhead+1:end,nrowhead+1:end) = num2cell(GFP_peaks_hand);
